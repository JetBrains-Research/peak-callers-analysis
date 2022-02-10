WORK_DIR=/mnt/stripe/shpynov/2021_chips
MODELS_DIR=$WORK_DIR/models
FASTA=$WORK_DIR/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

CHROMOSOME="chr15"
RANGE="chr15:1-101991190"
READS=500000

PEAKS1=500
PEAKS2=250

MULTS=(1.0 0.5 0.3 0.1)
N=5
THREADS=24

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

T=$'\t'
PF_NARROW=$(find $MODELS_DIR -name "H3K4me3*" | grep -v json)
echo "$PF_NARROW"
PF_BROAD=$(find $MODELS_DIR -name "H3K4me1*" | grep -v json)
echo "$PF_BROAD"

for I in $(seq 1 $N); do
  NAME="mixed_chr15_${I}"
  if [[ ! -f ${NAME}.bed ]]; then
    echo "Generating random peaks $NAME.bed"
    TF=$(mktemp)

    # Take top peaks by score, to avoid low scores only
    cat "$PF_NARROW" | grep "$CHROMOSOME$T" | sort -k5,5nr | head -n $PEAKS1 |\
      awk -v OFS='\t' '{print $1,$2,$3}' > $TF
    # Pick random peaks
    shuf -n $PEAKS2 $TF > ${NAME}_narrow.bed

    TF_BROAD=$(mktemp)
    bedtools subtract -a $PF_BROAD -b ${NAME}_narrow.bed > $TF_BROAD
    # Take top peaks by score, to avoid low scores only
    cat "$TF_BROAD" | grep "$CHROMOSOME$T" | sort -k5,5nr | head -n $PEAKS1 |\
      awk -v OFS='\t' '{print $1,$2,$3}' > $TF
    # Pick random peaks
    shuf -n $PEAKS2 $TF > ${NAME}_broad.bed

    cat ${NAME}_narrow.bed > $NAME.bedns
    cat ${NAME}_broad.bed >> $NAME.bedns
    cat $NAME.bedns | sort -k1,1 -k2,2n > $NAME.bed
    rm $NAME.bedns
  fi

  for MULT in "${MULTS[@]}"; do
    RESULT="${NAME}_${MULT}"
    # Use K4me1 as mixed mark
    MF="$MODELS_DIR/H3K4me1_${MULT}.json"

    if [[ ! -f ${RESULT}.fastq* ]]; then
      echo "Model $MF"
      echo "Random peaks $NAME.bed"
      echo "Processing $RESULT"
      chips simreads -p $NAME.bed -f $FASTA -o ${RESULT} -t bed -c 5 --numreads $READS \
        --model $MF --region $RANGE --scale-outliers --seed 42 --thread $THREADS
    fi
  done
done
