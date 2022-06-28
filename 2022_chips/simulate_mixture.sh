WORK_DIR=/data/2022_chips
MODELS_DIR=$WORK_DIR/models
FASTA=$WORK_DIR/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

CHROMOSOME="chr1"
RANGE="chr1:1-248956422"
READS=1000000
PEAKS=500 # Of each type

T=$'\t'
if [[ ! -f hg38.chrom.sizes.$CHROMOSOME ]]; then
    cat hg38.chrom.sizes | grep "$CHROMOSOME$T" > hg38.chrom.sizes.$CHROMOSOME
fi

MULTS=(1.0 0.5 0.2)
N=5
THREADS=16

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

PF_NARROW=$(find $MODELS_DIR -name "H3K4me3*" | grep -v json)
echo "$PF_NARROW"
PF_BROAD=$(find $MODELS_DIR -name "H3K4me1*" | grep -v json)
echo "$PF_BROAD"

function sample_peaks(){
  PF=$1
  PEAKS=$2
  GENOME=$3
  EXCLUDE=$4
  RESULT=$5
  echo "Original peaks $(wc -l "$PF")"
  echo "Generate $PEAKS random peaks"
  echo "Genome $GENOME"
  echo "Exclude $EXCLUDE"
  echo "Result $RESULT"
  TF=$(mktemp)
  echo "Pick $PEAKS random peaks"
  shuf -n $PEAKS $PF > $TF
  echo "Launching shuffle"
  bedtools shuffle -noOverlapping -i $TF -g $GENOME -excl $EXCLUDE | sort -k1,1 -k2,2n > $RESULT
  echo "Generated random peaks $(wc -l "$RESULT")"
  rm $TF
}

# Create empty exclude file
EF=$(mktemp)
touch $EF

for I in $(seq 1 $N); do
  NAME="mixed_${CHROMOSOME}_${I}"
  if [[ ! -f ${NAME}.bed ]]; then
    echo "Original peaks $(wc -l "$PF_NARROW") $(wc -l "$PF_BROAD")"
    echo "Generate random peaks file $NAME.bed"
    sample_peaks $PF_NARROW $PEAKS $WORK_DIR/hg38.chrom.sizes.$CHROMOSOME $EF ${NAME}_narrow.bed
    sample_peaks $PF_BROAD $PEAKS $WORK_DIR/hg38.chrom.sizes.$CHROMOSOME ${NAME}_narrow.bed ${NAME}_broad.bed
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
      echo "Generating $READS reads"
      chips simreads -p $NAME.bed -f $FASTA -o $RESULT -t bed -c 5 --numreads $READS --numcopies 1000 \
        --model $MF --region $RANGE --scale-outliers --seed 42 --thread $THREADS
    fi
  done
done
