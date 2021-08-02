FASTA=/mnt/stripe/shpynov/2021_chips/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
CHR15="chr15:1-101991190"
PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips
N=10
THREADS=24

OF=$WORK_DIR/fastq_mixed
mkdir -p $OF
cd $OF

T=$'\t'
PF_NARROW=$(find $PEAKS_DIR -name "H3K4me3*" | grep -v json)
echo "$PF_NARROW"
PF_BROAD=$(find $PEAKS_DIR -name "H3K36me3*" | grep -v json)
echo "$PF_BROAD"

for I in $(seq 1 $N); do
  NAME="mixed_chr15_${I}"
  if [[ ! -f ${NAME}.bed ]]; then
    TF=$(mktemp)
    TF_NARROW=$(mktemp)
    TF_BROAD=$(mktemp)
    bedtools subtract -a $PF_NARROW -b $PF_BROAD > $TF_NARROW
    bedtools subtract -b $PF_NARROW -a $PF_BROAD > $TF_BROAD

    # Take top 1000 peaks by score, to avoid low scores only
    cat "$TF_NARROW" | grep chr15 | sort -k5,5nr | head -n 1000 > $TF
    # Pick 250 random peaks
    shuf -n 250 $TF | sort -k1,1 -k2,2n > ${NAME}_narrow.bed

    # Take top 1000 peaks by score, to avoid low scores only
    cat "$TF_BROAD" | grep chr15 | sort -k5,5nr | head -n 1000 > $TF
    # Pick 250 random peaks
    shuf -n 250 $TF | sort -k1,1 -k2,2n > ${NAME}_broad.bed

    cat ${NAME}_narrow.bed > $NAME.bed
    cat ${NAME}_broad.bed >> $NAME.bed
    # Sort peaks
    sort -k1,1 -k2,2n -o $NAME.bed $NAME.bed

    echo "Random peaks $NAME.bed"

    MULTS=(1.0 0.5 0.2);
    for MULT in "${MULTS[@]}"; do
      RESULT="${NAME}_${MULT}"
      # Use K4me1 as mixed mark
      MF="$PEAKS_DIR/H3K4me1_${MULT}.json"
      echo "MODEL $MF"

      echo "Processing $RESULT 1mln"
      chips simreads -p $NAME.bed -f $FASTA -o ${RESULT}_1mln -t bed -c 5 --numreads 1000000 \
        --model $MF --region $CHR15 --scale-outliers --seed 42 --thread $THREADS

      echo "Processing $RESULT 500k"
      chips simreads -p $NAME.bed -f $FASTA -o ${RESULT}_500k -t bed -c 5 --numreads 500000 \
        --model $MF --region $CHR15 --scale-outliers --seed 42 --thread $THREADS

      echo "Processing $RESULT 200k"
      chips simreads -p $NAME.bed -f $FASTA -o ${RESULT}_200k -t bed -c 5 --numreads 200000 \
        --model $MF --region $CHR15 --scale-outliers --seed 42 --thread $THREADS

    done
  fi
done
