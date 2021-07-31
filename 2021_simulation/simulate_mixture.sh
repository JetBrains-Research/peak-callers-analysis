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
      MF_NARROW="$PEAKS_DIR/H3K4me3_${MULT}.json"
      echo "MODEL $MF_NARROW"
      MF_BROAD="$PEAKS_DIR/H3K36me3_${MULT}.json"
      echo "MODEL $MF_BROAD"

      echo "Processing $RESULT 1mln"
      chips simreads -p $TF_NARROW -f $FASTA -o narrow -t bed -c 5 --numreads 500000 \
        --model $MF_NARROW --region $CHR15 --scale-outliers --seed 42 --thread $THREADS
      chips simreads -p $TF_BROAD -f $FASTA -o broad -t bed -c 5 --numreads 500000 \
        --model $MF_BROAD --region $CHR15 --scale-outliers --seed 42 --thread $THREADS
      cat narrow.fastq > ${RESULT}_1mln.fastq
      cat broad.fastq >> ${RESULT}_1mln.fastq

      echo "Processing $RESULT 500k"
      chips simreads -p $TF_NARROW -f $FASTA -o narrow -t bed -c 5 --numreads 250000 \
        --model $MF_NARROW --region $CHR15 --scale-outliers --seed 42 --thread $THREADS
      chips simreads -p $TF_BROAD -f $FASTA -o broad -t bed -c 5 --numreads 250000 \
        --model $MF_BROAD --region $CHR15 --scale-outliers --seed 42 --thread $THREADS
      cat narrow.fastq > ${RESULT}_500k.fastq
      cat broad.fastq >> ${RESULT}_500k.fastq

      echo "Processing $RESULT 200k"
      chips simreads -p $TF_NARROW -f $FASTA -o narrow -t bed -c 5 --numreads 100000 \
        --model $MF_NARROW --region $CHR15 --scale-outliers --seed 42 --thread $THREADS
      chips simreads -p $TF_BROAD -f $FASTA -o broad -t bed -c 5 --numreads 100000 \
        --model $MF_BROAD --region $CHR15 --scale-outliers --seed 42 --thread $THREADS
      cat narrow.fastq > ${RESULT}_200k.fastq
      cat broad.fastq >> ${RESULT}_200k.fastq

      rm narrow.fastq
      rm broad.fastq
    done
  fi
done
