FASTA=/mnt/stripe/shpynov/2021_chips/fa/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
CHR15="chr15:1-101991190"
PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips
N=10

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

T=$'\t'
for PEAKS in macs2; do # sicer encode span; do
  echo $PEAKS
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    PF=$(find $PEAKS_DIR/$PEAKS -name "$M*" | grep -v json)
    echo "PEAKS $PF"

    for I in $(seq 1 $N); do
      NAME="${M}_${PEAKS}_chr15_${I}"
      TF=$(mktemp)
      cat "$PF" | grep chr15 >$TF
      # Pick random peaks
      shuf -n 500 $TF | sort -k1,1 -k2,2n >$NAME.bed
      echo "Random peaks $NAME.bed"

      MULTS=(1.0 0.5 0.2);
      for MULT in "${MULTS[@]}"; do
        MNAME="${NAME}_${MULT}"
        MF="$PEAKS_DIR/$PEAKS/${M}${MULT}.json"
        echo "MODEL $MF"

        echo "Processing $MNAME 1mln"
        if [[ ! -f ${MNAME}_1mln.fastq ]]; then
          chips simreads -p $NAME.bed -f $FASTA -o ${MNAME}_1mln -t bed -c 5 --numreads 1000000 \
            --model $MF --region $CHR15 --scale-outliers --seed 42 --thread 24
        fi

        echo "Processing $MNAME 500k"
        if [[ ! -f ${MNAME}_500k.fastq ]]; then
          chips simreads -p $NAME.bed -f $FASTA -o ${MNAME}_500k -t bed -c 5 --numreads 500000 \
            --model $MF --region $CHR15 --scale-outliers --seed 42 --thread 24
        fi

        echo "Processing $MNAME 200k"
        if [[ ! -f ${MNAME}_200k.fastq ]]; then
          chips simreads -p $NAME.bed -f $FASTA -o ${MNAME}_200k -t bed -c 5 --numreads 200000 \
            --model $MF --region $CHR15 --scale-outliers --seed 42 --thread 24
        fi

      done
    done
  done
done
