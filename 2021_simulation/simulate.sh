FASTA=/mnt/stripe/shpynov/2021_noise1/fa/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

T=$'\t'
for PEAKS in encode macs2 sicer span; do
  echo $PEAKS
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    MF=$PEAKS_DIR/$PEAKS/$M.json
    PF=$(find $PEAKS_DIR/$PEAKS -name "$M*" | grep -v json);
    for i in {1..10}; do
	NAME=${M}_${PEAKS}_chr15_$i
        TF=$(mktemp);
	cat $PF | grep chr15 > $TF
        # Pick random peaks
        shuf -n 500 $TF | sort -k1,1 -k2,2n > $NAME.bed
        rm $TF
        FA=$WORK_DIR/fasta/chr15.fasta
        echo "Simulating $i $MF chr15 $FA $PF_CHR"
        chips simreads -p $NAME.bed -f $FA -o ${NAME}_1mln -t bed -c 5 --numreads 1000000 \
          --model $MF --scale-outliers --seed 42 --thread 24;
        chips simreads -p $NAME.bed -f $FA -o ${NAME}_100k -t bed -c 5 --numreads 100000 \
          --model $MF --scale-outliers --seed 42 --thread 24;
        chips simreads -p $NAME.bed -f $FA -o ${NAME}_500k -t bed -c 5 --numreads 500000 \
          --model $MF --scale-outliers --seed 42 --thread 24;
    done;
  done;
done;
