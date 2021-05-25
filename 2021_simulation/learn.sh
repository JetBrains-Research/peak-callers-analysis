PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

# conda activate span_noise2a
for PEAKS in encode macs2 sicer span; do
  echo $PEAKS
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    BAM=$WORK_DIR/bams/$M.bam
    PF=$(find $PEAKS_DIR/$PEAKS -name "$M*" | grep -v chr);
    echo "Learning model for $BAM and $PF"
    chips learn -b $BAM -p $PF -t bed -c 5 -o $PEAKS_DIR/$PEAKS/$M
  done;
done;
