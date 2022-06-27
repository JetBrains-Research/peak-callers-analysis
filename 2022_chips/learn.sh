WORK_DIR=/data/2022_chips
PEAKS_DIR=$WORK_DIR/peaks

# conda activate chips
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
  BAM=$WORK_DIR/bams/$M.bam
  PF=$(find $PEAKS_DIR -name "$M*" | grep -v chr);
  echo "Learning model for $BAM and $PF"
  chips learn -b $BAM -p $PF -t bed -c 5 --scale-outliers -o $PEAKS_DIR/$M
done
