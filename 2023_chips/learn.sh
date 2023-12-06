WORK_DIR=~/data/2023_chips
PEAKS_DIR=$WORK_DIR/peaks
MODELS_DIR=$WORK_DIR/models

# ensure bams are indexed
for F in $WORK_DIR/bams/*.bam; do
    samtools index $F
done

# conda activate chips
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
  BAM=$WORK_DIR/bams/$M.bam
  PF=$(find $PEAKS_DIR -name "$M*" | grep -v chr);
  echo "Learning model for $BAM and $PF"
  chips learn -b $BAM -p $PF -t bed -c 5 --scale-outliers -o $PEAKS_DIR/$M
done

mkdir -p $MODELS_DIR
mv $PEAKS_DIR/*.json $MODELS_DIR