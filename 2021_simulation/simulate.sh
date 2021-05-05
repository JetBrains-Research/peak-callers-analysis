FASTA=/mnt/stripe/shpynov/2021_noise1/fa/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

# conda activate bio
for i in {1..20}; do
  echo "Processing chromosome $i"
  samtools faidx $FASTA chr$i > $WORK_DIR/fasta/chr$i.fasta;
  samtools faidx $WORK_DIR/fasta/chr$i.fasta;
done;

T=$'\t'
for PEAKS in encode macs2 sicer span; do
  echo $PEAKS
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    MF=$PEAKS_DIR/$PEAKS/$M.json
    PF=$(find $PEAKS_DIR/$PEAKS -name "$M*" | grep -v json | grep -v chr);
    for i in {15..16}; do
       PF_CHR=$PF.chr$i;
       cat $PF | grep "chr$i$T" | head -n 1000 > $PF_CHR;
#      for READS in 1000000 100000 10000000; do
#        for SPOT in 0.08 0.17;
        FA=$WORK_DIR/fasta/chr$i.fasta
        echo "Simulating $MF chr$i $FA $PF_CHR"
        OF=$WORK_DIR/fastq
        mkdir -p $OF
        cd $OF
        chips simreads -p $PF_CHR -f $FA -o ${M}_${PEAKS}_chr${i}_1mln -t bed -c 5 --numreads 1000000 \
          --model $MF --scale-outliers --seed 42 --thread 24;
    done;
  done;
done;

# Remove incorrect simulated files with suffixes _N
find $WORK_DIR/fastq/ -name "*1mln_*" | xargs rm