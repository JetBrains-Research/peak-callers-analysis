# conda activate span_noise2a
for PEAKS in encode macs2 sicer span; do
  echo $PEAKS
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    BAM=/mnt/stripe/shpynov/2021_noise1/bams/$M.bam
    PF=$(find /mnt/stripe/shpynov/2021_chips/peaks/$PEAKS -name "$M*" | grep -v chr);
    echo "Learning model for $BAM and $PF"
    chips learn -b $BAM -p $PF -t bed -c 5 -o /mnt/stripe/shpynov/2021_chips/peaks/$PEAKS/$M
  done;
done;

# Generate models with 0.1 0.2 and 0.05 FRIP multipliers
for MF in $(find peaks/ -name "*.json" | grep -v _); do
   echo $MF;
   FRIP=$(cat $MF | grep '"f": ' | sed -E 's/,|.*: //g');
   echo $FRIP;
   for M in 0.5 0.2 0.1; do
     echo $M;
     FRIPM=0$(echo "$FRIP * $M" | bc -l);
     echo $FRIPM;
     cat $MF | sed "s/$FRIP/$FRIPM/" > ${MF/.json/_$M.json};
  done;
done
