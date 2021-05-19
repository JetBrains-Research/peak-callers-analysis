PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

# conda activate span_noise2a

# Generate models with FRIP multipliers
for MF in $(find peaks/ -name "*.json" |  grep -v _0.2.json | grep -v _0.5.json); do
   echo $MF;
   FRIP=$(cat $MF | grep '"s": ' | sed -E 's/,|.*: //g');
   echo $FRIP;
   for M in 0.5 0.2; do
     echo $M;
     FRIPM=0$(echo "$FRIP * $M" | bc -l);
     echo $FRIPM;
     cat $MF | sed "s/$FRIP/$FRIPM/" > ${MF/.json/_$M.json};
  done;
done
