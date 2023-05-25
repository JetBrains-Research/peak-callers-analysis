WORK_DIR=~/data/2023_chips
MODELS_DIR=$WORK_DIR/models

MULTS=(0.5 0.2);

# Generate models with FRIP multipliers
for MF in $(find $MODELS_DIR -name "*.json"); do
   BN=$(basename $MF);
   if [[ -z "$(echo $BN | grep _)" ]]; then
      echo $MF;
      FRIP=$(cat $MF | grep '"s": ' | sed -E 's/,|.*: //g');
      echo $FRIP;
      for MULT in "${MULTS[@]}"; do
        echo $MULT;
        FRIPM=0$(echo "$FRIP * $MULT" | bc -l);
        echo $FRIPM;
        cat $MF | sed "s/$FRIP/$FRIPM/" > ${MF/.json/_$MULT.json};
      done;
      ln -sf $MF ${MF/.json/_1.0.json};
   fi;
done
