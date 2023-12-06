WORK_DIR=~/data/2023_chips
PEAKS_DIR=$WORK_DIR/peaks
MODELS_DIR=$WORK_DIR/models
FASTA=$WORK_DIR/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

CHROMOSOME="chr15"
RANGE="chr15:1-101991189"
READS=1000000
PEAKS=500 # Of each type
DISTANCE=5000

MULTS=(0.7 0.5 0.2 0.1);
N=5
THREADS=8

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

PF_NARROW=$(find $PEAKS_DIR -name "H3K4me3*" | grep -v json)
echo "$PF_NARROW"
PF_BROAD=$(find $PEAKS_DIR -name "H3K36me3*" | grep -v json)
echo "$PF_BROAD"

function sample_peaks(){
  PF=$1
  PEAKS=$2
  EXCLUDE=$3
  RESULT=$4
  echo "Original peaks $(wc -l "$PF")"
  echo "Generate $PEAKS random peaks"
  echo "Exclude $EXCLUDE"
  echo "Result $RESULT"
  TF=$(mktemp)
  echo "Excluding peaks"
  bedtools intersect -wa -v -a $PF -b $EXCLUDE > $TF.0
  TP=$(bc -l <<< "$PEAKS * 10")
  echo "Take top $TP peaks by score, to avoid low scores only"
  T=$'\t'
  cat $TF.0 | grep "$CHROMOSOME$T" | sort -k5,5nr | head -n $TP > $TF.1
  echo "Pick $PEAKS random peaks"
  shuf -n $PEAKS $TF.1 | sort -k1,1 -k2,2n > $TF.2
  echo "Make minimal distance $DISTANCE between peaks"
  cat $TF.2 | awk -v OFS='\t' -v D=$DISTANCE \
    'BEGIN {L=-D} {if ($2-L>D) {L=$3; print($1,$2,L,$4,$5);} \
     else {NL=L+D+$3-$2; print($1,L+D,NL,$4,$5); L=NL;};}' >\
     $RESULT
  rm $TF.*
  echo "Generated random peaks $(wc -l "$RESULT")"
}

# Create empty exclude file
EF=$(mktemp)
touch $EF

for I in $(seq 1 $N); do
  NAME="mixed_${CHROMOSOME}_${I}"
  if [[ ! -f ${NAME}.bed ]]; then
    echo "Original peaks $(wc -l "$PF_NARROW") $(wc -l "$PF_BROAD")"
    echo "Generate random peaks file $NAME.bed"
    sample_peaks $PF_NARROW $PEAKS $EF ${NAME}_narrow.bed
    sample_peaks $PF_BROAD $PEAKS ${NAME}_narrow.bed ${NAME}_broad.bed
    cat ${NAME}_narrow.bed > $NAME.bedns
    cat ${NAME}_broad.bed >> $NAME.bedns
    cat $NAME.bedns | sort -k1,1 -k2,2n > $NAME.bed
    rm $NAME.bedns
  fi

  for MULT in "${MULTS[@]}"; do
    RESULT="${NAME}_${MULT}"
    # Use K4me1 as mixed mark
    MF="$MODELS_DIR/H3K4me1_${MULT}.json"

    for NUMCOPIES in 100 1000 10000; do
      echo "Numcopies $NUMCOPIES"
      if [[ ! -f $RESULT.fastq ]]; then
        echo "Model $MF"
        echo "Random peaks $NAME.bed"
        echo "Processing $RESULT"
        echo "Generating $READS reads"
        chips simreads -p $NAME.bed -f $FASTA -o $RESULT -t bed -c 5 --numreads $READS --numcopies $NUMCOPIES \
          --model $MF --region $RANGE --scale-outliers --thread $THREADS
        rm ${RESULT}_*.fastq
      fi
    done
  done
done
