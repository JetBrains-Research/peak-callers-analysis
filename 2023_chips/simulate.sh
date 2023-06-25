WORK_DIR=~/data/2023_chips
MODELS_DIR=$WORK_DIR/models
FASTA=$WORK_DIR/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

CHROMOSOME="chr15"
RANGE="chr15:1-101991189"
READS=1000000
PEAKS=500
DISTANCE=5000

MULTS=(0.7 0.1)
N=5
THREADS=8

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

function sample_peaks(){
  PF=$1
  PEAKS=$2
  RESULT=$3
  echo "Original peaks $(wc -l "$PF")"
  echo "Generate $PEAKS random peaks"
  echo "Result $RESULT"
  TF=$(mktemp)
  TP=$(bc -l <<< "$PEAKS * 10")
  echo "Take top $TP peaks by score, to avoid low scores only"
  T=$'\t'
  cat "$PF" | grep "$CHROMOSOME$T" | sort -k5,5nr | head -n $TP > $TF.1
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

for M in H3K4me3 H3K4me1 H3K27ac H3K27me3 H3K36me3; do
  PF=$(find $MODELS_DIR -name "$M*" | grep -v json)
  echo "$PF"

  for I in $(seq 1 $N); do
    NAME="${M}_${CHROMOSOME}_${I}"
    if [[ ! -f ${NAME}.bed ]]; then
      sample_peaks $PF $PEAKS $NAME.bed
    fi

    for MULT in "${MULTS[@]}"; do
      RESULT="${NAME}_${MULT}"
      MF="$MODELS_DIR/${M}_${MULT}.json"
      # MACS2 fails to build model for smaller numcopies
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
done