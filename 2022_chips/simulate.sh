WORK_DIR=/data/2022_chips
MODELS_DIR=$WORK_DIR/models
FASTA=$WORK_DIR/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

CHROMOSOME="chr15"
RANGE="chr15:1-101991190"
READS=500000

PEAKS1=1000
PEAKS2=500

MULTS=(0.5 0.3 0.1)
N=5
THREADS=24

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

T=$'\t'
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
  PF=$(find $MODELS_DIR -name "$M*" | grep -v json)
  echo "$PF"

  for I in $(seq 1 $N); do
    NAME="${M}_${CHROMOSOME}_${I}"
    if [[ ! -f ${NAME}.bed ]]; then
      echo "Generate random peaks $NAME.bed"
      TF=$(mktemp)
      # Take top peaks by score, to avoid low scores only
      cat "$PF" | grep "$CHROMOSOME$T" | sort -k5,5nr | head -n $PEAKS1 > $TF
      # Pick random peaks
      shuf -n $PEAKS2 $TF | sort -k1,1 -k2,2n > $NAME.bed
    fi;

    for MULT in "${MULTS[@]}"; do
      RESULT="${NAME}_${MULT}"
      MF="$MODELS_DIR/${M}_${MULT}.json"

      if [[ ! -f ${RESULT}.fastq* ]]; then
        echo "Model $MF"
        echo "Random peaks $NAME.bed"
        echo "Processing $RESULT"
        chips simreads -p $NAME.bed -f $FASTA -o ${RESULT} -t bed -c 5 --numreads $READS \
          --model $MF --region $RANGE --scale-outliers --seed 42 --thread $THREADS
      fi
    done
  done
done