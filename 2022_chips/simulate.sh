WORK_DIR=/data/2022_chips
MODELS_DIR=$WORK_DIR/models
FASTA=$WORK_DIR/fasta/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta

CHROMOSOME="chr1"
RANGE="chr1:1-248956000" # Limit at 1kbp before chromosome end
READS=1000000
PEAKS=1000

T=$'\t'
if [[ ! -f hg38.chrom.sizes.$CHROMOSOME ]]; then
    cat hg38.chrom.sizes | grep "$CHROMOSOME$T" | awk -v OFS='\t' '{print $1, $2-10000}'> hg38.chrom.sizes.$CHROMOSOME
fi

MULTS=(1.0 0.5 0.2)
N=1
THREADS=8

OF=$WORK_DIR/fastq
mkdir -p $OF
cd $OF

function sample_peaks(){
  PF=$1
  PEAKS=$2
  GENOME=$3
  RESULT=$4
  echo "Original peaks $(wc -l "$PF")"
  echo "Generate $PEAKS random peaks"
  echo "Genome $GENOME"
  echo "Result $RESULT"
  TF=$(mktemp)
  echo "Pick $PEAKS random peaks"
  shuf -n $PEAKS $PF > $TF
  echo "Launching shuffle"
  bedtools shuffle -noOverlapping -i $TF -g $GENOME | sort -k1,1 -k2,2n > $RESULT
  echo "Generated random peaks $(wc -l "$RESULT")"
  rm $TF
}

for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
  PF=$(find $MODELS_DIR -name "$M*" | grep -v json)
  echo "$PF"

  for I in $(seq 1 $N); do
    NAME="${M}_${CHROMOSOME}_${I}"
    if [[ ! -f ${NAME}.bed ]]; then
      sample_peaks $PF $PEAKS $WORK_DIR/hg38.chrom.sizes.$CHROMOSOME $NAME.bed
    fi

    for MULT in "${MULTS[@]}"; do
      RESULT="${NAME}_${MULT}"
      MF="$MODELS_DIR/${M}_${MULT}.json"

      if [[ ! -f ${RESULT}.fastq* ]]; then
        echo "Model $MF"
        echo "Random peaks $NAME.bed"
        echo "Processing $RESULT"
        echo "Generating $READS reads"
        chips simreads -p $NAME.bed -f $FASTA -o $RESULT -t bed -c 5 --numreads $READS --numcopies 1000 \
          --model $MF --region $RANGE --scale-outliers --seed 42 --thread $THREADS
      fi
    done
  done
done