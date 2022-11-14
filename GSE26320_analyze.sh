echo "Analyze coverage of GSE262320 dataset"
WORK_DIR=$(pwd)
SHIFT=125
BLACKLIST=hg19-blacklist.v2.bed

echo "WORK_DIR $WORK_DIR"
echo "SHIFT $SHIFT"
echo "BLACKLIST $BLACKLIST"

mkdir -p analyze/tags

for BAM in bams/*.bam; do
   echo ${BAM};
   BN=$(basename $BAM);
   TAGS=analyze/tags/${BN/.bam/.tags}
   if [[ ! -f "$TAGS" ]]; then
      bedtools bamtobed -i ${BAM} |\
         awk -v OFS='\t' -v S=${SHIFT} \
         '{if ($6 != "-") {print($1, $2+S, $2+S+1)} else {if ($3-S>=1) {print($1, $3-S, $3-S+1)}}}' |\
         sort -u -k1,1 -k3,3n -k2,2n > ${TAGS}
   fi
done

echo "Compute bed4 files"
mkdir -p analyze/bed4
T=$(mktemp -d);
mkdir -p $T;
for F in $(find macs2 -name "*.narrowPeak") $(find macs2 -name "*.broadPeak") $(find span -name "*.peak") $(find sicer -name "*-summary-FDR*" | grep -v log); do
    echo $F;
    N=$(basename $F);
    cat $F | awk -v OFS='\t' '{print $1,$2,$3,$5}' | sort -k1,1 -k2,2n > $T/$N.unfiltered;
    bedtools intersect -v -a $T/$N.unfiltered -b $BLACKLIST > analyze/bed4/$N.bed4;
done;
rm -r $T;

echo "Add coverage information to peaks"
mkdir -p analyze/covs
CELLS=$(ls analyze/bed4/*.bed4 | sed 's#.*/##g' | sed -E 's/GSM[0-9]+_//g' | sed 's/_.*//g' | sort --unique);
echo $CELLS
for C in ${CELLS[@]}; do
   echo $C;
   for M in CTCF H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
       echo $M;
       for R in rep1 rep2; do
          echo $R;
          for F in analyze/bed4/*${C}_${M}_${R}*.bed4; do
             echo $F;
             echo $(ls analyze/tags/*${C}_${M}_${R}*.tags);
             echo $(ls analyze/tags/*${C}_WCE_${R}*.tags);
             N=$(basename $F);
             bedtools intersect -a $F -b analyze/tags/*${C}_${M}_${R}*.tags -wa -c > analyze/covs/$N.t;
             bedtools intersect -a analyze/covs/$N.t -b analyze/tags/*${C}_WCE_${R}*.tags -wa -c > analyze/covs/$N.tc;
             rm analyze/covs/$N.t;
          done
       done
   done
done