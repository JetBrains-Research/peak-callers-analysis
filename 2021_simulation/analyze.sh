PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

T=$'\t'
echo "Modification${T}PeaksSource${T}Mult${T}Library${T}I${T}TruePeaksFile${T}TruePeaks${T}Tool${T}PeaksFile${T}Fdr${T}Peaks${T}Precision${T}Recall" >report.tsv

for PEAKS in encode macs2 sicer span; do
  echo "Peaks $PEAKS"
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    echo "Modification $M"
    for LIB in 100k 500k 1mln; do
      for I in {1..10}; do
        echo "Chromosome chr15"
        TRUTH=$WORK_DIR/fastq/${M}_${PEAKS}_chr15_${I}.bed
        echo "True peaks file $TRUTH"
        TP=$(cat $TRUTH | wc -l)
        echo "True peaks $TP"

        MULTS=("" _0.1 _0.2 _0.5)
        for MULT in "${MULTS[@]}"; do
          echo MULT

          NAME="${M}_${PEAKS}_chr15_${I}${MULT}_${LIB}"
          echo "Processing $NAME"

          echo "MACS2 narrow"
          for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.narrowPeak"); do
            echo $F
            NP=$(cat $F | wc -l)
            FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g')
            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}Macs2$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

          echo "MACS2 --broad"
          for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.broadPeak"); do
            echo $F
            NP=$(cat $F | wc -l)
            FDR=$(echo $F | sed -E 's/.*_broad|_peaks.*//g')
            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}Macs2Broad$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

          echo "SPAN gap 5"
          for F in $(find $WORK_DIR/span/ -name "${NAME}_*_5.peak"); do
            echo $F
            NP=$(cat $F | wc -l)
            FDR=$(echo $F | sed -E 's/.*200_|_5.peak//g')
            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}SPAN-GAP5$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

#          echo "SPAN gap 0"
#          for F in $(find $WORK_DIR/span/ -name "${NAME}*_0.peak"); do
#            echo $F
#            NP=$(cat $F | wc -l)
#            FDR=$(echo $F | sed -E 's/.*200_|_0.peak//g')
#            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
#            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
#            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
#            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}SPAN-GAP0$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv
#          done
#          echo ""

          echo "SPAN NBHMM2NZ"
          for F in $(find $WORK_DIR/nbhmm2nz/ -name "${NAME}*.peak"); do
            echo $F
            NP=$(cat $F | wc -l)
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}SPAN-NBHMM2NZ$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""

          echo "SPAN NBHMM3NZ"
          for F in $(find $WORK_DIR/nbhmm3nz/ -name "${NAME}*.peak"); do
            echo $F
            NP=$(cat $F | wc -l)
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}SPAN-NBHMM3NZ$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""

          echo "SPAN Islands"
          for F in $(find $WORK_DIR/islands/ -name "${NAME}_*.peak"); do
            echo $F
            NP=$(cat $F | wc -l)
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l)
            echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TRUTH$T${TP}${T}SPAN-Islands$T$F${T}${FDR}$T${NP}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""
        done
      done
    done
  done
done
