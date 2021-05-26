PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

T=$'\t'
echo "Modification${T}PeaksSource${T}Mult${T}Library${T}I${T}TruePeaksFile${T}TruePeaks${T}TrueLength${T}Tool${T}PeaksFile${T}Fdr${T}Peaks${T}Length${T}Precision${T}Recall" >report.tsv

for PEAKS in macs2 sicer; do # encode span; do
  echo "Peaks $PEAKS"
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    echo "Modification $M"
    for LIB in 200k 500k 1mln; do
      for I in {1..5}; do
        echo "Chromosome chr15"
        TPF=$WORK_DIR/fastq/${M}_${PEAKS}_chr15_${I}.bed
        echo "True peaks file $TPF"
        TL=$(cat $TPF | awk '{L += $3-$2} END {print L}')
        TP=$(cat $TPF | wc -l)
        echo "True peaks $TP length $TL"

        MULTS=("" _0.2 _0.5)
        for MULT in "${MULTS[@]}"; do
          echo MULT

          NAME="${M}_${PEAKS}_chr15_${I}${MULT}_${LIB}"
          if [ -z "$MULT" ]; then
            MULT="_1.0"
          fi
          echo "Processing $NAME"

          echo "MACS2 narrow"
          for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.narrowPeak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}Macs2$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

          echo "MACS2 --broad"
          for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.broadPeak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_broad|_peaks.*//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}Macs2Broad$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

          echo "SICER"
          for F in $(find $WORK_DIR/sicer/ -name "${NAME}*-FDR*"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*-FDR//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SICER$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

          echo "SPAN gap 5"
          for F in $(find $WORK_DIR/span/ -name "${NAME}_*_5.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*200_|_5.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-GAP5$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""

          echo "SPAN gap 0"
          for F in $(find $WORK_DIR/span/ -name "${NAME}_*_0.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*200_|_0.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-GAP0$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >> report.tsv
          done
          echo ""


          echo "SPAN NBHMM2NZ"
          for F in $(find $WORK_DIR/nbhmm2nz/ -name "${NAME}*.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-NZ2$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""

          echo "SPAN NBHMM3NZ"
          for F in $(find $WORK_DIR/nbhmm3nz/ -name "${NAME}*.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-NZ3$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""

          echo "SPAN NBHMM4NZ"
          for F in $(find $WORK_DIR/nbhmm4nz/ -name "${NAME}*.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-NZ4$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""

          echo "SPAN NBHMM3Z"
          for F in $(find $WORK_DIR/nbhmm3z/ -name "${NAME}*.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-Z3$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""

          echo "SPAN NBHMM4Z"
          for F in $(find $WORK_DIR/nbhmm4z/ -name "${NAME}*.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-Z4$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""


          echo "SPAN Islands"
          for F in $(find $WORK_DIR/islands/ -name "${NAME}_*.peak"); do
            echo $F
            P=$(cat $F | wc -l)
            L=$(cat $F | awk '{L+=$3-$2} END {print L}')
            FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
            PR=$(bedtools intersect -a $TPF -b $F -wa -u | wc -l)
            RE=$(bedtools intersect -a $F -b $TPF -wa -u | wc -l)
            echo "Fdr $FDR Peaks $P Precision $PR Recall $RE"
            echo "${M}$T${PEAKS}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-Islands$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}" >>report.tsv
          done
          echo ""
        done
      done
    done
  done
done

