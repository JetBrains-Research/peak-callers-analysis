PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips
N=10

T=$'\t'
echo "Modification${T}Mult${T}Library${T}I${T}TruePeaksFile${T}TruePeaks${T}TrueLength${T}Tool${T}PeaksFile${T}Fdr${T}Peaks${T}Length${T}PrecisionP${T}RecallP${T}Intersection" > report.tsv

for M in mixed H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
  echo "Modification $M"
  for LIB in 200k 500k 1mln; do
    for I in $(seq 1 $N); do
      TPF=$WORK_DIR/fastq/${M}_chr15_${I}.bed
      TPF3=$(mktemp)
      cat $TPF | awk -v OFS=$'\t' '{print $1,$2,$3}' > $TPF3
      echo "True peaks file $TPF"
      TL=$(cat $TPF | awk '{L += $3-$2} END {print L}')
      TP=$(cat $TPF | wc -l)
      echo "True peaks $TP length $TL"

      MULTS=(1.0 0.5 0.2)
      for MULT in "${MULTS[@]}"; do
        NAME="${M}_chr15_${I}_${MULT}_${LIB}"
        echo "Processing $NAME"
        
        echo "MACS2 narrow"
        for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.narrowPeak"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g')
          PR=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          RE=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}Macs2$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""
        
        echo "MACS2 --broad"
        for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.broadPeak"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*_broad|_peaks.*//g')
          PR=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          RE=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}Macs2Broad$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""
        
        echo "SICER"
        for F in $(find $WORK_DIR/sicer/ -name "${NAME}*-FDR*"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*-FDR//g')
          PR=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          RE=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SICER$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""

        echo "SPAN gap 5"
        for F in $(find $WORK_DIR/span/ -name "${NAME}_*_5.peak"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*200_|_5.peak//g')
          PR=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          RE=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-GAP5$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""
          
        echo "SPAN NBHMM2NZ"
        for F in $(find $WORK_DIR/nbhmm2nz/ -name "${NAME}*.peak"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*_|\.peak//g')
          PR=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          RE=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN-NZ2$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""

      done
    done
  done
done


