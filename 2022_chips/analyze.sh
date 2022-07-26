WORK_DIR=$(pwd)
PEAKS_DIR=$WORK_DIR/peaks

echo "WORK_DIR $WORK_DIR"
echo "PEAKS_DIR $PEAKS_DIR"

MULTS=(1.0 0.5 0.2 0.1)
N=10

T=$'\t'
echo "Modification${T}Mult${T}Library${T}I${T}TruePeaksFile${T}TruePeaks${T}TrueLength${T}Tool${T}PeaksFile${T}Fdr${T}Peaks${T}Length${T}PrecisionP${T}RecallP${T}Intersection" > report.tsv

for M in mixed H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
  echo "Modification $M"
    for I in $(seq 1 $N); do
      TPF=$WORK_DIR/fastq/${M}_chr15_${I}.bed
      TPF3=$(mktemp)
      cat $TPF | awk -v OFS=$'\t' '{print $1,$2,$3}' > $TPF3
      echo "True peaks file $TPF"
      TL=$(cat $TPF | awk '{L += $3-$2} END {print L}')
      TP=$(cat $TPF | wc -l)
      echo "True peaks $TP length $TL"

      for MULT in "${MULTS[@]}"; do
        NAME="${M}_chr15_${I}_${MULT}"
        echo "Processing $NAME"
        
        echo "MACS2 narrow"
        for F in $(find $WORK_DIR/macs2/ -name "${NAME}*.narrowPeak"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g')
          PR=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          RE=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
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
          PR=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          RE=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
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
          PR=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          RE=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SICER$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""

        echo "SICER noinput"
        for F in $(find $WORK_DIR/sicer/ -name "${NAME}*.scoreisland"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=0.05
          PR=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          RE=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SICER$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""

        echo "SPAN"
        for F in $(find $WORK_DIR/span/ -name "${NAME}_*_3.peak"); do
          echo $F
          P=$(cat $F | wc -l)
          L=$(cat $F | awk '{L+=$3-$2} END {print L}')
          FDR=$(echo $F | sed -E 's/.*00_|_3.peak//g')
          PR=$(bedtools intersect -a $F -b $TPF3 -wa -u | wc -l)
          RE=$(bedtools intersect -a $TPF3 -b $F -wa -u | wc -l)
          INT=$(bedtools intersect -a $TPF3 -b $F | awk '{L+=$3-$2} END {print L}')
          ROW="${M}${T}$MULT${T}$LIB${T}$I${T}$TPF$T${TP}$T${TL}${T}SPAN$T$F${T}${FDR}$T${P}$T${L}$T${PR}$T${RE}$T${INT}"
          echo $ROW
          echo "$ROW" >> report.tsv
        done
        echo ""
      done
  done
done


