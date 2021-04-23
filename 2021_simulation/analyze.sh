PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

T=$'\t'
echo "Modification${T}PeaksSource${T}Chromosome${T}TruePeaks${T}Tool${T}Fdr${T}Peaks${T}Precision${T}Recall" > report.tsv

for PEAKS in encode macs2 sicer span; do
  echo $PEAKS
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    echo $M
    for i in {1..20}; do
      TRUTH=$(find $PEAKS_DIR/$PEAKS -name "$M*" | grep -E "\.chr${i}$");
      echo "True peaks file $TRUTH"; 
      TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; 

      echo "MACS2 narrow"
      for F in $(find $WORK_DIR/macs2/ -name "${M}_${PEAKS}_chr${i}_*.narrowPeak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}Macs2$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;  
      done;
      echo ""
      
      echo "MACS2 --broad"
      for F in $(find $WORK_DIR/macs2/ -name "${M}_${PEAKS}_chr${i}_*.broadPeak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*_broad|_peaks.*//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}Macs2Broad$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SICER is not supported, because it doesn't support FDR without control file"

      echo "SPAN gap 5"
      for F in $(find $WORK_DIR/span/ -name "${M}_${PEAKS}_chr${i}_*_5.peak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*200_|_5.peak//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-GAP5$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""            
      
      echo "SPAN gap 0"
      for F in $(find $WORK_DIR/span/ -name "${M}_${PEAKS}_chr${i}_*_0.peak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*200_|_0.peak//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-GAP0$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""            

      echo "SPAN NBHMM2NZ"
      for F in $(find $WORK_DIR/nbhmm2nz/ -name "${M}_${PEAKS}_chr${i}_*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-NBHMM2NZ$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN NBHMM3NZ"
      for F in $(find $WORK_DIR/nbhmm3nz/ -name "${M}_${PEAKS}_chr${i}_*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-NBHMM3NZ$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN NBHMM4NZ"
      for F in $(find $WORK_DIR/nbhmm4nz/ -name "${M}_${PEAKS}_chr${i}_*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-NBHMM4NZ$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN NBHMM4"
      for F in $(find $WORK_DIR/nbhmm4/ -name "${M}_${PEAKS}_chr${i}_*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-NBHMM4$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN Islands"
      for F in $(find $WORK_DIR/islands/ -name "${M}_${PEAKS}_chr${i}_*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}chr$i$T${TP}${T}SPAN-Islands$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

    done
  done 
done
