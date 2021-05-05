PEAKS_DIR=/mnt/stripe/shpynov/2021_chips/peaks
WORK_DIR=/mnt/stripe/shpynov/2021_chips

T=$'\t'
echo "Modification${T}PeaksSource${T}Library${T}I${T}TruePeaks${T}Tool${T}Fdr${T}Peaks${T}Precision${T}Recall" > report.tsv

for PEAKS in encode macs2 sicer span; do
  echo "Peaks $PEAKS"
  for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do
    echo "Modification $M"
    for lib in 500k 1mln; do
	echo $lib;
    for i in {1..20}; do
      echo "Chromosome chr15"
      TRUTH=$WORK_DIR/fastq/${M}_${PEAKS}_chr15_$i.bed
      echo "True peaks file $TRUTH"; 
      TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; 

      echo "MACS2 narrow"
      for F in $(find $WORK_DIR/macs2/ -name "${M}_${PEAKS}_chr15_${i}_${lib}_*.narrowPeak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}Macs2$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;  
      done;
      echo ""
      
      echo "MACS2 --broad"
      for F in $(find $WORK_DIR/macs2/ -name "${M}_${PEAKS}_chr15_${i}_${lib}*.broadPeak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*_broad|_peaks.*//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}Macs2Broad$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""
      
      echo "SPAN gap 5"
      for F in $(find $WORK_DIR/span/ -name "${M}_${PEAKS}_chr15_${i}_${lib}_*_5.peak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*200_|_5.peak//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-GAP5$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""            
      
      echo "SPAN gap 0"
      for F in $(find $WORK_DIR/span/ -name "${M}_${PEAKS}_chr15_${i}_${lib}*_0.peak"); do
        echo $F;
        NP=$(cat $F | wc -l); 
        FDR=$(echo $F | sed -E 's/.*200_|_0.peak//g'); 
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); 
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); 
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE"; 
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-GAP0$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""            

      echo "SPAN NBHMM2NZ"
      for F in $(find $WORK_DIR/nbhmm2nz/ -name "${M}_${PEAKS}_chr15_${i}_${lib}*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-NBHMM2NZ$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN NBHMM3NZ"
      for F in $(find $WORK_DIR/nbhmm3nz/ -name "${M}_${PEAKS}_chr15_${i}_${lib}*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-NBHMM3NZ$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN NBHMM4NZ"
      for F in $(find $WORK_DIR/nbhmm4nz/ -name "${M}_${PEAKS}_chr15_${i}_${lib}*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-NBHMM4NZ$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN NBHMM4"
      for F in $(find $WORK_DIR/nbhmm4/ -name "${M}_${PEAKS}_chr15_${i}_${lib}*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-NBHMM4$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

      echo "SPAN Islands"
      for F in $(find $WORK_DIR/islands/ -name "${M}_${PEAKS}_chr15_${i}_${lib}_*.peak"); do
        echo $F;
        NP=$(cat $F | wc -l);
        FDR=$(echo $F | sed -E 's/.*_|\.peak//g');
        PR=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l);
        RE=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l);
        echo "Fdr $FDR Peaks $NP Precision $PR Recall $RE";
        echo "${M}$T${PEAKS}${T}$lib${T}$i$T${TP}${T}SPAN-Islands$T${FDR}$T${NP}$T${PR}$T${RE}" >> report.tsv;
      done;
      echo ""

    done
done
  done 
done

