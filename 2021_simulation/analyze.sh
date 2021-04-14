T=$'\t'
echo "Modification${T}TruePeaks${T}Tool${T}Fdr${T}Peaks${T}Sensitivity${T}Specificity${T}" > report.tsv

echo "MACS2 narrow"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/macs2/ -name "$M*.narrowPeak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*_q|_peaks.*//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}Macs2$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done

echo "MACS2 --broad"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/macs2/ -name "$M*.broadPeak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*_broad|_peaks.*//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}Macs2Broad$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done

echo "SPAN"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/span/ -name "$M*.peak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*200_|_5.peak//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}SPAN$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done

echo "SPAN NBHMM2NZ"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/nbhmm2nz/ -name "${M}_*.peak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*_|\.peak//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}Nbhmm2nz$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done

echo "SPAN NBHMM3NZ"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/nbhmm3nz/ -name "${M}_*.peak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*_|\.peak//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}Nbhmm3nz$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done


echo "SPAN NBHMM4NZ"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/nbhmm4nz/ -name "${M}_*.peak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*_|\.peak//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}Nbhmm4nz$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done

echo "SPAN NBHMM4"
for M in H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me3; do echo $M; TRUTH=$(find /mnt/stripe/shpynov/2021_noise2/chr15/ -name "$M*"); echo $TRUTH; TP=$(cat $TRUTH | wc -l); echo "True peaks $TP"; for F in $(find /mnt/stripe/shpynov/2021_noise2/nbhmm4/ -name "${M}_*.peak"); do echo $F; PEAKS=$(cat $F | wc -l); FDR=$(echo $F | sed -E 's/.*_|\.peak//g'); SENSITIVITY=$(bedtools intersect -a $TRUTH -b $F -wa -u | wc -l); SPECIFICITY=$(bedtools intersect -a $F -b $TRUTH -wa -u | wc -l); echo "Fdr $FDR Peaks $PEAKS Sensitivity $SENSITIVITY Specificity $SPECIFICITY"; echo "$M$T$TP${T}Nbhmm4$T$FDR$T$PEAKS$T$SENSITIVITY$T$SPECIFICITY" >> report.tsv;  done;  echo ""; done




