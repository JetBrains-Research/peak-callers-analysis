!%% bash
# Bash commands to create markup by Immgen
WORKDIR=~/data/2023_Immgen
DHS_FILE=$WORKDIR/ENCFF754WCT_mm10_dhs_representative_sites.bed;

T=$'\t';
printf %s "chr${T}start${T}end" > $WORKDIR/intersect.tsv;
FILES=();
for F in $(find $WORKDIR/macs2/ -name "*.narrowPeak"); do
	FILES+=("$F");
	printf %s "${T}${F}" >> $WORKDIR/intersect.tsv;
done;
echo >> $WORKDIR/intersect.tsv;
bedtools multiinter -i "${FILES[@]}" |\
	bedtools merge -c $(seq -s, 6 $((${#FILES[@]} + 5))) -o max |\
	awk '{if (NR > 1) printf("\n"); printf("%s\t%s\t%s", $1, $2, $3); for (i=4; i<=NF; i++) printf("\t%d", int($i)); }'\
    >> $WORKDIR/intersect.tsv;

# Find out regions where all the peaks present
ALL="";
for _ in $(seq 1 ${#FILES[@]}); do
    ALL="${ALL}${T}1";
done;
cat $WORKDIR/intersect.tsv | grep "${ALL}" | awk -v OFS='\t' '{print $1,$2,$3}' > $WORKDIR/intersect_all.bed
wc -l $WORKDIR/intersect_all.bed

# Find regions interesting with at least 20% by DHS and report original DHS.
bedtools intersect -b $WORKDIR/intersect_all.bed -a $DHS_FILE -wa -u -F 0.5 > $WORKDIR/mm10_dhs_intersect_all.bed
wc -l $WORKDIR/mm10_dhs_intersect_all.bed


# Total 200 peaks
# 100 peaks
head -n 50 $WORKDIR/mm10_dhs_intersect_all.bed |\
    while read -r LINE; do echo "$LINE" | awk -v OFS='\t' '{print $1,$2,$3,"peaks"}'; done > $WORKDIR/markup.bed

#50 peakStart
head -n 100 $WORKDIR/mm10_dhs_intersect_all.bed | tail -n 50 |\
    while read -r LINE; do echo "$LINE" | awk '{ printf("%s\t%d\t%d\t%s\n", $1,$2,($2+$3)/2 - 1,"peakStart")}'; done\
    >> $WORKDIR/markup.bed

#50 peakEnd
head -n 150 $WORKDIR/mm10_dhs_intersect_all.bed | tail -n 50 |\
    while read -r LINE; do echo "$LINE" | awk '{printf("%s\t%d\t%d\t%s\n", $1,($2+$3)/2 + 1,$3,"peakEnd")}'; done\
    >> $WORKDIR/markup.bed
wc -l $WORKDIR/markup.bed

# extended markup
cat $WORKDIR/markup.bed | while read -r LINE; do echo "$LINE" | awk '{print($1,$2-2000,$3+2000)}'; done\
    > $WORKDIR/markup_ext.bed