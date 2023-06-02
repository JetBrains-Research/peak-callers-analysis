# 2023_Immune PEAK CALLING and NOISE

echo "Generating mixture"

# Noise
cd ~/data/2023_Immune

for READSMLN in 15 10 5 2 1; do
    CONTROLMLN=$(bc -l <<< "20-${READSMLN}")
    for M in H3K36me3 H3K4me3 H3K4me1 H3K27ac H3K27me3; do
        for CT in TCell BCell Monocyte; do
            C=$(ls bam/${CT}*_Control_*.bam | grep -v mln);
            T=$(ls bam/${CT}*_${M}_*.bam | grep -v mln);
            echo "$READSMLN $CONTROLMLN $M $CT $R $T $C";
            samtools index $T;
            samtools index $C;
            T_FRACTION=$(samtools idxstats $T | cut -f3 | awk -v ct=${READSMLN}000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
            T_FRACTION=$(echo $T_FRACTION | sed 's/,/./g');
            C_FRACTION=$(samtools idxstats $C | cut -f3 | awk -v ct=${CONTROLMLN}000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
            C_FRACTION=$(echo $C_FRACTION | sed 's/,/./g');
            echo "$T_FRACTION $C_FRACTION";
            OUT=${T/.bam/_${READSMLN}mln_noise.bam}
            echo $OUT
            samtools view -@ 4 -H $T > ${OUT}.sam
            samtools view -@ 4 -s ${T_FRACTION} $T >> ${OUT}.sam
            samtools view -@ 4 -s ${C_FRACTION} $C >> ${OUT}.sam
            samtools view -S -b ${OUT}.sam > ${OUT}.unsorted.bam
            samtools sort ${OUT}.unsorted.bam -o ${OUT}
            rm ${OUT}.sam ${OUT}.unsorted.bam
        done
    done
done

mv bam/*mln* ~/2023_Immune_noise/bam/
cd ~/data/2023_Immune_noise

echo "Peak calling"
# MACS2
# Narrow and broad
cd ~/data/2023_Immune
mkdir macs2
conda activate macs2
for M in H3K27ac H3K4me3 H3K4me1 H3K36me3 H3K27me3; do
    for CT in TCell BCell Monocyte; do
        C=$(ls bam/${CT}*_Control_*.bam);
        for T in $(ls bam/${CT}*_${M}_*.bam); do
            echo "$M $CT $T $C";
            TN=$(basename $T);
            FN=${TN/.bam/};
            if [[ -z $(find macs2/ -name "${FN}_*") ]]; then
                macs2 callpeak -t $T -c $C -n ${FN}_q0.05 -g hs -q 0.05;
                macs2 callpeak -t $T -c $C -n ${FN}_broad0.1 -g hs --broad --broad-cutoff 0.1;
                mv ${FN}_q0.05* macs2/;
                mv ${FN}_broad0.1* macs2/;
            fi;
        done;
    done;
done;

# SICER
conda activate sicer
D=$(pwd)
mkdir sicer
for M in H3K27ac H3K4me3 H3K4me1 H3K36me3 H3K27me3; do
    for CT in TCell BCell Monocyte; do
        C=$(ls bam/${CT}*_Control_*.bam);
        for T in $(ls bam/${CT}*_${M}_*.bam); do
            echo "$M $CT $R $T $C";
            TN=$(basename $T);
            TU=${TN/.bam/};
            CN=$(basename $C);
            CU=${CN/.bam/};
            PEAKS=sicer/${TU}-W200-G600-islands-summary-FDR0.01;
            echo $PEAKS;
            if [[ ! -f $PEAKS ]]; then
                TMP=$(mktemp -d);
                mkdir -p $TMP;
                echo "$TMP $TU $CU";
                bedtools bamtobed -i $T > $TMP/$TU.bed
                bedtools bamtobed -i $C > $TMP/$CU.bed
                cd $TMP;
                mkdir out;
                echo "SICER.sh $TMP $TU.bed $CU.bed $TMP/out hg38 1 200 150 0.74 600 0.01;"
                SICER.sh $TMP $TU.bed $CU.bed $TMP/out hg38 1 200 150 0.74 600 0.01;
                mv out/* $D/sicer;
                rm -rf $TMP;
                cd $D;
            fi;
        done;
    done;
done;

# SPAN
conda activate base
mkdir span
for M in H3K27ac H3K4me3 H3K4me1 H3K36me3 H3K27me3; do
    for CT in TCell BCell Monocyte; do
        C=$(ls bam/${CT}*_Control_*.bam);
        for T in $(ls bam/${CT}*_${M}_*.bam); do
            echo "$M $CT $R $T $C";
            TN=$(basename $T);
            PEAKS=span/${TN/.bam/}_100_q0.05.peak;
            echo $PEAKS
            if [[ ! -f $PEAKS ]]; then
                java -jar ~/span-1.1.5628.jar analyze -w span -t $T -c $C -b 100 -clip -p $PEAKS -cs hg38.chrom.sizes -fdr 0.05 -threads 4;
            fi;
        done;
    done;
done;
