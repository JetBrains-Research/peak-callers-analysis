# 2023_GSE26320 PEAK CALLING and NOISE

echo "Generating mixture"
cd ~/22023_GSE26320

for READSMLN in 7 5 2 1; do
    for M in H3K36me3 H3K4me3 H3K4me1 H3K27ac H3K27me3; do
        for CT in GM12878 K562 Huvec; do
            for R in rep1 rep2; do
                C=$(ls bams/*${CT}*_Input_${R}.bam);
                for T in $(ls bams/*${CT}*_${M}_${R}.bam | grep -v Input); do
                    echo "$M $CT $R $T $C";
                    OUT=${T/.bam/_${READSMLN}mln_noise.bam};
                    echo $OUT
                    if [[ ! -f "$OUT" ]]; then
                        CONTROLMLN=$(bc -l <<< "10-${READSMLN}")
                        echo "$READSMLN $CONTROLMLN";
                        samtools index $T;
                        samtools index $C;
                        T_FRACTION=$(samtools idxstats $T | cut -f3 | awk -v ct=${READSMLN}000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
                        C_FRACTION=$(samtools idxstats $C | cut -f3 | awk -v ct=${CONTROLMLN}000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
                        T_FRACTION=$(echo $T_FRACTION | sed 's/,/./g');
                        C_FRACTION=$(echo $C_FRACTION | sed 's/,/./g');
                        if [[ $(echo "$T_FRACTION>1" | bc -l) -eq 1 ]]; then T_FRACTION=1.0; fi
                        if [[ $(echo "$C_FRACTION>1" | bc -l) -eq 1 ]]; then C_FRACTION=1.0; fi
                        echo "$T_FRACTION $C_FRACTION";
                        samtools view -@ 4 -H $T > ${OUT}.sam
                        samtools view -@ 4 -s ${T_FRACTION} $T >> ${OUT}.sam
                        samtools view -@ 4 -s ${C_FRACTION} $C >> ${OUT}.sam
                        samtools view -S -b ${OUT}.sam > ${OUT}.unsorted.bam
                        samtools sort ${OUT}.unsorted.bam -o ${OUT}
                        rm ${OUT}.sam ${OUT}.unsorted.bam
                       fi;
                done;
            done;
        done;
    done;
done;

mv bams/*mln* ~/22023_GSE26320_noise/bams/
cd ~/2023_GSE26320_noise

echo "Peak calling"
# MACS2
# Narrow and broad

conda activate macs2
mkdir macs2

for CT in GM12878 K562 Huvec; do
    for M in H3K36me3 H3K4me3 H3K4me1 H3K27ac H3K27me3; do
        for R in rep1 rep2; do
            C=$(ls bams/*${CT}*_Input_${R}*.bam);
            for T in $(ls bams/*${CT}*_${M}_${R}*.bam | grep -v WCE); do
                echo "$M $CT $R $T $C";
                TN=$(basename $T);
                TU=${TN/.bam/};
                if [[ -z $(find macs2/ -name "${TU}*") ]]; then
	                macs2 callpeak -t $T -c $C -n ${TU}_q0.05 -g hs -q 0.05;
	                macs2 callpeak -t $T -c $C -n ${TU}_broad0.1 -g hs --broad --broad-cutoff 0.1;
                    mv ${TU}*q0.05* macs2/;
                    mv ${TU}*broad0.1* macs2/;
                fi;
            done;
        done;
    done;
done;

# SICER
conda activate sicer
mkdir sicer
D=$(pwd)
for CT in GM12878 K562 Huvec; do
    for M in H3K36me3 H3K4me3 H3K4me1 H3K27ac H3K27me3; do
        for R in rep1 rep2; do
            C=$(ls bams/*${CT}*_Input_${R}*.bam);
            for T in $(ls bams/*${CT}*_${M}_${R}*.bam | grep -v WCE); do
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
done;

# SPAN
mkdir span
conda activate base
for CT in GM12878 K562 Huvec; do
    for M in H3K36me3 H3K4me3 H3K4me1 H3K27ac H3K27me3; do
        for R in rep1 rep2; do
            C=$(ls bams/*${CT}*_Input_${R}*.bam);
            for T in $(ls bams/*${CT}*_${M}_${R}*.bam | grep -v WCE); do
                echo "$M $CT $R $T $C";
                TN=$(basename $T);
                PEAKS=span/${TN/.bam/}_100_q0.05.peak;
                echo $PEAKS
                if [[ ! -f $PEAKS ]]; then
                    java -jar  ~/span-1.1.5628.jar analyze -w span -t $T -c $C -b 100 -clip -p $PEAKS -cs hg38.chrom.sizes -fdr 0.05;
                fi;
            done;
        done;
    done;
done;
