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

mkdir -p ~/22023_GSE26320_noise/bams
mv bams/*mln* ~/22023_GSE26320_noise/bams/
cd ~/2023_GSE26320_noise

# Peak calling with Snakemake
conda activate snakemake

snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --cores all  --use-conda  --directory $(pwd) \
    --config fastq_dir=$(pwd)/bams  start_with_bams=True \
    macs2=True macs2_mode=narrow macs2_params="-q 0.05" macs2_suffix="q0.05"

snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --cores all  --use-conda  --directory $(pwd) \
    --config fastq_dir=$(pwd)/bams  start_with_bams=True \
    macs2=True macs2_mode=broad macs2_params="--broad --broad-cutoff=0.1" macs2_suffix="broad0.1"

snakemake -p -s ~/work/chipseq-smk-pipeline/Snakefile all --cores all  --use-conda  --directory $(pwd) \
    --config fastq_dir=$(pwd)/bams  start_with_bams=True \
    sicer=True span=True
