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
            OUT=${T/.bam/_${READSMLN}mln_noise.bam}
            echo $OUT
            if [[ ! -f $OUT ]]; then
                samtools index $T;
                samtools index $C;
                T_FRACTION=$(samtools idxstats $T | cut -f3 | awk -v ct=${READSMLN}000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
                C_FRACTION=$(samtools idxstats $C | cut -f3 | awk -v ct=${CONTROLMLN}000000 'BEGIN {total=0} {total += $1} END {print ct/total}')
                T_FRACTION=$(echo $T_FRACTION | sed 's/,/./g');
                C_FRACTION=$(echo $C_FRACTION | sed 's/,/./g');
                if [[ $(echo "$T_FRACTION>1" | bc -l) -eq 1 ]]; then echo $T; T_FRACTION=1.0; fi
                if [[ $(echo "$C_FRACTION>1" | bc -l) -eq 1 ]]; then echo $C; C_FRACTION=1.0; fi
                echo "$T_FRACTION $C_FRACTION"
                samtools view -@ 4 -H $T > ${OUT}.sam
                samtools view -@ 4 -s ${T_FRACTION} $T >> ${OUT}.sam
                samtools view -@ 4 -s ${C_FRACTION} $C >> ${OUT}.sam
                samtools view -S -b ${OUT}.sam > ${OUT}.unsorted.bam
                samtools sort ${OUT}.unsorted.bam -o ${OUT}
                rm ${OUT}.sam ${OUT}.unsorted.bam
            fi
        done
    done
done

mkdir -p ~/2023_Immune_noise/bam/
mv bam/*mln* ~/2023_Immune_noise/bam/
cd ~/data/2023_Immune_noise

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
