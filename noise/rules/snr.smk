####### Step: SNR estimation ##################
localrules: all_snr_results

rule all_snr_results:
    input:
         snr_max='snr/max.txt'

rule index_bams:
    input: '{anywhere}/{sample}.bam'
    output: '{anywhere}/{sample}.bam.bai'
    wrapper: '0.36.0/bio/samtools/index'

rule snr:
    input:
        bam=f"{config['bam_dir']}/{{sample}}.bam",
        bai=f"{config['bam_dir']}/{{sample}}.bam.bai"
    output: 'snr/{sample}.txt'
    resources:
        threads = 1,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    conda: '../envs/bio.env.yaml'
    shell:
        'python scripts/signal_to_noise_estimation.py {input.bam} -d 140 > {output}'

rule snr_max:
    input: expand('snr/{sample}.txt', sample=names_wo_ext(BAM_PATHS))
    output: 'snr/max.txt'
    shell:
        'SNR_MAX=0; F_MAX=""; '
         'for F in snr/*.txt; do SNR=$(cat $F); if (( $(echo "$SNR > $SNR_MAX" | bc -l) )); then SNR_MAX=$SNR; F_MAX=$F; fi; done; '
         'echo "$F_MAX\t$SNR_MAX" > {output}'


