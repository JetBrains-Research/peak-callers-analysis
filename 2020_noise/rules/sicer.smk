from glob import glob
import os
import re

######## Step: Peak Calling: SICER ##################
localrules: all_sicer_results

# !!!!!
# SICER doesn't support out of the box hg38 and mm10, needs to be tweaked a bit
# !!!!!

def sicer_all_peaks_input():
    files = []

    window_size=config['sicer_window']
    gap=config['sicer_gap']

    # XXX: change significance only via config, SICER rule takes the value from
    # config, not via wildcards

    for name in [re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam')]:
        # with control
        significance = config['sicer_fdr']
        files.append(f'sicer/{name}-W{window_size}-G{gap}-islands-summary-FDR{significance}')

    return files

def effective_genome_fraction(genome, chrom_sizes_path, pileup_bed):
    """From MACS2 documentation:
    The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
    Here are all precompiled parameters for effective genome size:
    hs: 2.7e9
    mm: 1.87e9
    ce: 9e7
    dm: 1.2e8"""

    # Get chr names covered with reads (e.g. if data is filtered by chromosome name
    # or some chrs excluded during alignment
    chromosomes = set()
    with open(str(pileup_bed)) as f:
        for line in f:
            chr = line.split()[0]
            chromosomes.add(chr)

    # Sized of chromosomes covered with reads
    chrom_sizes = {}
    with open(str(chrom_sizes_path)) as f:
        for line in f:
            chromosome, size = line.split()
            chrom_sizes[chromosome] = int(size)

    # Normalization if not all genome chromosomes are covered
    chromosomes_length = sum([chrom_sizes.get(c, 0) for c in chromosomes])
    genome_length = sum(chrom_sizes.values())

    if genome.startswith('mm'):
        size = 1.87e9
    elif genome.startswith('hg'):
        size = 2.7e9
    else:
        raise Exception('Unknown species {}'.format(genome))
    return (size / genome_length) * (1.0 * chromosomes_length / genome_length)

rule all_sicer_results:
    input: *sicer_all_peaks_input()


rule bam_to_pileup:
    input: 'mix/{name}.bam'
    output: 'pileup/{name}.bed'

    conda: '../envs/bio.env.yaml'
    shell: 'bedtools bamtobed -i {input} > {output}'


rule control_bam_to_pileup:
    output: 'pileup/control.bed'

    conda: '../envs/bio.env.yaml'
    shell: f'bedtools bamtobed -i {config["input2"]} > {{output}}'


rule download_chrom_sizes:
    output: f"{config['genome']}.chrom.sizes"
    log: f"logs/{config['genome']}.chrom.sizes.log"

    shell:
        'wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/{config[genome]}/bigZips/{config[genome]}.chrom.sizes &> {log}'

rule pileup_bed_effective_genome_fraction:
    input:
        pileup_bed=rules.bam_to_pileup.output,
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        temp(str(rules.bam_to_pileup.output) + ".egf")
    run:
        value = effective_genome_fraction(
            config['genome'], input.chrom_sizes, input.pileup_bed
        )
        shell("echo -n '{value}' > {output}")


rule call_peaks_sicer:
    input:
        signal_pileup='pileup/{name}.bed',
        control_pileup='pileup/control.bed',
        chrom_sizes=rules.download_chrom_sizes.output,
        effective_genome_fraction=rules.pileup_bed_effective_genome_fraction.output
    output: 'sicer/{name}-W{width}-G{gap, \d+}-{any_suffix}'
    log: 'logs/sicer/{name}-W{width}-G{gap}-{any_suffix}.log'

    conda: '../envs/py27.env.yaml'
    shadow: "shallow"
    params:
        signal_pileup_bed_fname=lambda wildcards, input: os.path.basename(input.signal_pileup),
        pileups_dir=lambda wildcards, input: os.path.split(str(input.signal_pileup))[0],
        peaks_file=lambda wildcards, output: os.path.basename(output[0]),
        fragment=config['sicer_fragment'],
        genome=config['genome'],
        significance=config['sicer_fdr']
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        # SICER.sh ["InputDir"] ["bed file"] ["control file"]
        #       ["OutputDir"] ["Species"] ["redundancy threshold"]
        #       ["window size (bp)"] ["fragment size"] ["effective genome fraction"]
        #       ["gap size (bp)"] [“FDR"]
        'echo "Significance threshold: {params.significance}" > {log}; '
        ' tmp_sicer=$(mktemp -d -p $(pwd) -t sicer-XXXXXXXXXX); mkdir -p $tmp_sicer; cd $tmp_sicer; '
        ' SICER.sh ../{params.pileups_dir} {params.signal_pileup_bed_fname} control.bed'
        '  $(pwd) {params.genome} 1 {wildcards.width}'
        '  {params.fragment} $(cat "../{input.effective_genome_fraction}")'
        '  {wildcards.gap} {params.significance} &>> ../{log}; '
        ' ls -lah  &>> ../{log}; mv {params.peaks_file} ../{output} &>> ../{log};'

