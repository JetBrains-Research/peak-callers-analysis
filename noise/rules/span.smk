import re
import os

localrules: all_span_results, download_span

######## Step: Peak Calling: SPAN ##################
rule all_span_results:
    input:
        span_peaks=expand(f'span/{{name}}_{config["span_bin"]}_{config["span_fdr"]}_{config["span_gap"]}.peak',
            name=[re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam')]
        )

rule all_span_tuned_results:
    input:
        span_peaks_tuned=expand(f'span/{{name}}_{config["span_bin"]}_tuned.peak',
            name=[re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam')]
        )

rule download_span:
    output: 'bin/span-0.11.0.jar'
    shell: 'wget -O {output} https://download.jetbrains.com/biolabs/span/span-0.11.0.4882.jar'


rule call_peaks_span:
    input:
        signal='mix/{name}.bam',
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        peaks=f'span/{{name}}_{{bin}}_{config["span_fdr"]}_{config["span_gap"]}.peak'
    log: f'logs/span/{{name}}_{{bin}}_{config["span_fdr"]}_{config["span_gap"]}.log'

    conda: '../envs/java8.env.yaml'
    threads: 4
    params:
        fdr=config["span_fdr"],
        gap=config["span_gap"],
        span_params=config['span_params']
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        'java -Xmx8G -jar {input.span} analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
        f'-c {config["input2"]} --peaks {{output.peaks}} --model span/fit/{{wildcards.name}}_{{wildcards.bin}}.span '
        '--workdir span --threads {threads} '
        '--bin {wildcards.bin} --fdr {params.fdr} --gap {params.gap} {params.span_params} &> {log}'


rule call_peaks_span_tuned:
    input:
        signal='mix/{name}.bam',
        span=rules.download_span.output,
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        peaks='span/{name}_{bin}_tuned.peak'
    log: 'logs/span/{name}_{bin}_tuned.log'

    conda: '../envs/java8.env.yaml'
    threads: 4
    params:
        span_markup=config['span_markup']
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        'java -Xmx8G -jar {input.span} analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
        f'-c {config["input2"]} --peaks {{output.peaks}} --model span/fit/{{wildcards.name}}_{{wildcards.bin}}.span '
        '--workdir span --threads {threads} '
        '--bin {wildcards.bin} --labels {params.span_markup} &> {log}'

