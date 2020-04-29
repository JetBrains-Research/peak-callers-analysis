import re
import os

localrules: all_span_results

######## Step: Peak Calling: SPAN ##################
rule all_span_results:
    input:
        span_peaks=expand(f'span/{{name}}_{config["span_bin"]}_{config["span_fdr"]}_{config["span_gap"]}.peak',
            name=[re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam')]
        )

rule all_span_tuned_results:
    input:
        span_peaks_tuned=expand(f'span/{{name}}_{config["span_bin"]}_tuned.peak',
            name=[re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam') if config['span_markup'] != '']
        )

rule all_span_replicated_results:
    input:
        span_rep_peaks=expand(f'span_rep/{{sample}}_{{noise}}_{config["span_bin"]}_{config["span_fdr"]}_{config["span_gap"]}.peak',
                          sample=[config['sample']], noise=[0, 1, 3, 5, 7, 9]
                          )


rule call_peaks_span:
    input:
        signal='mix/{name}.bam',
        chrom_sizes=rules.download_chrom_sizes.output
    output:
        peaks=f'span/{{name}}_{{bin}}_{{fdr}}_{{gap}}.peak'
    log: f'logs/span/{{name}}_{{bin}}_{{fdr}}_{{gap}}.log'

    conda: '../envs/java8.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        'java -Xmx8G -jar /mnt/stripe/shpynov/span-0.12.0.build.jar analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
        f'-c {config["input2"]} --peaks {{output.peaks}} --model span/fit/{{wildcards.name}}_{{wildcards.bin}}.span '
        '--workdir span --threads {threads} '
        '--bin {wildcards.bin} --fdr {wildcards.fdr} --gap {wildcards.gap} &> {log}'


rule call_peaks_span_tuned:
    input:
        signal='mix/{name}.bam',
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
        'java -Xmx8G -jar /mnt/stripe/shpynov/span-0.12.0.build.jar analyze -t {input.signal} --chrom.sizes {input.chrom_sizes} '
        f'-c {config["input2"]} --peaks {{output.peaks}} --model span/fit/{{wildcards.name}}_{{wildcards.bin}}.span '
        '--workdir span --threads {threads} '
        '--bin {wildcards.bin} --labels {params.span_markup} &> {log}'

def span_replicated_input_fun(wildcards):
    name, noise  = wildcards.name, wildcards.noise
    files = sorted([p for p in glob('mix/*.bam') if f'{name}_{noise}_' in p])
    print(name, noise, files)
    return dict(
        chrom_sizes=rules.download_chrom_sizes.output,
        **{f'signal{i}': f for i, f in enumerate(files)}  # Helpers to build DAG
    )


rule call_peaks_span_replicated:
    input: unpack(span_replicated_input_fun)
    output:
        peaks=f'span_rep/{{name}}_{{noise}}_{{bin}}_{{fdr}}_{{gap}}.peak'
    log: f'logs/span_rep/{{name}}_{{noise}}_{{bin}}_{{fdr}}_{{gap}}.log'
    conda: '../envs/java8.env.yaml'
    threads: 4
    params:
        signal=lambda wildcards, input: ','.join([v for k, v in input.items() if k.startswith('signal')])
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell:
        f'java -Xmx8G -jar /mnt/stripe/shpynov/span-0.12.0.build.jar analyze -t {{params.signal}} --chrom.sizes {{input.chrom_sizes}} '
        f'-c {config["input2"]} --peaks {{output.peaks}} '
        '--workdir span --threads {threads} '
        '--bin {wildcards.bin} --fdr {wildcards.fdr} --gap {wildcards.gap} &> {log}'
