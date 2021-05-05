import os
import re
from glob import glob

localrules: all_reads_coverage_results

######## Step: Visualization: Reads coverage ##################
rule all_reads_coverage_results:
    input:
         bws=expand('bw/{name}.bw', name=[re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam')]),

rule bam2bw:
    input:
         bam='mix/{filename}.bam',
         bai='mix/{filename}.bam.bai'
    output: 'bw/{filename, [^/]*}.bw'
    log: 'logs/bw/{filename}.log'

    conda: '../envs/deeptools.env.yaml'
    threads: 4
    resources:
        threads = 4,
        mem = 16, mem_ram = 12,
        time = 60 * 120
    shell: 'bamCoverage -b {input.bam} -p {threads} -o {output}'
