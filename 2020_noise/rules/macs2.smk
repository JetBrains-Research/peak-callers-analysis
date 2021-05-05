from glob import glob
import os
import re


######## Step: Call peaks MACS2 ##################

def macs_species(genome):
    """Convert genome to macs2 species encoding"""
    if re.match('^hg[0-9]+$', genome):
        return 'hs'
    elif re.match('^mm[0-9]+$', genome):
        return 'mm'
    raise Exception('Unknown species {}'.format(genome))

localrules: all_macs2_results

rule all_macs2_results:
    input:
        macs2_peaks=expand(
            f'macs2/{{name}}_{config["macs2_suffix"]}_peaks.{config["macs2_mode"]}Peak',
            name=[re.sub('\\.bam$', '', os.path.basename(p)) for p in glob('mix/*.bam')]
        )

rule call_peaks_macs2:
    input: signal='mix/{name}.bam'
    output: f'macs2/{{name}}_{config["macs2_suffix"]}_peaks.{config["macs2_mode"]}Peak'
    log: f'logs/macs2_{config["macs2_suffix"]}/{{name}}_{config["macs2_suffix"]}_{config["macs2_mode"]}.log'
    conda: '../envs/py27.env.yaml'
    params:
        macs2_params=config['macs2_params'],
        macs2_suffix=config['macs2_suffix'],
        species=macs_species(config['genome']),
        outdir=lambda wildcards, output: os.path.dirname(str(output[0])),
    shell:
        f'macs2 callpeak -t {{input.signal}} -c {config["input2"]} --outdir {{params.outdir}} ' +\
        '-n {wildcards.name}_{params.macs2_suffix} -g {params.species} ' 
        '{params.macs2_params} &> {log}'


