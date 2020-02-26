####### Step: Mixing Signal and Noise ##################
localrules: all_mix_results

rule all_mix_results:
    input: expand('mix/{sample}_{noise}_{n}.bam', sample=[config['sample']], noise=list(range(10))+[9.5, 9.9], n=range(int(config['n'])))

rule mix:
    input:
        track='{sample}.bam',
    output: 'mix/{sample}_{noise}_{n}.bam'
    resources:
        threads = 1,
        mem = 8, mem_ram = 4,
        time = 60 * 120
    conda: '../envs/bio.env.yaml'
    shell:
         f'bash scripts/join_files.sh {{input.track}} {config["input1"]} {{wildcards.noise}} {{output}};'

