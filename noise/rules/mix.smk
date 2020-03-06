####### Step: Mixing Signal and Noise ##################
localrules: all_mix_results

rule all_mix_results:
    input: expand('mix/{sample}_{noise}_{n}.bam', sample=[config['sample']], noise=[0, 1, 3, 5, 7, 9], n=range(int(config['n'])))

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

