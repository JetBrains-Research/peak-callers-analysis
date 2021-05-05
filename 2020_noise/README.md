This pipeline is based on https://github.com/JetBrains-Research/chipseq-smk-pipeline

1. BAM samples are in <work_dir>

2. Estimate signal to noise ratio for BAM files
```
snakemake -s <path>/Snakefile all --use-conda --directory <work_dir> --configfile <path>/noise/config.yaml --config bam_dir=<work_dir>
```

3. Find max signal to noise ratio among samples
```
find <work_dir> -name "max.txt" | xargs cat
```

4. Launch 10 noise rates x 3 replicates x 2 configs
```
cd <path_to_bams>
# CONFIG1
snakemake -s <path>/Snakefile all --use-conda --directory <work_dir> --configfile <path>/config.yaml --config bam_dir=<work_dir> sample=<sample_max_snr> input1=<input1.bam> n=3 input2=<input2.bam> span_markup=<labels.bed>

# CONFIG2
snakemake -s <path>/Snakefile all --use-conda --directory <work_dir> --configfile <path>/config2.yaml --config bam_dir=<work_dir> sample=<sample_max_snr> input1=<input1.bam> n=3 input2=<input2.bam> span_markup=<labels.bed>
``` 

