Pipeline to create mixed ChIP and control reads in given proportions.
This pipeline is based on https://github.com/JetBrains-Research/chipseq-smk-pipeline

1. BAM samples are in <work_dir>

2. Estimate signal to noise ratio for BAM files
```
snakemake -s Snakefile all --use-conda --directory <work_dir> --config bam_dir=<work_dir>
```

3. Find max signal-to-noise ratio among samples
```
find <work_dir> -name "max.txt" | xargs cat
```

4. Launch 10 noise rates x 3 replicates with chipseq-smk-pipeline
