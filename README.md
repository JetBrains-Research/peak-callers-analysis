Peak callers analysis
=====================

Analysis and comparison scripts for various peak callers, supported by the ChIP-seq analysis [pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline).

Peak calling
============

1. Fetch [chipseq-smk-pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) GitHub repository.
2. Prepare datasets by downloading files mentioned in [Datasets.xlsx](./Datasets.xlsx).
3. Alignment of ChIP-seq datasets to the reference genome.
```bash
echo "Alignment"
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --cores all --use-conda --directory $(pwd) --config genome=<genome> \
  fastq_dir=$(pwd)/fastq fastq_ext=fastq \
  --rerun-incomplete --rerun-trigger mtime;
```
4. Peak calling of ChIP-seq datasets.
```bash
echo "Peak calling with default settings (MACS2 narrow, HOMER factor)"
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --cores all --use-conda --directory $(pwd) --config genome=<genome> \
  start_with_bams=true \
  macs2=True sicer=True homer=True hotspot=True peakseq=True lanceotron=True omnipeak=True \
  --rerun-incomplete --rerun-trigger mtime;
  
echo "Peak calling other settings (MACS2 broad, HOMER histone)"
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --cores all --use-conda --directory $(pwd) --config genome=<genome> \
  start_with_bams=true \
  macs2=True macs2_mode=broad macs2_params="--broad --broad-cutoff 0.1" macs2_suffix=broad0.1 \
  homer=True homer_style=histone homer_suffix=regions.bed \
  --rerun-incomplete --rerun-trigger mtime;
```
5. Alignment and peak calling of Immgen ATAC-seq dataset.
```bash
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --use-conda --cores all --directory $(pwd) --config genome=mm10 \
  macs2=True bowtie2_params="-X 2000 --dovetail" span=True span_fragment=0 \
  --rerun-incomplete --rerun-trigger mtime;
```

Simulations
===========

See [Simulation instructions](./chips/Simulation.md) for details.

Requirements
============

Please ensure that you have the following Python packages installed:
* Jupyter
* Pandas
* PyRanges
* PyBigwig
* Seaborn
* Statannotations 
* Scipy

Please ensure that the following tools are available:
* bedtools
* samtools

