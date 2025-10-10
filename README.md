Peak callers analysis
=====================

Analysis and comparison scripts for various peak callers, supported by the ChIP-seq
analysis [pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline).

# Notebooks

* [1_benchmark.ipynb](./1_benchmark.ipynb) - Jupyter notebook with initial running time benchmarking of peak callers.
* [2_peaks.ipynb](2_peaks.ipynb) - Notebook with the number of peaks / length / jaccard between replicates, etc.
* [3_rnaseq.ipynb](./3_rnaseq.ipynb) - Notebook of peak calling comparison with RNA-seq data.
* [4_chips.ipynb](./4_chips.ipynb) - Benchmark of peak callers versus artificially simulated ChIP-seq data.
* [5_control.ipynb](./5_control.ipynb) - Assessment of the effect of control on peak calling.
* [6_immgen.ipynb](./6_immgen.ipynb) - Benchmark of peak callers on ATAC-seq data.
* [7_omnipeak.ipynb](7_omnipeak.ipynb) - Deep model analysis of Omnipeak models.
* [8_hyperparameters.ipynb](./8_hyperparameters.ipynb) - Analysis of hyperparameters of Omnipeak.

# Datasets

Prepare datasets by downloading files mentioned in [Datasets.xlsx](./Datasets.xlsx).

1. Download `fastq` files from the tab `GSE26320` into `~/data/2023_GSE26320` folder.
2. Download `bam` files from the tab `RoadmapEpigenomics` into `~/data/2023_Immune` folder.
3. Download `bed.gz` files from the tab `ABF` into `~/data/2018_chipseq_y20o20` folder.
   Convert them to `bam` format using `samtools`.
4. Download `bam` files from the tab `CTCF` into `~/data/2025_TFs` folder.
5. Download `fastq` files from the tab `Immgen` into `~/data/2025_Immgen` folder.
6. Download `tsv` files with transcription counts from the tab `RNAseq` into `~/data/2025_transcription` folder.
7. Download `bam` files from the tab `Chips` into `~/data/2025_chips` folder.

**Files layout** - please place `fastq` datasets into `fastq` subfolder, and `bam` datasets into `bam` subfolder.<br>
Datasets **without control** should be prepared by copying all the raw data without control files into the corresponding folders with `_no_control` suffix.<br>
Please ensure to use **a correct genome** version for the datasets -  `mm10` for `Immgen`, `hg19` for `ABF` and `hg38` for
the rest.

# Peak calling

1. Fetch [chipseq-smk-pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) GitHub repository into `~/work/chipseq-smk-pipeline`.
2. Navigate to the dataset folder.
3. Launch alignment of datasets to the reference genome (optional).

```bash
echo "Alignment"
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --cores all --use-conda --directory $(pwd) --config genome=<genome> \
  fastq_dir=$(pwd)/fastq fastq_ext=fastq \
  --rerun-incomplete --rerun-trigger mtime;
```
Use additional `bowtie2_params="-X 2000 --dovetail"` parameter for ATAC-seq alignment. 

3. Peak calling of ChIP-seq / ATAC-seq datasets.

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

# Simulations

See [Simulation instructions](./chips/Simulation.md) for details.

# Requirements

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

