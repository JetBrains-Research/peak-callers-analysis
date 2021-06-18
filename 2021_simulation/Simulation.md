Simulation
==========

**Chips** is available from bioconda and on [GitHub](https://github.com/gymreklab/chips).
Paper in [Bioinformatics](https://link.springer.com/article/10.1186/s12859-021-04097-5).

# Data
Download data for human CD14-positive monocyte cells.

ChIP-seq
* [H3K4me1](https://www.encodeproject.org/files/ENCFF076WOE/)
* [H3K4me3](https://www.encodeproject.org/files/ENCFF001FYS/)
* [H3K27ac](https://www.encodeproject.org/files/ENCFF000CEN/)
* [H3K27me3](https://www.encodeproject.org/files/ENCFF001FYR/)
* [H3K36me3](https://www.encodeproject.org/files/ENCFF000CFB/)

Control
* [ENCFF825XKT](https://www.encodeproject.org/files/ENCFF825XKT/) (for H3K4me1)
* [ENCFF001HUV](https://www.encodeproject.org/files/ENCFF001HUV/) (for H3K4me3 and H3K27me3)
* [ENCFF692GVG](https://www.encodeproject.org/files/ENCFF692GVG/) (for H3K27ac and H3K36me3)

Reference
* [hg38 reference fa](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/)

# Download data
```
WORK_DIR=/mnt/stripe/shpynov/2021_chips
cd $WORK_DIR

# Download precomputed encode peaks
wget https://www.encodeproject.org/files/ENCFF366GZW/@@download/ENCFF366GZW.bed.gz -O H3K4me1.bed.gz
wget https://www.encodeproject.org/files/ENCFF651GXK/@@download/ENCFF651GXK.bed.gz -O H3K4me3.bed.gz  
wget https://www.encodeproject.org/files/ENCFF039XWV/@@download/ENCFF039XWV.bed.gz -O H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF666NYB/@@download/ENCFF666NYB.bed.gz -O H3K27me3.bed.gz      
wget https://www.encodeproject.org/files/ENCFF213IBM/@@download/ENCFF213IBM.bed.gz -O H3K36me3.bed.gz

# Download fastq reads data 
wget https://www.encodeproject.org/files/ENCFF076WOE/@@download/ENCFF076WOE.fastq.gz - O H3K4me1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001FYS/@@download/ENCFF001FYS.fastq.gz -O H3K4me3.fastq.gz
wget https://www.encodeproject.org/files/ENCFF000CEN/@@download/ENCFF000CEN.fastq.gz -O H3K27ac.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001FYR/@@download/ENCFF001FYR.fastq.gz -O H3K27me3.fastq.gz      
wget https://www.encodeproject.org/files/ENCFF000CFB/@@download/ENCFF000CFB.fastq.gz -O H3K36me3.fastq.gz
 
wget https://www.encodeproject.org/files/ENCFF825XKT/@@download/ENCFF825XKT.fastq.gz -O input_H3K4me1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001HUV/@@download/ENCFF001HUV.fastq.gz -O input_H3K4me3_H3K27me3.fastq.gz
wget https://www.encodeproject.org/files/ENCFF692GVG/@@download/ENCFF692GVG.fastq.gz -O input_H3K27ac_H3K36me3.fastq.gz

# Dowload fasta reference
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz 

gunzip *.gz
```

For MACS2, SICER, SPAN peaks launch [ChIP-seq snakemake pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) from fastq files.
Assuming that it is installed to `/mnt/stripe/shpynov/chipseq-smk-pipeline`.


# Learn models and create modified models with tweaked FRIP

```
# Ensure that chips is available!
bash learn.sh
bash frip.sh
```

# Simulate reads

```
# Ensure that chips is available!
bash simulate.sh
mkdir fastq
mv *.fastq fastq/
```

# Prepare control for peak calling 

1. Launch chipseq pipeline on input files only to obtain bam files. 
2. Filter to chromosome 15
```
for F in $(ls input*.bam | grep -v chr15); do 
    echo $F; 
    samtools view $F chr15 -b > ${F/.bam/_chr15.bam}; 
done
 
for F in input*chr15.bam; do 
    echo $F; 
    bedtools bamtofastq -i $F -fq ${F/.bam/.fastq}; 
done
```
3. Copy resulting input fastq files into the `/fastq` folder.

# Launch peak callers
```
cd /mnt/stripe/shpynov/chipseq-smk-pipeline
conda activate snakemake

WORK_DIR=/mnt/stripe/shpynov/2021_chips

# Perform peak calling using chipseq snakemake pipeline
for FDR in 0.05 0.01 1e-3 1e-4 1e-5 1e-6; do
  echo "FDR $FDR"
  
  echo "MACS2 narrow"
  snakemake all --cores 24 --use-conda --directory $WORK_DIR \--config genome=hg38 \
    fastq_dir=$WORK_DIR/fastq fastq_ext=fastq macs2_params="-q $FDR" macs2_mode=narrow macs2_suffix=q$FDR \
    span_fdr=$FDR sicer_fdr=$FDR;
  
  echo "MACS2 broad"
  snakemake all --cores 24 --use-conda --directory $WORK_DIR --config genome=hg38 \
    fastq_dir=$WORK_DIR/fastq fastq_ext=fastq macs2_params="--broad --broad-cutoff $FDR" macs2_suffix=broad$FDR \
    span_fdr=$FDR sicer_fdr=$FDR;
  
done
```


# Launch SPAN modifications

```
cd $WORK_DIR
for FDR in 0.05 0.01 1e-3 1e-4 1e-5 1e-6; do 
    snakemake -f span-modifications-smk/Snakefile all --cores 24 --config fdr=$FDR; 
done
```


#Report

Prepare data report by collecting all the peak calling files and overlap with ground truth. 
```
bash analyze.sh
```

#Visualize results

Launch `analysis.ipynb` jupyter notebook for analysis and visualization of results.
