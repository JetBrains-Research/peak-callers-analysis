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

mkdir fastq
cd fastq

# Download fastq reads data 
wget https://www.encodeproject.org/files/ENCFF076WOE/@@download/ENCFF076WOE.fastq.gz - O H3K4me1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001FYS/@@download/ENCFF001FYS.fastq.gz -O H3K4me3.fastq.gz
wget https://www.encodeproject.org/files/ENCFF000CEN/@@download/ENCFF000CEN.fastq.gz -O H3K27ac.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001FYR/@@download/ENCFF001FYR.fastq.gz -O H3K27me3.fastq.gz      
wget https://www.encodeproject.org/files/ENCFF000CFB/@@download/ENCFF000CFB.fastq.gz -O H3K36me3.fastq.gz
 
wget https://www.encodeproject.org/files/ENCFF825XKT/@@download/ENCFF825XKT.fastq.gz -O input_H3K4me1.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001HUV/@@download/ENCFF001HUV.fastq.gz -O input_H3K4me3_H3K27me3.fastq.gz
wget https://www.encodeproject.org/files/ENCFF692GVG/@@download/ENCFF692GVG.fastq.gz -O input_H3K27ac_H3K36me3.fastq.gz

cd ..
mkdir fasta
cd fasta

# Dowload fasta reference
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz 
gunzip *.gz
# Create index
samtools faidx GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta 

cd ..
```

For MACS2, SICER, SPAN peaks launch [ChIP-seq snakemake pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) from fastq files.
Assuming that it is installed to `~/work/chipseq-smk-pipeline`.
```
FDR=0.05; snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
   all --use-conda --cores all --directory $(pwd) --config fastq_ext=fastq.gz \
   fastq_dir=$(pwd)/fastq genome=hg38 \
   macs2_mode=narrow macs2_params="-q $FDR" macs2_suffix="q$FDR" --rerun-incomplete
```

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


# Simulate mixed_reads

```
# Ensure that chips is available!
bash simulate_mixture.sh
mv *.fastq fastq/
```


# Prepare control for peak calling 

1. Launch chipseq pipeline on input files only to obtain bam files. 
2. Filter to chromosome 15
```
for F in $(ls *input*.bam | grep -v chr15); do 
    echo $F;
    samtools index $F; 
    samtools view $F chr15 -b > ${F/.bam/_chr15.bam}; 
done
 
for F in *input*chr15.bam; do 
    echo $F; 
    bedtools bamtofastq -i $F -fq ${F/.bam/.fastq}; 
done
```
3. Copy resulting input fastq files into the `/fastq` folder.
   This step is important, otherwise peak calling will be performed without input.

# Launch peak callers
```
# Perform peak calling using chipseq snakemake pipeline

conda activate snakemake

WORK_DIR=/mnt/stripe/shpynov/2021_chips
GENOME=hg38
  
echo "MACS2 narrow"
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --cores 24 --use-conda --directory $WORK_DIR --config genome=$GENOME \
  fastq_dir=$WORK_DIR/fastq fastq_ext=fastq macs2_mode=narrow macs2_params="-q 0.05" macs2_suffix=q0.05 \
  --rerun-incomplete;
  
echo "MACS2 broad"
snakemake --printshellcmds -s ~/work/chipseq-smk-pipeline/Snakefile \
  all --cores 24 --use-conda --directory $WORK_DIR --config genome=$GENOME \
  fastq_dir=$WORK_DIR/fastq fastq_ext=fastq macs2_mode=broad macs2_params="--broad --broad-cutoff 0.1" macs2_suffix=broad0.1 \
  --rerun-incomplete;
```

# Launch SPAN modifications

```
cd $WORK_DIR
for FDR in 0.05; do 
    snakemake -f SpanModificationsSnakefile all --cores 24 --config fdr=$FDR; 
done
```


#Report

Prepare data report by collecting all the peak calling files and overlap with ground truth. 
```
bash analyze.sh
```

#Visualize results

Launch `analysis.ipynb` jupyter notebook for analysis and visualization of results.
