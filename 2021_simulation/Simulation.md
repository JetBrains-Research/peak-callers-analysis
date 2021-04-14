Simulation
==========

**Chips** (prev. Tulip) is available from bioconda and on [GitHub](https://github.com/gymreklab/chips).

Works only with fasta and peaks on chromosome 15, otherwise fails.

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

# Download peaks data
```
cd /mnt/stripe/shpynov/2021_noise2
wget https://www.encodeproject.org/files/ENCFF366GZW/@@download/ENCFF366GZW.bed.gz -O H3K4me1.bed.gz
wget https://www.encodeproject.org/files/ENCFF651GXK/@@download/ENCFF651GXK.bed.gz -O H3K4me3.bed.gz  
wget https://www.encodeproject.org/files/ENCFF039XWV/@@download/ENCFF039XWV.bed.gz -O H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF666NYB/@@download/ENCFF666NYB.bed.gz -O H3K27me3.bed.gz      
wget https://www.encodeproject.org/files/ENCFF213IBM/@@download/ENCFF213IBM.bed.gz -O H3K36me3.bed.gz
gunzip *.gz
```

# Filter fasta
`samtools faidx /mnt/stripe/shpynov/2021_noise1/fa/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta chr$i  >> GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta;`

# Learn models

### Pattern
`chips learn -b ${bam} -p ${peaks} -t bed -c 5 -o ${name}`

`bash learn.sh`

# Simulate reads

### Pattern
```
chips simreads -p ${peaks} \                                                                                                                                                                   
-f ${ref} \                                                                                                                                                                                   
-o fastqs/${name}_noize${i/./} \                                                                                                                                                              
-t bed -c 5 --numreads ${n_reads} \                                                                                                                                                           
--model ${model} \                                                                                                                                                                            
--scale-outliers --seed 12 --thread 24 --numcopies 100000;
```

`bash simulate.sh`

# Launch peak callers
```
cd /mnt/stripe/shpynov/chipseq-smk-pipeline
for FDR in 0.1 0.05 0.01 1e-3 1e-4 1e-5 1e-6; do
snakemake all --conda-frontend conda --cores 24 --use-conda --directory /mnt/stripe/shpynov/2021_noise2 --config genome=hg38 fastq_dir=/mnt/stripe/shpynov/2021_noise2/fastq fastq_ext=fastq macs2_params="-q $FDR" macs2_mode=narrow macs2_suffix=$FDR
snakemake all --conda-frontend conda --cores 24 --use-conda --directory /mnt/stripe/shpynov/2021_noise2 --config genome=hg38 fastq_dir=/mnt/stripe/shpynov/2021_noise2/fastq fastq_ext=fastq macs2_params="--broad --broad-cutoff $FDR" macs2_suffix=broad$FDR;
snakemake all --conda-frontend conda --cores 24 --use-conda --directory /mnt/stripe/shpynov/2021_noise2 --config genome=hg38 fastq_dir=/mnt/stripe/shpynov/2021_noise2/fastq fastq_ext=fastq span_fdr=$FDR
done
```


# Launch SPAN modifications
SPAN modification `span234.jar` is built from the branch span234.

```
cd /mnt/stripe/shpynov/2021_noise2
for FDR in 0.1 0.05 0.01 1e-3 1e-4 1e-5 1e-6; do snakemake all --cores 12 --config fdr=$FDR; done
```


#Analyze 
Prepare data for visualization and analysis in jupyter notebook.
`bash analyze.sh`
