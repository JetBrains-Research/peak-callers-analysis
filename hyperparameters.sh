#!/usr/bin/env bash

mkdir -p ~/data/2025_hyperparameters
cd ~/data/2025_hyperparameters

# Download SPAN
wget https://download.jetbrains.com/biolabs/span/span-2.0.6610.jar -O span-2.0.jar

# Download chromosome sizes
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes

# Download data from ENCODE
wget https://www.encodeproject.org/files/ENCFF278QPY/@@download/ENCFF278QPY.bam -O GM12878_H3K4me3_rep1.bam
wget https://www.encodeproject.org/files/ENCFF278QPY/@@download/ENCFF278QPY.bam -O GM12878_H3K4me3_rep2.bam
wget https://www.encodeproject.org/files/ENCFF500GXC/@@download/ENCFF500GXC.bam -O GM12878_Input_rep1.bam
wget https://www.encodeproject.org/files/ENCFF402BOP/@@download/ENCFF402BOP.bam -O GM12878_Input_rep2.bam

for SNR in 0 0.001 0.01 0.05 0.1 0.3; do
 for LOW in 0 0.1 0.2 0.3 0.5; do

  echo "$SNR $LOW"

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t GM12878_H3K4me3_rep1.bam -c GM12878_Input_rep1.bam \
    --keep-cache --hmm-snr $SNR --hmm-low $LOW --model GM12878_H3K4me3_rep1_${SNR}_${LOW}.span --debug \
    --peaks GM12878_H3K4me3_rep1_${SNR}_${LOW}.peak --chr chr1 --sensitivity -0.6 --gap 0;

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t GM12878_H3K4me3_rep2.bam -c GM12878_Input_rep2.bam \
    --keep-cache --hmm-snr $SNR --hmm-low $LOW --model GM12878_H3K4me3_rep2_${SNR}_${LOW}.span --debug \
    --peaks GM12878_H3K4me3_rep2_${SNR}_${LOW}.peak --chr chr1 --sensitivity -0.6 --gap 0;

 done;
done | tee k4me3_gm12878.txt;

wget https://www.encodeproject.org/files/ENCFF910QDY/@@download/ENCFF910QDY.bam -O GM12878_H3K36me3_rep1.bam
wget https://www.encodeproject.org/files/ENCFF245LKG/@@download/ENCFF245LKG.bam -O GM12878_H3K36me3_rep2.bam

for SNR in 0 0.001 0.01 0.05 0.1 0.3; do
 for LOW in 0 0.1 0.2 0.3 0.5; do

  echo "$SNR $LOW"

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t GM12878_H3K36me3_rep1.bam -c GM12878_Input_rep1.bam \
    --keep-cache --hmm-snr $SNR --hmm-low $LOW --model GM12878_H3K36me3_rep1_${SNR}_${LOW}.span --debug \
    --peaks GM12878_H3K36me3_rep1_${SNR}_${LOW}.peak --chr chr1 --sensitivity -0.6 --gap 0;

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t GM12878_H3K36me3_rep2.bam -c GM12878_Input_rep2.bam \
    --keep-cache --hmm-snr $SNR --hmm-low $LOW --model GM12878_H3K36me3_rep2_${SNR}_${LOW}.span --debug \
    --peaks GM12878_H3K36me3_rep2_${SNR}_${LOW}.peak --chr chr1 --sensitivity -0.6 --gap 0;

 done;
done | tee k36me3_gm12878.txt;

wget https://www.encodeproject.org/files/ENCFF813OSH/@@download/ENCFF813OSH.bam -O HSMM_H3K36me3_rep1.bam
wget https://www.encodeproject.org/files/ENCFF345NRK/@@download/ENCFF345NRK.bam -O HSMM_H3K36me3_rep2.bam
wget https://www.encodeproject.org/files/ENCFF049ADX/@@download/ENCFF049ADX.bam -O HSMM_Input_rep1.bam
wget https://www.encodeproject.org/files/ENCFF003IHT/@@download/ENCFF003IHT.bam -O HSMM_Input_rep2.bam

for SNR in 0 0.001 0.01 0.05 0.1 0.3; do
 for LOW in 0 0.1 0.2 0.3 0.5; do

  echo "$SNR $LOW"

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t HSMM_H3K36me3_rep1.bam -c HSMM_Input_rep1.bam \
    --keep-cache --hmm-snr $SNR --hmm-low $LOW --model HSMM_H3K36me3_rep1_${SNR}_${LOW}.span --debug \
    --peaks HSMM_H3K36me3_rep1_${SNR}_${LOW}.peak --chr chr1 --sensitivity -0.6 --gap 0;

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t HSMM_H3K36me3_rep2.bam -c HSMM_Input_rep2.bam \
  --keep-cache --hmm-snr $SNR --hmm-low $LOW --model HSMM_H3K36me3_rep2_${SNR}_${LOW}.span --debug \
  --peaks HSMM_H3K36me3_rep2_${SNR}_${LOW}.peak --chr chr1 --sensitivity -0.6 --gap 0;

 done;
done | tee k36me3_hsmm.txt;


wget https://www.encodeproject.org/files/ENCFF862NDZ/@@download/ENCFF862NDZ.bam -O HepG2_H3K27ac_rep1.bam
wget https://www.encodeproject.org/files/ENCFF926NHE/@@download/ENCFF926NHE.bam -O HepG2_H3K27ac_rep2.bam
wget https://www.encodeproject.org/files/ENCFF085CKF/@@download/ENCFF085CKF.bam -O HepG2_Input_rep1.bam
wget https://www.encodeproject.org/files/ENCFF100YPY/@@download/ENCFF100YPY.bam -O HepG2_Input_rep2.bam

for SNR in 0 0.001 0.01 0.05 0.1 0.3; do
 for LOW in 0 0.1 0.2 0.3 0.5; do

  echo "$SNR $LOW"

    java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t HepG2_H3K27ac_rep1.bam -c HepG2_Input_rep1.bam \
   --keep-cache --model HepG2_H3K27ac_rep1_${SNR}_${LOW}.span --debug \
   --peaks HepG2_H3K27ac_rep1_${SNR}_${LOW}.peak --chr chr1 --hmm-snr $SNR --hmm-low $LOW --sensitivity -0.6 --gap 0;

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t HepG2_H3K27ac_rep2.bam -c HepG2_Input_rep2.bam \
   --keep-cache --model HepG2_H3K27ac_rep2_${SNR}_${LOW}.span --debug \
   --peaks HepG2_H3K27ac_rep2_${SNR}_${LOW}.peak --chr chr1 --hmm-snr $SNR --hmm-low $LOW --sensitivity -0.6 --gap 0;

 done;
done | tee k27ac_hepg2.txt;


wget https://www.encodeproject.org/files/ENCFF026YFQ/@@download/ENCFF026YFQ.bam -O CD34_H3K27ac_hg38_ENCFF026YFQ.bam
wget https://www.encodeproject.org/files/ENCFF112IIL/@@download/ENCFF112IIL.bam -O CD34_Control_hg38_ENCFF112IIL.bam


for SNR in 0 0.001 0.01 0.05 0.1 0.3; do
 for LOW in 0 0.1 0.2 0.3 0.5; do

  echo "$SNR $LOW"

  java -jar span-2.0.jar analyze -cs hg38.chrom.sizes -t CD34_H3K27ac_hg38_ENCFF026YFQ.bam \
    -c CD34_Control_hg38_ENCFF112IIL.bam --keep-cache \
    --model CD34_H3K27ac_hg38_ENCFF026YFQ_${SNR}_${LOW}.span --debug \
    --peaks CD34_H3K27ac_hg38_ENCFF026YFQ_${SNR}_${LOW}.peak --chr chr1 \
    --hmm-snr $SNR --hmm-low $LOW --sensitivity -0.6 --gap 0;

 done;
done | tee k27ac_cd34.txt;


wget https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/H3K27ac/OD1_k27ac_hg19.bed.gz
wget https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/H3K27ac/YD11_k27ac_hg19.bed.gz
wget https://artyomovlab.wustl.edu/publications/supp_materials/aging/chipseq/Y20O20/bedgz/input/input.bed.gz

for SNR in 0 0.001 0.01 0.05 0.1 0.3; do
 for LOW in 0 0.1 0.2 0.3 0.5; do

  echo "$SNR $LOW"

  java -jar span-2.0.jar analyze -t OD1_k27ac_hg19.bed.gz --chrom.sizes hg19.chrom.sizes -c input.bed.gz \
    --peaks OD1_k27ac_hg19_${SNR}_${LOW}.peak --model OD1_k27ac_hg19_${SNR}_${LOW}.span --fdr 1e-4 \
    --bin 100 --keep-cache --debug --chr chr1  --hmm-snr $SNR --hmm-low $LOW --sensitivity -0.6 --gap 0;

  java -jar span-2.0.jar analyze -t YD11_k27ac_hg19.bed.gz --chrom.sizes hg19.chrom.sizes -c input.bed.gz \
    --peaks YD11_k27ac_hg19_${SNR}_${LOW}.peak --model YD11_k27ac_hg19_${SNR}_${LOW}.span --fdr 1e-4 \
    --bin 100 --keep-cache --debug --chr chr1 --hmm-snr $SNR --hmm-low $LOW --sensitivity -0.6 --gap 0;

 done;
done | tee k27ac_abf.txt;

