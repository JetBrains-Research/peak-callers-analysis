#!/usr/bin/env bash

mkdir -p ~/data/2025_peps
cd ~/data/2025_peps

# Download OMNIPEAK
OMNIPEAK_JAR=omnipeak.jar
wget https://download.jetbrains.com/biolabs/omnipeak/omnipeak-1.0.6679.jar -O $OMNIPEAK_JAR

# Download chromosome sizes
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes

# Download data from ENCODE
wget https://www.encodeproject.org/files/ENCFF278QPY/@@download/ENCFF278QPY.bam -O GM12878_H3K4me3_rep1.bam
wget https://www.encodeproject.org/files/ENCFF910QDY/@@download/ENCFF910QDY.bam -O GM12878_H3K36me3_rep1.bam
wget https://www.encodeproject.org/files/ENCFF500GXC/@@download/ENCFF500GXC.bam -O GM12878_Input_rep1.bam

# Launch peak calling
# H3K4me3 - Default
java --add-modules=jdk.incubator.vector -jar $OMNIPEAK_JAR analyze -cs hg38.chrom.sizes \
  -t GM12878_H3K4me3_rep1.bam -c GM12878_Input_rep1.bam \
  --keep-cache --model GM12878_H3K4me3_rep1.omni --debug \
  --peaks GM12878_H3K4me3_rep1.peak --chr chr1;
# H3K4me3 - Stringent PEP sensitivity threshold
java --add-modules=jdk.incubator.vector -jar $OMNIPEAK_JAR analyze -cs hg38.chrom.sizes \
  -t GM12878_H3K4me3_rep1.bam -c GM12878_Input_rep1.bam \
  --keep-cache --model GM12878_H3K4me3_rep1.omni --debug \
  --peaks GM12878_H3K4me3_rep1_pep-50.peak --chr chr1 --sensitivity -50;
# H3K4me3 - Relaxed PEP sensitivity threshold
java --add-modules=jdk.incubator.vector -jar $OMNIPEAK_JAR analyze -cs hg38.chrom.sizes \
  -t GM12878_H3K4me3_rep1.bam -c GM12878_Input_rep1.bam \
  --keep-cache --model GM12878_H3K4me3_rep1.omni --debug \
  --peaks GM12878_H3K4me3_rep1_pep-1e-5.peak --chr chr1 --sensitivity -1e-5;

# H3K36me3 - Default
java --add-modules=jdk.incubator.vector -jar $OMNIPEAK_JAR analyze -cs hg38.chrom.sizes \
  -t GM12878_H3K36me3_rep1.bam -c GM12878_Input_rep1.bam \
  --keep-cache --model GM12878_H3K36me3_rep1.omni --debug \
  --peaks GM12878_H3K36me3_rep1.peak --chr chr1;
# H3K36me3 - Stringent PEP sensitivity threshold
java --add-modules=jdk.incubator.vector -jar $OMNIPEAK_JAR analyze -cs hg38.chrom.sizes \
  -t GM12878_H3K36me3_rep1.bam -c GM12878_Input_rep1.bam \
  --keep-cache --model GM12878_H3K36me3_rep1.omni --debug \
  --peaks GM12878_H3K36me3_rep1_pep-50.peak --chr chr1 --sensitivity -50;
# H3K36me3 - Relaxed PEP sensitivity threshold
java --add-modules=jdk.incubator.vector -jar $OMNIPEAK_JAR analyze -cs hg38.chrom.sizes \
  -t GM12878_H3K36me3_rep1.bam -c GM12878_Input_rep1.bam \
  --keep-cache --model GM12878_H3K36me3_rep1.omni --debug \
  --peaks GM12878_H3K36me3_rep1_pep-1e-5.peak --chr chr1 --sensitivity -1e-5;

echo "DONE";

