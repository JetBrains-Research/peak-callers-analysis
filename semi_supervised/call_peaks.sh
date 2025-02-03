#!/bin/bash

# This bash script runs SPAN for every signal file in the experiment with appropriate input and provided parameters.
# Usage: call_peaks.sh [<FDR> [<GAP>]]
# <FDR> and <GAP> are passed to SPAN calls as --fdr and --gap parameters, respectively.
# The default values are 0.05 and 0

# Locate SPAN jar file.

SPAN=$(ls | grep "span.*\.jar$" | head -n 1)

if [[ -z ${SPAN} ]]; then
	>&2 echo "Error: SPAN jar file missing. Please check that SPAN jar file is in the current directory."
	exit 1;
else
	echo "Found SPAN jar file: $SPAN";
fi;

# Parse FDR and GAP arguments.

if [[ $# -ge 1 ]]; then
	FDR=$1;
else
	FDR=0.05;
fi;

if [[ $# -ge 2 ]]; then
	GAP=$2;
else
	GAP=0;
fi;

# Launch SPAN with given FDR and GAP for every .bed.gz signal file with appropriate input.

for ID in $(find ./chip-seq-benchmark/H3K36me3/ -maxdepth 1 -name 'H3K36me3_*.bed.gz' | sed 's#\./##g; s#^.*H3K36me3_##g; s#\.bed\.gz$##g')
do :
	if [[ -f "./chip-seq-benchmark/H3K36me3/peaks${ID}_${GAP}_$FDR.bed" ]];
	then echo "Skipping existing file peaks${ID}_${GAP}_$FDR.bed" & continue;
	fi
	java -jar ${SPAN} analyze --bed ./chip-seq-benchmark/H3K36me3/peaks${ID}_${GAP}_${FDR}.bed \
	    --cs ./hg19.chrom.sizes --fdr ${FDR} -t ./chip-seq-benchmark/H3K36me3/H3K36me3_${ID}.bed.gz \
	    -c ./chip-seq-benchmark/Input/Input_${ID}.bed.gz --gap ${GAP} --keep-cache
done

for ID in $(find ./chip-seq-benchmark/H3K4me3/ -maxdepth 1 -name 'H3K4me3_*.bed.gz' | sed 's#\./##g; s#^.*H3K4me3_##g; s#\.bed\.gz$##g')
do :
	if [[ -f "./chip-seq-benchmark/H3K4me3/peaks${ID}_${GAP}_$FDR.bed" ]];
	then echo "Skipping existing file peaks${ID}_${GAP}_$FDR.bed" & continue;
	fi
	java -jar ${SPAN} analyze --bed ./chip-seq-benchmark/H3K4me3/peaks${ID}_${GAP}_${FDR}.bed \
	    --cs ./hg19.chrom.sizes --fdr ${FDR} -t ./chip-seq-benchmark/H3K4me3/H3K4me3_${ID}.bed.gz \
	    -c ./chip-seq-benchmark/Input/Input_${ID}.bed.gz --gap ${GAP} --keep-cache
done
