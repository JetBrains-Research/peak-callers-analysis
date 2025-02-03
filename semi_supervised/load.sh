#!/bin/bash

# Load and parse the index of .bed.gz files.

if [[ ! -s "./index_wget.txt" ]]
then
	echo "Getting index file..."
	wget "http://hubs.hpc.mcgill.ca/~thocking/bed/index.txt"
	# We only need the first column, since the rest is .bedGraph files and file sizes;
	# we also have to convert relative URLs to absolute ones.
	tail -n +2 ./index.txt | sed 's# .*$##g;s#^#http://hubs.hpc.mcgill.ca/~thocking/bed/#g;' > ./index_wget.txt
fi

while read -r line
do
	# Extract sample ID and ChIP-seq target description from the URL.
	SAMPLE_ID=$(echo ${line} | sed 's#^.*/bed/McGill##g;s#/[^/]*$##g;')
	TARGET=$(echo ${line} | sed 's#^.*/##g;s#.bed.gz$##g;')
	echo "Getting $line..."
	# We have to rename the file from "McGillXXXX" to "TARGET_XXXX", since SPAN caching relies on unique file names.
	wget -c -O ./chip-seq-benchmark/${TARGET}/${TARGET}_${SAMPLE_ID}.bed.gz ${line}
done < "./index_wget.txt"

# Download chrom.sizes file.

wget -c "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes"

# Download 4-fold test data.

wget -c "http://members.cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/4foldcv-test-folds.csv"
