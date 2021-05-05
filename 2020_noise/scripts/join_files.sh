#!/bin/bash


if [[ $1 == "-h" || $1 == "--help" ]]; then
  echo "Program to mix ChIP-seq track and control aligned reads.
  Usage: [track file (.bam)] [control file (.bam)] [noise rate (0..9)] [output_file (.bam)]"
  exit 0
fi

if [[ "$#" -ne 4 ]]; then
  echo "Incorrect number of arguments (use -h for help)"
  exit 1
fi

TRACK=$1
CONTROL=$2
NOISE=$3
OUTPUT=$4

TRACK_LINES=$(sambamba view -t 4 -c "$TRACK")
CONTROL_LINES=$(sambamba view -t 4 -c "$CONTROL")

TMP_TRACK=$(mktemp -t track-XXXXXXXXXX.bam)
TRACK_SAMPLE_FRACTION=$(bc -l <<< "1-$NOISE/10")
sambamba view -f bam -o "$TMP_TRACK" -t 4 -s $TRACK_SAMPLE_FRACTION "$TRACK"

TMP_CONTROL=$(mktemp -t control-XXXXXXXXXX.bam)
CONTROL_SAMPLE_FRACTION=$(bc -l <<< "$TRACK_LINES*$NOISE/$CONTROL_LINES/10")
sambamba view -f bam -o "$TMP_CONTROL" -t 4 -s $CONTROL_SAMPLE_FRACTION "$CONTROL"

sambamba merge -t 4 "$OUTPUT" "$TMP_TRACK" "$TMP_CONTROL"
rm "$TMP_TRACK" "$TMP_CONTROL"