#!/bin/bash

#bed coverage stat: bash bed_coverage_stat.sh A.bed B.bam 检验B.bam在A.bed中的覆盖率

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 A.bed B.bam"
    exit 1
fi

A_bed=$1
B_bam=$2

if [ ! -f "$A_bed" ]; then
    echo "Error: $1 file not found."
    exit 1
fi

if [ ! -f "$B_bam" ]; then
    echo "Error: $2 file not found."
    exit 1
fi

A_length=$(awk '{sum += $3 - $2} END {print sum}' "$A_bed")

if [ "$A_length" -eq 0 ]; then
    echo "Error: A.bed has zero length."
    exit 1
fi

overlap_length=$(bedtools intersect -a "$A_bed" -b "$B_bam" -v | awk '{sum += $3 - $2} END {print sum}')
overlap_ratio=$(echo "scale=4; 1 - $overlap_length / $A_length" | bc)
overlap_ratio_percentage=$(echo "$overlap_ratio * 100" | bc)

echo "Overlap Ratio: $overlap_ratio_percentage%"
echo "A_length: $A_length"
echo "Overlap Length: $overlap_length"

