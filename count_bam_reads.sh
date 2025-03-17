#!/bin/bash
for i in *.bam; do
    echo "${i}" >> reads_counts.txt
    samtools view -c ${i} >> reads_counts.txt
done