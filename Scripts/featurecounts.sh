#!/bin/bash

# Give path to your aligned reads (.bam files)
cd /home/vaibhav/portfolio/RNA-seq/Data/7_aligned_reads

# Loop over all BAM files
for bam in *.bam; doi
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a /home/vaibhav/portfolio/RNA-seq/Data/Homo_sapiens.GRCh38.114.gtf \
        -o /home/vaibhav/RNA-seq/portfolio/RNA-seq/Data/8_Read_count/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done
