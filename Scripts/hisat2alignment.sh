#!/bin/bash

FASTQ_DIR="/home/vaibhav/portfolio/RNA-seq/3_Fastq_samplewise"
GENOME_INDEX="/home/vaibhav/portfolio/RNA-seq/Data/grch38/genome"
OUTDIR="/home/vaibhav/portfolio/RNA-seq/6_Aligned_reads/"
LOGFILE="/home/vaibhav/portfolio/RNA-seq/Data/alignment.log"


> $LOGFILE

FILES=(
    "LNCAP_Hypoxia_S1.fastq.gz"
    "LNCAP_Hypoxia_S2.fastq.gz"
    "LNCAP_Normoxia_S1.fastq.gz"
    "LNCAP_Normoxia_S2.fastq.gz"
    "PC3_Hypoxia_S1.fastq.gz"
    "PC3_Hypoxia_S2.fastq.gz"
    "PC3_Normoxia_S1.fastq.gz"
    "PC3_Normoxia_S2.fastq.gz"
)

for f in "${FILES[@]}"; do
    SAMPLE=$(basename "$f" .fastq.gz)

    echo "Processing $SAMPLE at $(date)" | tee -a $LOGFILE
    START=$(date +%s)

    hisat2 -p 8 -q -x $GENOME_INDEX -U $FASTQ_DIR/$f | \
    samtools sort -@ 8 -o $OUTDIR/${SAMPLE}.bam

    samtools index $OUTDIR/${SAMPLE}.bam

    END=$(date +%s)
    echo "Finished $SAMPLE in $((END - START)) sec" | tee -a $LOGFILE
    echo "------------------------------------" | tee -a $LOGFILE
done

