#!/bin/bash

# Paths
QUALIMAP="./Software/qualimap_v2.3/qualimap"
BAM_DIR="./6_Aligned_reads2"
GTF="./Data/Homo_sapiens.GRCh38.114.gtf"
OUTDIR_BASE="./7_Aligned_reads_quality"
MEM="12G"

mkdir -p "$OUTDIR_BASE"

for bamfile in $BAM_DIR/*.bam; do
    sample=$(basename "$bamfile" .bam)
    sample_outdir="$OUTDIR_BASE/${sample}_qc"

    echo "=== Running Qualimap RNA-seq for $sample ==="
    $QUALIMAP rnaseq \
        -bam "$bamfile" \
        -gtf "$GTF" \
        --outdir "$sample_outdir" \
        --java-mem-size=$MEM

    echo ">>> Finished $sample. Results saved in $sample_outdir"
    echo
done

