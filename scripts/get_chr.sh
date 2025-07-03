#!/bin/bash

BAM_DIR="$1"
OUTPUT_DIR="$2"
CHROM="$3"
START="$4"
END="$5"

[ -z "$CHROM" ] && CHROM="chr8"

mkdir -p "$OUTPUT_DIR"

for bam in "${BAM_DIR}"/*.bam; do
    bam_name=$(basename "$bam" .bam)
    temp_bam="${OUTPUT_DIR}/${bam_name}_temp.bam"
    temp_sorted_bam="${OUTPUT_DIR}/${bam_name}_sorted.bam"
    clean_name=$(echo "$bam_name" | sed 's/\.rna_seq\.genomic\.gdc_realn//g' | sed 's/-/_/g')
    final_bam="${OUTPUT_DIR}/${clean_name}.bam"
    final_bai="${OUTPUT_DIR}/${clean_name}.bam.bai"

    if [ -n "$START" ] && [ -n "$END" ]; then
        samtools view -b "$bam" "${CHROM}:${START}-${END}" > "$temp_bam"
    else
        samtools view -b "$bam" "$CHROM" > "$temp_bam"
    fi

    if [[ ! -s "$temp_bam" ]]; then
        echo "ERROR: Extraction failed for $bam, skipping..."
        rm -f "$temp_bam"
        continue
    fi
    samtools sort -o "$temp_sorted_bam" "$temp_bam"
    samtools index "$temp_sorted_bam"
    mv "$temp_sorted_bam" "$final_bam"
    mv "${temp_sorted_bam}.bai" "$final_bai"
    rm "$temp_bam"
done
