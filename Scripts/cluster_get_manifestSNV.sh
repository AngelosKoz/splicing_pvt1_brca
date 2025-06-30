#!/bin/bash

#SBATCH --job-name=get_manifestSNV  # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./get_manifestSNV.out  # standard output file
#SBATCH --error=./get_manifestSNV.err   # standard error file
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks=4                      # number of cores to allocate
#SBATCH --time=10:00:00                 # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=16G                       # memory to allocate

VENV_PATH=/data/spliceenv
source $VENV_PATH/bin/activate
export R_LIBS_USER=/home/R/x86_64-pc-linux-gnu-library/4.1

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam_manifest) BAM_MANIFEST="$2"; shift ;;
        --json_bam) JSON_BAM="$2"; shift ;;
        --json_snv) JSON_SNV="$2"; shift ;;
        --id) ID="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --wxs_file) WXS_FILE="$2"; shift ;;
        --wgs_file) WGS_FILE="$2"; shift ;;
        --mode) MODE="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$BAM_MANIFEST" || -z "$JSON_BAM" || -z "$ID" || -z "$OUTPUT_DIR" || -z "$MODE" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --bam_manifest <path> --json_bam <path> --id <ID> --output_dir <path> --mode <wxs|wgs|both|simple> [--json_snv <path>] [--wxs_file <path>] [--wgs_file <path>]"
    exit 1
fi

if [[ "$MODE" == "wxs" && -z "$WXS_FILE" ]]; then
    echo "Error: --wxs_file is required when --mode is 'wxs'."
    exit 1
fi
if [[ "$MODE" == "wgs" && -z "$WGS_FILE" ]]; then
    echo "Error: --wgs_file is required when --mode is 'wgs'."
    exit 1
fi
if [[ "$MODE" == "both" && ( -z "$WXS_FILE" || -z "$WGS_FILE" ) ]]; then
    echo "Error: Both --wxs_file and --wgs_file are required when --mode is 'both'."
    exit 1
fi

# Preprocess WGS file if mode is wgs or both
if [[ "$MODE" == "wgs" || "$MODE" == "both" ]]; then
    FILTERED_WGS_FILE="$OUTPUT_DIR/filt_$(basename "$WGS_FILE")"
    FILTERED_WGS_DIR=$(dirname "$FILTERED_WGS_FILE")
    
    # Ensure the output directory for the filtered WGS file exists
    mkdir -p "$FILTERED_WGS_DIR"
    
    echo "Preprocessing WGS file: $WGS_FILE"
    awk -v rs_counter=1 'BEGIN {FS=OFS="\t"} NR == 1 {print "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "RS", "FORMAT", "NORMAL", "TUMOR", "file_name_wgs"; next}{
        if ($7 == "PASS") {
            rs = "rsNA" rs_counter
            if ($8 ~ /rs[0-9]+/) {
                match($8, /rs[0-9]+/, id); rs = id[0]
            } else {
                rs_counter++
            }
            print $1, $2, $3, $4, $5, $6, $7, rs, $9, $10, $11, $12
        }
    }' "$WGS_FILE" > "$FILTERED_WGS_FILE"
    
    echo "Filtered WGS file saved at: $FILTERED_WGS_FILE"
    WGS_FILE="$FILTERED_WGS_FILE"  # Update WGS_FILE to point to the filtered file
fi


Rscript refine_manifestSNV.R \
    --bam_manifest "$BAM_MANIFEST" \
    --json_bam "$JSON_BAM" \
    ${JSON_SNV:+--json_snv "$JSON_SNV"} \
    --id "$ID" \
    --output_dir "$OUTPUT_DIR" \
    --mode "$MODE" \
    ${WXS_FILE:+--wxs_file "$WXS_FILE"} \
    ${WGS_FILE:+--wgs_file "$WGS_FILE"}

deactivate



### ------------------------- ### 
########   Runs   ########

: <<COMMENT
## Prostate 
sbatch cluster_get_manifestSNV.sh \
    --bam_manifest /data/tissues/tcga_metadata/Prostate/bam/bam_prostate_tumor_gdc_manifest.2025-02-26.213511.txt \
    --json_bam /data/tissues/tcga_metadata/Prostate/bam/bam_prostate_metadata.repository.2025-02-26.json \
    --id Prostate_tumor \
    --output_dir . \
    --mode simple
COMMENT

: <<COMMENT
## Ovary
sbatch cluster_get_manifestSNV.sh \
    --bam_manifest /data/tissues/tcga_metadata/Ovary/bam/bam_ovary_tumor_gdc_manifest.2025-02-25.txt \
    --json_bam /data/tissues/tcga_metadata/Ovary/bam/bam_ovary_metadata.repository.2025-02-25.json \
    --id Ovary_tumor \
    --output_dir . \
    --mode simple
COMMENT

: <<COMMENT
## Uterus
sbatch cluster_scripts/cluster_get_manifestSNV.sh \
    --bam_manifest /data/tissues/tcga_metadata/Uterus/bam/bam_uterus_tumor_gdc_manifest.2025-02-26.230133.txt \
    --json_bam /data/tissues/tcga_metadata/Uterus/bam/bam_uterus_metadata.repository.2025-02-26.json \
    --id Uterus_tumor \
    --output_dir . \
    --mode simple
COMMENT

: <<COMMENT
## Testis
sbatch cluster_scripts/cluster_get_manifestSNV.sh \
    --bam_manifest /data/tissues/tcga_metadata/Testis/bam/bam_testis_tumor_gdc_manifest.2025-02-26.225629.txt \
    --json_bam /data/tissues/tcga_metadata/Testis/bam/bam_testis_tumor_metadata.repository.2025-02-26.json \
    --id Testis_tumor \
    --output_dir . \
    --mode simple
COMMENT

: <<COMMENT
## Adrenal Gland
sbatch cluster_scripts/cluster_get_manifestSNV.sh \
    --bam_manifest /data/tissues/tcga_metadata/Adrenal_gland/bam/bam_adrenal_gdc_manifest.2025-02-26.222333.txt \
    --json_bam /data/tissues/tcga_metadata/Adrenal_gland/bam/bam_adrenal_tumor_metadata.repository.2025-02-26.json \
    --id Adrenal_gland_tumor \
    --output_dir . \
    --mode simple
COMMENT
