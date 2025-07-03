#!/bin/bash

#SBATCH --job-name=SSselect        # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./SSselect.out    # standard output file
#SBATCH --error=./SSselect.err     # standard error file
#SBATCH --nodes=1                   # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=5                 # number of cores to allocate
#SBATCH --time=100-00:00:00         # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=20G                   # memory to allocate

VENV_PATH=/data/spliceenv
source $VENV_PATH/bin/activate


while [[ "$#" -gt 0 ]]; do
    case $1 in
        --manifest) MANIFEST="$2"; shift ;;
        --directory) DIRECTORY="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        --condition_column) CONDITION_COLUMN="$2"; shift ;;
        --subtype_column) SUBTYPE_COLUMN="$2"; shift ;;
        --condition_identifier) CONDITION_IDENTIFIER="$2"; shift ;;
        --modes) MODES="$2"; shift ;;
        --ss_tup) SS_TUP="$2"; shift ;;
        --cutoff) CUTOFF="$2"; shift ;;
        --drop_threshold) DROP_THRESHOLD="$2"; shift ;; 
        --individual_drop) INDIVIDUAL_DROP="$2"; shift ;;
        --drop_cut) DROP_CUT="$2"; shift ;;
        --single_selected_file) SINGLE_SELECTED_FILE="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


if [[ -z "$MANIFEST" || -z "$DIRECTORY" || -z "$OUTDIR" || -z "$CONDITION_COLUMN" || -z "$MODES" || -z "$SS_TUP" || -z "$CUTOFF" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --manifest <path> --directory <path> --outdir <path> --condition_column <column> --modes <all|binary|single> --ss_tup <tuple> --cutoff <value> [--subtype_column <column>] [--condition_identifier <value>] [--drop_threshold <value>] [--individual_drop <value>] [--drop_cut <value>] [--single_selected_file <path>]"
    exit 1
fi

if [[ "$MODES" == "both" ]]; then
    MODES="all binary"
fi

if [[ -z "$SUBTYPE_COLUMN" ]]; then
    SUBTYPE_COLUMN=""
fi

if [[ -z "$CONDITION_IDENTIFIER" ]]; then
    CONDITION_IDENTIFIER=""
fi

if [[ -z "$DROP_THRESHOLD" ]]; then
    DROP_THRESHOLD=10
fi

if [[ -z "$INDIVIDUAL_DROP" ]]; then
    INDIVIDUAL_DROP=2
fi

if [[ -z "$DROP_CUT" ]]; then
    DROP_CUT=0
fi


: << COMMENT
sbatch cluster_get_common_SS.sh \
    --manifest <path_to_manifest> \                   <-- Manifest file containing metadata, along with the sample information (file_name) and a condition column
    --directory <path_to_directory> \                 <-- Directory containing the files processed by the tss_altSS.sh script (.altSS, .TSS)  
    --outdir <path_to_output_directory> \             <-- Output directory for the results
    --condition_column <condition_column_name> \      <-- Column in the manifest used to group the samples. It can be any column in the manifest, (e.g., "condition", "subtype", "tissue", "pam50_subtype", etc.).
    --modes <all|binary|single> \                     <-- Using "both" will run [all binary]. Binary is for pairwise comparison of the splice sites, while all is for identifying across all different conditions (condition_column). Use single if only 1 condition (can select with condition_identifier)
    --ss_tup <splice_identifier,length> \             <-- file identifier (3.altSS, 5.altSS, 3.TSS, 5.TSS) and number of splice sites in each sample (0 for median) --  (e.g., 3.altSS,50, 5.TSS,0)
    --cutoff <cutoff_value> \                         <-- Sanity check for current cutoff usage (total sum of reads)
    [--subtype_column <subtype_column_name>] \        <-- Optional: Extra information column appended to final df
    [--condition_identifier <specific_condition>] \   <-- Optional: Use with single mode. Specific condition to filter the samples (e.g., Tumor, Normal, Control, BRCA/PRAD/etc -- if tissues). 
    [--drop_threshold <drop_threshold_value>] \       <-- Optional for single mode: Allowed drop threshold for the splice sites (default: 10) (e.g if 50 total splice sites, after 40 the individual_drop will be used). Initially merges greedy
    [--individual_drop <individual_drop_value>]       <-- Optional for single mode: Allowed drop threshold for the splice sites (default: 2) (e.g if a sample would remove more than 2 splice sites, it is ignored)
    [--drop_cut <drop_cut_value>]                     <-- Optional for single mode: Allowed drop threshold for the splice sites (default: 0) (e.g if a sample would reduce the current splice sites below this threshold it is ignored.)
    [--single_selected_file <path>]                   <-- Optional: Path to a single selected file for processing.

COMMENT


python3 parse_transform_multi.py \
    --manifest "$MANIFEST" \
    --directory "$DIRECTORY" \
    --outdir "$OUTDIR" \
    --condition_column "$CONDITION_COLUMN" \
    ${SUBTYPE_COLUMN:+--subtype_column "$SUBTYPE_COLUMN"} \
    ${CONDITION_IDENTIFIER:+--condition_identifier "$CONDITION_IDENTIFIER"} \
    --modes $MODES \
    --ss_tup "$SS_TUP" \
    --cutoff "$CUTOFF" \
    --drop_threshold "$DROP_THRESHOLD" \
    ${DROP_CUT:+--drop_cut "$DROP_CUT"} \
    --individual_drop "$INDIVIDUAL_DROP" \
    ${SINGLE_SELECTED_FILE:+--single_selected_file "$SINGLE_SELECTED_FILE"}
deactivate



### ------------------------- ### 
########   Runs   ########
# ss_final_out

: << COMMENT
## PVT1 ss3 validation (Control + Tumor) -- data used was generated in the same manner as candidates
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/pvt1/ss_final_out_PVT1/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/pvt1/ss_final_out_PVT1/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype \
    --modes binary \
    --ss_tup 3.altSS,50 --cutoff 10
COMMENT

: << COMMENT
# -- PVT1 Tumor Only -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/pvt1/ss_final_out_PVT1/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/pvt1/ss_final_out_PVT1/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,50 --cutoff 10 \
    --modes single --drop_threshold 10 --individual_drop 1
COMMENT

: << COMMENT
## PVT1 ss5 -- data used was generated in the same manner as candidates
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/pvt1/ss_final_out_PVT1/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/pvt1/ss_final_out_PVT1/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype \
    --modes binary \
    --ss_tup 5.altSS,0 --cutoff 10
COMMENT

: << COMMENT
# -- ERBB4 -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/erbb4/ss_final_out_ERBB4/cutoff_5/altSS \
    --outdir /data/chromatin_associated_genes/erbb4/ss_final_out_ERBB4/cutoff_5/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 5 \
    --modes single --drop_threshold 1 --individual_drop 0
COMMENT

: << COMMENT
# -- FTX -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/ftx/ss_final_out_FTX/cutoff_5/altSS \
    --outdir /data/chromatin_associated_genes/ftx/ss_final_out_FTX/cutoff_5/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,20 --cutoff 5 \
    --modes single --drop_threshold 15 --individual_drop 1
COMMENT

: << COMMENT
# -- RAD51B -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/rad51b/ss_final_out_RAD51B/cutoff_5/altSS \
    --outdir /data/chromatin_associated_genes/rad51b/ss_final_out_RAD51B/cutoff_5/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,10 --cutoff 5 \
    --modes single --drop_threshold 12 --individual_drop 0
COMMENT

: << COMMENT
# -- SLC30A8 -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/slc30a8/ss_final_out_SLC30A8/cutoff_5/altSS \
    --outdir /data/chromatin_associated_genes/slc30a8/ss_final_out_SLC30A8/cutoff_5/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,8 --cutoff 5 \
    --modes single --drop_threshold 6 --individual_drop 0
COMMENT

: << COMMENT
# -- DLEU2 -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/dleu2/ss_final_out_DLEU2/cutoff_5/altSS \
    --outdir /data/chromatin_associated_genes/dleu2/ss_final_out_DLEU2/cutoff_5/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,10 --cutoff 5 \
    --modes single --drop_threshold 8 --individual_drop 1 --drop_cut 8
COMMENT

: << COMMENT
# -- TPRG1 -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/tprg1/ss_final_out_TPRG1/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/tprg1/ss_final_out_TPRG1/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 10 \
    --modes single --drop_threshold 5 --individual_drop 1 --drop_cut 5
COMMENT

: << COMMENT
# -- CASC15 -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/casc15/ss_final_out_CASC15/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/casc15/ss_final_out_CASC15/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 10 \
    --modes single --drop_threshold 15 --individual_drop 1 --drop_cut 10
COMMENT

: << COMMENT
# -- PTPRT -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/ptprt/ss_final_out_PTPRT/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/ptprt/ss_final_out_PTPRT/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 10 \
    --modes single --drop_threshold 10 --individual_drop 1 --drop_cut 10
COMMENT

: << COMMENT
# -- SMYD3 -- #
sbatch cluster_get_common_SS.sh \
    --manifest /data/concat_manifest_with_details_V4.tsv \
    --directory /data/chromatin_associated_genes/smyd3/ss_final_out_SMYD3/cutoff_10/altSS \
    --outdir /data/chromatin_associated_genes/smyd3/ss_final_out_SMYD3/cutoff_10/altSS \
    --condition_column condition --subtype_column pam50_subtype --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 10 \
    --modes single --drop_threshold 5 --individual_drop 1 --drop_cut 6
COMMENT



#===============================================================================================================================================================================#


: << COMMENT
## Prostate PVT1 ss3 ##
sbatch cluster_get_common_SS.sh \
    --manifest /data/tissues/prostate/Prostate_tumor_bam_manifest.txt \
    --directory /data/tissues/prostate/ss_final_out_PVT1/cutoff_10/altSS \
    --outdir /data/tissues/prostate/ss_final_out_PVT1/cutoff_10/altSS \
    --condition_column condition --condition_identifier Tumor \
    --ss_tup 3.altSS,50 --cutoff 10 \
    --modes single --drop_threshold 5 --individual_drop 1 --drop_cut 30
COMMENT

: << COMMENT
## Ovary PVT1 using selected BRCA 34 SS ss3 ##
sbatch cluster_get_common_SS.sh \
    --manifest /data/tissues/ovary/Ovary_tumor_bam_manifest.txt \
    --directory /data/tissues/ovary/ss_final_out_PVT1/cutoff_5/altSS \
    --outdir /data/tissues/ovary/ss_final_out_PVT1/cutoff_5/altSS \
    --condition_column condition --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 10 \
    --modes single --drop_threshold 10 --individual_drop 1
COMMENT

: << COMMENT
## Uterus PVT1 using selected BRCA 34 SS ss3 ##
sbatch cluster_get_common_SS.sh \
    --manifest /data/tissues/uterus/Uterus_tumor_bam_manifest.txt \
    --directory /data/tissues/uterus/ss_final_out_PVT1_combined/cutoff_5/altSS \
    --outdir /data/tissues/uterus/ss_final_out_PVT1_combined/cutoff_5/altSS \
    --condition_column condition --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 5 \
    --modes single_selected --single_selected_file /data/chromatin_associated_genes/pvt1/pvt134ss.list
COMMENT

: << COMMENT
## Testis PVT1 using selected BRCA 34 SS ss3 ##
sbatch cluster_get_common_SS.sh \
    --manifest /data/tissues/testis/Testis_tumor_bam_manifest.txt \
    --directory /data/tissues/testis/ss_final_out_PVT1/cutoff_5/altSS \
    --outdir /data/tissues/testis/ss_final_out_PVT1/cutoff_5/altSS \
    --condition_column condition --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 5 \
    --modes single_selected --single_selected_file /data/chromatin_associated_genes/pvt1/pvt134ss.list
COMMENT

: << COMMENT
## Adrenal gland PVT1 using selected BRCA 34 SS ss3 ##
sbatch cluster_get_common_SS.sh \
    --manifest /data/tissues/adrenal/Adrenal_tumor_bam_manifest.txt \
    --directory /data/tissues/adrenal/ss_final_out_PVT1/cutoff_5/altSS \
    --outdir /data/tissues/adrenal/ss_final_out_PVT1/cutoff_5/altSS \
    --condition_column condition --condition_identifier Tumor \
    --ss_tup 3.altSS,0 --cutoff 5 \
    --modes single_selected --single_selected_file /data/chromatin_associated_genes/pvt1/pvt134ss.list
COMMENT

