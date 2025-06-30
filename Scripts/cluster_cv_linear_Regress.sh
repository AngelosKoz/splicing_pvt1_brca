#!/bin/bash

#SBATCH --job-name=cvLR          # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./cvLR.out      # standard output file
#SBATCH --error=./cvLR.err       # standard error file
#SBATCH --nodes=1                # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=20              # number of cores to allocate
#SBATCH --time=100-00:00:00      # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=120G               # memory to allocate

VENV_PATH=/data/spliceenv
source $VENV_PATH/bin/activate
export R_LIBS_USER=/home/R/x86_64-pc-linux-gnu-library/4.1
## Check installed R packages in case of error
##Rscript -e '.libPaths(); installed.packages()'

: << COMMENT
sbatch /home/akozonakis/cluster_scripts/cluster_cv_linear_Regress.sh \
    --manifest_file <path_to_manifest> \              <-- Manifest file containing metadata, along with the sample information (file_name) and a condition column
    --counts_file <path_to_counts_file> \             <-- File containing gene or exon counts
    --norm <yes|no> \                                 <-- Normalize the counts. Input should contain "Geneid" as first column, then samples (genes x samples) (RPM --> ceiling at 0.99 --> natural log with +10^-9 peusdocount) (yes or no -- depends on count input). 
                                                          If you prefer different normalization, perform it, then use that as input with --norm no. In this case, df needs to contain "file_name" as first column, then genes (samples x genes)
    --se_file <path_to_splicing_efficiency_file> \    <-- File containing splicing efficiency data (output 3.*SS and 5.*SS from get_tss_altSS.sh)
    --umrs_file <path_to_umrs_file> \                 <-- File containing uniquely mapped reads (columns expected: file_name Total_Reads RPM_NormFactor). Required only when selecting --norm yes
    --output_dir <path_to_output_directory> \         <-- Output directory for the results
    --candidate <candidate_name> \                    <-- Candidate name for splice sites
    --mode <lm | glm | both> \                        <-- Specify the model type (lm or glm)
    [--selected_genes_file <path_to_selected_genes>] \<-- Optional: File containing selected genes to filter
    [--group_column <column_name>] \                  <-- Optional: Column in the manifest used to group the samples. If not used, performs "agnostic" run
    [--group_values <values>] \                       <-- Optional: Values in the group column to filter samples
    [--create_plots <y|yes|off>] \                    <-- Optional: Create plots (default: off)
    [--plot_mode <best|all>] \                        <-- Optional: For plotting observed vs predicted. Best will use best folds samples, while all will use all the observed. This is purely for plotting purposes. (default: best)
    [--breaks <value>] ]                              <-- Optional: The amount of samples to be considered as minimum for a 10-fold cross validation. If that value is unmet, 5-fold cross validation is performed. (default: 500)
    [--custom_lm <path_to_predictors>] \              <-- Optional: Path to a file containing custom predictors for the linear model. This file should contain row separated predictors, one per line. If not provided, the script will use the default predictors based on the input data. Overrides the default predictors.

COMMENT

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --manifest_file) MANIFEST_FILE="$2"; shift ;;
        --counts_file) COUNTS_FILE="$2"; shift ;;
        --norm) NORM="$2"; shift ;;
        --se_file) SE_FILE="$2"; shift ;;
        --umrs_file) UMRS_FILE="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --candidate) CANDIDATE="$2"; shift ;;
        --mode) MODE="$2"; shift ;;
        --selected_genes_file) SELECTED_GENES_FILE="$2"; shift ;;
        --group_column) GROUP_COLUMN="$2"; shift ;;
        --group_values) GROUP_VALUES="$2"; shift ;;
        --create_plots) CREATE_PLOTS="$2"; shift ;;
        --plot_mode) PLOT_MODE="$2"; shift ;;
        --breaks) BREAKS="$2"; shift ;;
        --custom_lm) CUSTOM_LM="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$MANIFEST_FILE" || -z "$COUNTS_FILE" || -z "$NORM" || -z "$SE_FILE" || -z "$OUTPUT_DIR" || -z "$CANDIDATE" || -z "$MODE" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --manifest_file <path> --counts_file <path> --norm <yes|no> --se_file <path> --umrs_file <path> --output_dir <path> --candidate <name> --mode <lm|glm|both> [--selected_genes_file <path>] [--group_column <column>] [--group_values <values>] [--create_plots <y|yes|off>] [--plot_mode <best|all>] [--breaks <value>]"
    exit 1
fi

if [[ -z "$CREATE_PLOTS" ]]; then
    CREATE_PLOTS="off"
fi

if [[ -z "$PLOT_MODE" ]]; then
    PLOT_MODE="best"
fi

if [[ -z "$BREAKS" ]]; then
    BREAKS=500
fi


Rscript refine_cv_linear.r \
    --manifest_file "$MANIFEST_FILE" \
    --counts_file "$COUNTS_FILE" \
    --norm "$NORM" \
    --se_file "$SE_FILE" \
    --umrs_file "$UMRS_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --candidate "$CANDIDATE" \
    --mode "$MODE" \
    ${SELECTED_GENES_FILE:+--selected_genes_file "$SELECTED_GENES_FILE"} \
    ${GROUP_COLUMN:+--group_column "$GROUP_COLUMN"} \
    ${GROUP_VALUES:+--group_values "$GROUP_VALUES"} \
    --create_plots "$CREATE_PLOTS" \
    --plot_mode "$PLOT_MODE" \
    --breaks "$BREAKS" \
    ${CUSTOM_LM:+--custom_lm "$CUSTOM_LM"}

deactivate
echo "Job completed successfully."
echo "Output files are located in $OUTPUT_DIR"


### ------------------------- ### 
########   Runs   ########

# -- PVT1 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/pvt1/condition_altSS_concat_Control_Tumor_SE_SS3.altSS_L50.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/pvt1/run_cvlm \
    --candidate PVT1 --group_column condition --mode both --group_values Tumor
COMMENT

# -- ERBB4 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/erbb4/Tumor_SE_SS3.altSS_L25_ad1id0dc0c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/erbb4/run_cvlm \
    --candidate ERBB4 --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- FTX -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/ftx/Tumor_SE_SS3.altSS_L20_ad15id1dc0c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/ftx/run_cvlm \
    --candidate FTX --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- RAD51B -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/rad51b/Tumor_SE_SS3.altSS_L10_ad12id0dc0c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/rad51b/run_cvlm \
    --candidate RAD51B --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- SLC30A8 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/slc30a8/Tumor_SE_SS3.altSS_L8_ad6id0c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/slc30a8/run_cvlm \
    --candidate SLC30A8 --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- PTPRT -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/PVT1_RPMlog_counts.tsv \
    --se_file /data/chromatin_associated_genes/ptprt/Tumor_SE_SS3.altSS_L30_ad10id1dc10c10.tsv \
    --norm no \
    --output_dir /data/chromatin_associated_genes/ptprt/run_cvlm \
    --candidate PTPRT --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- SMYD3 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/PVT1_RPMlog_counts.tsv \
    --se_file /data/chromatin_associated_genes/smyd3/Tumor_SE_SS3.altSS_L18_ad5id1dc6c10.tsv \
    --norm no \
    --output_dir /data/chromatin_associated_genes/smyd3/run_cvlm \
    --candidate SMYD3 --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- DLEU2 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/dleu2/DLEU2_Tumor_SE_SS3.altSS_L10_ad8id1dc8c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/dleu2/run_cvlm \
    --candidate DLEU2 --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- CASC15 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/casc15/CASC15_Tumor_SE_SS3.altSS_L11_ad15id1dc14c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/casc15/run_cvlm \
    --candidate CASC15 --group_column condition --mode glmnet --group_values Tumor
COMMENT

# -- TPRG1 -- #
: << COMMENT
sbatch cluster_cv_linear_Regress.sh \
    --manifest_file /data/concat_manifest_with_details_V4.tsv \
    --counts_file /data/counts/exonLevel_allCounts_filtbyexpr.tsv \
    --se_file /data/chromatin_associated_genes/tprg1/TPRG1_Tumor_SE_SS3.altSS_L20_ad5id1dc5c5.tsv \
    --norm yes \
    --umrs_file /data/umrs/brca_uniquely_mapped_reads_all.txt \
    --output_dir /data/chromatin_associated_genes/tprg1/run_cvlm \
    --candidate TPRG1 --group_column condition --mode glmnet --group_values Tumor
COMMENT