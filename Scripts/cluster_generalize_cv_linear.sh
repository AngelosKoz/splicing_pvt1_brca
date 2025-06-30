#!/bin/bash

#SBATCH --job-name=cvGenLR          # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./cvGenLR.out      # standard output file
#SBATCH --error=./cvGenLR.err       # standard error file
#SBATCH --nodes=1                # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=30              # number of cores to allocate
#SBATCH --time=100-00:00:00      # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=160G               # memory to allocate

VENV_PATH=/data/spliceenv
source $VENV_PATH/bin/activate
export R_LIBS_USER=/home/R/x86_64-pc-linux-gnu-library/4.1


: << COMMENT
sbatch /home/akozonakis/cluster_scripts/cluster_generalize_cv_linear.sh \
    --manifest_file <path_to_manifest> \              <-- Manifest file containing metadata, along with the sample information (file_name) and a condition column
    --coefficients_file <path_to_coefficients_file> \ <-- File containing the coefficients from the linear regression model (output from cv lm and cv glmnet models)
    --model_obj <path_to_model_obj> \                 <-- File containing pre-trained models (optional, for 'object' input_mode)
    --counts_file <path_to_counts_file> \             <-- File containing gene or exon counts
    --norm <yes|no> \                                 <-- Normalize the counts. Input should contain "Geneid" as first column, then samples (genes x samples) (RPM --> ceiling at 0.99 --> natural log with +10^-9 peusdocount) (yes or no -- depends on count input). 
                                                          If you prefer different normalization, perform it, then use that as input with --norm no. In this case, df needs to contain "file_name" as first column, then genes (samples x genes)
    --se_file <path_to_splicing_efficiency_file> \    <-- File containing splicing efficiency data (output 3.*SS and 5.*SS from get_tss_altSS.sh)
    --umrs_file <path_to_umrs_file> \                 <-- File containing uniquely mapped reads (columns expected: file_name Total_Reads RPM_NormFactor). Required only when selecting --norm yes
    --output_dir <path_to_output_directory> \         <-- Output directory for the results
    --candidate <candidate_name> \                    <-- Candidate name for splice sites
    --input_mode <coeff|object> \                     <-- Specify the input mode ('coeff' for coefficients file or 'object' for pre-trained models)
    --mode <lm | glm | both> \                        <-- Specify the model type (lm or glm)
    [--selected_genes_file <path_to_selected_genes>] \<-- Optional: File containing selected genes to filter
    [--create_plots <y|yes|no>] \                    <-- Optional: Create plots (default: no)
    [--gen_plot <rmse|rsq|no>] \                     <-- Optional: Generate a generalization plot ('rmse', 'rsq', or 'no' -- default: no)
    [--gen_metrics <path_to_metrics_file>] \         <-- Optional: Path to the metrics file for comparison (required if --gen_plot is selected)
    [--gen_col <column_name>] \                      <-- Optional: Column name in the metrics file to use for comparison (required if --gen_plot is selected)

COMMENT

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --manifest_file) MANIFEST_FILE="$2"; shift ;;
        --coefficients_file) COEFFICIENTS_FILE="$2"; shift ;;
        --model_obj) MODEL_OBJ="$2"; shift ;;
        --counts_file) COUNTS_FILE="$2"; shift ;;
        --norm) NORM="$2"; shift ;;
        --se_file) SE_FILE="$2"; shift ;;
        --umrs_file) UMRS_FILE="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --candidate) CANDIDATE="$2"; shift ;;
        --input_mode) INPUT_MODE="$2"; shift ;;
        --mode) MODE="$2"; shift ;;
        --selected_genes_file) SELECTED_GENES_FILE="$2"; shift ;;
        --create_plots) CREATE_PLOTS="$2"; shift ;;
        --ss_table) SS_TABLE="$2"; shift ;;
        --gen_plot) GEN_PLOT="$2"; shift ;;
        --gen_metrics) GEN_METRICS="$2"; shift ;;
        --gen_col) GEN_COL="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$MANIFEST_FILE" || -z "$COUNTS_FILE" || -z "$NORM" || -z "$SE_FILE" || -z "$OUTPUT_DIR" || -z "$CANDIDATE" || -z "$INPUT_MODE" || -z "$MODE" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --manifest_file <path> --coefficients_file <path> --model_obj <path> --counts_file <path> --norm <yes|no> --se_file <path> --umrs_file <path> --output_dir <path> --candidate <name> --input_mode <coeff|object> --mode <lm|glm|both> [--selected_genes_file <path>] [--create_plots <y|yes|no>] [--gen_plot <rmse|rsq|no>] [--gen_metrics <path>] [--gen_col <column_name>]"
    exit 1
fi

if [[ "$INPUT_MODE" == "coeff" && -z "$COEFFICIENTS_FILE" ]]; then
    echo "Error: --coefficients_file must be provided when --input_mode is 'coeff'."
    exit 1
fi

if [[ "$INPUT_MODE" == "object" && -z "$MODEL_OBJ" ]]; then
    echo "Error: --model_obj must be provided when --input_mode is 'object'."
    exit 1
fi

if [[ "$GEN_PLOT" != "no" && ( -z "$GEN_METRICS" || -z "$GEN_COL" ) ]]; then
    echo "Error: --gen_metrics and --gen_col must be provided when --gen_plot is selected."
    exit 1
fi

if [[ -z "$CREATE_PLOTS" ]]; then
    CREATE_PLOTS="no"
fi

mkdir -p "$OUTPUT_DIR"
echo "Arguments passed check"
Rscript refine_generalize_linearCV.r \
    --manifest_file "$MANIFEST_FILE" \
    ${COEFFICIENTS_FILE:+--coefficients_file "$COEFFICIENTS_FILE"} \
    ${MODEL_OBJ:+--model_obj "$MODEL_OBJ"} \
    --counts_file "$COUNTS_FILE" \
    --norm "$NORM" \
    --se_file "$SE_FILE" \
    ${UMRS_FILE:+--umrs_file "$UMRS_FILE"} \
    --output_dir "$OUTPUT_DIR" \
    --candidate "$CANDIDATE" \
    --input_mode "$INPUT_MODE" \
    --mode "$MODE" \
    ${SELECTED_GENES_FILE:+--selected_genes_file "$SELECTED_GENES_FILE"} \
    ${SS_TABLE:+--ss_table "$SS_TABLE"} \
    --create_plots "$CREATE_PLOTS" \
    ${GEN_PLOT:+--gen_plot "$GEN_PLOT"} \
    ${GEN_METRICS:+--gen_metrics "$GEN_METRICS"} \
    ${GEN_COL:+--gen_col "$GEN_COL"}

deactivate

echo "Job completed successfully."
echo "Output files are located in $OUTPUT_DIR"



### ------------------------- ### 
########   Runs   ########



# Pre-normalized (glmnet) -- Ovary
: << COMMENT
sbatch cluster_generalize_cv_linear.sh \
    --manifest_file /data/tissues/ovary/Ovary_tumor_manifest.txt \
    --coefficients_file /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_coefficients_brcatumor.txt \
    --counts_file /data/tissues/ovary/Ovary_RPMlog_counts.tsv \
    --norm no \
    --se_file /data/tissues/ovary/Ovary_tumor_SE_SS3.altSS_L53_ad10id1c10.tsv \
    --output_dir /data/tissues/ovary \
    --candidate PVT1_ovary --input_mode coeff --mode glmnet \
    --gen_plot rmse \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse
COMMENT

# Pre-normalized (glmnet) -- Prostate
: << COMMENT
sbatch cluster_generalize_cv_linear.sh \
    --manifest_file /data/tissues/prostate/Prostate_tumor_bam_manifest.txt \
    --coefficients_file /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_coefficients_brcatumor.txt \
    --counts_file /data/tissues/prostate/Prostate_RPMlog_counts.tsv \
    --norm no \
    --se_file /data/tissues/prostate/Tumor_SE_SS3.altSS_L50_ad5id1dc30c10.tsv \
    --output_dir /data/tissues/prostate \
    --candidate PVT1_prostate --input_mode coeff --mode glmnet \
    --gen_plot rmse \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse
COMMENT

# Pre-normalized (glmnet) -- Uterus
: << COMMENT
sbatch cluster_generalize_cv_linear.sh \
    --manifest_file /data/tissues/uterus/Uterus_tumor_bam_manifest.txt \
    --coefficients_file /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_coefficients_brcatumor.txt \
    --counts_file /data/tissues/uterus/Uterus_RPMlog_counts.tsv \
    --norm no \
    --se_file /data/tissues/uterus/single_selected_SE_SS3.altSS_L27_single_selected.tsv \
    --output_dir /data/tissues/uterus \
    --candidate PVT1_uterus --input_mode coeff --mode glmnet \
    --gen_plot rmse \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse
COMMENT

# Pre-normalized (glmnet) -- Testis
: << COMMENT
sbatch cluster_generalize_cv_linear.sh \
    --manifest_file /data/tissues/testis/Testis_tumor_bam_manifest.txt \
    --coefficients_file /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_coefficients_brcatumor.txt \
    --counts_file /data/tissues/testis/Testis_RPMlog_counts.tsv \
    --norm no \
    --se_file /data/tissues/testis/single_selected_SE_SS3.altSS_L36_single_selected.tsv \
    --output_dir /data/tissues/testis \
    --candidate PVT1_testis --input_mode coeff --mode glmnet \
    --gen_plot rmse \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse
COMMENT

# Pre-normalized (glmnet) -- Adrenal gland
: << COMMENT
sbatch cluster_generalize_cv_linear.sh \
    --manifest_file /data/tissues/adrenal/Adrenal_gland_tumor_bam_manifest.txt \
    --coefficients_file /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_coefficients_brcatumor.txt \
    --counts_file /data/tissues/adrenal/Adrenal_RPMlog_counts.tsv \
    --norm no \
    --se_file /data/tissues/adrenal/single_selected_SE_SS3.altSS_L43_single_selected.tsv \
    --output_dir /data/tissues/adrenal \
    --candidate PVT1_adrenal --input_mode coeff --mode glmnet \
    --gen_plot rmse \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse
COMMENT

