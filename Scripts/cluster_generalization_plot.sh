#!/bin/bash

#SBATCH --job-name=genPlot          # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./genPlot.out      # standard output file
#SBATCH --error=./genPlot.err       # standard error file
#SBATCH --nodes=1                   # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=5                 # number of cores to allocate
#SBATCH --time=100-00:00:00         # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=20G                   # memory to allocate


VENV_PATH=/data/spliceenv
source $VENV_PATH/bin/activate
export R_LIBS_USER=/home/R/x86_64-pc-linux-gnu-library/4.1

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --calculated_metrics) CALCULATED_METRICS="$2"; shift ;;
        --gen_metrics) GEN_METRICS="$2"; shift ;;
        --gen_col) GEN_COL="$2"; shift ;;
        --gen_plot) GEN_PLOT="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --candidate) CANDIDATE="$2"; shift ;;
        --mode) MODE="$2"; shift ;;
        --input_mode) INPUT_MODE="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


if [[ -z "$CALCULATED_METRICS" || -z "$GEN_METRICS" || -z "$GEN_COL" || -z "$GEN_PLOT" || -z "$OUTPUT_DIR" || -z "$CANDIDATE" || -z "$MODE" || -z "$INPUT_MODE" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --calculated_metrics <path> --gen_metrics <path> --gen_col <column_name> --gen_plot <rmse|rsq> --output_dir <path> --candidate <name> --mode <mode> --input_mode <input_mode>"
    exit 1
fi
mkdir -p "$OUTPUT_DIR"

Rscript generalization_plot.R \
    --calculated_metrics "$CALCULATED_METRICS" \
    --gen_metrics "$GEN_METRICS" \
    --gen_col "$GEN_COL" \
    --gen_plot "$GEN_PLOT" \
    --output_dir "$OUTPUT_DIR" \
    --candidate "$CANDIDATE" \
    --mode "$MODE" \
    --input_mode "$INPUT_MODE"

deactivate

echo "Generalization plot job completed successfully."
echo "Output files are located in $OUTPUT_DIR"




### ------------------------- ### 
########   Runs   ########

## Columns in BRCA GLMNET results: mean_raw_r_sq   best_raw_r_sq   mean_rmse   min_rmse

## RMSE -- Ovary
: << COMMENT
mkdir -p /data/tissues/ovary/gen_ovaryRes/full && cd /data/tissues/ovary/gen_ovaryRes/full
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/ovary/gen_ovaryRes/generalized_results/PVT1_ovary_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_ovary --mode glmnet --input_mode coeff
COMMENT

## RMSE -- Ovary -- 365
: << COMMENT
mkdir -p /data/tissues/ovary/gen_ovaryRes/brca_365 && cd /data/tissues/ovary/gen_ovaryRes/brca_365
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/ovary/gen_ovaryRes/generalized_results/PVT1_ovary_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/brca_glmnet_365_metrics.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_ovary --mode glmnet --input_mode coeff
COMMENT

#===============================================================================================================================================================================#


## RMSE -- Prostate
: << COMMENT
mkdir -p /data/tissues/prostate/gen_prostateRes/full && cd /data/tissues/prostate/gen_prostateRes/full
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/prostate/gen_prostateRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_prostate --mode glmnet --input_mode coeff
COMMENT

## RMSE -- Prostate -- 365
: << COMMENT
mkdir -p /data/tissues/prostate/gen_prostateRes/brca_365 && cd /data/tissues/prostate/gen_prostateRes/brca_365
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/prostate/gen_prostateRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/brca_glmnet_365_metrics.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_prostate --mode glmnet --input_mode coeff
COMMENT

#===============================================================================================================================================================================#


## RMSE -- Uterus
: << COMMENT
mkdir -p /data/tissues/uterus/gen_uterusRes/full && cd /data/tissues/uterus/gen_uterusRes/full
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/uterus/gen_uterusRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_uterus --mode glmnet --input_mode coeff
COMMENT

## RMSE -- Uterus -- 365
: << COMMENT
mkdir -p /data/tissues/uterus/gen_uterusRes/brca_365 && cd /data/tissues/uterus/gen_uterusRes/brca_365
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/uterus/gen_uterusRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/brca_glmnet_365_metrics.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_uterus --mode glmnet --input_mode coeff
COMMENT

#===============================================================================================================================================================================#


## RMSE -- Testis
: << COMMENT
mkdir -p /data/tissues/testis/gen_testisRes/full && cd /data/tissues/testis/gen_testisRes/full
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/testis/gen_testisRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_testis --mode glmnet --input_mode coeff
COMMENT

## RMSE -- Testis -- 365
: << COMMENT
mkdir -p /data/tissues/testis/gen_testisRes/brca_365 && cd /data/tissues/testis/gen_testisRes/brca_365
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/testis/gen_testisRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/brca_glmnet_365_metrics.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_testis --mode glmnet --input_mode coeff
COMMENT

#===============================================================================================================================================================================#


## RMSE -- Adrenal gland
: << COMMENT
mkdir -p /data/tissues/adrenal/gen_adrenalRes/full && cd /data/tissues/adrenal/gen_adrenalRes/full
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/adrenal/gen_adrenalRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_adrenal --mode glmnet --input_mode coeff
COMMENT

## RMSE -- Adrenal gland -- 365
: << COMMENT
mkdir -p /data/tissues/adrenal/gen_adrenalRes/brca_365 && cd /data/tissues/adrenal/gen_adrenalRes/brca_365
sbatch cluster_generalization_plot.sh \
    --calculated_metrics /data/tissues/adrenal/gen_adrenalRes/generalized_results/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --gen_metrics /data/chromatin_associated_genes/pvt1/brca_glmnet_365_metrics.txt \
    --gen_col min_rmse --gen_plot rmse \
    --output_dir . \
    --candidate PVT1_adrenal --mode glmnet --input_mode coeff
COMMENT

