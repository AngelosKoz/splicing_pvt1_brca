#!/bin/bash

#SBATCH --job-name=enrich          # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./enrich.out      # standard output file
#SBATCH --error=./enrich.err       # standard error file
#SBATCH --nodes=1                  # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=4                 # number of cores to allocate
#SBATCH --time=10:00:00            # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=16G                  # memory to allocate

VENV_PATH=/data/spliceenv # Python environment path, adjust as needed
source $VENV_PATH/bin/activate
export R_LIBS_USER=/home/user/R/x86_64-pc-linux-gnu-library/4.1 # R environment path, adjust as needed

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_file) INPUT_FILE="$2"; shift ;;
        --target_file) TARGET_FILE="$2"; shift ;;
        --output_dir) OUTPUT_DIR="$2"; shift ;;
        --filter1) FILTER1="$2"; shift ;;
        --filter2) FILTER2="$2"; shift ;;
        --filter3) FILTER3="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$INPUT_FILE" || -z "$TARGET_FILE" || -z "$OUTPUT_DIR" || -z "$FILTER1" ]]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --input_file <path> --target_file <path> --output_dir <path> --filter1 <filter> [--filter2 <filter>] [--filter3 <filter>]"
    exit 1
fi

Rscript check_enrich.R \
    --input_file "$INPUT_FILE" \
    --target_file "$TARGET_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --filter1 "$FILTER1" \
    ${FILTER2:+--filter2 "$FILTER2"} \
    ${FILTER3:+--filter3 "$FILTER3"}

deactivate


### ------------------------- ### 
########   Runs   ########

# PVT1 1e-9 -- GLMNET
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/pvt1/PVT1_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# PVT1 1e-9 -- cv lm
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/pvt1/PVT1_cv_lm_metrics_tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/pvt1/PVT1_enrich_res/lm \
    --filter1 max_adj_r_sq,0.2 \
    --filter2 mean_cor,0.2
COMMENT


# ERBB4 -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/erbb4/ERBB4_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_uniq_geneId \
    --output_dir /data/chromatin_associated_genes/erbb4/ERBB4_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# FTX -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/ftx/FTX_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/ftx/FTX_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# RAD51B -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/rad51b/RAD51B_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/rad51b/RAD51B_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# SLC30A8 -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/slc30a8/SLC30A8_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/slc30a8/SLC30A8_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# CASC15 -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/casc15/CASC15_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/casc15/CASC15_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# DLEU2 -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/dleu2/DLEU2_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/dleu2/DLEU2_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# TPRG1 -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/tprg1/TPRG1_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/tprg1/TPRG1_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# SMYD3 -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/smyd3/SMYD3_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/smyd3/SMYD3_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT

# PTPRT -- GLMNET #
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/chromatin_associated_genes/ptprt/PTPRT_cv_glmnet_metrics_Tumor.txt \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/chromatin_associated_genes/ptprt/PTPRT_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT


#===============================================================================================================================================================================#


# GLMNET Generalization -- Prostate
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/tissues/prostate/gen_prostateRes/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/tissues/prostate/prostate_enrich_res \
    --filter1 R_squared,0.2
COMMENT

# GLMNET Generalization -- Ovary
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/tissues/ovary/gen_ovaryRes/PVT1_ovary_generalized_glmnet_coeff_metrics.tsv \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/tissues/ovary/ovary_enrich_res \
    --filter1 R_squared,0.2
COMMENT

# GLMNET Generalization -- Uterus
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/tissues/uterus/gen_uterusRes/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/tissues/uterus/uterus_enrich_res \
    --filter1 R_squared,0.2
COMMENT

# GLMNET Generalization -- Adrenal
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/tissues/adrenal/gen_adrenalRes/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/tissues/adrenal/adrenal_enrich_res \
    --filter1 R_squared,0.2
COMMENT

# GLMNET Generalization -- Testis
: << COMMENT
sbatch cluster_enrich.sh \
    --input_file /data/tissues/testis/gen_testisRes/PVT1_generalized_glmnet_coeff_metrics.tsv \
    --target_file /data/mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir /data/tissues/testis/testis_enrich_res \
    --filter1 R_squared,0.2
COMMENT


# SalmonTPM
: << COMMENT
mkdir -p enrich_results && cd enrich_results
sbatch cluster_enrich.sh \
    --input_file PVT1salmonTPM_cv_glmnet_metrics_Tumor.txt \
    --target_file mir200abc141_429_205_mirPath_noNA_uniq \
    --output_dir ./salmonTPM_enrich_res/glmnet \
    --filter1 best_raw_r_sq,0.2 \
    --filter2 cor_value_R2_RMSE,0 \
    --filter3 mean_rmse,1.5
COMMENT
