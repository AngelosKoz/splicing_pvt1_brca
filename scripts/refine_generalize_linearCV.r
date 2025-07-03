set.seed(1234)
library(data.table)
library(dplyr)
library(pls)
library(glm2)
library(glmnet)
library(ROCR)
library(Rcpp)
library(ggplot2)
library(corrplot)
library(argparse)

parser <- ArgumentParser(description = "Generalized Cross-Validation Linear Regression (LM + GLMNET)")
parser$add_argument("-mani", "--manifest_file", required = TRUE, help = "Manifest file containing metadata, along with the sample information (file_name) and a condition column")
parser$add_argument("-coeff", "--coefficients_file", required = FALSE, help = "Path to the original model's linear model coefficients file (for 'coeff' input_mode)")
parser$add_argument("-model_obj", "--model_obj", required = FALSE, help = "Path to the .rds file containing pre-trained models (for 'object' input_mode)")
parser$add_argument("-input_mode", "--input_mode", required = TRUE, choices = c("coeff", "object"), help = "Input mode: 'coeff' for coefficients file or 'object' for pre-trained models")
parser$add_argument("-counts", "--counts_file", required = TRUE, help = "Path to the counts file")
parser$add_argument("-norm", "--norm", required = TRUE, choices = c("yes", "no"), help = "Perform normalization ('yes' or 'no')")
parser$add_argument("-se", "--se_file", required = TRUE, help = "Path to the splicing efficiency file")
parser$add_argument("-umrs", "--umrs_file", required = FALSE, help = "Path to the UMRs file (mandatory if --norm yes is selected)")
parser$add_argument("-out", "--output_dir", required = TRUE, help = "Path to the output directory")
parser$add_argument("-candidate", "--candidate", required = TRUE, help = "Candidate name for splice sites")
parser$add_argument("-sel_gene", "--selected_genes_file", required = FALSE, help = "Optional: Path to the selected genes file")
parser$add_argument("-plots", "--create_plots", required = FALSE, default = "no", choices = c("y", "yes", "no"), help = "Optional: Create plots ('y', 'yes', or 'no') -- (default: no)")
parser$add_argument("-mode", "--mode", required = TRUE, choices = c("lm", "glmnet", "both"), help = "Mode of regression: 'lm' for linear model, 'glmnet' for GLMNET, or 'both'")
parser$add_argument("-ss_table", "--ss_table", required = FALSE, help = "Optional: Path to the splice site table file")
parser$add_argument("-gen_plot", "--gen_plot", required = FALSE, default = "no", choices = c("rmse", "rsq", "no"), help = "Optional: Generate a generalization plot ('rmse', 'rsq', or 'no' -- default: no)")
parser$add_argument("-gen_metrics", "--gen_metrics", required = FALSE, help = "Optional: Path to the metrics file for comparison (required if --gen_plot is selected)")
parser$add_argument("-gen_col", "--gen_col", required = FALSE, help = "Optional: Column name in the metrics file to use for comparison (required if --gen_plot is selected)")

args <- parser$parse_args()
manifest_file <- args$manifest_file
coeff_file <- args$coefficients_file
model_obj <- args$model_obj
input_mode <- args$input_mode
counts_file <- args$counts_file
norm <- args$norm
se_file <- args$se_file
umrs_file <- args$umrs_file
output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
candidate <- args$candidate
selected_genes_file <- args$selected_genes_file
create_plots <- args$create_plots
mode <- args$mode
ss_table <- args$ss_table
gen_plot <- args$gen_plot
gen_metrics <- args$gen_metrics
gen_col <- args$gen_col

if (input_mode == "coeff" && is.null(coeff_file)) {
    stop("Error: --coefficients_file must be provided when input_mode is 'coeff'.")
}
if (input_mode == "object" && is.null(model_obj)) {
    stop("Error: --model_obj must be provided when input_mode is 'object'.")
}
if (!is.null(coeff_file) && !is.null(model_obj)) {
    stop("Error: Only one of --coefficients_file or --model_obj should be provided.")
}
if (gen_plot != "no" && (is.null(gen_metrics) || is.null(gen_col))) {
    stop("Error: --gen_metrics and --gen_col must be provided when --gen_plot is selected.")
}

cat("Arguments parsed successfully.\n")
cat("Input mode:", input_mode, "\n")
cat("Manifest file:", manifest_file, "\n")
cat("Counts file:", counts_file, "\n")
cat("Splicing efficiency file:", se_file, "\n")
cat("Normalization:", norm, "\n")
cat("Output directory:", output_dir, "\n")
cat("Candidate:", candidate, "\n")
cat("Selected genes file:", selected_genes_file, "\n")
cat("Create plots:", create_plots, "\n")
cat("Splice site table:", ss_table, "\n")
cat("Generalization plot:", gen_plot, "\n")
if (input_mode == "coeff") cat("Coefficients file:", coeff_file, "\n")
if (input_mode == "object") cat("Model object:", model_obj, "\n")
if (gen_plot != "no") {
    cat("Metrics file for generalization plot:", gen_metrics, "\n")
    cat("Column for generalization plot:", gen_col, "\n")
}

# Ensure required files exist
if (!file.exists(manifest_file)) stop("Error: Manifest file does not exist.")
if (!file.exists(counts_file)) stop("Error: Counts file does not exist.")
if (!file.exists(se_file)) stop("Error: Splicing efficiency file does not exist.")
if (norm == "yes" && !file.exists(umrs_file)) stop("Error: UMRs file does not exist.")
if (input_mode == "coeff" && !file.exists(coeff_file)) stop("Error: Coefficients file does not exist.")
if (input_mode == "object" && !file.exists(model_obj)) stop("Error: Model object file does not exist.")

cat("Loading input files...\n")
manifest <- fread(manifest_file)
cat("Manifest loaded. Rows:", nrow(manifest), "\n")

counts_in <- fread(counts_file)
cat("Counts file loaded. Rows:", nrow(counts_in), "Columns:", ncol(counts_in), "\n")

se_file <- fread(se_file)
cat("Splicing efficiency file loaded. Rows:", nrow(se_file), "Columns:", ncol(se_file), "\n")

# RPM + loge() counts
if (norm == "yes") {
    umrs <- fread(umrs_file)
    cat("UMRs file loaded. Rows:", nrow(umrs), "Columns:", ncol(umrs), "\n")
    cat("Performing normalization on counts...\n")
    counts <- as.data.frame(counts_in[,-1])  # Removing the first column (Geneid)
    rownames(counts) <- counts_in$Geneid
    counts_nonzero <- counts[rowSums(counts) != 0,]  # Remove rows with all zero counts
    counts_filtered <- counts_nonzero[(rowSums(counts_nonzero >= 10) >= 10), ]  # At least 10 counts in at least 10 samples
    counts_transposed <- as.data.frame(t(counts_filtered))
    counts_transposed$file_name <- rownames(counts_transposed)
    counts_merged <- merge(counts_transposed, umrs, by="file_name")
    gene_cols <- grep("^ENS", names(counts_merged), value = TRUE)
    counts_rpm <- counts_merged %>%
        mutate(across(all_of(gene_cols), ~ . * RPM_NormFactor)) %>%
        mutate(across(all_of(gene_cols), ~ pmin(.x, quantile(.x, 0.99, na.rm = TRUE)))) %>%
        mutate(across(all_of(gene_cols), ~ log(.x + 1e-9)))
    fwrite(counts_rpm, file.path(output_dir, paste0(candidate, "_RPMlog_counts.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Normalization completed. Normalized counts saved to:", file.path(output_dir, paste0(candidate, "_RPMlog_counts.tsv")), "\n")
} else {
    cat("Normalization skipped. Counts considered normalized.\n")
    counts_rpm <- counts_in
}

cat("Processing splice site file...\n")
colnames(se_file)[1] <- "file_name"
se_df <- se_file %>%
    rename(file_name = 1) %>%
    select(file_name, matches(":"))
splice_sites <- grep(":", names(se_df), value = TRUE)
nsplicesites <- length(splice_sites)
geneSS <- merge(counts_rpm, se_df, by = "file_name")
cat("Splice site file processed. Number of splice sites:", nsplicesites, "\n")
cat("Genes RPM Shape:", dim(counts_rpm), "\n")
cat("Splice sites Shape:", dim(se_df), "\n")
cat("geneSS Shape:", dim(geneSS), "\n")
cat("Merging manifest with splice site data...\n")
samp_df <- se_df[, .(file_name)]
manifest <- merge(manifest, samp_df, by = "file_name", all.x = FALSE)
if (nrow(manifest) == 0) {
    stop("Error: No matching rows found between manifest and splice site data.")
}
cat("Manifest merged successfully. Rows retained:", nrow(manifest), "\n")

# Filter genes if selected_genes_file is provided
if (!is.null(selected_genes_file)) {
    cat("Filtering genes using selected genes file...\n")
    selected_genes <- fread(selected_genes_file, header = FALSE)$V1
    subset_geneSS <- geneSS[, c("file_name", intersect(names(geneSS), selected_genes), splice_sites), with = FALSE] # removed ,with = FALSE
} else {
    cat("No gene filtering applied. Using all genes.\n")
    subset_geneSS <- geneSS  # Use all genes if no file is provided
}

subset_geneSS$Intercept <- 1
metrics_list <- list()

print(subset_geneSS[1:5,1:3])
print(dim(subset_geneSS))

base_dir <- file.path(output_dir, "generalized_results")
if (!dir.exists(base_dir)) dir.create(base_dir)

if (!is.null(ss_table)) {
    cat("Loading splice site table...\n")
    ss_mapping <- fread(ss_table, header = FALSE)
    colnames(ss_mapping) <- c("splice_site", "label")
    cat("Splice site table loaded. Rows:", nrow(ss_mapping), "\n")
} else {
    ss_mapping <- NULL
}

if (input_mode == "coeff") {
    cat("Running in 'coeff' input_mode...\n")
    cat("Loading coefficients file...\n")
    coef_data <- fread(coeff_file)
    coef_gene <- unique(coef_data$Gene)
    coef_ss <- unique(coef_data$Term)
    counts_rpm <- counts_rpm %>%
        relocate(file_name) %>% 
        select(file_name, matches("^ENSG") & any_of(coef_gene))
    ##
    ##
    for (gene in unique(coef_data$Gene)) {
        if (!gene %in% names(subset_geneSS)) {
            warning(paste("Gene column not found in subset_geneSS. Skipping:", gene))
            next
        }
        cat("Processing gene:", gene, "\n")
        
        gene_coefs <- coef_data[Gene == gene]
        predictors <- gene_coefs$Term[gene_coefs$Term != "Intercept"]
        formula_str <- paste("`", gene, "` ~ ", paste(predictors, collapse = " + "), sep = "")
        formula <- as.formula(formula_str)
        cat("Formula for ", gene, "\t", formula_str, "\n")
        ##
        gene_data <- subset_geneSS[, c("file_name", "Intercept", predictors, gene), with = FALSE] # removed ,with = FALSE
        fitted_values <- as.matrix(gene_data[, c("Intercept", predictors), with = FALSE]) %*% gene_coefs$Estimate
        observed_values <- gene_data[[gene]]
        ##
        residuals <- observed_values - fitted_values
        mse <- mean(residuals^2, na.rm = TRUE)
        rmse <- sqrt(mse)
        r_squared <- 1 - (sum(residuals^2, na.rm = TRUE) / sum((observed_values - mean(observed_values, na.rm = TRUE))^2, na.rm = TRUE))
        metrics_list[[gene]] <- data.frame(Gene = gene, RMSE = rmse, R_squared = r_squared)
        ##
        ##
        if (create_plots %in% c("y", "yes")) {
            cat("Creating plots for gene:", gene, "\n")

            plot_dir <- file.path(base_dir, "plots")
            if (!dir.exists(plot_dir)) dir.create(plot_dir)
            if (!dir.exists(file.path(plot_dir, "coeff"))) dir.create(file.path(plot_dir, "coeff"))
            if (!dir.exists(file.path(plot_dir, "fit"))) dir.create(file.path(plot_dir, "fit"))
            ##
            ##
            # Plot Observed vs. Predicted
            plot_data <- data.frame(Observed = observed_values, Predicted = as.numeric(fitted_values))
            cor_val <- cor.test(plot_data$Observed, plot_data$Predicted)$estimate
            plot <- ggplot(plot_data, aes(x = Observed, y = Predicted)) +
                geom_point(size = 0.5) +
                geom_smooth(method = 'lm', color = "grey30", lwd = 1.1) +
                theme_bw() +
                labs(title = paste(candidate, gene, "~ Generalized Linear Model (All Observed)"),
                     subtitle = paste("Pearson's correlation:", round(cor_val, 3), "| R-squared:", round(r_squared, 3)),
                     x = paste(candidate, "Observed Values (loge(RPM))"),
                     y = "Predicted Values")
            ggsave(filename = file.path(plot_dir, "fit", paste0(candidate, "_", gene, "_", mode, "_observed_vs_predicted_all.pdf")),
                plot = plot, width = 6, height = 5, device = "pdf")
        }
    }
} else if (input_mode == "object") {
    cat("Running in 'object' input_mode...\n")
    # Load pre-trained model object
    cat("Loading pre-trained model object from:", model_obj, "\n")
    model_list <- readRDS(model_obj)
    cat("Model object loaded successfully. Number of models:", length(model_list), "\n")
    print(names(model_list))

    for (gene in names(model_list)) {
        cat("Processing gene:", gene, "\n")
        model <- model_list[[gene]]

        cat("Structure of model for ", gene, "\n")
        print(str(model))
        ##
        ## GLMNET models
        ##
        if ("glmnet" %in% class(model$model)) {
            cat("Detected glmnet model for gene:", gene, "\n")
            lambda_min <- model$lambda_min 
            model <- model$model
        } else if ("lm" %in% class(model)) {
            cat("Detected lm model for gene:", gene, "\n")
        } else {
            warning("Unsupported model type for gene:", gene, ". Skipping...")
            next
        }
        cat("Class of extracted model for gene:", gene, "is:", class(model), "\n")
        if (!gene %in% names(subset_geneSS)) {
            warning("Gene not found in subset_geneSS. Skipping:", gene)
            next
        }
        if (inherits(model, "glmnet")) {
            cat("Processing glmnet model for gene:", gene, "\n")
            predictors <- rownames(coef(model, s = lambda_min))[-1]  # Exclude intercept
            if (length(predictors) == 0) {
                warning("No valid predictors for gene:", gene, ". Skipping.")
                next
            }
            gene_data <- subset_geneSS[, c("file_name", predictors, gene), with = FALSE]
            x <- as.matrix(gene_data[, predictors, with = FALSE])
            y <- gene_data[[gene]]
            if (var(y, na.rm = TRUE) == 0) {
                warning("Skipping gene:", gene, "- constant observed values.")
                next
            }
            cat("Attempting to predict for gene:", gene, "\n")
            set.seed(1234)
            pred_vals <- predict(model, newx = x, s = lambda_min)
            if (is.null(pred_vals)) {
                warning("Prediction failed for gene:", gene, ". Skipping.")
                next
            }
        ##
        ## Linear models lm()
        ##
        } else if (inherits(model, "lm")) {
            cat("Processing lm model for gene:", gene, "\n")
            predictors <- names(coef(model))[-1]
            if (length(predictors) == 0) {
                warning("No valid predictors for gene:", gene, ". Skipping.")
                next
            }
            gene_data <- subset_geneSS[, c("file_name", predictors, gene), with = FALSE]
            formula <- as.formula(paste(gene, "~", paste(predictors, collapse = " + ")))
            set.seed(1234)
            pred_vals <- predict(model, newdata = gene_data)
            y <- gene_data[[gene]]
        }
        residuals <- y - pred_vals
        mse <- mean(residuals^2, na.rm = TRUE)
        rmse <- sqrt(mse)
        r_squared <- 1 - (sum(residuals^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE))
        metrics_list[[gene]] <- data.frame(Gene = gene, RMSE = rmse, R_squared = r_squared)
        cat("Metrics calculated for gene:", gene, "\n")
    }
}


if (length(metrics_list) == 0) {
    stop("Error: metrics_list is empty. No valid metrics were calculated.")
}


metrics_list <- metrics_list[!sapply(metrics_list, is.null)]
if (length(metrics_list) == 0) {
    stop("Error: All elements in metrics_list are NULL or invalid.")
}


metrics_df <- do.call(rbind, metrics_list)
if (is.null(metrics_df) || nrow(metrics_df) == 0) {
    stop("Error: metrics_df is empty or invalid. No valid metrics were calculated.")
}


if (!"R_squared" %in% colnames(metrics_df)) {
    stop("Error: R_squared column is missing in metrics_df.")
}
if (!is.numeric(metrics_df$R_squared)) {
    stop("Error: R_squared column in metrics_df is not numeric.")
}


metrics_df <- metrics_df[order(-metrics_df$R_squared), ]
output_file <- file.path(base_dir, paste0(candidate, "_generalized_", mode, "_", input_mode, "_metrics.tsv"))
fwrite(metrics_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Metrics saved to:", output_file, "\n")



if (gen_plot != "no") {
    cat("Generating generalization plot...\n")
    

    if (!file.exists(gen_metrics)) {
        stop("Error: The specified --gen_metrics file does not exist.")
    }
    original_metrics <- fread(gen_metrics)
    
    if (!(gen_col %in% colnames(original_metrics))) {
        stop(paste("Error: Column", gen_col, "not found in the provided --gen_metrics file."))
    }
    

    original_metrics <- original_metrics[, .(Gene, Original_metric = get(gen_col))]
    merged_metrics <- merge(metrics_df, original_metrics, by = "Gene", all = FALSE)
    if (nrow(merged_metrics) == 0) {
        stop("Error: No common genes found between the generated metrics and the external metrics.")
    }


    y_metric <- if (gen_plot == "rmse") "RMSE" else "R_squared"
    x_metric <- "Original_metric"
    ##
    ## Generalization plot
    ##
    plot <- ggplot(merged_metrics, aes(x = .data[[x_metric]], y = .data[[y_metric]])) +
        geom_point(size = 0.6, alpha = 0.4) +
        geom_smooth(method = "lm", color = "red", linetype = "dashed") +
        geom_abline(slope = 1, intercept = 0, linetype = "dotdash", color = "blue") +  # identity line
        theme_minimal() +
        labs(
            title = paste("Generalization Plot:", toupper(gen_plot)),
            x = paste("External Metric (", gen_col, ")", sep = ""),
            y = paste("Generated Metric (", y_metric, ")", sep = "")
        )
    plot_dir <- file.path(base_dir, "plots", "generalization")
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    plot_file <- file.path(plot_dir, paste0(candidate, "_", mode, "_", input_mode, "_generalization_", gen_plot, ".pdf"))
    ggsave(plot_file, plot, width = 12, height = 10, device = "pdf")
    cat("Generalization plot saved to:", plot_file, "\n")
}
