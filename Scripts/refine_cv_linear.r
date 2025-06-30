set.seed(1234)

library(data.table)
library(glmnet)
library(ggplot2)
library(argparse)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggrepel)
library(ggplotify)

parser <- ArgumentParser(description = "Cross-Validation Linear Regression (LM + GLMNET)")
parser$add_argument("-mani", "--manifest_file", required = TRUE, help = "Manifest file containing metadata, along with the sample information (file_name) and a condition column")
parser$add_argument("-counts", "--counts_file", required = TRUE, help = "Path to the counts file")
parser$add_argument("-norm", "--norm", required = TRUE, choices = c("yes", "no"), help = "Perform normalization ('yes' or 'no')")
parser$add_argument("-se", "--se_file", required = TRUE, help = "Path to the splicing efficiency file")
parser$add_argument("-umrs", "--umrs_file", required = FALSE, help = "Path to the UMRs file (mandatory if --norm yes is selected)")
parser$add_argument("-out", "--output_dir", required = TRUE, help = "Path to the output directory")
parser$add_argument("-candidate", "--candidate", required = TRUE, help = "Candidate name for splice sites (e.g., 'PVT1')")
parser$add_argument("-sel_gene", "--selected_genes_file", required = FALSE, help = "Optional: Path to the selected genes file")
parser$add_argument("-group", "--group_column", required = FALSE, help = " Optional: Column in the manifest used to group the samples. If not used, performs 'agnostic' run")
parser$add_argument("-plots", "--create_plots", required = FALSE, default = "off", choices = c("y", "yes", "off"), help = "Optional: Create plots ('y', 'yes', or 'off') -- (default: off)")
parser$add_argument("-pmode", "--plot_mode", required = FALSE, default = "best", choices = c("best", "all"), help = "Optional: For plotting observed vs predicted. Best will use best folds samples, while all will use all the observed. This is purely for plotting purposes. (default: best)")
parser$add_argument("-brk", "--breaks", required = FALSE, default = 500, type = "integer", help = "<-- Optional: The amount of samples to be considered as minimum for a 10-fold cross validation. If that value is unmet, 5-fold cross validation is performed. (default: 500)")
parser$add_argument("-mode", "--mode", required = TRUE, choices = c("lm", "glmnet", "both"), help = "Mode of regression: 'lm' for linear model, 'glmnet' for GLMNET, or 'both'")
parser$add_argument("-group_values", "--group_values", required = FALSE, help = "Optional: Specific values in the group_column to filter (comma-separated). If not provided, all unique values in group_column will be used.")

args <- parser$parse_args()
manifest_file <- args$manifest_file
counts_file <- args$counts_file
norm <- args$norm
se_file <- args$se_file
umrs_file <- args$umrs_file
output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
candidate <- args$candidate
selected_genes_file <- args$selected_genes_file
group_column <- args$group_column
create_plots <- args$create_plots
plot_mode <- args$plot_mode
breaks <- args$breaks
mode <- args$mode
group_values <- args$group_values

if (norm == "yes" && is.null(umrs_file)) {
    stop("Error: --umrs_file is mandatory when --norm yes is selected.")
}

cat("Arguments parsed successfully.\n")
cat("Manifest file:", manifest_file, "\n")
cat("Counts file:", counts_file, "\n")
cat("Normalization:", norm, "\n")
cat("Output directory:", output_dir, "\n")
cat("Breaks for cross-validation:", breaks, "\n")

manifest <- fread(manifest_file)
counts_in <- fread(counts_file)
se_file <- fread(se_file)
cat("Input files loaded successfully.\n")



# RPM + loge() counts
if (norm == "yes") {
    umrs <- fread(umrs_file)
    cat("Uniquely mapped reads succesfully loaded:", nrow(umrs), "\n")
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
    fwrite(counts_rpm, file.path(output_dir, paste0(candidate, "_RPMlog_pseudo_counts1e-4.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Normalization completed. Normalized counts saved to:", file.path(output_dir, paste0(candidate, "_RPMlog_pseudo_counts.tsv")), "\n")
} else {
    cat("Normalization skipped. Counts considered normalized.\n")
    counts_rpm <- as.data.frame(counts_in)
    counts_rpm$file_name <- as.character(counts_rpm$file_name)
    counts_rpm <- counts_rpm %>%
        relocate(file_name) %>%
        select(file_name, where(~ is.numeric(.x) && sum(.x, na.rm = TRUE) != 0))
    cat("Counts file shape:", dim(counts_rpm), "\n")
}


cat("Processing splice site file...\n")
colnames(se_file)[1] <- "file_name"
se_df <- se_file %>% select(file_name, matches(":"))
splice_sites <- grep(":", names(se_df), value = TRUE)
nsplicesites <- length(splice_sites)
cat("Splice site DF Shape:", dim(se_df), "\n")
print(se_df[1:5, 1:5]) # debug
cat("Splice sites before merging:\n")
print(splice_sites)





geneSS <- merge(counts_rpm, se_df, by = "file_name")
cat("Splice site file processed. Number of splice sites:", nsplicesites, "\n")
cat("geneSS Shape after merging:", dim(geneSS), "\n")
cat("Splice sites after merging:\n")
print(intersect(names(geneSS), splice_sites))

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
    subset_geneSS <- geneSS[, c("file_name", intersect(names(geneSS), selected_genes), splice_sites)]  # Removed 'with = FALSE'
    cat("Gene filtering completed. Number of genes after filtering:", ncol(subset_geneSS) - 1 - nsplicesites, "\n")
} else {
    cat("No gene filtering applied. Using all genes.\n")
    subset_geneSS <- geneSS  # Use all genes if no file is provided
}

# Process based on group_column value
if (!is.null(group_column)) {
    cat("Grouping data by column:", group_column, "\n")
    unique_groups <- unique(manifest[[group_column]])
    
    # Filter specific group values if provided
    if (!is.null(group_values)) {
        group_values <- unlist(strsplit(group_values, ","))
        unique_groups <- unique_groups[unique_groups %in% group_values]
        if (length(unique_groups) == 0) {
            stop("Error: No matching group values found in the manifest.")
        }
    }
    
    subset_list <- list()
    for (group in unique_groups) {
        group_id <- manifest[get(group_column) == group, .(file_name)]
        group_data <- subset_geneSS[subset_geneSS$file_name %in% group_id$file_name,]
        if (nrow(group_data) > 0) {
            subset_list[[group]] <- group_data
        } else {
            cat("Skipping group:", group, "- no data available.\n")
        }
    }
    if (length(subset_list) == 0) {
        stop("Error: No non-empty groups found. Please check your data.")
    }
    cat("Final subset list names and values:\n")
    for (name in names(subset_list)) {
        cat("Subset name:", name, "- Number of rows:", nrow(subset_list[[name]]), "\n")
        cat("Subset shape:", dim(subset_list[[name]]), "\n")  # Add this line to print the shape of each subset
    }
    cat("Data grouped into", length(subset_list), "non-empty groups.\n")
} else {
    cat("No grouping applied. Using full dataset.\n")
    if (nrow(subset_geneSS) > 0) {
        subset_list <- list(agnostic = subset_geneSS)  # full dataset if no group_column is provided
    } else {
        stop("Error: Full dataset is empty. Please check your data.")
    }
}

cat("subset_geneSS Shape:", dim(subset_geneSS), "\n")
print(subset_geneSS[1:5, 1:5])  # Print first 5 rows and columns for debugging


##
##
predictor_columns <- grep(":", names(subset_geneSS), value = TRUE)  # For Splice Sites
predictor_columns_bqt <- sapply(predictor_columns, function(x) paste0("`", x, "`"))
response_column <- grep("ENS", names(subset_geneSS), value = TRUE)
##
##
## Model loops -- Linear Model and GLMNET
##
##
##### GLMNET Regression with Cross Validation #####
if (mode == "glmnet" || mode == "both") {
    cat("GLMNET regression cross-validation selected.\n")
    cat("Starting GLMNET regression cross-validation...\n")
    base_dir <- file.path(output_dir, "cv_glmnet_results")
    glmnet_models <- list()  # Initialize a list to store all glmnet models
    for (sub_Name in names(subset_list)) {
        cat("Processing subset:", sub_Name, "\n")
        sub <- subset_list[[sub_Name]]
        if (!dir.exists(base_dir)) dir.create(base_dir)
        if (!dir.exists(file.path(base_dir, sub_Name))) dir.create(file.path(base_dir, sub_Name))
        if (create_plots %in% c("y", "yes")) {
            if (!dir.exists(file.path(base_dir, sub_Name, "plots"))) dir.create(file.path(base_dir, sub_Name, "plots"))
            if (!dir.exists(file.path(base_dir, sub_Name, "plots", "coeff"))) dir.create(file.path(base_dir, sub_Name, "plots", "coeff"))
            if (!dir.exists(file.path(base_dir, sub_Name, "plots", "fit"))) dir.create(file.path(base_dir, sub_Name, "plots", "fit"))
        }
        ##
        combined_rsq <- data.frame()
        glm_coeffs_list <- list()
        for (gene in response_column) {
            # Check for constant y across all samples before any folds
            y_all <- sub[[gene]]
            if (length(unique(y_all)) <= 1) {
                cat("Skipping gene:", gene, "- constant response across all samples.\n")
                next
            }
            raw_r_sq <- c()
            mean_rmse_error <- c()
            set.seed(1234)
            data <- sub[sample(nrow(sub)),]
            brk <- if (nrow(sub) > breaks) 10 else 5
            set.seed(1234)
            folds <- cut(seq(1, nrow(data)), breaks = brk, labels = FALSE)
            cat("Processing subset:", sub_Name, ", Gene:", gene, "\n")
            formula_str <- paste(gene, "~", paste(predictor_columns_bqt, collapse = "+"))
            cat("Formula for glmnet:", formula_str, "\n")
            for (i in 1:brk) {
                testIndexes <- which(folds==i, arr.ind=TRUE)
                testData <- data[testIndexes, ]
                trainData <- data[-testIndexes, ]
                y_train <- trainData[[gene]]
                y_test <- testData[[gene]]
                ##
                uniq_train <- unique(y_train)
                tab_train <- table(y_train)
                fail_train <- (length(uniq_train) <= 1) ||
                              (length(uniq_train) == 2 && any(tab_train < 2))
                uniq_test <- unique(y_test)
                tab_test <- table(y_test)
                fail_test <- (length(uniq_test) <= 1) ||
                             (length(uniq_test) == 2 && any(tab_test < 2))
                if (fail_train || fail_test) {
                    cat("Skipping fold", i, "for gene", gene, ": not enough replicates for each unique y value in train or test fold\n")
                    cat("  y_train unique:", uniq_train, " counts:", as.vector(tab_train), "\n")
                    cat("  y_test unique:", uniq_test, " counts:", as.vector(tab_test), "\n")
                    mean_rmse_error[i] <- NA_real_
                    raw_r_sq[i] <- NA_real_
                    next
                }
                x_train <- as.matrix(trainData[, predictor_columns, drop = FALSE])
                x_test <- as.matrix(testData[, predictor_columns, drop = FALSE])
                set.seed(1234)
                cv_fit <- tryCatch(
                    cv.glmnet(x_train, y_train, alpha=0.5, type.measure="mse", family="gaussian"),
                    error = function(e) {
                        cat("cv.glmnet error for gene", gene, "in fold", i, ":", conditionMessage(e), "\n")
                        return(NULL)
                    }
                )
                if (is.null(cv_fit)) {
                    mean_rmse_error[i] <- NA_real_
                    raw_r_sq[i] <- NA_real_
                    next
                }
                lam <- cv_fit$lambda.min
                glm_fit <- glmnet(x_train ,y_train ,alpha=0.5, lambda=lam)
                pred_vals <- predict(glm_fit, newx=x_test, type='response', s=lam)
                mean_rmse_error[i] <- sqrt(mean((y_test - pred_vals)^2))
                SSE <- sum((y_test - pred_vals)^2)
                SST <- sum((y_test - mean(y_test))^2)
                raw_r_sq[i] <- 1 - (SSE / SST)
            }
            # Add check: skip gene if all folds were skipped (all mean_rmse_error are NA)
            if (all(is.na(mean_rmse_error))) {
                cat("Skipping gene:", gene, "- all folds skipped due to constant y in train or test fold.\n")
                next
            }
            best_fold <- which.min(mean_rmse_error)
            mean_raw_r_sq <- mean(raw_r_sq)
            if (length(raw_r_sq[best_fold]) == 1){
                best_raw_r_sq <- raw_r_sq[best_fold]
            } else {
                best_raw_r_sq <- NA_real_
            }
            mean_rmse <- mean(mean_rmse_error)
            min_rmse <- min(mean_rmse_error)
            v_raw_r_sq <- unlist(raw_r_sq)
            v_rmse <- unlist(mean_rmse_error)
            v_index <- complete.cases(v_raw_r_sq, v_rmse)
            filt_raw_r_sq <- v_raw_r_sq[v_index]
            filt_rmse <- v_rmse[v_index]
            if (length(unique(filt_raw_r_sq)) > 2 && length(unique(filt_rmse)) > 2) {
                test_cor <- cor.test(filt_raw_r_sq, filt_rmse) # Expect a negative correlation (lower RMSE ~ higher R²) 
            } else {
                test_cor <- list(estimate = NA_real_, p.value = NA_real_)
                print("Skipping correlation test: insufficient observations")
            }
            combined_rsq <- rbind(combined_rsq, data.frame( Gene = gene, mean_raw_r_sq = mean_raw_r_sq, best_raw_r_sq = best_raw_r_sq, mean_rmse = mean_rmse, min_rmse = min_rmse, 
                                                            cor_value_R2_RMSE = test_cor$estimate, cor_pvalue = test_cor$p.value, row.names = NULL))  # Extra addittion to solve issue. Remove if problematic
            ##
            ##
            ## refit the model using best fold
            testIndexes <- which(folds == best_fold, arr.ind=TRUE)
            testData <- data[testIndexes, ]
            trainData <- data[-testIndexes, ]
            y_train <- trainData[[gene]]
            x_train <- as.matrix(trainData[, predictor_columns, drop = FALSE])
            y_test <- testData[[gene]]
            x_test <- as.matrix(testData[, predictor_columns, drop = FALSE]) 
            y_full <- data[[gene]]
            x_full <- as.matrix(data[, predictor_columns, drop = FALSE])
            if (length(unique(y_train)) > 1 && length(unique(y_test)) > 1) {
                set.seed(1234)
                cv_fit_best <- tryCatch(
                    cv.glmnet(x_train, y_train, alpha = 0.5, type.measure = "mse", family = "gaussian"),
                    error = function(e) {
                        cat("cv.glmnet error (best fold) for gene", gene, ":", conditionMessage(e), "\n")
                        return(NULL)
                    }
                )
                if (is.null(cv_fit_best)) {
                    glmnet_models[[gene]] <- NULL
                    glm_coeffs_list[[gene]] <- data.frame(Gene = gene, Term = NA_real_, GLMNET_coefficients = NA_real_, lambda_min = NA_real_)
                    next
                }
                lam_best <- cv_fit_best$lambda.min
                glm_fit_best <- glmnet(x_train, y_train, alpha = 0.5, lambda = lam_best)
                pred_vals_best <- predict(glm_fit_best, newx = x_test, s = lam_best)
                pred_vals_full <- predict(glm_fit_best, newx = x_full, s = lam_best)
                coefs_mat <- as.matrix(coef(glm_fit_best, s = lam_best))
                coefs_df <- data.frame(Term = rownames(coefs_mat), Estimate = coefs_mat[, 1])
                glm_coeffs_list[[gene]] <- data.frame(Gene = gene, Term = rownames(coefs_mat), GLMNET_coefficients = coefs_mat[, 1], lambda_min = lam_best)
                glmnet_models[[gene]] <- list(model = glm_fit_best, lambda_min = lam_best)
            } else {
                glmnet_models[[gene]] <- NULL
                glm_coeffs_list[[gene]] <- data.frame(Gene = gene, Term = NA_real_, GLMNET_coefficients = NA_real_, lambda_min = NA_real_)
            }
            ##
            ##
            ## Plot Observed vs. Predicted
            ##
            ##
            if (create_plots %in% c("y", "yes")) {
                pred_vals <- if (plot_mode == "best") pred_vals_best else pred_vals_full
                observed_vals <- if (plot_mode == "best") y_test else y_full

                if (length(pred_vals) == length(observed_vals)) {
                    z <- data.frame(Observed = observed_vals, Predicted = as.numeric(pred_vals))
                    plot2 <- ggplot(z, aes(x = Observed, y = Predicted)) +
                        geom_point(size = 0.5) +
                        geom_smooth(method = 'lm', color = "grey30", lwd = 1.1) +
                        theme_bw() +
                        labs(
                            title = paste(candidate, gene, "~ Elastic Net (", if (plot_mode == "best") "Best Fold" else "All Observed", ")"),
                            subtitle = paste("RMSE (Best Fold):", round(min_rmse, 3), "| R² (Best Fold):", round(best_raw_r_sq, 3)),
                            x = paste(candidate, "Observed Values (loge(RPM))"),
                            y = "Predicted Values")
                    ggsave(filename = file.path(base_dir, sub_Name, "plots", "fit", paste0(candidate, "_", gene, "_", sub_Name, "_glmnet_Fit_", if (plot_mode == "best") "bestFold" else "fullDF", ".pdf")),
                           plot = plot2, device = "pdf", width = 6, height = 5)
                } else {
                    warning("Length mismatch between predicted and observed values. Skipping observed vs. predicted plot.")
                }
            }
        }
        ##
        ## GLMNET coefficients
        ##
        glm_coeffs_df <- do.call(rbind, glm_coeffs_list)
        write.table(glm_coeffs_df, file.path(base_dir, sub_Name, paste0(candidate, "_cv_glmnet_coefficients_", sub_Name, ".txt")), sep = "\t", row.names = FALSE, quote = FALSE)
        ##
        combined_rsq$cor_padj_BH <- p.adjust(combined_rsq$cor_pvalue, method = "BH")
        combined_rsq <- combined_rsq[order(combined_rsq$min_rmse), ]
        write.table(combined_rsq, file.path(base_dir, sub_Name, paste0(candidate, "_cv_glmnet_metrics_", sub_Name, ".txt")), sep = "\t", row.names = FALSE, quote = FALSE)
        ##
        saveRDS(glmnet_models, file = file.path(base_dir, sub_Name, paste0(candidate, "_cv_glmnet_models_", sub_Name, ".rds")))
        cat("All GLMNET models saved to:", file.path(base_dir, sub_Name, paste0(candidate, "_cv_glmnet_models_", sub_Name, ".rds")), "\n")
        cat("Subset processing completed:", sub_Name, "\n")
    }
    
    cat("GLMNET regression cross-validation completed.\n")
}


##
##
# Linear Model with Cross Validation
if (mode == "lm" || mode == "both") {
    cat("Linear model cross-validation selected.\n")
    cat("Starting linear model cross-validation...\n")
    base_dir <- file.path(output_dir, "cv_lm_results")
    lm_models <- list()
    for (sub_Name in names(subset_list)) {
        cat("Processing subset:", sub_Name, "\n")
        sub <- subset_list[[sub_Name]]
        if (!dir.exists(base_dir)) dir.create(base_dir)
        if (!dir.exists(file.path(base_dir, sub_Name))) dir.create(file.path(base_dir, sub_Name))
        if (create_plots %in% c("y", "yes")) {
            if (!dir.exists(file.path(base_dir, sub_Name, "plots"))) dir.create(file.path(base_dir, sub_Name, "plots"))
            if (!dir.exists(file.path(base_dir, sub_Name, "plots", "coeff"))) dir.create(file.path(base_dir, sub_Name, "plots", "coeff"))
            if (!dir.exists(file.path(base_dir, sub_Name, "plots", "fit"))) dir.create(file.path(base_dir, sub_Name, "plots", "fit"))
        }
        combined_rsq <- data.frame()
        lm_coeffs_list <- list()
        for (gene in response_column) {
            cat("Processing gene:", gene, "in subset:", sub_Name, "\n")
            corel <- c()
            adj_r_sq <- c()
            set.seed(1234)
            data <- sub[sample(nrow(sub)),]
            formula_str <- paste(gene, "~", paste(predictor_columns_bqt, collapse="+"))
            formula <- as.formula(formula_str)
            ##
            ## Create folds
            brk <- if (dim(sub)[1] > breaks) 10 else 5
            set.seed(1234)
            folds <- cut(seq(1, nrow(data)), breaks = brk, labels = FALSE)
            cat("Processing subset:", sub_Name, ", Gene:", gene,"\n")
            ##
            df_coeff <- data.frame(matrix(ncol = length(predictor_columns_bqt), nrow = 0)) # added for boxplot
            colnames(df_coeff) <- predictor_columns_bqt # added for boxplot
            for (i in 1:brk) {
                testIndexes <- which(folds==i, arr.ind=TRUE)
                testData <- data[testIndexes, ]
                trainData <- data[-testIndexes, ]
                # Check for enough non-NA cases in trainData and testData
                y_train <- trainData[[gene]]
                y_test <- testData[[gene]]
                uniq_train <- unique(y_train)
                tab_train <- table(y_train)
                fail_train <- (length(uniq_train) <= 1) ||
                              (length(uniq_train) == 2 && any(tab_train < 2))
                uniq_test <- unique(y_test)
                tab_test <- table(y_test)
                fail_test <- (length(uniq_test) <= 1) ||
                             (length(uniq_test) == 2 && any(tab_test < 2))
                if (fail_train || fail_test) {
                    cat("Skipping fold", i, "for gene", gene, ": not enough replicates for each unique y value in train or test fold\n")
                    cat("  y_train unique:", uniq_train, " counts:", as.vector(tab_train), "\n")
                    cat("  y_test unique:", uniq_test, " counts:", as.vector(tab_test), "\n")
                    corel[i] <- NA_real_
                    adj_r_sq[i] <- NA_real_
                    next
                }
                ##
                ##
                lm_fit <- tryCatch(
                    lm(formula, data=trainData),
                    error = function(e) {
                        cat("lm error for gene", gene, "in fold", i, ":", conditionMessage(e), "\n")
                        return(NULL)
                    }
                )
                if (is.null(lm_fit)) {
                    corel[i] <- NA_real_
                    adj_r_sq[i] <- NA_real_
                    next
                }
                c <- coef(lm_fit)[-1]  # added for boxplot, exclude intercept
                df_coeff <- rbind(df_coeff, c) # added for boxplot
                pred_vals <- predict(lm_fit, testData, type='response')
                if (sum(is.finite(pred_vals) & is.finite(testData[[gene]])) >= 2) {
                    corel[i] <- cor.test(pred_vals, testData[[gene]])$estimate
                } else {
                    warning(paste("Not enough valid values for correlation.", gene))
                    corel[i] <- NA_real_
                }
                adj_r_sq[i] <- summary(lm_fit)$adj.r.squared
            }
            best_fold <- which.max(adj_r_sq)
            mean_adj_r_sq <- mean(adj_r_sq)
            max_adj_r_sq <- max(adj_r_sq)
            mean_cor <- mean(unlist(corel))
            ## check if enough values to perform the correlation test
            v_adj_r_sq <- unlist(adj_r_sq)
            v_corel <- unlist(corel)
            v_index <- complete.cases(v_adj_r_sq, v_corel)
            filt_adj_r_sq <- v_adj_r_sq[v_index]
            filt_corel <- v_corel[v_index]
            if (length(filt_corel) > 2 && length(unique(filt_corel)) > 2) {
                test_cor <- cor.test(filt_adj_r_sq, filt_corel)
            } else {
                test_cor <- list(estimate = NA, p.value = NA)
                print("Skipping correlation test: insufficient observations")
            }
            combined_rsq <- rbind(combined_rsq, data.frame(Gene = gene, mean_adj_r_sq = mean_adj_r_sq, max_adj_r_sq = max_adj_r_sq,
                                                           mean_cor = mean_cor, cor_value = test_cor$estimate, cor_pvalue = test_cor$p.value))
            ##
            ##
            ## refit the model using best fold
            testIndexes <- which(folds == best_fold, arr.ind=TRUE)
            testData <- data[testIndexes, ]
            trainData <- data[-testIndexes, ]
            y_train <- trainData[[gene]]
            # Add check for best_fold: if nothing to fit, fill with NA_real_
            if (length(y_train) == 0 || all(is.na(y_train)) || length(unique(na.omit(y_train))) < 2) {
                lm_models[[gene]] <- NULL
                coefs <- data.frame(lm_Predictor = NA_character_, Estimate = NA_real_, StdError = NA_real_,
                    tvalue = NA_real_, PrGTt = NA_real_)
                sig_stars <- NA_character_
            } else {
                lm_fit <- lm(formula, data=trainData)
                lm_models[[gene]] <- lm_fit  # Store the model
                coefs <- data.frame(Term = rownames(summary(lm_fit)$coefficients), summary(lm_fit)$coefficients)
                colnames(coefs) <- c("lm_Predictor", "Estimate", "StdError", "tvalue", "PrGTt")
                sig_stars <- ifelse(coefs$PrGTt < 0.001, "***", ifelse(coefs$PrGTt < 0.01, "**", ifelse(coefs$PrGTt < 0.05, "*", "")))
            }
            lm_coeffs_list[[gene]] <- data.frame(Gene = gene, Term = coefs$lm_Predictor, Estimate = coefs$Estimate, StdError = coefs$StdError,
                                                 tvalue = coefs$tvalue, PrGTt = coefs$PrGTt, Significance = sig_stars)
            ##
            ##
            ## Plot  Observed vs. Predicte
            if (create_plots %in% c("y", "yes")) {
                pred_vals <- if (plot_mode == "best") {
                    predict(lm_fit, testData, type = 'response')  # Best fold only
                } else {
                    predict(lm_fit, data, type = 'response')  # All observed
                }
                ##
                ##
                z <- data.frame(if (plot_mode == "best") testData[[gene]] else data[[gene]], pred_vals)
                names(z) <- c(gene, 'Predicted')
                plot2 <- ggplot(z, aes(y = Predicted, x = if (plot_mode == "best") testData[[gene]] else data[[gene]])) +
                    geom_point(size = 0.5) +
                    theme_bw() +
                    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 12)) +
                    labs(title = paste(plot_mode, gene, nsplicesites,"~*Splice Sites"),
                         subtitle = paste("Pearson's correlation:", round(cor_val, 3), "| Adjusted R-squared:", round(max_adj_r_sq, 3)),
                         y = "Predicted Values",
                         x = paste(candidate, "Observed Values (loge(RPM))")) +
                    geom_smooth(method = lm, col = "grey30", lwd = 1.1)
                ggsave(file.path(base_dir, sub_Name, "plots", "fit", paste0(candidate, "_", gene, "_", sub_Name, "_exonLevel_CorrR2_", if (plot_mode == "best") "bestFold" else "fullDF", ".pdf")), plot2, device = "pdf")
            }
        }
        ## Linear model coefficients
        lm_coeffs_df <- do.call(rbind, lm_coeffs_list)
        write.table(lm_coeffs_df, file.path(base_dir, sub_Name, paste0(candidate, "_cv_lm_coefficients_", sub_Name, ".txt")), sep = "\t", row.names = FALSE, quote = FALSE)
        ##
        combined_rsq$cor_padj_BH <- p.adjust(combined_rsq$cor_pvalue, method = "BH")
        combined_rsq <- combined_rsq[order(-combined_rsq$mean_adj_r_sq), ]
        write.table(combined_rsq, file.path(base_dir, sub_Name, paste0(candidate, "_cv_lm_metrics_", sub_Name, ".txt")), sep="\t", row.names=FALSE, quote=FALSE)
        ##
        saveRDS(lm_models, file = file.path(base_dir, sub_Name, paste0(candidate, "_cv_lm_models_", sub_Name, ".rds")))
        cat("All LM models saved to:", file.path(base_dir, sub_Name, paste0(candidate, "_cv_lm_models_", sub_Name, ".rds")), "\n")
        cat("Subset processing completed:", sub_Name, "\n")
    }
    cat("Linear model cross-validation completed.\n")
}+

if (!(mode %in% c("lm", "glmnet", "both"))) {
    stop("Error: Invalid mode selected. Please choose 'lm', 'glmnet', or 'both'.")
}
cat("Script execution completed successfully.\n")

