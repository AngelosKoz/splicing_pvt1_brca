library(data.table)
library(argparse)
library(ggplot2)
set.seed(1234)


parser <- ArgumentParser(description = "Generate Generalization Plot")
parser$add_argument("--calculated_metrics", required = TRUE, help = "Path to the calculated metrics file (metrics_df)")
parser$add_argument("--gen_metrics", required = TRUE, help = "Path to the external metrics file for comparison")
parser$add_argument("--gen_col", required = TRUE, help = "Column name in the external metrics file to use for comparison")
parser$add_argument("--gen_plot", required = TRUE, choices = c("rmse", "rsq"), help = "Metric for the generalization plot ('rmse' or 'rsq')")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to save the generalization plot")
parser$add_argument("--candidate", required = TRUE, help = "Candidate name for the plot")
parser$add_argument("--mode", required = TRUE, help = "Mode of regression (e.g., 'lm', 'glmnet')")
parser$add_argument("--input_mode", required = TRUE, help = "Input mode (e.g., 'coeff', 'object')")
args <- parser$parse_args()


# Generalization calculated metrics
if (!file.exists(args$calculated_metrics)) {
    stop("Error: The specified --calculated_metrics file does not exist.")
}
calculated_metrics <- fread(args$calculated_metrics)
cat("Calculated metrics loaded. Rows:", nrow(calculated_metrics), "\n")


if (!file.exists(args$gen_metrics)) {
    stop("Error: The specified --gen_metrics file does not exist.")
}
external_metrics <- fread(args$gen_metrics)
cat("External metrics loaded. Rows:", nrow(external_metrics), "\n")


# Check if the specified column exists in the external metrics file
if (!(args$gen_col %in% colnames(external_metrics))) {
    stop(paste("Error: Column", args$gen_col, "not found in the provided --gen_metrics file."))
}


external_metrics <- external_metrics[, .(Gene, Original_metric = get(args$gen_col))]
merged_metrics <- merge(calculated_metrics, external_metrics, by = "Gene", all = FALSE)
if (nrow(merged_metrics) == 0) {
    stop("Error: No common genes found between the generalization metrics and the external metrics.")
}
cat("Merged metrics. Rows:", nrow(merged_metrics), "\n")


y_metric <- if (args$gen_plot == "rmse") "RMSE" else "R_squared"
x_metric <- "Original_metric"
plot <- ggplot(merged_metrics, aes(x = .data[[x_metric]], y = .data[[y_metric]])) +
    geom_point(size = 0.6, alpha = 0.4) +
    geom_smooth(method = "lm", color = "red", linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotdash", color = "blue") +  # Add identity line
    theme_minimal() +
    labs(
        title = paste("Generalization Plot:", toupper(args$gen_plot)),
        x = paste("External Metric (", args$gen_col, ")", sep = ""),
        y = paste("Generated Metric (", y_metric, ")", sep = "")
    )


# Pearson
cor_result <- cor.test(merged_metrics[[x_metric]], merged_metrics[[y_metric]], method = "pearson")
print(cor_result)
corr_value <- round(cor_result$estimate, 3)
cor_pval <- formatC(cor_result$p.value, format = "e")  # Format p-value in scientific notation
# Spearman
spear_cor_result <- cor.test(merged_metrics[[x_metric]], merged_metrics[[y_metric]], method = "spearman")
print(spear_cor_result)
spear_corr_value <- round(spear_cor_result$estimate, 3)
spear_cor_pval <- formatC(spear_cor_result$p.value, format = "e")  # Format p-value in scientific notation


plot <- plot +
    annotate("text", x = min(merged_metrics[[x_metric]], na.rm = TRUE), 
             y = max(merged_metrics[[y_metric]], na.rm = TRUE), 
             label = paste("Pearson cor value =", corr_value, "\n", "Pearson cor pval =", cor_pval, "\n",
                           "Spearman rho =", spear_corr_value, "\n", "Spearman rho pval =", spear_cor_pval), 
             hjust = 0, vjust = 1, size = 4, color = "blue")

plot_dir <- file.path(args$output_dir, "plots", "generalization")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
plot_file <- file.path(plot_dir, paste0(args$candidate, "_", args$mode, "_", args$input_mode, "_generalization_", args$gen_plot, ".pdf"))
ggsave(plot_file, plot, width = 12, height = 10, device = "pdf")
cat("Generalization comparison plot saved to:", plot_file, "\n")
