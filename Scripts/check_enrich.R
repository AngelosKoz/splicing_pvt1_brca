library(data.table)
library(argparse)
set.seed(1234)

parse_column_threshold <- function(arg, default_threshold) {
    parts <- strsplit(arg, ",")[[1]]
    column <- parts[1]
    threshold <- ifelse(length(parts) > 1, as.numeric(parts[2]), default_threshold)
    return(list(column = column, threshold = threshold))
}


perform_fisher_test <- function(lm_genes, filtered_genes, mirnaTG_genes) {
    total_genes <- unique(lm_genes)
    cat("Total unique genes:", length(total_genes), "\n")
    # Genes that are miRNA targets
    mirna_targets <- intersect(total_genes, mirnaTG_genes)
    non_mirna_targets <- setdiff(total_genes, mirnaTG_genes)
    cat("Number of miRNA target genes:", length(mirna_targets), "\n")
    cat("Number of non-miRNA target genes:", length(non_mirna_targets), "\n")
    pass_cutoff <- filtered_genes
    dont_pass_cutoff <- setdiff(total_genes, filtered_genes)
    cat("Number of genes passing the cutoff:", length(pass_cutoff), "\n")
    # Contigency table
    in_mirna_and_pass <- length(intersect(mirna_targets, pass_cutoff))
    in_mirna_and_dont_pass <- length(intersect(mirna_targets, dont_pass_cutoff))
    not_in_mirna_and_pass <- length(intersect(non_mirna_targets, pass_cutoff))
    not_in_mirna_and_dont_pass <- length(intersect(non_mirna_targets, dont_pass_cutoff))
    cat("Fisher's test inputs:\n")
    cat("  miRNA targets passing cutoff:", in_mirna_and_pass, "\n")
    cat("  miRNA targets not passing cutoff:", in_mirna_and_dont_pass, "\n")
    cat("  Non-miRNA targets passing cutoff:", not_in_mirna_and_pass, "\n")
    cat("  Non-miRNA targets not passing cutoff:", not_in_mirna_and_dont_pass, "\n")
    contingency_table <- matrix(c(in_mirna_and_pass, in_mirna_and_dont_pass, not_in_mirna_and_pass, not_in_mirna_and_dont_pass), nrow = 2, byrow = TRUE)
    dimnames(contingency_table) <- list(
        'miRNA Target' = c('Yes', 'No'),
        'Pass Cutoff' = c('Yes', 'No'))
    fisher_test_result <- fisher.test(contingency_table)
    fisher_result_df <- data.frame(
        OddsRatio = fisher_test_result$estimate,
        PValue = fisher_test_result$p.value,
        ConfIntLow = fisher_test_result$conf.int[1],
        ConfIntHigh = fisher_test_result$conf.int[2]
    )
    return(list(table = contingency_table, result = fisher_result_df))
}


parser <- ArgumentParser(description = "Check for enrichment using Fisher's test.")
parser$add_argument("--input_file", required = TRUE, help = "Path to the linear model results file.")
parser$add_argument("--target_file", required = TRUE, help = "Path to the miRNA target genes file.")
parser$add_argument("--output_dir", required = TRUE, help = "Path to save the enrichment results.")
parser$add_argument("--filter1", required = TRUE, help = "First filter column and optional threshold (e.g., 'max_adj_r_sq,0.2').")
parser$add_argument("--filter2", required = FALSE, help = "Second filter column and optional threshold (e.g., 'cor_value,0').")
parser$add_argument("--filter3", required = FALSE, help = "Third filter column and optional threshold, used for upper bound calculation (e.g., 'mean_rmse,1.5').")
args <- parser$parse_args()


outdir <- args$output_dir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
cat("-")
filter1 <- parse_column_threshold(args$filter1, 0.2)
print(filter1)
filter2 <- if (!is.null(args$filter2)) parse_column_threshold(args$filter2, 0) else NULL
filter3 <- if (!is.null(args$filter3)) parse_column_threshold(args$filter3, 1.5) else NULL


if (!file.exists(args$input_file)) {
    stop("Error: Input file does not exist: ", args$input_file)
}
if (!file.exists(args$target_file)) {
    stop("Error: Target file does not exist: ", args$target_file)
}


lm_results <- fread(args$input_file)
mirnaTG <- fread(args$target_file, header = TRUE)
mirnaTG_genes <- na.omit(mirnaTG$Geneid)
cat("Number of miRNA target genes loaded:", length(mirnaTG_genes), "\n")


# Calculate IQR-based upper bound (filter 3)
if (!is.null(filter3)) {
    iqr <- IQR(lm_results[[filter3$column]], na.rm = TRUE)
    q3 <- quantile(lm_results[[filter3$column]], 0.75, na.rm = TRUE)
    upper_bound <- q3 + (filter3$threshold * iqr)
    cat("Filter3 upper bound calculated:", upper_bound, "\n")
}


if (!(filter1$column %in% colnames(lm_results))) {
    stop("Error: Column specified in filter1 does not exist in the input file: ", filter1$column)
}
if (!is.null(filter2) && !(filter2$column %in% colnames(lm_results))) {
    stop("Error: Column specified in filter2 does not exist in the input file: ", filter2$column)
}
if (!is.null(filter3) && !(filter3$column %in% colnames(lm_results))) {
    stop("Error: Column specified in filter3 does not exist in the input file: ", filter3$column)
}


cat("Filtering results...\n")
print(head(lm_results))
cond1 <- !is.na(lm_results[[filter1$column]]) & lm_results[[filter1$column]] > filter1$threshold
cond2 <- if (!is.null(filter2)) (!is.na(lm_results[[filter2$column]]) & lm_results[[filter2$column]] < filter2$threshold) else TRUE # for glmnet
#cond2 <- if (!is.null(filter2)) (!is.na(lm_results[[filter2$column]]) & lm_results[[filter2$column]] > filter2$threshold) else TRUE # for lm 
cond3 <- if (!is.null(filter3)) (!is.na(lm_results[[filter3$column]]) & lm_results[[filter3$column]] < upper_bound) else TRUE
filtered_results <- lm_results[cond1 & cond2 & cond3]


filtered_genes <- unique(filtered_results$Gene)
lm_genes <- unique(lm_results$Gene)
cat("Number of genes in linear model results:", length(lm_genes), "\n")
cat("Number of genes after filtering:", length(filtered_genes), "\n")
filtered_genes_df <- data.table(Gene = filtered_genes, mirnaTG = ifelse(filtered_genes %in% mirnaTG_genes, "yes", "no"))
write.table(filtered_genes_df, file = file.path(outdir, "pass_genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


cat("Performing Fisher's test...\n")
fisher_results <- perform_fisher_test(lm_genes, filtered_genes, mirnaTG_genes)
cat("Fisher's test results:\n")
print(fisher_results$result)
cat("Contingency table:\n")
print(fisher_results$table)

cat("Saving results...\n")
write.table(fisher_results$table, file = file.path(outdir, "contingency_table.tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
write.table(fisher_results$result, file = file.path(outdir, "fisher_test_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Enrichment analysis completed. Results saved to:", outdir, "\n")