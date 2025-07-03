library(data.table)
library(tidyverse)
library(ggplot2)
library(ggcorrplot)



se34_subtypes <- fread("/data/chromatin_associated_genes/pvt1/pvt134SE_tumorsub_varsortIR.tsv")
se34_all <- fread("/data/chromatin_associated_genes/pvt1/pvt134SE_varsortIR.tsv")
sal19 <- fread("/data/salmon/salmon_predictors/salmonTPM_top19_varsortIR.tsv")
se_sal <- fread("/data/salmon/salmon_predictors/pvt134SE_salmonTPM_top19_varsortIR.tsv")


plot_corr_heatmap <- function(data, plot_title = "Correlation Heatmap", output_path = NULL, annotate = TRUE) {
    if ("label" %in% colnames(data)) {
        data <- data %>% select(-c("label", "V1"))
    }
    numeric_data <- data %>% mutate(across(everything(), ~ as.numeric(as.character(.))))
    cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
    p <- ggcorrplot(
        cor_matrix,
        lab = annotate,
        lab_size = 3,
        method = "square",
        colors = c("blue", "white", "red"),
        title = plot_title,
        ggtheme = theme_minimal()
    )
    print(p)
    if (!is.null(output_path)) {
        dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
        ggsave(output_path, p, width = 30, height = 20)
    }
}



plot_corr_heatmap(sal19, plot_title = "Salmon TPM – Top 19 Isoforms", 
                  output_path = "pvt1_heatmaps/salmonTPM_669_top19_corr_heatmap.pdf")
plot_corr_heatmap(se34_subtypes, plot_title = "PVT1 SE Top 34 Sites", 
                  output_path = "pvt1_heatmaps/se34_592_corr_heatmap.pdf")
plot_corr_heatmap(se_sal, plot_title = "PVT1 SE + Salmon TPM – Top 19", 
                  output_path = "pvt1_heatmaps/se_sal_592_corr_heatmap.pdf")


#####
get_var_r <- function(input, directory, top_var = 10, top_ind = 0.1, gene_naming = "", plots = TRUE) {
    plots_dir <- file.path(directory, paste0(gene_naming, "_top_", top_var), "var_plots")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    label_col <- NULL
    if ("label" %in% names(input)) {
        label_col <- input$label
        input <- input %>% select(-label)
    }
    if ("V1" %in% names(input)) input <- input %>% select(-V1)
    if ("IntronRet" %in% names(input)) input <- input %>% select(-IntronRet)
    numeric_data <- input %>%
        mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
        select(where(is.numeric))
    ## variances
    variances <- apply(numeric_data, 2, var)
    sorted_var <- sort(variances, decreasing = TRUE)
    total_var <- sum(sorted_var)
    var_explained <- (sorted_var / total_var) * 100
    cumulative_var <- cumsum(var_explained)
    top_indices <- order(var_explained, decreasing = TRUE)[1:top_var]
    top_features <- names(var_explained)[top_indices]
    top_feature_table <- tibble(Splice_Site = top_features, Variance = var_explained[top_features])
    write_tsv(top_feature_table, file.path(directory, paste0(gene_naming, "_top_", top_var), paste0(gene_naming, "_Top_", top_var, "_Feature_Variances.tsv")))
    low_features <- names(var_explained)[1:10]
    low_feature_table <- tibble(Splice_Site = low_features, Variance = var_explained[low_features])
    write_tsv(low_feature_table, file.path(directory, paste0(gene_naming, "_top_", top_var), paste0(gene_naming, "_Lowest_10_Feature_Variances.tsv")))
    if (plots != FALSE) {
        ##
        ## Plot 1: Raw Variances
        ##
        df1 <- tibble(Splice_Site = names(sorted_var), Variance = sorted_var) %>%
            mutate(Percent = round((Variance / sum(Variance)) * 100, 1))
        p1 <- ggplot(df1, aes(x = reorder(Splice_Site, -Variance), y = Variance)) +
            geom_bar(stat = "identity", fill = "skyblue") +
            geom_text(aes(label = paste0(Percent, "%")), vjust = -0.5, size = 3) +
            geom_hline(yintercept = min(sorted_var[top_indices]), color = "red") +
            labs(title = "Raw Variance of Splice Sites", x = "Splice Sites", y = "Variance") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(file.path(plots_dir, paste0(gene_naming, "_RawVarBar.pdf")), p1, width = 16, height = 10)
        ##
        ## Plot 2: Percentage of Total Variance Explained
        ##
        df2 <- tibble(Splice_Site = names(var_explained), VarExplained = var_explained)
        p2 <- ggplot(df2, aes(x = 1:length(VarExplained), y = VarExplained)) +
            geom_line(color = "blue") +
            geom_point(color = "darkred", size = 2) +
            geom_hline(yintercept = min(var_explained[top_indices]), color = "red") +
            labs(title = "Percentage of Total Variance Explained", x = "Features", y = "% Variance Explained")
        ggsave(file.path(plots_dir, paste0(gene_naming, "_PercentVar.pdf")), p2, width = 14, height = 9)
        ##
        ##
        if (!is.null(label_col)) {
            input_full <- bind_cols(tibble(label = label_col), input[, ..top_features])
            long_df <- pivot_longer(input_full, -label, names_to = "splice_site", values_to = "value") %>%
                mutate(splice_site = factor(splice_site, levels = top_features))
            ##
            ## Box plot
            ##
            p_box <- ggplot(long_df, aes(x = splice_site, y = value, fill = label)) +
                geom_boxplot(outlier.size = 0.5, show.legend = TRUE) +
                #stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
                stat_summary(fun = mean, geom = "point",
                             aes(group = label),  # group-wise mean
                             shape = 23, size = 2, fill = "white", position = position_dodge(width = 0.75)) +
                labs(title = paste("Splicing Efficiency", gene_naming), y = "Splicing Efficiency") +
                theme(axis.text.x = element_text(angle = 90))
            ggsave(file.path(directory, paste0(gene_naming, "_top_", top_var), paste0(gene_naming, "_EDA_boxplot.pdf")),
                   p_box, width = 14, height = 9)
        }
    }   
    return(list(top = top_features, low = low_features))
}
