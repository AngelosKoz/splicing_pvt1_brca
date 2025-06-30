set.seed(1234)
library(data.table)
library(ggplot2)


df <- fread("/data/chromatin_associated_genes/pvt1/condition_altSS_concat_Control_Tumor_SE_SS3.altSS_L50")
exondf_lm <- fread("/data/counts/allExon_allSS_RPM.txt")

variance <- sapply(df[, .SD, .SDcols = 2:35], var)
top_sets <- list(
  top6 = names(sort(variance, decreasing = TRUE)[1:6]),
  top12 = names(sort(variance, decreasing = TRUE)[1:12]),
  top20 = names(sort(variance, decreasing = TRUE)[1:20]),
  top34 = names(sort(variance, decreasing = TRUE)[1:34])
)


top_sets[1]
x <- grep("128100979:128100981", names(exondf_lm[1:5, 30888:30927]))
names(exondf_lm[1:5, 30888:30927])[x]

response_column <- grep("ENSG", names(exondf_lm), value = TRUE)
response_column <- c("ENSG00000251562.11", "ENSG00000015479.20")

results_df <- data.frame()
correlation_df <- data.frame()
for (set_name in names(top_sets)) {
    predictors <- top_sets[[set_name]]
    keepcolumns <- c("V1", predictors)
    df_lm <- merge(df[, ..keepcolumns], exondf_lm, by.x = "V1", by.y = "V1")
    predictors_bqt <- sapply(predictors, function(x) paste0("`", x, "`"))
    predictors_corF <- predictors
    message("Running ", set_name, " with ", length(predictors), " predictors")
    print(predictors_bqt)
    break
    for (gene in response_column) {
        formula_str <- paste(gene, "~", paste(predictors_bqt, collapse = "+"))
        formula <- as.formula(formula_str)
        lm_fit <- lm(formula, data = df_lm)
        r2 <- summary(lm_fit)$r.squared
        results_df <- rbind(results_df, data.frame(Gene = gene, Model = "LM", Setting = set_name, Statistic = r2))
        for (pred in predictors_corF) {
            cor_test <- cor.test(df_lm[[gene]], df_lm[[pred]], method = "pearson")
            correlation_df <- rbind(correlation_df, data.frame(Gene = gene, Predictor = pred, Setting = set_name, P_Value = cor_test$p.value))
        }
        fitted_values <- predict(lm_fit, newdata = df_lm)
        observed_values <- df_lm[[gene]]
        cor_val <- cor.test(observed_values, fitted_values)$estimate
        r_squared <- summary(lm_fit)$r.squared   
        plot_data <- data.frame(Observed = observed_values, Predicted = as.numeric(fitted_values))
        plot <- ggplot(plot_data, aes(x = Observed, y = Predicted)) +
            geom_point(size = 0.5) +
            geom_smooth(method = 'lm', color = "grey30", lwd = 1.1) +
            theme_bw() +
            labs(
                title = paste(gene, "~ Linear Model (Top", gsub("top", "", set_name), "Splice Sites)"),
                subtitle = paste("Pearson's correlation:", round(cor_val, 3), "| R-squared:", round(r_squared, 3)),
                x = paste(gene, "Observed Values (log(RPM))"),
                y = "Predicted Values"
            )
        ggsave(filename = file.path("fix_noCV_plots", paste0(gene, "_Top", gsub("top", "", set_name), "_corr_OvP.pdf")),
               plot = plot, width = 6, height = 5, device = "pdf")
  }
}
