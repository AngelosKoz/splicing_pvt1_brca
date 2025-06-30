set.seed(1234)
library(dplyr)
library(data.table)

df <- fread("/data/chromatin_associated_genes/pvt1/condition_altSS_concat_Control_Tumor_SE_SS3.altSS_L50")
variance <- sapply(df[, .SD, .SDcols = 2:35], var)
top6_var_cols <- names(sort(variance, decreasing = TRUE)[1:6])
top12_var_cols <- names(sort(variance, decreasing = TRUE)[1:12])
top20_var_cols <- names(sort(variance, decreasing = TRUE)[1:20])
top34_var_cols <- names(sort(variance, decreasing = TRUE)[1:34])

keepcolumns <- c("V1", top6_var_cols)
df <- df[, ..keepcolumns]

exondf_lm <- fread("/data/counts/allExon_allSS_RPM.txt")

#predictor_columns <- top6_var_cols
#predictor_columns <- top12_var_cols
#predictor_columns <- top20_var_cols
predictor_columns <- top34_var_cols # results might slightly differ since initially, set.seed was not set
predictor_columns_corF <- predictor_columns
predictor_columns_bqt <- sapply(predictor_columns, function(x) paste0("`", x, "`"))
response_column <- grep("ENSG", names(exondf_lm), value = TRUE)
##
results_df <- data.frame(Gene = character(), Model = character(), Statistic = numeric())
correlation_df <- data.frame(Gene = character(), Predictor = character(), P_Value = numeric())
##
##
n=0
for (gene in response_column) {
    n = n+ 1
    print(gene)
    print(n)
    formula_str <- paste(gene, "~", paste(predictor_columns_bqt, collapse = "+"))
    formula <- as.formula(formula_str)
    x <- model.matrix(formula, data = df_lm)[, -1]
    y <- df_lm[[gene]]
    ##
    lm_fit <- lm(formula, data = df_lm)
    results_df <- rbind(results_df, data.frame(Gene = gene, Model = "LM", Statistic = summary(lm_fit)$r.squared))
    ##
    for (pred in predictor_columns_corF) {
        cor_test <- cor.test(df_lm[[gene]], df_lm[[pred]], method = "pearson")
        correlation_df <- rbind(correlation_df, data.frame(Gene = gene, Predictor = pred, P_Value = cor_test$p.value))
    }
}




