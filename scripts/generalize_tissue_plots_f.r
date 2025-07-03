library(data.table)
library(dplyr)
library(ggplot2)

brca <- fread("/data/chromatin_associated_genes/pvt1/PVT1_cv_glmnet_metrics_brcatumor.txt")
brca_pass <- fread("/data/chromatin_associated_genes/pvt1/PVT1_enrich_res/glmnet/pass_genes.txt")
tissues <- c("Prostate", "Ovary", "Uterus", "Testis", "Adrenal")

for (tissue in tissues) {
  path <- paste0("/generalize/gen_", tolower(tissue), "Res/PVT1_generalized_glmnet_coeff_metrics.tsv")
  assign(tolower(tissue), fread(path))
}

for (tissue in tissues) {
  obj <- get(tolower(tissue))
  r2_02 <- obj %>% filter(R_squared > 0.2)
  assign(paste0(tolower(tissue), "_r2_02"), r2_02)
  assign(paste0(tolower(tissue), "_brca_02"), merge(brca, r2_02, by = "Gene"))
}

for (tissue in tissues) {
  obj <- get(tolower(tissue))
  r2_01 <- obj %>% filter(R_squared > 0.1)
  assign(paste0(tolower(tissue), "_r2_01"), r2_01)
  assign(paste0(tolower(tissue), "_brca_01"), merge(brca, r2_01, by = "Gene"))
}


rmse_list <- list(
  "BRCA Mean RMSE" = brca$mean_rmse,
  "BRCA Min RMSE" = brca$min_rmse
)
for (tissue in tissues) {
  rmse_list[[paste0(tissue, " RMSE")]] <- get(tolower(tissue))$RMSE
}

group_labels_base <- names(rmse_list)
n_models <- sapply(rmse_list, length)
group_labels_full <- paste0(group_labels_base, " (n=", n_models, ")")

rmse_df <- data.frame(
  value = unlist(rmse_list),
  group = factor(rep(group_labels_full, n_models), levels = group_labels_full),
  color_group = factor(rep(group_labels_base, n_models), levels = group_labels_base)
)

tissue_colors <- c(
  "BRCA Mean RMSE" = "#5e121c",
  "BRCA Min RMSE" = "#E41A1C",
  "Prostate RMSE" = "#4DAF4A",
  "Ovary RMSE" = "#377EB8",
  "Uterus RMSE" = "#984EA3",
  "Testis RMSE" = "#A65628",
  "Adrenal RMSE" = "#FF7F00"
)


p <- ggplot(rmse_df, aes(x = group, y = value, fill = color_group)) +
  geom_boxplot(alpha = 0.8, notch = TRUE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_colors, guide = "none") +
  labs(title = "RMSE Distributions Across Tissues", x = "Tissue", y = "RMSE") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("/generalize/rmse_boxplot.pdf", p, width = 10, height = 8)


common_genes <- brca_pass$Gene
rmse_common_list <- list(
  "BRCA Mean RMSE" = brca$mean_rmse[brca$Gene %in% common_genes],
  "BRCA Min RMSE" = brca$min_rmse[brca$Gene %in% common_genes]
)
for (tissue in tissues) {
  obj <- get(tolower(tissue))
  rmse_common_list[[paste0(tissue, " RMSE")]] <- obj$RMSE[obj$Gene %in% common_genes]
}

n_common_models <- sapply(rmse_common_list, length)
group_labels_common <- paste0(names(rmse_common_list), " (n=", n_common_models, ")")
rmse_common_df <- data.frame(
  value = unlist(rmse_common_list),
  group = factor(rep(group_labels_common, n_common_models), levels = group_labels_common),
  color_group = factor(rep(group_labels_base, n_common_models), levels = group_labels_base)
)


p_common <- ggplot(rmse_common_df, aes(x = group, y = value, fill = color_group)) +
  geom_boxplot(alpha = 0.8, notch = TRUE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 1, fill = "white", color = "black") +
  scale_fill_manual(values = tissue_colors, guide = "none") +
  labs(title = "RMSE Distributions Across Tissues (Common Genes with BRCA_365)", x = "Tissue", y = "RMSE") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("/generalize/rmse_boxplot_common_genes_brca365.pdf", p_common, width = 10, height = 8)
