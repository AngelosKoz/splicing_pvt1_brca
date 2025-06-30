set.seed(1234)
library(data.table)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

exon_exp <- fread("/data/counts/exonLevel_allCounts_filtbyexpr.tsv")
colnames(exon_exp)[1] <- "gene_id"

mani <- fread("/data/concat_manifest_with_details_V4.tsv")
salmonTPM <- fread("/data/salmon/salmon_predictors/salmon_pvt1tpm_fix.tsv")
colnames(salmonTPM)[1] <- "file_name"

maniInf <- mani[, c(3, 4)]
maniInf <- maniInf[!duplicated(maniInf$file_name), ]
salm <- merge(salmonTPM, maniInf, by = "file_name")
# write.table(salm, "/data/salmon/salmon_predictors/salmon_pvt1tpm_Subfix.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

salm_df <- fread("/data/salmon/salmon_predictors/salmon_pvt1tpm_Subfix.tsv", sep = "\t")
samples <- data.frame(file_name = salm_df$file_name)

identifier <- fread("/data/survival/brca_metadata/concat_manifest_with_details.tsv")
mani_entt <- identifier[, c(1, 2, 5)]

maniInf <- merge(maniInf, mani_entt, by = "file_name")
survival_data <- fread("/data/survival/brca_metadata/curated_survival_BRCA.txt")

surv_df <- merge(maniInf, samples, by = "file_name")
surv_df$sample <- substr(surv_df$entity_submitter_id, 1, 15)

final_df <- merge(surv_df, survival_data, by = "sample")
final_df <- final_df[, .(file_id, pam50_subtype, OS.time, OS, DFI.time, DFI, PFI.time, PFI, DSS.time, DSS)]
final_df$pam50_subtype <- factor(final_df$pam50_subtype, levels = c("LumA", "LumB", "Her2", "Basal"))

survival_metrics <- list(
    "OS" = list(time = "OS.time", event = "OS"),
    "DFI" = list(time = "DFI.time", event = "DFI"),
    "PFI" = list(time = "PFI.time", event = "PFI"),
    "DSS" = list(time = "DSS.time", event = "DSS")
  )

# All subtype & pairwise survival plots
subtype_pairs <- combn(levels(final_df$pam50_subtype), 2, simplify = FALSE)
for (metric in names(survival_metrics)) {
    metric_info <- survival_metrics[[metric]]
    dir.create(metric, showWarnings = FALSE)

    # All subtypes survival plot
    surv_object <- Surv(time = final_df[[metric_info$time]], event = final_df[[metric_info$event]])
    km_fit_all <- survfit(surv_object ~ pam50_subtype, data = final_df)
    subtype_counts <- table(final_df$pam50_subtype)
    legend_labels <- paste(names(subtype_counts), "(n =", subtype_counts, ")")

    plot_all <- ggsurvplot(
        km_fit_all, 
        data = final_df, 
        pval = TRUE, 
        censor = TRUE, 
        legend.title = "Subtype",
        legend.labs = legend_labels,
        xlab = "Time (days)",
        ylab = "Survival Probability",
        ggtheme = theme_minimal(),
        risk.table = TRUE,
        risk.table.title = "Number at risk by time"
      )

  ggsave(
    filename = file.path(metric, paste0(metric, "_All_Subtypes_SalmonTPM_KaplanMeier.pdf")),
    plot = plot_all$plot,
    width = 12, height = 12, device = "pdf"
  )

  # Pairwise survival plots
  for (pair in subtype_pairs) {
      pair_df <- final_df[final_df$pam50_subtype %in% pair, ]
      surv_object_pair <- Surv(time = pair_df[[metric_info$time]], event = pair_df[[metric_info$event]])
      km_fit_pair <- survfit(surv_object_pair ~ pam50_subtype, data = pair_df)
      pair_counts <- table(pair_df$pam50_subtype)
      legend_labels_pair <- paste(names(pair_counts), "(n =", pair_counts, ")")

      plot_pair <- ggsurvplot(
          km_fit_pair,
          data = pair_df,
          pval = TRUE,
          censor = TRUE,
          legend.title = "Subtype",
          legend.labs = legend_labels_pair,
          xlab = "Time (days)",
          ylab = "Survival Probability",
          ggtheme = theme_minimal(),
          risk.table = TRUE,
          risk.table.title = "Number at risk by time"
        )

    comparison_name <- paste(pair, collapse = "_vs_")
    ggsave(
        filename = file.path(metric, paste0(metric, "_", comparison_name, "_SalmonTPM_KaplanMeier.pdf")),
        plot = plot_pair$plot,
        width = 12, height = 12, device = "pdf"
      )
  }
}
