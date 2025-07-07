set.seed(1234)
library(data.table)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)

SE_df <- fread("/data/chromatin_associated_genes/pvt1/All_Subtypes_top_37_FeatLabels_SS3.altSS_L50")
samples <- data.frame(file_name = SE_df$V1)

identifier <- fread("/data/survival/brca_metadata/concat_manifest_with_details.tsv")
survival_data <- fread("/data/survival/brca_metadata/curated_survival_BRCA.txt")
surv_df <- merge(identifier, samples, by = "file_name")
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

# All subtype and pairwise survival plots
subtype_pairs <- combn(levels(final_df$pam50_subtype), 2, simplify = FALSE)
for (metric in names(survival_metrics)) {
    metric_info <- survival_metrics[[metric]]
    dir.create(metric, showWarnings = FALSE)
  
    surv_object <- Surv(final_df[[metric_info$time]], final_df[[metric_info$event]])
    km_fit_all <- survfit(surv_object ~ pam50_subtype, data = final_df)
  
    plot_all <- ggsurvplot(
        km_fit_all,
        data = final_df,
        pval = TRUE,
        censor = TRUE,
        legend.title = "Subtype",
        legend.labs = paste0(levels(final_df$pam50_subtype), " (n = ", table(final_df$pam50_subtype), ")"),
        xlab = "Time (days)",
        ylab = "Survival Probability",
        ggtheme = theme_minimal(),
        risk.table = TRUE,
        risk.table.title = "Number at risk by time"
    )
  
    ggsave(file.path(metric, paste0(metric, "_All_Subtypes_KaplanMeier.pdf")), plot = plot_all$plot, width = 12, height = 12)
  
    # Pairwise survival plots
    for (pair in subtype_pairs) {
        pair_df <- final_df[pam50_subtype %in% pair, ]
    surv_obj_pair <- Surv(pair_df[[metric_info$time]], pair_df[[metric_info$event]])
    km_fit_pair <- survfit(surv_obj_pair ~ pam50_subtype, data = pair_df)
    
    plot_pair <- ggsurvplot(
        km_fit_pair,
        data = pair_df,
        pval = TRUE,
        censor = TRUE,
        legend.title = "Subtype",
        legend.labs = paste0(names(table(pair_df$pam50_subtype)), " (n = ", table(pair_df$pam50_subtype), ")"),
        xlab = "Time (days)",
        ylab = "Survival Probability",
        ggtheme = theme_minimal(),
        risk.table = TRUE,
        risk.table.title = "Number at risk by time"
        )
        
    ggsave(file.path(metric, paste0(metric, "_", paste(pair, collapse = "_vs_"), "_KaplanMeier.pdf")),
           plot = plot_pair$plot, width = 12, height = 12)
  }
}
