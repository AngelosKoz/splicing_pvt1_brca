library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

ss6 <- fread("/data/chromatin_associated_genes/pvt1/noCV_linear_model/linear_model_results_allGenes_noIR_top6SS.txt")
ss6 <- ss6[ss6$Model=="LM",]
ss12 <- fread("/data/chromatin_associated_genes/pvt1/noCV_linear_model/linear_model_results_allGenes_noIR_top12SS.txt")
ss12 <- ss12[ss12$Model=="LM",]
ss20 <- fread("/data/chromatin_associated_genes/pvt1/noCV_linear_model/linear_model_results_allGenes_noIR_top20SS.txt")
ss20 <- ss20[ss20$Model=="LM",]
ss34 <- fread("/data/chromatin_associated_genes/pvt1/noCV_linear_model/linear_model_results_SE_allGenes_noIR_allSS.txt")
ss34 <- ss34[ss34$Model=="LM",]
pvt1exp <- fread("/data/chromatin_associated_genes/pvt1/noCV_linear_model/linear_model_results_exonPVT1_allGenes_noIR_allSS.txt")
pvt1exp <- pvt1exp[pvt1exp$Model=="LM",]
pvt1expIR <- fread("/data/chromatin_associated_genes/pvt1/noCV_linear_model/linear_model_results_exonPVT1_allGenes.txt")
pvt1expIR <- pvt1expIR[pvt1expIR$Model=="LM",]


ss6$Dataset <- "PVT1 6 SS"
ss12$Dataset <- "PVT1 12 SS"
ss20$Dataset <- "PVT1 20 SS"
ss34$Dataset <- "PVT1 34 SS"

pvt1exp$Dataset <- "PVT1 loge(exonRPM)"
pvt1expIR$Dataset <- "PVT1 loge(exonRPM) + IR"

threshold <- 0.2
#all_data <- rbindlist(list(ss6, ss12, ss20, ss34), use.names = TRUE, fill = TRUE)
all_data <- rbindlist(list(pvt1exp, pvt1expIR, ss34), use.names = TRUE, fill = TRUE)
all_data <- rbindlist(list(pvt1exp, ss34), use.names = TRUE, fill = TRUE)
thresh <- all_data[, .(Count = sum(Statistic > threshold)), by = Dataset]
all_data <- merge(all_data, thresh, by = "Dataset")
all_data$plabel <- paste0(all_data$Dataset, " (n>", threshold, " = ", all_data$Count, ")")
all_data <- unique(all_data[, .(Statistic, plabel)])


porder <- c("PVT1 6 SS", "PVT1 12 SS", "PVT1 20 SS", "PVT1 34 SS")
porder <- c("PVT1 loge(exonRPM)", "PVT1 loge(exonRPM)", "PVT1 loge(exonRPM) + IR", "PVT1 34 SS")
porder <- c("PVT1 loge(exonRPM)", "PVT1 loge(exonRPM)", "PVT1 34 SS")

ordered_labels <- grep("PVT1", levels(factor(all_data$plabel)), value = TRUE)
ordered_labels <- ordered_labels[order(match(sub(" \\(.*", "", ordered_labels), porder))]
ordered_labels
all_data$plabel <- factor(all_data$plabel, levels = ordered_labels)

custom_colors <- c(
  "PVT1 loge(exonRPM)" = "red",
#  "PVT1 loge(exonRPM) + IR" = "steelblue",
  "PVT1 34 SS (n>0.2 = 770)" = "#7B3294"  # purple
)

p <- ggplot(all_data, aes(x = Statistic, color = plabel)) +
    stat_ecdf(geom = "step", size = 1) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "black") +
#    scale_color_manual(values = custom_colors) + # comment out for default
    labs(
        title = "ECDF of R² (Gene ~ PVT1 Splice Sites Splicing Efficiency)",
        subtitle = paste("Vertical line at R² =", threshold),
        x = expression(R^2),
        y = expression(F(R^2)),
        color = "Dataset (Genes with R² > 0.2)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
    )
print(p)
##ggsave("ECDF_pvt1exp_noIR.pdf", plot = p, width = 10, height = 10)

box_data <- all_data[plabel %in% c("PVT1 loge(exonRPM) (n>0.2 = 53)", "PVT1 34 SS (n>0.2 = 770)")]
box_data$plabel <- factor(box_data$plabel, levels = c("PVT1 loge(exonRPM) (n>0.2 = 53)", "PVT1 34 SS (n>0.2 = 770)"))

comparison_list <- list(c("PVT1 loge(exonRPM) (n>0.2 = 53)", "PVT1 34 SS (n>0.2 = 770)"))
box_data <- box_data[-30888]

p_box <- ggplot(box_data, aes(x = plabel, y = Statistic, fill = plabel)) +
    geom_boxplot(width = 0.5, outlier.size = 1, alpha = 0.7) +
    scale_y_continuous(breaks = seq(0, 0.6, 0.2), limits = c(0, .65)) +
#    stat_compare_means(
#        method = "t.test",
#        label = "p.signif", # swap for star/pval
#        comparisons = comparison_list,
#        label.y = 0.62 ) +
    stat_compare_means(
        method = "t.test",
        label = "p.format", # swap for star/pval
        comparisons = comparison_list,
        label.y = 0.57) +
    labs(
        title = "R² Distributions",
        x = NULL,
        y = expression(R^2)) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank())
print(p_box)
##ggsave("boxplot_pvt1exp_vs_ss34_pval_no_pvt1.pdf", plot = p_box, width = 10, height = 8)









######## dotted ecdf ################
ecdf_lines <- all_data[, {
  stat_sorted <- sort(Statistic)
  y_value <- ecdf(stat_sorted)(threshold)
  .(y_value_at_threshold = y_value)
}, by = plabel]
ecdf_lines


p <- ggplot(all_data, aes(x = Statistic, color = plabel)) +
  stat_ecdf(geom = "step", size = 1) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "black", size = 0.4) +
  geom_hline(data = ecdf_lines, aes(yintercept = y_value_at_threshold, color = plabel),
             linetype = "dotted", size = 0.6, show.legend = FALSE) +
  labs(
    title = "ECDF of R² (Gene ~ PVT1 Splicing Models)",
    subtitle = paste("Dashed: R² =", threshold, "| Dotted: % genes ≤ 0.2"),
    x = expression(R^2),
    y = expression(F(R^2)),
    color = "Dataset (Genes with R² > 0.2)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
print(p)
############ dotted ecdf #############



ggsave("ECDF_pvt1_allSS.pdf", plot = p, width = 10, height = 10)
