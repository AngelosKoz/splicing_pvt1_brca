library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(parallel)
library(data.table)
library(ggplot2)


outdir <- "GO_KEGG"
tissue <- "prostate" # ovary, uterus, testis and adrenal
rsq_n="rsq0_2" # or rsq0_1
prefix <- paste0(tissue, "_", rsq_n)


#results_df <- fread("/data/chromatin_associated_genes/pvt1/PVT1_enrich_res/pass_genes.txt") # for BRCA 365
results_df <- fread(paste0("/data/tissues/", tissue, "/gen_", tissue, "Res/", rsq_n, "/merged_brca.tsv"))
geneid <- as.data.frame(results_df$Gene)
names(geneid) <- "gene"
geneid$gene <- sub("\\..*", "", geneid$gene)
nrow(geneid)
##
##
### Biological Processes ###
GO_BP <- enrichGO(gene = geneid$gene,
                  keyType       = 'ENSEMBL',
                  OrgDb         = org.Hs.eg.db, 
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  minGSSize     = 5,
                  readable      = TRUE)
gobp <- as.data.frame(summary(GO_BP))
write.table(gobp, file.path(outdir, paste0(prefix, "_GO_BP.txt")), sep="\t", quote=F)
##head(gobp)
nrow(gobp)
##
##
## Dot plot for Biological Processes
pdf(file.path(outdir, paste0(prefix, "_GO_BP_DOT.pdf")), width = 10, height = 12)
dotplot(GO_BP) + 
  ggtitle("Biological Processes")
dev.off()
## Net plot for Biological Processes
pdf(file.path(outdir, paste0(prefix, "_GO_BP_NET.pdf")), width = 15, height = 10)
cnetplot(GO_BP, categorySize = 10) + 
  ggtitle("Biological Processes")
dev.off()
##
##
### Molecular Functions ###
GO_MF <- enrichGO(gene = geneid$gene,
                  keyType       = 'ENSEMBL',
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  minGSSize     = 5,
                  readable      = TRUE)
gomf <- as.data.frame(summary(GO_MF))
write.table(gomf, file.path(outdir, paste0(prefix, "_GO_MF.txt")), sep="\t", quote=F)
##head(gomf)
nrow(gomf)

##
##
## Dot plot for Molecular Functions
pdf(file.path(outdir, paste0(prefix, "_GO_MF_DOT.pdf")), width = 10, height = 12)
dotplot(GO_MF) + 
  ggtitle("Molecular Functions")
dev.off()
## Net plot for Molecular Functions
pdf(file.path(outdir, paste0(prefix, "_GO_MF_NET.pdf")), width = 15, height = 10)
cnetplot(GO_MF, categorySize = 10) + 
  ggtitle("Molecular Functions")
dev.off()
##
##
### Cellular Component ###
GO_CC <- enrichGO(gene = geneid$gene,
                  keyType       = 'ENSEMBL',
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  minGSSize     = 5,
                  readable      = TRUE)
gocc <- as.data.frame(summary(GO_CC))
write.table(gocc, file.path(outdir, paste0(prefix, "_GO_CC.txt")), sep="\t", quote=F)
##head(gocc)
nrow(gocc)
##
##
## Dot plot for Cellular Component
pdf(file.path(outdir, paste0(prefix, "_GO_CC_DOT.pdf")), width = 10, height = 12)
dotplot(GO_CC) + 
  ggtitle("Cellular Compartment")
dev.off()
## Net plot for Cellular Component
pdf(file.path(outdir, paste0(prefix, "_GO_CC_NET.pdf")), width = 15, height = 10)
cnetplot(GO_CC, categorySize = 10) + 
  ggtitle("Cellular Compartment")
dev.off()
