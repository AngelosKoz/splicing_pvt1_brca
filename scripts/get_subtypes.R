library(data.table)
library(genefu)
library(edgeR)
library(biomaRt)
library(jsonlite)
library(dplyr)
library(tidyr)

all_counts <- fread("/data/counts/BRCA_geneLevel_allCounts_noOverlapRaw.tsv")
tumor_ids <- fread("/data/tissues/tcga_metadata/BRCA/bam/tumor_sub_gdc_manifest.txt")
tumor_sample_ids <- sub("\\.rna_seq.*", "", tumor_ids$filename)

all_counts <- all_counts[, c("Geneid", colnames(all_counts)[colnames(all_counts) %in% tumor_sample_ids]), with = FALSE]
gene_ids <- all_counts$Geneid
counts_matrix <- as.matrix(all_counts[, -1, with = FALSE])
rownames(counts_matrix) <- gene_ids
counts_matrix <- counts_matrix[rowSums(counts_matrix) != 0,]  # Remove rows with all zero counts


## Use edgeR for normalisation
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge)
normalized_counts <- cpm(dge, log = TRUE)
dim(normalized_counts)
gene_ids <- rownames(normalized_counts)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

## Map Ensembl gene IDs to gene symbols
gene_mapping <- getBM(filters = "ensembl_gene_id_version", 
                      attributes = c("ensembl_gene_id_version", "hgnc_symbol"),
                      values = gene_ids,
                      mart = mart)

gene_mapping <- gene_mapping[!is.na(gene_mapping$hgnc_symbol), ] # remove NA
gene_mapping <- gene_mapping[gene_mapping$hgnc_symbol != "", ] #remove non-mapped

norm_counts_map <- normalized_counts[rownames(normalized_counts) %in% gene_mapping$ensembl_gene_id_version, ]
rownames(norm_counts_map) <- gene_mapping$hgnc_symbol[match(rownames(norm_counts_map), gene_mapping$ensembl_gene_id_version)]
duplicates <- duplicated(rownames(norm_counts_map))
sum(duplicates)  # How many duplicates there are
norm_counts_map <- norm_counts_map[!duplicates, ] # Remove duplicates, keep first occurrence

## load pam50 signature data
data(pam50.robust)
pam50_genes <- rownames(pam50.robust$centroids)
##
## check common between pam50 signature and our data
common_genes <- intersect(rownames(norm_counts_map), pam50_genes)
length(common_genes)
missing_genes <- setdiff(pam50_genes, common_genes)
missing_genes # genes from pam50 not in our set

## Keep only pam50 signature genes
norm_counts_map_filt <- t(norm_counts_map[rownames(norm_counts_map) %in% common_genes, ])
dim(norm_counts_map_filt)
annot_data <- data.frame(Gene.Symbol = rownames(norm_counts_map_filt))


## expression data:  samples as columns and genes as rows (HGNC symbols)
## Molecular Subtype classification using PAM50 signature
subtype_results <- molecular.subtyping(sbt.model = "pam50", 
                                       data = norm_counts_map_filt, 
                                       annot = annot_data, # reference vignette (https://www.bioconductor.org/packages/release/bioc/manuals/genefu/man/genefu.pdf)
                                       do.mapping = FALSE) # needs entrezid to map

saveRDS(subtype_results, file = "subtype_results.rds")
write.table(subtype_results$subtype, file = "subtype_name.csv", quote=FALSE, row.names=TRUE, col.names=NA, sep="\t")
write.table(subtype_results$subtype.proba, file = "subtype_probs.csv", quote=FALSE, row.names=TRUE, col.names=NA, sep="\t")
write.table(subtype_results$subtype.crisp, file = "subtype_class.csv", quote=FALSE, row.names=TRUE, col.names=NA, sep="\t")

subtype_results <- readRDS("Subtypes/subtype_results.rds")

sum(subtype_results$subtype.crisp[,1]) #155 Basal
sum(subtype_results$subtype.crisp[,2]) #87 Her2
sum(subtype_results$subtype.crisp[,3]) #263 LumA
sum(subtype_results$subtype.crisp[,4]) #182 LumB
sum(subtype_results$subtype.crisp[,5]) #9 Normal