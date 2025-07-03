set.seed(1234)
library(data.table)
library(edgeR)
library(data.table)


##
##
## Read a single file and format it (output from featureCounts)
read_count_file <- function(file_path) {
    count_data <- fread(file_path, skip = 1)  # skip the header line
    patient_id <- tools::file_path_sans_ext(basename(file_path))
    patient_id <- sub("\\.rna_seq\\.genomic\\.gdc_realn_counts", "", patient_id)
    setnames(count_data, c("Geneid", "Chr", "Start", "End", "Strand", "Length", patient_id))
    count_data[, (patient_id) := as.numeric(get(patient_id))]
    keep_columns <- c("Geneid", patient_id)
    return(count_data[, ..keep_columns])
}
##
##
##
count_dir <- "allCounts_Control_noOver" #80 samples
count_files <- list.files(count_dir, pattern="rna_seq\\.genomic\\.gdc_realn_counts\\.txt$", full.names = TRUE)
count_dir <- "allCounts_Tumor_solid_noOver" #113 samples
count_files <- c(count_files, list.files(count_dir, pattern="rna_seq\\.genomic\\.gdc_realn_counts\\.txt$", full.names = TRUE))
ount_dir <- "allCounts_Tumor_noOver" #573 samples
count_files <- c(count_files, list.files(count_dir, pattern="rna_seq\\.genomic\\.gdc_realn_counts\\.txt$", full.names = TRUE))
##
##
count_list <- lapply(count_files, read_count_file)
length(count_list)
all_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)



##!##!## Start from here if already done the previous part -- preferable if you have Raw counts ##!##!##
all_counts <- fread("BRCA_exonLevel_allCounts_noOverlapRaw.tsv")
all_countsF <- as.data.frame(all_counts[,-1])  # removing the first column (Geneid)
rownames(all_countsF) <- all_counts$Geneid
all_countsF0 <- all_countsF[rowSums(all_countsF) != 0,]  # Remove rows with all zero counts
all_countsF <- all_countsF0[(rowSums(all_countsF0 >= 10) >= 3), ] # Atleast 10 counts, in atleast 3 samples
all_countsF$id <- rownames(all_countsF)
rownames(all_countsF1) <- all_countsF1$id
all_countsF <- all_countsF1[, 2]
##
## edgeR filterByExpr
##
dge <- DGEList(counts = all_countsF0, group = group)
keep <- filterByExpr(dge) # This is the filtering option
dge <- dge[keep,, keep.lib.sizes=FALSE]
all_countsF_filtered <- all_countsF0[keep,]
write.table(all_countsF_filtered, file = "exonLevel_allCounts_filtbyexpr.tsv", sep="\t",row.names=T, col.names=NA, quote=F) # our 30,887 list
