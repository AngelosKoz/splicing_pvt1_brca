library(data.table)
library(argparse)
set.seed(1234)

parser <- ArgumentParser(description = "Merge count files")
parser$add_argument("--input_dir", required = TRUE, help = "Directory containing count files")
parser$add_argument("--output_file", required = TRUE, help = "Output file path")
parser$add_argument("--bam_type", required = TRUE, choices = c("full", "sliced"), help = "Type of BAM file: 'full' or 'sliced'")
args <- parser$parse_args()

file_pattern <- if (args$bam_type == "full") {
    "rna_seq\\.genomic\\.gdc_realn_counts\\.txt$"
} else {
    "_counts_genePVT1\\.txt$"
}

count_files <- list.files(args$input_dir, pattern = file_pattern, full.names = TRUE)
cat("Files to be processed:\n")
print(count_files)

# Read and format a single count file
read_count_file <- function(file_path) {
    cat("Processing file:", file_path, "\n")  
    count_data <- fread(file_path, skip = 1)
    patient_id <- tools::file_path_sans_ext(basename(file_path))
    patient_id <- if (args$bam_type == "full") {
        sub("\\.rna_seq\\.genomic\\.gdc_realn_counts", "", patient_id)  # original
    } else {
        sub("\\_counts_genePVT1", "", patient_id)  # for sliced gene PVT1 (gene level expression)
        patient_id <- gsub("_", "-", patient_id) # for sliced gene PVT1
    }
    setnames(count_data, c("Geneid", "Chr", "Start", "End", "Strand", "Length", patient_id))
    count_data[, (patient_id) := as.numeric(get(patient_id))]
    keep_columns <- c("Geneid", patient_id)
    return(count_data[, ..keep_columns])
}

count_list <- lapply(count_files, read_count_file)
all_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)
write.table(all_counts, args$output_file, row.names = FALSE, sep = "\t", quote = FALSE)
