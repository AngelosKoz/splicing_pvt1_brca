library(data.table)
library(jsonlite)
library(dplyr)
library(tidyr)
library(argparse)

parser <- ArgumentParser(description = "Refine manifest script")
parser$add_argument("--bam_manifest", required = TRUE, help = "Path to the BAM manifest")
parser$add_argument("--json_bam", required = TRUE, help = "Path to the BAM JSON")
parser$add_argument("--json_snv", required = FALSE, help = "Path to the SNV JSON (optional)")
parser$add_argument("--id", required = TRUE, help = "ID prefix for the output file (e.g., Ovary, BRCA_Tumor)")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to save the output file")
parser$add_argument("--wxs_file", required = FALSE, help = "Path to the WXS file")
parser$add_argument("--wgs_file", required = FALSE, help = "Path to the WGS file")
parser$add_argument("--mode", required = TRUE, choices = c("wxs", "wgs", "both", "simple"), help = "Mode of operation: wxs, wgs, both, or simple")
args <- parser$parse_args()

if (args$mode == "wxs" && is.null(args$wxs_file)) {
    stop("Error: --wxs_file is required when --mode is 'wxs'")
}
if (args$mode == "wgs" && is.null(args$wgs_file)) {
    stop("Error: --wgs_file is required when --mode is 'wgs'")
}
if (args$mode == "both" && (is.null(args$wxs_file) || is.null(args$wgs_file))) {
    stop("Error: Both --wxs_file and --wgs_file are required when --mode is 'both'")
}

bam_files <- fread(args$bam_manifest)
json_data <- fromJSON(args$json_bam, flatten = TRUE)

outdir <- args$output_dir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

json_info <- data.frame(
    case_id = sapply(json_data$associated_entities, function(x) x$case_id),
    file_name = json_data$file_name
)

cat("Initial dimensions of JSON info:", dim(json_info), "\n")
cat("Unique case IDs:", length(unique(json_info$case_id)), "\n")
cat("Unique file names:", length(unique(json_info$file_name)), "\n")

json_info <- json_info[!duplicated(json_info$case_id), ]
cat("Dimensions after removing duplicates:", dim(json_info), "\n")

json_info$file_name <- sub("\\.rna_seq\\.genomic\\.gdc_realn\\.bam$", "", json_info$file_name)
mani <- json_info

output_file <- file.path(outdir, paste0(args$id, "_bam_manifest.txt"))
write.table(json_info, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Case ids saved at:", output_file, "\n")

if (args$mode == "simple") {
    cat("Mode is 'simple'. Only initial manifest created. Exiting...\n")
    quit(save = "no", status = 0)
}

######################
### Add WXS Data   ###
######################

if (args$mode %in% c("wxs", "both")) {
    wxs <- fread(args$wxs_file)
    cat("Initial dimensions of WXS data:", dim(wxs), "\n")

    if (!is.null(args$json_snv)) {
        json_data_wxs <- fromJSON(args$json_snv, flatten = TRUE)
        json_info_wxs <- json_data_wxs %>%
            unnest(associated_entities) %>%
            filter(grepl("wxs.*.maf", file_name, ignore.case = TRUE)) %>%
            select(entity_id, case_id, file_name) %>%
            distinct()

        wxs_filt <- wxs[, c(5, 6, 7, 8, 10, 14, 11, 12, 13, 18, 19, 143, 144, 33, 34)]
        print(head(wxs_filt), "\n")
        wxs_filt <- wxs_filt %>% left_join(json_info_wxs, by = c("Tumor_Sample_UUID" = "entity_id"))
        print(head(wxs_filt), "\n")
        cat("Dimensions after filtering WXS data:", dim(wxs_filt), "\n")

        wxs_add <- wxs_filt[, c(16, 17)]
        colnames(wxs_add)[2] <- "file_name_wxs"

        wxs_add_unique <- wxs_add %>%
            group_by(case_id) %>%
            summarise(file_name_wxs = first(file_name_wxs))  # Keeps the first entry per case_id
        mani_wxs <- mani %>%
            left_join(wxs_add_unique, by = "case_id")
        cat("Dimensions after joining WXS data with manifest:", dim(mani_wxs), "\n")

        output_manifest <- file.path(outdir, paste0(args$id, "_manifest_with_detailsWXS.tsv"))
        write.table(mani_wxs, output_manifest, quote = FALSE, sep = "\t", row.names = FALSE)
        cat("Manifest with WXS details saved at:", output_manifest, "\n")

        missing_files <- setdiff(wxs_add$file_name_wxs, mani_wxs$file_name_wxs)
        cat("\n\nMissing files:\n", missing_files, "\n\n\n")

        wxs_final <- wxs_filt[wxs_filt$case_id %in% unique(na.omit(mani_wxs$case_id)), ]
        colnames(wxs_final)[17] <- "file_name_wxs"

        # Handle empty mutation (rs) names
        mask <- wxs_final$dbSNP_RS == ""
        counter <- seq_len(sum(mask))
        wxs_final$dbSNP_RS[mask] <- paste0("rsNA", counter)

        # Separate mutations 
        wxs_rs <- wxs_final %>%
            filter(dbSNP_RS != "novel") %>%
            distinct(dbSNP_RS, .keep_all = TRUE)
        cat("Dimensions of wxs_rs:", dim(wxs_rs), "\n")

        wxs_novel <- wxs_final %>%
        filter(dbSNP_RS == "novel") %>%
        select(Chromosome, Start_Position, End_Position, Variant_Type, Strand) %>%
        distinct() %>%
        arrange(Chromosome, Start_Position, End_Position) %>%
        mutate(dbSNP_RS = paste0("novel", row_number()))

        keepcol <- colnames(wxs_novel)
        wxs_rs1 <- wxs_rs[, ..keepcol]
        wxs_uniq <- rbind(wxs_rs1, wxs_novel)

        print(head(wxs_final))
        cat("wxs_final \n\n\n")


        wxs_final_fix <- wxs_final %>%
            left_join(wxs_novel %>%
                select(Chromosome, Start_Position, End_Position, Strand, Variant_Type, newID = dbSNP_RS),
                by = c("Chromosome","Start_Position","End_Position","Strand","Variant_Type")) %>%
            mutate(dbSNP_RS = if_else(dbSNP_RS == "novel", newID, dbSNP_RS)) %>%
            select(-newID)
        print(head(wxs_final_fix)) 
        cat("wxs_final_FIX \n\n\n")

        ##
        ## We convert 0/0 0/1 and 1/1 genotype to 0, 1, 2 respectively
        wxs_parsed <- wxs_final_fix %>%
            mutate(tumor_gt_parsed = sub("^(.*?):.*", "\\1", vcf_tumor_gt)) %>%
            mutate(tumor_gt_code = case_when(
                    tumor_gt_parsed == "0/0" ~ 0,
                    tumor_gt_parsed == "0/1" ~ 1,
                    tumor_gt_parsed == "1/1" ~ 2,
                    TRUE                     ~ NA_real_))

        print(sum(wxs_parsed$tumor_gt_code))

        df_joined <- wxs_parsed %>%
            left_join(mani, by = "case_id", relationship = "many-to-many") %>%
            mutate(genotype_code = tumor_gt_code) %>%
            filter(!is.na(file_name))

        df_joined <- df_joined %>%
            distinct(file_name, dbSNP_RS, genotype_code) %>%
            complete(file_name = mani$file_name,
                     dbSNP_RS = unique(wxs_parsed$dbSNP_RS),
                     fill = list(genotype_code = NA))

        final <- df_joined %>%
            pivot_wider(id_cols     = file_name,
                        names_from  = dbSNP_RS,
                        values_from = genotype_code,
                        values_fill = NA)
        df_final <- as.data.frame(final)
        ##

        # Fill missing values with 0 for rows corresponding to non-NA file_name_wxs in mani_wxs. NA values in other rows means no WXS file was available
        non_na_files <- mani_wxs$file_name[!is.na(mani_wxs$file_name_wxs)]
        df_final_fill <- df_final %>%
            mutate(across(-file_name, ~ if_else(file_name %in% non_na_files & is.na(.x), 0, .x)))
        num_full_na_rows <- sum(rowSums(is.na(df_final_fill[, -1])) == ncol(df_final_fill[, -1]))
        num_rows_with_one <- sum(apply(df_final_fill[, -1], 1, function(row) any(row == 1, na.rm = TRUE)))
        total_ones <- sum(df_final_fill[, -1] == 1, na.rm = TRUE)

        cat("Number of rows fully NA in the final dataframe:", num_full_na_rows, "\n")
        cat("Number of rows with at least one '1' in the final dataframe:", num_rows_with_one, "\n")
        cat("Total number of '1's in the final dataframe:", total_ones, "\n")

        output_wxs <- file.path(outdir, paste0(args$id, "_wxs_mutRes_fill.tsv"))
        write.table(df_final_fill, output_wxs, quote = FALSE, sep = "\t", row.names = FALSE)
        cat("Final WXS data saved at:", output_wxs, "\n")
    }
}

######################
### Add WGS Data   ###
######################
## WGS file requires pre-processing:

#cat WGS_unfiltered.txt | awk -v rs_counter=1 'BEGIN {FS=OFS="\t"} NR == 1 {print "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "RS", "FORMAT", "NORMAL", "TUMOR", "file_name_wgs"; next}{
#    rs = "rsNA" rs_counter
#    if ($8 ~ /rs[0-9]+/) {
#        match($8, /rs[0-9]+/, id); rs = id[0]
#    } else {
#        rs_counter++
#    }
#    print $1, $2, $3, $4, $5, $6, $7, rs, $9, $10, $11, $12}' > out


if (args$mode %in% c("wgs", "both")) {
    wgs <- fread(args$wgs_file)
    cat("Initial dimensions of WGS data:", dim(wgs), "\n")

    if (!is.null(args$json_snv)) {
        json_data_wgs <- fromJSON(args$json_snv, flatten = TRUE)
        json_info_wgs <- json_data_wgs %>%
            unnest(associated_entities) %>%
            filter(grepl("wgs.*.vcf", file_name, ignore.case = TRUE)) %>%
            select(case_id, file_name) %>%
            distinct()
        colnames(json_info_wgs)[2] <- "file_name_wgs"

        mani_wgs <- mani %>%
            left_join(json_info_wgs, by = "case_id")
        cat("Dimensions after joining WGS data with manifest:", dim(mani_wgs), "\n")

        ##output_manifest <- file.path(outdir, paste0(args$id, "_manifest_with_detailsWGS.tsv"))
        ##write.table(mani_wgs, output_manifest, quote = FALSE, sep = "\t", row.names = FALSE)

        save_manifest <- mani_wgs %>%
            group_by(case_id, file_name) %>%
            summarise(file_name_wgs = paste(na.omit(file_name_wgs), collapse = ", "), .groups = "drop") %>%
            complete(case_id = mani$case_id, file_name = mani$file_name, fill = list(file_name_wgs = NA))

        output_save_manifest <- file.path(outdir, paste0(args$id, "_concat_manifest_with_detailsWGS.tsv"))
        write.table(save_manifest, output_save_manifest, quote = FALSE, sep = "\t", row.names = FALSE)
        cat("Saved manifest with combined WGS details at:", output_save_manifest, "\n")

        wgs_filt <- wgs %>%
            filter(file_name_wgs %in% mani_wgs$file_name_wgs)
        wgs_case <- wgs_filt %>%
            left_join(mani_wgs, by = "file_name_wgs")
        cat("Dimensions of filtered WGS data:", dim(wgs_filt), "\n")
        wgs_uniq <- wgs_case %>% distinct(CHROM, POS, REF, ALT, case_id, .keep_all = TRUE)

        print(head(wgs_uniq), "\n\n\n")


        wgs_rsNA <- wgs_uniq %>%
            filter(grepl("^rsNA", RS))
        wgs_rsID <- wgs_uniq %>%
            filter(!grepl("^rsNA", RS))
        cat("Dimensions of RS NA:", dim(wgs_rsNA), "\n")
        cat("Dimensions of RS ID:", dim(wgs_rsID), "\n")

        wgs_parsed <- wgs_uniq %>%
            mutate(tumor_gt_parsed = sub("^(.*?):.*", "\\1", TUMOR),
                   normal_gt_parsed = sub("^(.*?):.*", "\\1", NORMAL)) %>%
            mutate(tumor_gt_code = case_when(
                       tumor_gt_parsed %in% c("0/0", "0|0") ~ 0,
                       tumor_gt_parsed %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
                       tumor_gt_parsed %in% c("1/1", "1|1") ~ 2,
                       TRUE ~ NA_real_),
                   normal_gt_code = case_when(
                       normal_gt_parsed %in% c("0/0", "0|0") ~ 0,
                       normal_gt_parsed %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
                       normal_gt_parsed %in% c("1/1", "1|1") ~ 2,
                       TRUE ~ NA_real_))
        cat("Sum of tumor_gt_code:", sum(wgs_parsed$tumor_gt_code, na.rm = TRUE), "\n")
        cat("Sum of normal_gt_code:", sum(wgs_parsed$normal_gt_code, na.rm = TRUE), "\n")

        print(head(wgs_parsed))

        wgs_final_l <- wgs_parsed %>%
            select(file_name, RS, tumor_gt_code) %>%
            complete(file_name = mani_wgs$file_name,
                     RS = unique(wgs_parsed$RS),
                     fill = list(tumor_gt_code = NA))

        wgs_final <- wgs_final_l %>%
            pivot_wider(id_cols = file_name,
                        names_from = RS,
                        values_from = tumor_gt_code,
                        values_fill = NA)
        wgs_final <- as.data.frame(wgs_final)

        # Fill missing values with 0 for rows corresponding to non-NA file_name in mani_wgs
        non_na_files_wgs <- mani_wgs$file_name[!is.na(mani_wgs$file_name_wgs)]
        wgs_final_fill <- wgs_final %>%
            mutate(across(-file_name, ~ if_else(file_name %in% non_na_files_wgs & is.na(.x), 0, .x)))

        num_full_na_rows_wgs <- sum(rowSums(is.na(wgs_final_fill[, -1])) == ncol(wgs_final_fill[, -1]))
        num_rows_with_one_wgs <- sum(apply(wgs_final_fill[, -1], 1, function(row) any(row == 1, na.rm = TRUE)))
        total_ones_wgs <- sum(wgs_final_fill[, -1] == 1, na.rm = TRUE)

        cat("Number of rows fully NA in the WGS dataframe:", num_full_na_rows_wgs, "\n")
        cat("Number of rows with at least one '1' in the WGS dataframe:", num_rows_with_one_wgs, "\n")
        cat("Total number of '1's in the WGS dataframe:", total_ones_wgs, "\n")

        output_wgs <- file.path(outdir, paste0(args$id, "_wgs_mutRes_fill.tsv"))
        write.table(wgs_final_fill, output_wgs, quote = FALSE, sep = "\t", row.names = FALSE)
        cat("Final WGS data saved at:", output_wgs, "\n")
    }
}