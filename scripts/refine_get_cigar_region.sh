#!/bin/bash
set -ue
sortbam="" 
paired_end="" 
gene="" 
stranded="" 
CHROM="chr8"  # PVT1
START=$((127794541 - 10000))  # PVT1 - 10k
END=$((128187101 + 10000))    # PVT1 + 10k

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sortbam) sortbam="$2"; shift ;;
        --paired_end) paired_end="$2"; shift ;;
        --stranded) stranded="$2"; shift ;;
        --gene) gene="$2"; shift ;;
        --chrom) CHROM="$2"; shift ;;
        --start) START="$2"; shift ;;
        --end) END="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done


## Check for .bam files
bam_count=$(ls *.bam 2>/dev/null | wc -l)
if [[ "${bam_count}" -eq 0 ]]; then
    echo "No BAM files found in the current directory."
    exit 1
fi

### Below are argument checks ###
## Check for inclusion of a gene name.
#The -z tests if the string is non zero, since this is a mandatory argument.
if [[ -z "$gene" ]]; then
    echo "Error: No gene name provided. Please specify a gene using --gene option."
    exit 1
fi
##
## Only accepts sort and no as input. In any other case, this will stop the script
if [[ -n "$sortbam" && "$sortbam" != "sort" && "$sortbam" != "no" ]]; then
    echo "Invalid value for --sortbam. Only 'sort' or 'no' are accepted."
    exit 1
fi

## Only accepts y, yes and no as input. In any other case, this will stop the script
if [[ -n "$paired_end" && "$paired_end" != "y" && "$paired_end" != "yes" && "$paired_end" != "no" ]]; then
    echo "Invalid value for --paired_end. Please use 'y', 'yes', or 'no'."
    exit 1
fi

## Only accepts y, yes and no as input. In any other case, this will stop the script
if [[ -n "$stranded" && "$stranded" != "y" && "$stranded" != "yes" && "$stranded" != "no" ]]; then
    echo "Invalid value for --stranded. Please use 'y', 'yes', or 'no'."
    exit 1
fi

## Checks if the file exists and if not downloads it.
if [[ ! -f hg38.chrom.sizes ]]; then
    echo "Downloading chromosome sizes hg38."
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
fi


### First we will be doing some pre-processing from the gtf file, to extract only the gene information for later usage during the splicing extraction.
## Make sure only one gtf file is in the same directory
##
awk -v gene="\"$gene\";" -v chrom="$CHROM" '$0 ~ chrom && $3 == "exon" && $0 ~ "gene_name " gene' gencode*.gtf > "${gene}_exon.gtf"

## BED has chr start end ENSE_TranscrName "score" strand ##
##
# Get 5 prime splice sites (accounting for strand)
awk '{if ($7 == "+") print $1"\t"$5-1"\t"$5+1"\t"$24"\t"$20"\t"".""\t"$7; else print $1"\t"$4-1"\t"$4+1"\t"$24"\t"$20"\t"".""\t"$7}' "${gene}_exon.gtf" | awk '{ gsub(/[";]/, ""); print $1"\t"$2"\t"$3"\t"$4"_"$5"\t"$6"\t"$7}' > ss5.bed

# Get 3 prime splice sites (accounting for strand)
awk '{if ($7 == "+") print $1"\t"$4-1"\t"$4+1"\t"$24"\t"$20"\t"".""\t"$7; else print $1"\t"$5-1"\t"$5+1"\t"$24"\t"$20"\t"".""\t"$7}' "${gene}_exon.gtf" | awk '{ gsub(/[";]/, ""); print $1"\t"$2"\t"$3"\t"$4"_"$5"\t"$6"\t"$7}' > ss3.bed

# Sort the splice site BED files before grouping
sort -k1,1 -k2,2n -k3,3n ss3.bed > ss3.sorted.bed
sort -k1,1 -k2,2n -k3,3n ss5.bed > ss5.sorted.bed

## Addittion to run subtract
##
##
awk -v gene="\"$gene\";" -v chrom="$CHROM" '$0 ~ chrom && $3 == "gene" && $0 ~ "gene_name " gene' gencode*.gtf > "${gene}_gene.gtf"

##
# Correct for 5 prime (TES)
awk '{if ($7 == "+") print $1"\t"$5-1"\t"$5+1"\t"$14"\t"".""\t"$7; else print $1"\t"$4-1"\t"$4+1"\t"$14"\t"".""\t"$7}' "${gene}_gene.gtf" | awk '{ gsub(/[";]/, ""); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ss5.tes
sort -k1,1 -k2,2n -k3,3n ss5.tes > ss5.sorted.tes
bedtools subtract -s -A -a ss5.sorted.bed -b ss5.sorted.tes > ss5.sub.bed
##
# Correct for 3 prime (TSS)
awk '{if ($7 == "+") print $1"\t"$4-1"\t"$4+1"\t"$14"\t"".""\t"$7; else print $1"\t"$5-1"\t"$5+1"\t"$14"\t"".""\t"$7}' "${gene}_gene.gtf" | awk '{ gsub(/[";]/, ""); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > ss3.tss
sort -k1,1 -k2,2n -k3,3n ss3.tss > ss3.sorted.tss
bedtools subtract -s -A -a ss3.sorted.bed -b ss3.sorted.tss > ss3.sub.bed

# Collapse duplicates
bedtools groupby -i ss5.sub.bed -g 1,2,3 -c 4,5,6 -o first > "ss5uniq${gene}.bed"
bedtools groupby -i ss3.sub.bed -g 1,2,3 -c 4,5,6 -o first > "ss3uniq${gene}.bed"

rm ss3.bed ss5.bed ss3.tss ss5.tes ss*.sorted.* ss*.sub.bed


mkdir -p processed_bam_files
mkdir -p ss_final_out_"${gene}"/out_unfiltered
mkdir -p ss_final_out_"${gene}"/cutoff_0
mkdir -p ss_final_out_"${gene}"/cutoff_5
mkdir -p ss_final_out_"${gene}"/cutoff_10



for file in *.bam; do
    if [[ -f "$file" ]]; then
        naming="${file%.bam}"
        echo "Correcting for $file"
        if [[ ${sortbam} == "sort" ]]; then
            echo "Sorting bam option enabled. Now sorting ${naming}"
            samtools view -b "$file" "${CHROM}:${START}-${END}" | samtools sort -o - - | bedtools bamtobed -cigar -i - > "${naming}.bed"
        else
            echo "Skipping Sort"
            samtools view -b "$file" "${CHROM}:${START}-${END}" | bedtools bamtobed -cigar -i - > "${naming}.bed"
        fi
        ## Keep reads with mapping quality = 255
        ## Keep only Chromosomes 1-22,X,Y
        ## Remove any chromosomes that might contain "_"
        awk 'BEGIN {FS=OFS="\t"} ($5 == 255) && ($1 ~ /^chr([0-9]+|X|Y)$/) && ($1 !~ /_/) {print $1, $2, $3, $4, $6, $7}' "${naming}.bed" > "${naming}.qf.bed"
        rm "${naming}.bed"
        ## Since we discarded the mapping quality column above, we want to end up with 7 columns.
        ## For now we have 6 columns (1-4, 5-strand(+/-) and 6-split_information). Paired end will include on the 5th column the pair number. For single end we add a "0" column
        if [[ ${paired_end} == "y" || ${paired_end} == "yes" ]]; then
            echo "Paired end option selected. Correcting for ${naming}"
            if [[ ${stranded} == "y" || ${stranded} == "yes" ]]; then
                echo "Stranded option selected. Getting strand information"
                sed -i'' -e 's/\/1\t-/\t1\t+/g' "${naming}.qf.bed"
                sed -i'' -e 's/\/1\t+/\t1\t-/g' "${naming}.qf.bed"
                sed -i'' -e 's/\/2/\t2/g' "${naming}.qf.bed"
            else # If paired end and unstranded
                echo "Unstranded option is set as default. Ignoring strand info."
                sed -i'' -e 's/\//\t/g' "${naming}.qf.bed"
            fi
        else # If single end
            echo "Single end option is default. Correcting for ${naming}"
            if [[ ${stranded} == "y" || ${stranded} == "yes" ]]; then
                echo "Stranded option selected. Getting strand information"
                awk '{if ($5~/+/) print $1"\t"$2"\t"$3"\t"$4"\t""0""\t""-""\t"$7; else print $1"\t"$2"\t"$3"\t"$4"\t""0""\t""+""\t"$7}' "${naming}.qf.bed" > "${naming}.bed.cor"
                rm "${naming}.qf.bed"
                mv "${naming}.bed.cor" "${naming}.qf.bed"
            else #If single end and unstranded
                echo "Unstranded option is set as default. Ignoring strand info."
                awk '{print $1"\t"$2"\t"$3"\t"$4"\t0\t.\t"$6}' "${naming}.qf.bed" > "${naming}.bed.cor" # Added extra
                rm "${naming}.qf.bed" # Added extra
                mv "${naming}.bed.cor" "${naming}.qf.bed" # Added extra
            fi
        fi
        awk '$7 ~/N/' "${naming}.qf.bed" > "${naming}.split"
        awk '$7 !~/N/' "${naming}.qf.bed" > "${naming}.non_split"
        rm "${naming}.qf.bed"

        if [[ ${stranded} == "y" || ${stranded} == "yes" ]]; then
        ## bedtools coverage output uses 0-based: chrStart includes first base, chrEnd excludes last base: for 5' :  first 2 bp of the intron, for 3' : last 2 bp of the intron
            echo "Extracting coverage from 5' and 3' prime."
            ## Coverage for 5' prime if stranded (-s means we take into account strand info). First for split and then for non split reads
            bedtools coverage -s -a "ss5uniq${gene}.bed" -b "${naming}.split" > "${naming}_split_ss5"
            bedtools coverage -s -a "ss5uniq${gene}.bed" -b "${naming}.non_split" > "${naming}_non_split_ss5"
            ## Coverage for 3' prime
            bedtools coverage -s -a "ss3uniq${gene}.bed" -b "${naming}.split" > "${naming}_split_ss3"
            bedtools coverage -s -a "ss3uniq${gene}.bed" -b "${naming}.non_split" > "${naming}_non_split_ss3"
        else
            echo "Extracting coverage from 5' and 3' prime."
            ## Coverage for 5' prime if unstranded (we just remove -s)
            bedtools coverage -a "ss5uniq${gene}.bed" -b "${naming}.split" > "${naming}_split_ss5"
            bedtools coverage -a "ss5uniq${gene}.bed" -b "${naming}.non_split" > "${naming}_non_split_ss5"
            ## Coverage for 3' prime
            bedtools coverage -a "ss3uniq${gene}.bed" -b "${naming}.split" > "${naming}_split_ss3"
            bedtools coverage -a "ss3uniq${gene}.bed" -b "${naming}.non_split" > "${naming}_non_split_ss3"
        fi

        # Combine the non split and split splice sites for 5 prime
        paste "${naming}_split_ss5" "${naming}_non_split_ss5" > "${naming}.out_ss5"
        # Combine the non split and split splice sites for 3 prime
        paste "${naming}_split_ss3" "${naming}_non_split_ss3" > "${naming}.out_ss3"

        awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17}' "${naming}.out_ss5" > "${naming}.ss5"
        awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17}' "${naming}.out_ss3" > "${naming}.ss3"
    fi
    rm ${naming}_*split_ss* ${naming}.out_ss* "${naming}.split" "${naming}.non_split"

    ## Add splicing efficiency = split / (split + non-split)
    echo "Extracting splicing efficiency for ${naming}, almost there."

    # Will run by default these three cutoff options, adapt as needed
    # Cutoff 0
    python3 /data/akozonakis/get_splicing_eff.py "${naming}.ss3" "${naming}_filt_c0.ss3" --freedom 0
    python3 /data/akozonakis/get_splicing_eff.py "${naming}.ss5" "${naming}_filt_c0.ss5" --freedom 0
    # Cutoff 5
    python3 /data/akozonakis/get_splicing_eff.py "${naming}.ss3" "${naming}_filt_c5.ss3" --freedom 5
    python3 /data/akozonakis/get_splicing_eff.py "${naming}.ss5" "${naming}_filt_c5.ss5" --freedom 5
    # Cutoff 10
    python3 /data/akozonakis/get_splicing_eff.py "${naming}.ss3" "${naming}_filt_c10.ss3" --freedom 10
    python3 /data/akozonakis/get_splicing_eff.py "${naming}.ss5" "${naming}_filt_c10.ss5" --freedom 10

    mv "${naming}.ss"* ss_final_out_"${gene}"/out_unfiltered/
    mv "${naming}_filt_c0.ss"* ss_final_out_"${gene}"/cutoff_0/
    mv "${naming}_filt_c5.ss"* ss_final_out_"${gene}"/cutoff_5/
    mv "${naming}_filt_c10.ss"* ss_final_out_"${gene}"/cutoff_10/

    mv ${file}* processed_bam_files/ # Move the processed bam and bai files to a new directory so as to not re-do in case of script failure

done


