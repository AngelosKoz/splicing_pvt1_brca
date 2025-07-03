#!/bin/bash

#SBATCH --job-name=customSE           # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./customSE.out          # standard output file
#SBATCH --error=./customSE.err           # standard error file
#SBATCH --nodes=1                       # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=20                      # number of cores to allocate
#SBATCH --time=100-00:00:00                 # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=60G                        # memory to allocate

VENV_PATH=/data/spliceenv
source $VENV_PATH/bin/activate

if [[ $# -lt 9 ]]; then
    echo "Usage: sbatch $0 <WORKDIR> <SUFFIX> <SORTBAM> <STRANDED> <PAIRED_END> <GENE> <CHROM> <START> <END>"
    exit 1
fi

WORKDIR="$1"
SUFFIX="$2"
SORTBAM="$3"
STRANDED="$4"
PAIRED_END="$5"
GENE="$6"
CHROM="$7"
START="$8"
END="$9"


cd "$WORKDIR" || { echo "Error: Cannot change to directory ${WORKDIR}"; exit 1; }
echo "Current working directory: ${WORKDIR}"
echo "Start time and date: $(date)"

mkdir -p bam_files_"${SUFFIX}"
find . -type f \( -name "*.bam" -o -name "*.bam.bai" \) -exec mv {} bam_files_"${SUFFIX}"/ \;
cd "bam_files_${SUFFIX}" || { echo "Error: Cannot change to directory bam_files_${SUFFIX}"; exit 1; }
ln -s /data/gtfs/gencode.v44.annotation.gtf gencode.v44.annotation.gtf # this makes sure a gtf file is present in the current directory


## Check if the files are sorted. 1 means sorted, 0 means unsorted (based on coordinates)
#: << COMMENT
output_file="bam_sortRes_${SUFFIX}.txt"
for bam_file in ./*.bam; do
    if [[ -f "$bam_file" ]]; then
        is_sorted=$(samtools view -H "$bam_file" | grep -c "SO:coordinate")
        if [[ "$is_sorted" -eq 0 ]]; then
            echo "$(basename "$bam_file") is not sorted. Sorting now..."
            sorted_bam="${bam_file%.bam}.sorted.bam"
            samtools sort -o "$sorted_bam" "$bam_file" && mv "$sorted_bam" "$bam_file"
        else
            echo "$(basename "$bam_file") is already sorted."
        fi
        echo "$(basename "$bam_file"): $is_sorted" >> "$output_file"
    fi
done

## Index bam files that do not contain a .bai
echo "Starting indexing"
for bam in ./*.bam; do [ -f "$bam.bai" ] || (echo "Indexing $(basename "$bam")" && samtools index "$bam"); done
echo "Indexing completed."
#COMMENT




### ------------------------- ### 
########   Runs for BRCA  ########
#sbatch cluster_SE.sh --WORKDIR /bam_files_dir --SUFFIX suffix_name --SORTBAM y||no --STRANDED y||no --PAIRED_END y||no --GENE gene_name --CHROM chrN --START pos_start --END pos_end



#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y PVT1 chr8 127784541 128197101  # ENSG00000249859.14
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y ERBB4 chr2 211365717 212548841 # ENSG00000178568.16
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y FTX chrX 73930435 74303574     # ENSG00000230590.13
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y RAD51B chr14 67809779 68740218
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y SLC30A8 chr8 116940273 117186714
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y DLEU2 chr13 49956670 50125720  # ENSG00000231607.15
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y DRAIC chr15 69462921 69843120  # ENSG00000245750.11
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y CASC15 chr6 21664184 22654455  #  ENSG00000272168.10
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y TPRG1 chr3 188947214 189325304 # ENSG00000188001.10
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y PTPRT chr20 42072752 43189970  # ENSG00000196090.14
#sbatch cluster_SE.sh tcga_data/BRCA_all_soft/ BRCA_all no no y SMYD3 chr1 245749342 246507312 # ENSG00000185420.19

### ------------------------- ### 

echo "Extracting Splicing Efficiency"
/data/akozonakis/refine_get_cigar_region.sh \
    --sortbam "$SORTBAM" \
    --stranded "$STRANDED" \
    --paired_end "$PAIRED_END" \
    --gene "$GENE" \
    --chrom "$CHROM" \
    --start "$START" \
    --end "$END"
echo "Splicing Efficiency extraction completed."


echo "Separating TSS and SS"
bash /data/akozonakis/get_tss_altSS.sh ${WORKDIR}/bam_files_${SUFFIX}/ss_final_out_${GENE}


###### count uniquely-mapped reads ######
# count uniquely-mapped reads in parallel (~6-8h for 400 files)
#: << COMMENT

shopt -s nullglob
bam_list=( *.bam )
(( ${#bam_list[@]} )) || { echo "No BAM files – aborting"; exit 1; }

printf '%s\n' "${bam_list[@]}" > bam_list.txt        # one path per line

threads_bam=1
jobs_parallel=${SLURM_NTASKS:-30}
chroms=$(printf "chr%s " {1..22} X Y)

umrs="${SUFFIX}_uniquely_mapped_reads.txt"
printf "file_name\tTotal_Reads\tRPM_NormFactor\n" > "$umrs"


parallel --no-notice --bar --jobs "$jobs_parallel" -a bam_list.txt '
    echo "### counting {}" >&2                       # progress -> stderr

    reads=$(samtools view -@ '"$threads_bam"' -c -q 255 -F 4 {} '"$chroms"')
    if [ "$reads" -eq 0 ]; then
        printf "%s\t0\t0\n" "{/.}"
    else
        rpm=$(awk -v n="$reads" "BEGIN{printf \"%.6f\", 1e6/n}")
        printf "%s\t%d\t%s\n" "{/.}" "$reads" "$rpm"
    fi
' >> "$umrs"

sed -i 's/\.rna_seq\.genomic\.gdc_realn//g' "$umrs"
echo "Finished uniquely-mapped-read counting at $(date)"
#COMMENT

deactivate

echo "Finish time and date: ${date}"
