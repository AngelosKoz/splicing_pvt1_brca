#!/bin/bash

#SBATCH --job-name=RfeatCount             # name of script, this is just for reporting/accounting purposes
#SBATCH --output=./RFC.out                # standard output file
#SBATCH --error=./RFC.err                 # standard error file
#SBATCH --nodes=1                         # number of nodes to allocate, if your application does not run in parallel (MPI) mode set this to 1
#SBATCH --ntasks=40                       # number of cores to allocate
#SBATCH --time=100-00:00:00               # set a limit on the total run time, hrs:min:sec
#SBATCH --mem=120G                        # memory to allocate
set -euo pipefail

#source /miniconda3/bin/activate
source /miniconda3/etc/profile.d/conda.sh
conda activate base

command -v parallel >/dev/null 2>&1 || {
   echo "ERROR: GNU Parallel not found in \$PATH"; exit 1; }

if [[ $# -lt 5 ]]; then
    echo "Usage: sbatch $0 <WORKDIR> <GTF_FILE> <OUTPUT_DIR> <COUNTS_LEVEL> <COUNT_ID>"
    echo "Example: sbatch $0 /path/workdir /path/gtf_file /path/output_dir exon gene_id"
    exit 1
fi

WORKDIR="$1"
GTF_FILE="$2"
OUTPUT_DIR="$3"
COUNTS_LEVEL="$4"  # e.g., exon, gene
COUNT_ID="$5"      # e.g., gene_id, exon_id

mkdir -p "${OUTPUT_DIR}/summaries"

echo -e "WORKDIR       : $WORKDIR
GTF           : $GTF_FILE
OUTPUT_DIR    : $OUTPUT_DIR
COUNTS_LEVEL  : $COUNTS_LEVEL
COUNT_ID      : $COUNT_ID
Started       : $(date)\n"

MAX_JOBS=4
THREADS_PER_JOB=10

export THREADS_PER_JOB COUNTS_LEVEL COUNT_ID GTF_FILE OUTPUT_DIR

# For paired-end reads, use -p
: << COMMENT
find "$WORKDIR" -name '*.bam' -print0 |
parallel -0 -j "$MAX_JOBS" --no-notice --bar --progress --plain \
    featureCounts -T "$THREADS_PER_JOB" -p -Q 255 \
                  -t "$COUNTS_LEVEL" -F GTF -g "$COUNT_ID" \
                  -a "$GTF_FILE" \
                  -o "$OUTPUT_DIR"/{/.}_counts.txt {} 
COMMENT

# For single-end reads, use -s 1
: << COMMENT
find "$WORKDIR" -name '*.bam' -print0 |
parallel -0 -j "$MAX_JOBS" --no-notice --bar --progress --plain \
    featureCounts -T "$THREADS_PER_JOB" -s 1 -Q 255 \
                  -t "$COUNTS_LEVEL" -F GTF -g "$COUNT_ID" \
                  -a "$GTF_FILE" \
                  -o "$OUTPUT_DIR"/{/.}_counts.txt {} 
COMMENT

mv "${OUTPUT_DIR}"/*.summary "${OUTPUT_DIR}/summaries/"

conda deactivate

echo "Finish time and date: $(date)"
