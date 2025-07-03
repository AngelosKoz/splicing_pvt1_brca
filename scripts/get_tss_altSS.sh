#!/bin/bash

if [[ $# -eq 1 ]]; then
    TARGET_FOLDER="$1"
    if [[ ! -d "$TARGET_FOLDER" ]]; then
        echo "Error: Specified folder '$TARGET_FOLDER' does not exist."
        exit 1
    fi
    echo "Processing specified folder: $TARGET_FOLDER"
    FOLDERS=("$TARGET_FOLDER")
else
    echo "No folder specified. Processing all 'ss_final_out_*' folders."
    FOLDERS=($(pwd)/ss_final_out_*)
fi

if [[ ! -f hg38.fa ]]; then
    echo "Downloading fasta file hg38 -- v44 used"
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz
    mv GRCh38.primary_assembly.genome.fa hg38.fa
fi

for main_folder in "${FOLDERS[@]}"; do
    echo "Processing folder: $main_folder"

    for cutoff_folder in "$main_folder"/cutoff_*; do
        echo "Inside cutoff folder: $cutoff_folder"

        mkdir -p "$cutoff_folder/noncanonSS"
        mkdir -p "$cutoff_folder/altSS"
        mkdir -p "$cutoff_folder/fasta_seq"

        for ss_file in "$cutoff_folder"/*.rna_seq.genomic.gdc_realn_filt*.ss*; do
            [[ -e "$ss_file" ]] || continue
            echo "      – splice-site table: $(basename "$ss_file")"

            #
            case "$ss_file" in
                *.ss3) cls="ss3"; motif="AG" ;;
                *.ss5) cls="ss5"; motif="GT" ;;
                *)     echo "        ! Unknown suffix, skip"; continue ;;
            esac

            tail -n +2 "$ss_file" | awk -v type="$cls" '
            BEGIN { FS=OFS="\t" }
            {
                if (type=="ss3") {                          # 3′ splice site
                    if ($6~/+/)        { $2 -= 2 }              # + strand: Start-2…End  -- $2 = start, $3 = end
                    else               { $2 -= 1; $3 += 1 }     # − strand: Start-1…End+1
                } else {                                    # 5′ splice site
                    if ($6~/+/)        { $2 -= 1; $3 += 1 }     # + strand: Start-1…End+1
                    else               { $2 -= 2 }              # − strand: Start-2…End
                }
                print
            }' > fasta_fix.bed


            # Get fasta sequences
            bedtools getfasta -tab -s -fi ./hg38.fa -bed fasta_fix.bed -fo splice_site_sequences
            awk 'BEGIN {FS="\t"; OFS="\t"} {match($1, /^(.+):([0-9]+)-([0-9]+)\(([+-])\)$/, m); $2 = toupper($2); print m[1], m[2], m[3], m[4], $2}' splice_site_sequences > splice_site_sequences_fix
            out_altSS="${ss_file%%.*}.${cls}.altSS"
            printf "chr\tchrStart\tchrEnd\tunique_name\tform1\tstrand\tsplitReads\tnonSplitReads\treadSum\tsplice_eff\n" > "$out_altSS"
            out_noncanonSS="${ss_file%%.*}.${cls}.noncanonSS"            
            printf "chr\tchrStart\tchrEnd\tunique_name\tform1\tstrand\tsplitReads\tnonSplitReads\treadSum\tsplice_eff\n" > "$out_noncanonSS"


            # Revert coordinates to the original ones (3-bp)
            awk -v cls="$cls" -v motif="$motif" '
                BEGIN { FS=OFS="\t" }
                NR==FNR { seq[$1 FS $2 FS $3 FS $4]=$5; next }
                {
                    key=$1 FS $2 FS $3 FS $6
                    canon = (cls=="ss3" ? substr(seq[key],1,2)==toupper(motif) : substr(seq[key],length(seq[key])-1,2)==toupper(motif))
                    if (cls=="ss3") {
                        if ($6~/\+/) { $2 += 2              }
                        else         { $2 += 1; $3 -= 1     }
                    } else {          # ss5
                        if ($6~/\+/) { $2 += 1; $3 -= 1     }
                        else         { $2 += 2              }
                    }
                    dest = canon ? "'"$out_altSS"'" : "'"$out_noncanonSS"'"
                    print >> dest
                }
            ' splice_site_sequences_fix fasta_fix.bed


            base_name=$(basename "${ss_file%%.*}")
            mv splice_site_sequences_fix "$cutoff_folder/fasta_seq/${base_name}.${cls}_fix.fasta"
            mv splice_site_sequences "$cutoff_folder/fasta_seq/${base_name}.${cls}.fasta"
            rm fasta_fix.bed #splice_site_sequences
        done

        mv "$cutoff_folder"/*.noncanonSS "$cutoff_folder/noncanonSS/"
        mv "$cutoff_folder"/*.altSS "$cutoff_folder/altSS/"
    done
done