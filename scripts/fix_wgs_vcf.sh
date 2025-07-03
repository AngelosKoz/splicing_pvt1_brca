#!/bin/bash

INPUT=$1
OUTPUT=$2

RS_COUNTER=1

cat "$INPUT" | awk -v rs_counter=1 'BEGIN {FS=OFS="\t"} NR == 1 {print "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "RS", "FORMAT", "NORMAL", "TUMOR", "file_name_wgs"; next}{
    rs = "rsNA" rs_counter
    if ($8 ~ /rs[0-9]+/) {
        match($8, /rs[0-9]+/, id); rs = id[0]
    } else {
        rs_counter++
    }
    print $1, $2, $3, $4, $5, $6, $7, rs, $9, $10, $11, $12}' > "$OUTPUT"
