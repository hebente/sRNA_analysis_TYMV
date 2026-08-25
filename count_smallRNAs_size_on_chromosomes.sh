#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 input.sam|input.bam long_output.tsv matrix_output.tsv"
    exit 1
fi

infile="$1"
out_long="$2"
out_matrix="$3"

if [[ ! -f "$infile" ]]; then
    echo "Error: input file not found: $infile" >&2
    exit 1
fi

# Decide how to read the file
if [[ "$infile" =~ \.bam$ ]]; then
    if ! command -v samtools >/dev/null 2>&1; then
        echo "Error: samtools is required to read BAM files." >&2
        exit 1
    fi
    input_cmd=(samtools view "$infile")
elif [[ "$infile" =~ \.sam$ ]]; then
    input_cmd=(cat "$infile")
else
    echo "Error: input must end with .sam or .bam" >&2
    exit 1
fi

"${input_cmd[@]}" | awk -v out_long="$out_long" -v out_matrix="$out_matrix" '
BEGIN {
    OFS = "\t"
}
{
    flag = $2
    chr  = $3
    seq  = $10

    # skip unmapped
    if (and(flag, 4)) next
    if (chr == "*" || chr == "") next

    len = length(seq)

    if (len == 21 || len == 22 || len == 23 || len == 24) {
        count[chr, len]++
        chr_seen[chr] = 1
    }
}
END {
    # long format
    print "chromosome", "read_length", "read_count" > out_long
    for (k in count) {
        split(k, a, SUBSEP)
        print a[1], a[2], count[k] >> out_long
    }
    close(out_long)

    # collect chromosome names
    nchr = 0
    for (c in chr_seen) {
        chr_list[++nchr] = c
    }

    # simple sort
    for (i = 1; i <= nchr; i++) {
        for (j = i + 1; j <= nchr; j++) {
            if (chr_list[i] > chr_list[j]) {
                tmp = chr_list[i]
                chr_list[i] = chr_list[j]
                chr_list[j] = tmp
            }
        }
    }

    # matrix header
    printf "len" > out_matrix
    for (i = 1; i <= nchr; i++) {
        printf "\t%s", chr_list[i] >> out_matrix
    }
    printf "\n" >> out_matrix

    # rows 21..24
    for (l = 21; l <= 24; l++) {
        printf "%d", l >> out_matrix
        for (i = 1; i <= nchr; i++) {
            c = chr_list[i]
            val = count[c, l]
            if (val == "") val = 0
            printf "\t%d", val >> out_matrix
        }
        printf "\n" >> out_matrix
    }
    close(out_matrix)
}
'

echo "Done:"
echo "  Long table:   $out_long"
echo "  Matrix table: $out_matrix"