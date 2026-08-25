#!/usr/bin/env bash

module add biotools/Samtools-1.19.2
module add biotools/Bedtools-2.30.0



set -euo pipefail

usage() {
    cat <<EOF
Usage:
  bash srna_5prime_by_chr_region.sh -b reads.bam -o outprefix [options]

Required:
  -b BAM              Input aligned BAM
  -o PREFIX           Output prefix

Optional:
  -s SIZES            Comma-separated read sizes to keep [default: 21,22,23,24]
  -q MAPQ             Minimum MAPQ [default: 0]
  -r BED              Optional BED file of regions
  -m MODE             Region mode: overlap or 5p [default: overlap]
  -g FILE             Optional chromosome-group mapping file
  -k                  Keep temporary files
  -h                  Show help

Notes:
  - Requires: samtools, gawk
  - If -r is used, also requires: bedtools
  - Region BED:
      * columns 1-3 required
      * column 4 optional region ID
  - Chromosome-group file:
      * two columns: chromosome<TAB>group
      * example:
            chr1    nuclear
            chr2    nuclear
            chrM    mito
            chrC    plastid
  - The script reports the ORIGINAL small RNA 5' nucleotide:
      * forward-strand read  -> first base of read
      * reverse-strand read  -> complement of last base of read
  - Output rows are sorted by:
      * Group alphabetically
      * Size numerically
EOF
}

BAM=""
OUT=""
SIZES="21,22,23,24"
MAPQ=0
REGIONS=""
MODE="overlap"
CHR_GROUPS=""
KEEP_TEMP=0

while getopts ":b:o:s:q:r:m:g:kh" opt; do
    case "$opt" in
        b) BAM="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        s) SIZES="$OPTARG" ;;
        q) MAPQ="$OPTARG" ;;
        r) REGIONS="$OPTARG" ;;
        m) MODE="$OPTARG" ;;
        g) CHR_GROUPS="$OPTARG" ;;
        k) KEEP_TEMP=1 ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "Error: option -$OPTARG requires an argument" >&2
            usage
            exit 1
            ;;
        \?)
            echo "Error: invalid option -$OPTARG" >&2
            usage
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))

if [[ -z "$BAM" || -z "$OUT" ]]; then
    echo "Error: -b and -o are required" >&2
    usage
    exit 1
fi

if [[ ! -f "$BAM" ]]; then
    echo "Error: BAM not found: $BAM" >&2
    exit 1
fi

if [[ "$MODE" != "overlap" && "$MODE" != "5p" ]]; then
    echo "Error: -m must be 'overlap' or '5p'" >&2
    exit 1
fi

command -v samtools >/dev/null 2>&1 || { echo "Error: samtools not found" >&2; exit 1; }
command -v gawk >/dev/null 2>&1 || { echo "Error: gawk not found" >&2; exit 1; }

if [[ -n "$REGIONS" ]]; then
    [[ -f "$REGIONS" ]] || { echo "Error: BED not found: $REGIONS" >&2; exit 1; }
    command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools not found" >&2; exit 1; }
fi

if [[ -n "$CHR_GROUPS" ]]; then
    [[ -f "$CHR_GROUPS" ]] || { echo "Error: chromosome-group file not found: $CHR_GROUPS" >&2; exit 1; }
fi

TMPDIR=$(mktemp -d "${OUT}.tmp.XXXXXX")
if [[ "$KEEP_TEMP" -eq 0 ]]; then
    trap 'rm -rf "$TMPDIR"' EXIT
fi

READS_BED="${TMPDIR}/reads.filtered.bed"
READS_5P_BED="${TMPDIR}/reads.filtered.5p.bed"
REGIONS_NAMED="${TMPDIR}/regions.named.bed"
REGIONS_INTERSECT="${TMPDIR}/regions.intersect.tsv"

CHR_COUNTS="${OUT}.chromosome_5prime_counts.tsv"
CHR_FREQ="${OUT}.chromosome_5prime_freq.tsv"
REG_COUNTS="${OUT}.region_5prime_counts.tsv"
REG_FREQ="${OUT}.region_5prime_freq.tsv"

CHR_COUNTS_TMP="${TMPDIR}/chromosome_5prime_counts.unsorted.tsv"
CHR_FREQ_TMP="${TMPDIR}/chromosome_5prime_freq.unsorted.tsv"
REG_COUNTS_TMP="${TMPDIR}/region_5prime_counts.unsorted.tsv"
REG_FREQ_TMP="${TMPDIR}/region_5prime_freq.unsorted.tsv"

echo "[INFO] Extracting mapped reads and computing 5' nucleotide..." >&2

samtools view -F 2308 -q "$MAPQ" "$BAM" | \
gawk -v sizes="$SIZES" '
BEGIN {
    OFS = "\t"
    n = split(sizes, a, ",")
    for (i = 1; i <= n; i++) keep[a[i]] = 1
}
function ref_len_from_cigar(cigar,   tok, num, op, len) {
    len = 0
    while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
        tok = substr(cigar, RSTART, RLENGTH)
        num = tok
        sub(/[A-Z=]$/, "", num)
        op = substr(tok, length(tok), 1)
        if (op ~ /[MDN=X]/) len += num
        cigar = substr(cigar, RSTART + RLENGTH)
    }
    return len
}
function comp(base) {
    base = toupper(base)
    if (base == "A") return "T"
    if (base == "C") return "G"
    if (base == "G") return "C"
    if (base == "T") return "A"
    return "N"
}
{
    flag = $2
    chrom = $3
    pos1 = $4
    cigar = $6
    seq = toupper($10)
    qname = $1

    if (seq == "*" || cigar == "*") next

    read_len = length(seq)
    if (!(read_len in keep)) next

    start0 = pos1 - 1
    span = ref_len_from_cigar(cigar)
    end0 = start0 + span

    rev = and(flag, 16)
    if (rev) {
        strand = "-"
        nt5 = comp(substr(seq, read_len, 1))
        fivep = end0 - 1
    } else {
        strand = "+"
        nt5 = substr(seq, 1, 1)
        fivep = start0
    }

    nt5 = toupper(nt5)
    if (nt5 !~ /^[ACGTN]$/) nt5 = "N"

    print chrom, start0, end0, qname, read_len, strand, nt5, fivep
}
' > "$READS_BED"

if [[ ! -s "$READS_BED" ]]; then
    echo "Error: no reads passed filters" >&2
    exit 1
fi

echo "[INFO] Writing chromosome/group summaries..." >&2

gawk -v mapfile="$CHR_GROUPS" -v counts_out="$CHR_COUNTS_TMP" -v freq_out="$CHR_FREQ_TMP" '
BEGIN {
    OFS = "\t"
    if (mapfile != "") {
        while ((getline < mapfile) > 0) {
            if ($0 ~ /^#/ || NF < 2) continue
            chr_map[$1] = $2
        }
        close(mapfile)
    }
}
{
    chr  = $1
    size = $5
    nt   = $7

    grp = (chr in chr_map ? chr_map[chr] : chr)

    key = grp SUBSEP size
    counts[key, nt]++
    total[key]++
    seen[key] = 1
}
END {
    print "Group","Size","A","C","G","T","N","Total" > counts_out
    print "Group","Size","A","C","G","T","N","Total" > freq_out

    for (key in seen) {
        split(key, x, SUBSEP)
        grp = x[1]
        sz  = x[2]

        A = counts[key, "A"] + 0
        C = counts[key, "C"] + 0
        G = counts[key, "G"] + 0
        T = counts[key, "T"] + 0
        N = counts[key, "N"] + 0
        Tot = total[key] + 0

        print grp, sz, A, C, G, T, N, Tot >> counts_out

        if (Tot > 0) {
            print grp, sz, A/Tot, C/Tot, G/Tot, T/Tot, N/Tot, Tot >> freq_out
        }
    }
}
' "$READS_BED"

{
    head -n 1 "$CHR_COUNTS_TMP"
    tail -n +2 "$CHR_COUNTS_TMP" | sort -t $'\t' -k1,1 -k2,2n
} > "$CHR_COUNTS"

{
    head -n 1 "$CHR_FREQ_TMP"
    tail -n +2 "$CHR_FREQ_TMP" | sort -t $'\t' -k1,1 -k2,2n
} > "$CHR_FREQ"

if [[ -n "$REGIONS" ]]; then
    echo "[INFO] Preparing region BED..." >&2

    gawk '
    BEGIN { OFS = "\t" }
    {
        if ($0 ~ /^#/ || NF < 3) next
        id = (NF >= 4 && $4 != "" ? $4 : $1 ":" $2 "-" $3)
        print $1, $2, $3, id
    }
    ' "$REGIONS" > "$REGIONS_NAMED"

    if [[ "$MODE" == "overlap" ]]; then
        echo "[INFO] Intersecting reads with regions using overlap mode..." >&2
        bedtools intersect -wa -wb -a "$READS_BED" -b "$REGIONS_NAMED" > "$REGIONS_INTERSECT"
    else
        echo "[INFO] Intersecting reads with regions using 5p mode..." >&2
        gawk 'BEGIN{OFS="\t"} {print $1, $8, $8+1, $4, $5, $6, $7, $8}' "$READS_BED" > "$READS_5P_BED"
        bedtools intersect -wa -wb -a "$READS_5P_BED" -b "$REGIONS_NAMED" > "$REGIONS_INTERSECT"
    fi

    echo "[INFO] Writing region summaries..." >&2

    gawk -v counts_out="$REG_COUNTS_TMP" -v freq_out="$REG_FREQ_TMP" '
    BEGIN {
        OFS = "\t"
    }
    {
        region = $12
        size   = $5
        nt     = $7

        key = region SUBSEP size
        counts[key, nt]++
        total[key]++
        seen[key] = 1
    }
    END {
        print "Group","Size","A","C","G","T","N","Total" > counts_out
        print "Group","Size","A","C","G","T","N","Total" > freq_out

        for (key in seen) {
            split(key, x, SUBSEP)
            grp = x[1]
            sz  = x[2]

            A = counts[key, "A"] + 0
            C = counts[key, "C"] + 0
            G = counts[key, "G"] + 0
            T = counts[key, "T"] + 0
            N = counts[key, "N"] + 0
            Tot = total[key] + 0

            print grp, sz, A, C, G, T, N, Tot >> counts_out

            if (Tot > 0) {
                print grp, sz, A/Tot, C/Tot, G/Tot, T/Tot, N/Tot, Tot >> freq_out
            }
        }
    }
    ' "$REGIONS_INTERSECT"

    {
        head -n 1 "$REG_COUNTS_TMP"
        tail -n +2 "$REG_COUNTS_TMP" | sort -t $'\t' -k1,1 -k2,2n
    } > "$REG_COUNTS"

    {
        head -n 1 "$REG_FREQ_TMP"
        tail -n +2 "$REG_FREQ_TMP" | sort -t $'\t' -k1,1 -k2,2n
    } > "$REG_FREQ"
fi

echo "[INFO] Done." >&2
echo "[INFO] Outputs:" >&2
echo "  $CHR_COUNTS" >&2
echo "  $CHR_FREQ" >&2
if [[ -n "$REGIONS" ]]; then
    echo "  $REG_COUNTS" >&2
    echo "  $REG_FREQ" >&2
fi