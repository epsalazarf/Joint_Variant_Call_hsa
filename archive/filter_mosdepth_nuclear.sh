#!/bin/bash
# ============================================================
# filter_mosdepth_nuclear.sh
# Filters a mosdepth summary file to only include canonical
# human chromosomes (chr1–22, chrX, chrY), recalculates totals,
# and outputs a subset table like the input.
#
# Usage:
#   ./filter_mosdepth_nuclear.sh input.mosdepth.summary.txt > output.filtered.txt
# ============================================================

# Check input
if [ $# -lt 1 ]; then
    echo "Usage: $0 <mosdepth_summary_file>" >&2
    exit 1
fi

FILE="$1"

awk '
BEGIN {
    FS=OFS="\t";
    print "chrom", "length", "bases", "mean", "min", "max"
}
NR==1 && $1 == "chrom" { next }  # Skip header
# Match canonical chromosomes: chr1–22, chrX, chrY only
/^chr([1-9]|1[0-9]|2[0-2])\t/ {
    print $0;
    total_length += $2;
    total_bases  += $3;
    if (min_total == "") min_total = $5;
    if ($5 < min_total) min_total = $5;
    if ($6 > max_total) max_total = $6;
}
END {
    if (total_length > 0) {
        mean_total = total_bases / total_length;
        printf "%s\t%d\t%d\t%.2f\t%d\t%d\n",
            "total", total_length, total_bases, mean_total, min_total, max_total;
    }
}' "$FILE"
