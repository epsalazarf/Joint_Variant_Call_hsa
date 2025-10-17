#!/usr/bin/env bash
# rerun-bamqc-step3d.sh
# Rerun incomplete BQSR AnalyzeCovariates steps on existing outputs
# Assumes tables already exist in $OUTPUT_PATH

set -euo pipefail

# === CONFIGURATION ===
OUTPUT_PATH="${1:-$PWD}"   # path to processed BAM outputs


# Detect environment
if [ -n "$SSH_CLIENT" ] || [ -n "$SSH_TTY" ] || [ -n "$SSH_CONNECTION" ]; then
    module load gatk
    module load r
fi

# === FUNCTION ===
step3_bqsr() {
    local bam_base="$1"
    local table="${OUTPUT_PATH}/${bam_base}.rmdup.mqfilt.bqsr_table.txt"
    local table_recal="${OUTPUT_PATH}/${bam_base}.rmdup.mqfilt.bqsr_table_recal.txt"
    local output="${OUTPUT_PATH}/${bam_base}.rmdup.mqfilt.bqsr.bqsr-cov"

    echo; echo ">>> Processing ${bam_base} <<<"
    echo "> Step 3d: Analyze Covariates"

    # skip if already complete
    if [ -s "$output".pdf ]; then
        echo " [SKIP] Output exists: $output.pdf"
        return 0
    fi

    # sanity check
    if [[ ! -s "$table" || ! -s "$table_recal" ]]; then
        echo " [!] Missing required tables for $bam_base — skipping"
        return 1
    fi

    # run command
    echo " [RUN] gatk AnalyzeCovariates ..."
    gatk AnalyzeCovariates \
        -before "$table" \
        -after "$table_recal" \
        --verbosity ERROR \
        -csv "$output".csv \
        -plots "$output".pdf 
}

# === MAIN LOOP ===
echo "=== Rerunning incomplete BQSR AnalyzeCovariates steps ==="
find "$OUTPUT_PATH" -type f -name "*.rmdup.mqfilt.bqsr_table.txt" | while read -r table; do
    bam_base=$(basename "$table" .rmdup.mqfilt.bqsr_table.txt)
    step3_bqsr "$bam_base"
done
