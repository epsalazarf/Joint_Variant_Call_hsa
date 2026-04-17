#!/usr/bin/env bash

# =============================================================================
# Title       : GATK4 BQSR Evaluate [FENIX]
# Description : Post-recalibration evaluation for analysis-ready BAMs. Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-04-17
# Version     : 1.0
# Usage       : 04_bqsr_evaluate.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

BAM_FILE="${1:?Usage: $(basename "$0") <sample.rmdup.mqfilt.bqsr.bam> [output_path]}"
OUTPUT_PATH="${2:-$PWD}"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GATK4 BQSR Evaluate [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking input files..."

if [ -f "$BAM_FILE" ]; then
  echo "[<]  $BAM_FILE"
  echo "[i]  Output: ${OUTPUT_PATH}"
  BAM_prefix=${BAM_FILE%%.*bam}
  FINAL_FILE="${OUTPUT_PATH}/${BAM_prefix}.bqsr_table_recal.txt"
else
  echo "[X]  CANCELLED. File not found: ${BAM_FILE}"
  exit 1
fi

# Options
BQSR_COV=true   # Run Step 2: AnalyzeCovariates (only reason to run this script)

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
else
  env_type="local"
fi

echo "[i]  Environment: $env_type"

# Parse YAML config into Bash variables
eval "$(
  awk -v env="$env_type" '
    BEGIN { in_env=0 }
    $1 ~ env":" { in_env=1; next }
    in_env && /^[^[:space:]]/ { in_env=0 }
    in_env && /^[[:space:]]+[a-zA-Z0-9_]+:/ {
      gsub(":", "=", $1)
      sub(/^[[:space:]]+/, "", $1)
      gsub(/^"/, "", $2); gsub(/"$/, "", $2)
      print $1 $2
    }
  ' "$CONFIG_FILE"
)"

# Load modules on remote (work-around due to faulty parser [ARC02])
if [[ "$env_type" == "remote" ]]; then
  echo "[i]  Loading modules..."
  module load gatk
  module load r
fi

# Guard: reference file paths
if [ -z "$ref_gnm" ] || [ -z "$ref_vars" ]; then
  echo "[X]  Missing required reference paths. Check config: $CONFIG_FILE"
  exit 1
fi

# Guard: reference file availability
echo "[i]  References:"
[ -f "${ref_gnm}" ]  || { echo "[X]  Reference genome not found: ${ref_gnm}"; exit 1; }
echo "[i]    Genome  : ${ref_gnm}"
[ -f "${ref_vars}" ] || { echo "[X]  Reference variants not found: ${ref_vars}"; exit 1; }
echo "[i]    Variants: ${ref_vars}"

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## Step 1: Post-BQSR BaseRecalibrator
step1_post_bqsr_baserecalibrator() {
  local step_name="Step 1: Post-BQSR BaseRecalibrator"
  local infile="$BAM_FILE"
  local outfile="${OUTPUT_PATH}/${BAM_prefix}.bqsr_table_recal.txt"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  gatk BaseRecalibrator \
    --input "$infile" \
    --reference "$ref_gnm" \
    --known-sites "$ref_vars" \
    --verbosity ERROR \
    --native-pair-hmm-threads 4 \
    --output "${outfile}.tmp"

  mv "${outfile}.tmp" "$outfile"

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 2: AnalyzeCovariates (optional)
step2_analyze_covariates() {
  local step_name="Step 2: Analyze Covariates"
  local table="${OUTPUT_PATH}/${BAM_prefix}.bqsr_table.txt"
  local table_recal="${OUTPUT_PATH}/${BAM_prefix}.bqsr_table_recal.txt"
  local cov_report="${OUTPUT_PATH}/${BAM_prefix}.bqsr_covariates"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ "$BQSR_COV" = true ] || { echo "[i]  Skipped (BQSR_COV=false)"; return 0; }

  [ -f "$table" ]       || { echo "[X]  Missing pre-BQSR table: $table"; exit 1; }
  [ -f "$table_recal" ] || { echo "[X]  Missing post-BQSR table: $table_recal"; exit 1; }
  [ -s "${cov_report}.pdf" ] && { echo "[i]  Already completed (${cov_report}.pdf exists)"; return 0; }

  gatk AnalyzeCovariates \
    --before-report-file "$table" \
    --after-report-file "$table_recal" \
    --verbosity ERROR \
    --intermediate-csv-file "${cov_report}.csv" \
    --plots-report-file "${cov_report}.pdf"

  [ -s "${cov_report}.pdf" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: ${cov_report}.pdf"; exit 1; }
  echo "[>]  ${cov_report}.pdf"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Finisher: check final output and report
finisher() {
  if [ -f "$FINAL_FILE" ]; then
    echo
    echo "[>]  $FINAL_FILE"
    echo "[$]  GATK4 BQSR Evaluate [FENIX] completed successfully!"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 0
  else
    echo
    echo "[X]  GATK4 BQSR Evaluate [FENIX] INCOMPLETE. Final output not found: ${FINAL_FILE}"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 1
  fi
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

main() {
  step1_post_bqsr_baserecalibrator
  step2_analyze_covariates
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
