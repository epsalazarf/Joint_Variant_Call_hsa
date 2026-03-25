#!/usr/bin/env bash

# =============================================================================
# Title       : GATK4 BAM QC [FENIX]
# Description : Quality control for mapped BAM files. Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-03-25
# Version     : 1.2
# Usage       : 02_gatk_bam_qc_workflow.sh <bam_files> [output_path] [rg_string] [add_rg]
#             : bam_files — comma-separated BAM path(s); quote when passing multiple
#             : rg_string — optional legacy @RG string; triggers step 0 when add_rg=auto
#             : add_rg    — auto (default) | true | false  (controls step 0 execution)
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

BAM_ARG="${1:?Usage: $(basename "$0") <bam_files> [output_path] [rg_string] [add_rg]}"
OUTPUT_PATH="${2:-$PWD}"
rg_string="${3:-}"
ADD_RG="${4:-auto}"   # auto: run step 0 only if rg_string provided | true: always | false: never

# Parse comma-separated BAM list into array
IFS=',' read -r -a BAM_INPUTS <<< "$BAM_ARG"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GATK4 BAM QC [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking input files..."

# Validate all input BAMs
for bam in "${BAM_INPUTS[@]}"; do
  if [ -f "$bam" ]; then
    echo "[<]  $bam"
  else
    echo "[X]  CANCELLED. File not found: ${bam}"
    exit 1
  fi
done
echo "[i]  Output: ${OUTPUT_PATH}"

# Derive sample name from first BAM filename (basename only; fixes path-prefix bug)
first_bam=$(basename "${BAM_INPUTS[0]}")
BAM_sample="${first_bam%.sort.bam}"
BAM_sample="${BAM_sample%.bam}"
# For multiple BAMs (multiplexed): strip _rg suffix — SAMPLE_NAME must not contain underscores
[[ ${#BAM_INPUTS[@]} -gt 1 ]] && BAM_sample="${BAM_sample%%_*}"
echo "[i]  Sample: ${BAM_sample}"

FINAL_FILE="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr.bam"

# Omni-skip if final output already exists
[ -s "$FINAL_FILE" ] && { echo "[i]  Already completed ($FINAL_FILE). Skipping."; exit 0; }

# Determine step 0 execution based on ADD_RG flag
run_step0=false
if [[ "$ADD_RG" == "true" ]]; then
  run_step0=true
elif [[ "$ADD_RG" == "auto" && -n "$rg_string" ]]; then
  run_step0=true
fi

# Guard: AddOrReplaceReadGroups accepts only one input BAM
if [[ "$run_step0" == "true" && ${#BAM_INPUTS[@]} -gt 1 ]]; then
  echo "[X]  CANCELLED. Step 0 (AddOrReplaceReadGroups) requires a single input BAM."
  echo "[X]  Multiple BAMs + RG correction is unsupported. Merge first, or set add_rg=false."
  exit 1
fi

# Resolve step 1 input sources (original BAMs, or RG-corrected BAM after step 0)
if [[ "$run_step0" == "true" ]]; then
  STEP1_INPUTS=("${OUTPUT_PATH}/${BAM_sample}.rg.bam")
else
  STEP1_INPUTS=("${BAM_INPUTS[@]}")
fi

# Options
njobs=4
BQSR_COV=true    # Run Step 3d: AnalyzeCovariates (auto-off if remote)
RUN_METRICS=true # Run Step 5: Alignment & Insert Size Metrics
MQ_FILTER=false  # Run Step 2: MQ filtering (enable for aDNA; may reduce BQSR model quality on low-coverage data)
REMOVE_DUPS=false # false = mark duplicates (GATK best practice); true = remove (saves storage)
HOUSEKEEP=true

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
  module load samtools
  module load mosdepth
  module load r
  BQSR_COV=false
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

## Step 0: Assign Read Groups (legacy only — skip for BAMs from 01_bwa_map_fastq_reads.sh)
step0_add_readgroup() {
  local step_name="Step 0: ReadGroup Assignation"
  local infile="${BAM_INPUTS[0]}"
  local outfile="${OUTPUT_PATH}/${BAM_sample}.rg.bam"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  if [[ "$rg_string" == @RG* ]]; then
    local rgID="$(sed -E 's/.*ID:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
    local rgPL="$(sed -E 's/.*PL:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
    local rgPU="$(sed -E 's/.*PU:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
    local rgLB="$(sed -E 's/.*LB:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
    local rgSM="$(sed -E 's/.*SM:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
  else
    local rgID="$BAM_sample"
    local rgPL="ILLUMINA"
    local rgPU="$BAM_sample"
    local rgLB="$BAM_sample"
    local rgSM="$BAM_sample"
  fi

  gatk AddOrReplaceReadGroups \
    --INPUT "$infile" \
    --RGID "$rgID" \
    --RGPL "$rgPL" \
    --RGPU "$rgPU" \
    --RGLB "$rgLB" \
    --RGSM "$rgSM" \
    --VERBOSITY ERROR \
    --OUTPUT "${outfile}.tmp"

  mv "${outfile}.tmp" "$outfile"

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 1: Mark/Remove Duplicates (Spark)
step1_mark_duplicates_spark() {
  local step_name="Step 1: $([ "$REMOVE_DUPS" = true ] && echo Remove || echo Mark) Duplicates (Spark)"
  local step_timestamp=$(date +%s)
  local outfile="${OUTPUT_PATH}/${BAM_sample}.rmdup.bam"
  local metrics="${OUTPUT_PATH}/${BAM_sample}-dups.txt"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  # Validate all inputs
  for bam in "${STEP1_INPUTS[@]}"; do
    [ -f "$bam" ] || { echo "[X]  Missing input: $bam"; exit 1; }
  done
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  # Build --input flags for each BAM (supports merge of multiple files)
  local input_flags=()
  for bam in "${STEP1_INPUTS[@]}"; do
    input_flags+=(--input "$bam")
  done

  # Mark or remove duplicates
  local dup_args=()
  [[ "$REMOVE_DUPS" == "true" ]] && dup_args=(--remove-all-duplicates)

  local shard_tmp=".${RANDOM}.parts"
  gatk --java-options "-Xmx24G" MarkDuplicatesSpark \
    "${input_flags[@]}" \
    --spark-runner LOCAL \
    --spark-master local["${njobs}"] \
    "${dup_args[@]}" \
    --verbosity ERROR \
    --output-shard-tmp-dir "$shard_tmp" \
    --metrics-file "$metrics" \
    --output "${outfile}.tmp"
    #--conf "spark.executor.cores=${njobs}" \

  rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*
  rm -rf "$shard_tmp"

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 1 (alt): Mark/Remove Duplicates (Picard)
step1_mark_duplicates_picard() {
  local step_name="Step 1: $([ "$REMOVE_DUPS" = true ] && echo Remove || echo Mark) Duplicates (Picard)"
  local step_timestamp=$(date +%s)
  local outfile="${OUTPUT_PATH}/${BAM_sample}.rmdup.bam"
  local metrics="${OUTPUT_PATH}/${BAM_sample}-dups.txt"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  # Validate all inputs
  for bam in "${STEP1_INPUTS[@]}"; do
    [ -f "$bam" ] || { echo "[X]  Missing input: $bam"; exit 1; }
  done
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  # Build --INPUT flags for each BAM (supports merge of multiple files)
  local input_flags=()
  for bam in "${STEP1_INPUTS[@]}"; do
    input_flags+=(--INPUT "$bam")
  done

  # Mark or remove duplicates
  local dup_args=()
  [[ "$REMOVE_DUPS" == "true" ]] && dup_args=(--REMOVE_DUPLICATES)

  gatk MarkDuplicates \
    "${input_flags[@]}" \
    "${dup_args[@]}" \
    --VERBOSITY ERROR \
    --METRICS_FILE "$metrics" \
    --OUTPUT "${outfile}.tmp"

  rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 2: Mapping Quality Filter
step2_mapping_quality_filter() {
  local step_name="Step 2: Mapping Quality Filter"
  local step_timestamp=$(date +%s)
  local infile="${OUTPUT_PATH}/${BAM_sample}.rmdup.bam"
  local outfile="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bam"
  local counts="${OUTPUT_PATH}/${BAM_sample}.mqfilt-counts.txt"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ "$MQ_FILTER" = true ] || { echo "[i]  Skipped (MQ_FILTER=false)"; return 0; }

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  set -o xtrace
  samtools view "$infile" \
    -F 4 -q 30 \
    --threads "$njobs" \
    --write-index \
    --save-counts "$counts" \
    --bam \
    --output "${outfile}.tmp"
  set +o xtrace

  rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 3: Base Quality Score Recalibration (BQSR)
step3_bqsr() {
  local step_name="Step 3: Base Quality Recalibration"
  local step_timestamp=$(date +%s)
  local infile
  if [ "$MQ_FILTER" = true ]; then
    infile="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bam"
  else
    infile="${OUTPUT_PATH}/${BAM_sample}.rmdup.bam"
  fi
  local table="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr_table.txt"
  local bam_out="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr.bam"
  local table_recal="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr_table_recal.txt"
  local cov_report="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr.cov"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$bam_out" ] && { echo "[i]  Already completed ($bam_out exists)"; return 0; }

  ## 3a: Build recalibration model
  echo
  echo "[*]  Step 3a: Model Building"
  if [ -s "$table" ]; then
    echo "[i]   Model table exists: $table"
  else
    gatk BaseRecalibrator \
      --input "$infile" \
      --reference "$ref_gnm" \
      --known-sites "$ref_vars" \
      --verbosity ERROR \
      --output "${table}.tmp"
    mv "${table}.tmp" "$table"
    echo "[!]   Step 3a: Model Building"
  fi

  ## 3b: Apply recalibration
  echo
  echo "[*]  Step 3b: Apply Recalibration"
  if [ -s "$bam_out" ]; then
    echo "[i]   Recalibrated BAM exists: $bam_out"
  else
    gatk ApplyBQSR \
      --input "$infile" \
      --reference "$ref_gnm" \
      --bqsr-recal-file "$table" \
      --create-output-bam-index \
      --verbosity ERROR \
      --output "${bam_out}.tmp"
    rename -v "${bam_out}.tmp" "${bam_out}" "${bam_out}.tmp"*
    echo "[!]   Step 3b: Apply Recalibration"
  fi

  echo "[&]  Step time (3a+3b): $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"

  ## 3c: Post-recalibration evaluation table
  echo
  echo "[*]  Step 3c: Evaluate Recalibration"
  if [ -s "$table_recal" ]; then
    echo "[i]   Post-recalibration table exists: $table_recal"
  else
    gatk BaseRecalibrator \
      --input "$bam_out" \
      --reference "$ref_gnm" \
      --known-sites "$ref_vars" \
      --verbosity ERROR \
      --output "${table_recal}.tmp"
    mv "${table_recal}.tmp" "$table_recal"
    echo "[!]   Step 3c: Evaluate Recalibration"
  fi

  ## 3d: Analyze Covariates (optional, local only; requires R)
  echo
  echo "[*]  Step 3d: Analyze Covariates"
  [ "$BQSR_COV" = true ] || { echo "[i]   Skipped (BQSR_COV=false)"; return 0; }

  if [ -s "${cov_report}.pdf" ]; then
    echo "[i]   Covariate plots exist: ${cov_report}.pdf"
  else
    gatk AnalyzeCovariates \
      --before-report-file "$table" \
      --after-report-file "$table_recal" \
      --verbosity ERROR \
      --intermediate-csv-file "${cov_report}.csv" \
      --plots-report-file "${cov_report}.pdf"
  fi

  echo "[>]  $bam_out"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 4: Coverage Calculation
step4_mosdepth() {
  local step_name="Step 4: Calculate Coverage (mosdepth)"
  local step_timestamp=$(date +%s)
  local infile="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr.bam"
  local prefix="${OUTPUT_PATH}/${BAM_sample}.rmb"
  local outfile="${prefix}.mosdepth.summary.txt"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"
  echo "[&]  mosdepth --fast-mode --no-per-base $prefix $infile"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  set -o xtrace
  mosdepth \
    --threads "$njobs" \
    --fast-mode \
    --no-per-base \
    --fasta "$ref_gnm" \
    "$prefix" \
    "$infile"
  set +o xtrace

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 5: Alignment & Insert Size Metrics (optional)
step5_metrics() {
  local step_name="Step 5: Collect Alignment & Insert Size Metrics"
  local step_timestamp=$(date +%s)
  local infile="${OUTPUT_PATH}/${BAM_sample}.rmdup.mqfilt.bqsr.bam"
  local align_metrics="${OUTPUT_PATH}/${BAM_sample}.alignment_metrics.txt"
  local insert_metrics="${OUTPUT_PATH}/${BAM_sample}.insert_size_metrics.txt"
  local insert_hist="${OUTPUT_PATH}/${BAM_sample}.insert_size_histogram.pdf"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ "$RUN_METRICS" = true ] || { echo "[i]  Skipped (RUN_METRICS=false)"; return 0; }

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }

  ## 5a: Alignment metrics
  echo "[*]  Step 5a: Alignment Summary Metrics"
  if [ -s "$align_metrics" ]; then
    echo "[i]   Alignment metrics exist: $align_metrics"
  else
    gatk CollectAlignmentSummaryMetrics \
      --INPUT "$infile" \
      --REFERENCE_SEQUENCE "$ref_gnm" \
      --VERBOSITY ERROR \
      --OUTPUT "${align_metrics}.tmp"
    mv "${align_metrics}.tmp" "$align_metrics"
    echo "[!]   Step 5a: Alignment Summary Metrics"
  fi

  ## 5b: Insert size metrics
  echo "[*]  Step 5b: Insert Size Metrics"
  if [ -s "$insert_metrics" ] && [ -s "$insert_hist" ]; then
    echo "[i]   Insert size metrics exist"
  else
    gatk CollectInsertSizeMetrics \
      --INPUT "$infile" \
      --VERBOSITY ERROR \
      --OUTPUT "${insert_metrics}.tmp" \
      --Histogram_FILE "${insert_hist}.tmp"
    mv "${insert_metrics}.tmp" "$insert_metrics"
    mv "${insert_hist}.tmp" "$insert_hist"
    echo "[!]   Step 5b: Insert Size Metrics"
  fi

  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Housekeeping: Remove intermediate files (optional)
housekeeping() {
  echo
  echo "[*]  Housekeeping: Remove intermediate files"

  [ "$HOUSEKEEP" = true ] || { echo "[i]  Skipped (HOUSEKEEP=false)"; return 0; }

  if [ -f "$FINAL_FILE" ]; then
    # Remove RG-corrected BAM only if step 0 ran (file may not exist for new-workflow runs)
    [ -f "${OUTPUT_PATH}/${BAM_sample}.rg.bam" ] && rm -v "${OUTPUT_PATH}/${BAM_sample}.rg.bam"
    rm -vf "${OUTPUT_PATH}/${BAM_sample}".rmdup.ba*
    [ "$MQ_FILTER" = true ] && rm -vf "${OUTPUT_PATH}/${BAM_sample}".rmdup.mqfilt.ba*
    echo "[W]  Intermediate files removed."
  else
    echo "[X]  Final file not found. Intermediate files not removed."
  fi
}

## Finisher: check final output and report
finisher() {
  if [ -f "$FINAL_FILE" ]; then
    echo
    echo "[>]  $FINAL_FILE"
    echo "[$]  GATK4 BAM QC [FENIX] completed successfully!"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 0
  else
    echo
    echo "[X]  GATK4 BAM QC [FENIX] INCOMPLETE. Final output not found: ${FINAL_FILE}"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 1
  fi
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

main() {
  [[ "$run_step0" == "true" ]] && step0_add_readgroup
  step1_mark_duplicates_spark
  #step1_mark_duplicates_picard
  step2_mapping_quality_filter
  step3_bqsr
  step4_mosdepth
  step5_metrics
  housekeeping
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
