#!/usr/bin/env bash

# =============================================================================
# Title       : GATK HaplotypeCaller [FENIX]
# Description : Per-sample GVCF generation and chromosome splitting. Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-03-24
# Version     : 1.1
# Usage       : 03_gatk_haplotype_caller.sh <bqsr.bam> [output_path]
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

BAM_FILE="${1:?Usage: $(basename "$0") <bqsr.bam> [output_path]}"
OUTPUT_PATH="${2:-$PWD}"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GATK HaplotypeCaller [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking input files..."

if [ -f "$BAM_FILE" ]; then
  echo "[<]  $BAM_FILE"
  echo "[i]  Output: ${OUTPUT_PATH}"
  BAM_name=$(basename "$BAM_FILE")
  BAM_prefix=${BAM_FILE%%.*bam}
  FINAL_FILE="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.canon_chr.g.vcf.gz"
else
  echo "[X]  CANCELLED. File not found: ${BAM_FILE}"
  exit 1
fi

# Options
njobs=2
SPLIT_VAR_TYPE=true
HOUSEKEEP=false

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
  MEM="20G"
else
  env_type="local"
  MEM="12G"
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
  module unload oracle-java
  module load oracle-java/25.0.2
  module load gatk
  module load samtools
  module load bcftools
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

## Step 1: HaplotypeCaller (GVCF mode)
step1_run_haplotype_caller() {
  local step_name="Step 1: HaplotypeCaller"
  local infile="$BAM_FILE"
  local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  gatk --java-options "-Xms$MEM -Xmx$MEM -XX:ParallelGCThreads=2" HaplotypeCaller \
    --input "$infile" \
    --reference "$ref_gnm" \
    --dbsnp "$ref_vars" \
    --emit-ref-confidence GVCF \
    --verbosity ERROR \
    --create-output-variant-index \
    --output "${outfile}"

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Step 2: Index Raw GVCF
step2_index_raw_gvcf() {
  local step_name="Step 2: Index GVCF (tbi)"
  local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"
  local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz.tbi"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  bcftools index "$infile" \
    --tbi \
    --threads "$njobs" \
    --force

  [ -s "$outfile" ] && { echo "[>]  $outfile"; echo "[!]  $step_name"; }
}

## Step 3: Extract Canonical Chromosomes (chr1-22, X, Y, M)
step3_extract_canon_chroms() {
  local step_name="Step 3: Extract Canonical Chromosomes"
  local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"
  local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.canon_chr.g.vcf.gz"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  set -o xtrace
  bcftools view "$infile" \
    --regions "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM" \
    --threads "$njobs" \
    --write-index=tbi \
    --output-type b \
    --output "${outfile}"
  set +o xtrace

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
}

## Step 4: Split GVCFs by Chromosome
step4_split_chroms_gvcf() {
  local step_name="Step 4: Split Chromosomes GVCFs"
  local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.canon_chr.g.vcf.gz"
  local outdir="${OUTPUT_PATH}/chrom_gvcf"
  local outfile="${outdir}/${BAM_prefix}.raw_vars"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ ! -f "$infile" ]]; then
    echo "[X]  Missing input: $infile" >&2
    return 1
  fi

  echo "[i]  Creating chromosome directory..."
  mkdir -pv "$outdir" || { echo "[X]  Cannot create output directory: $outdir" >&2; return 1; }

  echo "[&]  Splitting ${infile}..."

  local chroms
  if ! chroms=$(bcftools index -s "$infile" | cut -f1); then
    echo "[X]  Failed to extract chromosome list from index" >&2
    return 1
  fi

  if [[ -z "$chroms" ]]; then
    echo "[X]  No chromosomes found in $infile" >&2
    return 1
  fi

  local count=0
  while read -r C; do
    [[ -z "$C" ]] && continue

    echo "[&]  Extracting chromosome $C..."

    if ! bcftools view "$infile" \
      --regions "$C" \
      --threads "$njobs" \
      --write-index=tbi \
      --output-type z \
      --output "${outfile}.${C}.g.vcf.gz"
    then
      echo "[X]  bcftools view failed for chromosome: $C" >&2
      return 1
    fi

    echo "[>]  ${outfile}.${C}.g.vcf.gz"
    (( ++count ))

  done <<< "$chroms"

  if (( count == 0 )); then
    echo "[X]  No chromosome files were produced" >&2
    return 1
  fi

  echo "[!]  $step_name ($count chromosomes)"
}

## Housekeeping: Remove intermediate files (optional)
housekeeping() {
  echo
  echo "[*]  Housekeeping: Remove intermediate files"

  [ "$HOUSEKEEP" = true ] || { echo "[i]  Skipped (HOUSEKEEP=false)"; return 0; }

  if [ -f "$FINAL_FILE" ]; then
    rm -v \
      "${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz" \
      "${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz.tbi"
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
    echo "[$] GATK HaplotypeCaller [FENIX] completed successfully!"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 0
  else
    echo
    echo "[X]  GATK HaplotypeCaller [FENIX] INCOMPLETE. Final output not found: ${FINAL_FILE}"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 1
  fi
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

main() {
  #step0_map_reads_per_sample
  step1_run_haplotype_caller
  step2_index_raw_gvcf
  step3_extract_canon_chroms
  step4_split_chroms_gvcf
  #housekeeping
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
