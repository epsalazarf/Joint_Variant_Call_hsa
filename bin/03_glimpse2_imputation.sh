#!/usr/bin/env bash

# =============================================================================
# Title       : GLIMPSE2 Imputation [FENIX]
# Description : Per-sample low-coverage WGS imputation and phasing. Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-03-27
# Version     : 1.0
# Usage       : 03_glimpse2_imputation.sh <bqsr.bam> [output_path]
# Source      : GLIMPSE2 — https://odelaneau.github.io/GLIMPSE/
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

BAM_FILE="${1:?Usage: $(basename "$0") <bqsr.bam> [output_path]}"
OUTPUT_PATH="${2:-$PWD}"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GLIMPSE2 Imputation [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking input files..."

if [ -f "$BAM_FILE" ]; then
  echo "[<]  $BAM_FILE"
  echo "[i]  Output: ${OUTPUT_PATH}"
  BAM_name=$(basename "$BAM_FILE")
  SAMPLE_NAME="${BAM_name%.rmdup.mqfilt.bqsr.bam}"
  SAMPLE_NAME="${SAMPLE_NAME%.bam}"
  FINAL_FILE="${OUTPUT_PATH}/${SAMPLE_NAME}.imputed.canon_chr.vcf.gz"
else
  echo "[X]  CANCELLED. File not found: ${BAM_FILE}"
  exit 1
fi

# Options
njobs=4
HOUSEKEEP=false

# Canonical chromosomes
CANON_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

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
  module load glimpse2 2>/dev/null \
    || echo "[W]  glimpse2 module unavailable — checking PATH"
  module load bcftools
  module load samtools
fi

# Guard: GLIMPSE2 binaries in PATH
for bin in GLIMPSE2_chunk GLIMPSE2_split_reference GLIMPSE2_phase GLIMPSE2_ligate; do
  command -v "$bin" > /dev/null 2>&1 \
    || { echo "[X]  $bin not found in PATH"; exit 1; }
done

# Guard: reference paths present in config
if [ -z "${ref_panel:-}" ] || [ -z "${ref_gmap:-}" ]; then
  echo "[X]  Missing required reference paths (ref_panel, ref_gmap). Check config: $CONFIG_FILE"
  exit 1
fi

# Guard: reference directory availability
echo "[i]  References:"
[ -d "${ref_panel}" ] || { echo "[X]  Reference panel directory not found: ${ref_panel}"; exit 1; }
echo "[i]    Panel   : ${ref_panel}"
[ -d "${ref_gmap}" ]  || { echo "[X]  Genetic map directory not found: ${ref_gmap}"; exit 1; }
echo "[i]    Gmap    : ${ref_gmap}"

# Derived paths
CACHE_DIR="${ref_panel}/../glimpse2_cache"
CHUNKS_WORKDIR="${OUTPUT_PATH}/imputed_chunks"
CHROM_DIR="${OUTPUT_PATH}/chrom_gvcf"

mkdir -p "$CACHE_DIR" "$CHUNKS_WORKDIR" "$CHROM_DIR"

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## Step 1: GLIMPSE2_chunk — define imputation windows per chromosome (cached)
step1_chunk_chromosomes() {
  local step_name="Step 1: Chunk Chromosomes"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  local n_done=0 n_skip=0

  for CHR in $CANON_CHROMS; do
    local ref_panel_chr="${ref_panel}/reference_panel.${CHR}.bcf"
    local gmap_chr="${ref_gmap}/${CHR}.b38.gmap.gz"
    local chunks_file="${CACHE_DIR}/chunks.${CHR}.txt"

    if [ ! -f "$ref_panel_chr" ]; then
      echo "[i]  Skipping ${CHR}: reference panel not found ($(basename "$ref_panel_chr"))"
      continue
    fi
    if [ ! -f "$gmap_chr" ]; then
      echo "[i]  Skipping ${CHR}: genetic map not found ($(basename "$gmap_chr"))"
      continue
    fi

    if [ -s "$chunks_file" ]; then
      echo "[i]  Already chunked: ${CHR}"
      n_skip=$(( n_skip + 1 ))
      continue
    fi

    echo "[&]  Chunking: ${CHR}"
    set -o xtrace
    GLIMPSE2_chunk \
      --input "$ref_panel_chr" \
      --region "$CHR" \
      --map "$gmap_chr" \
      --sequential \
      --output "${chunks_file}.tmp"
    set +o xtrace

    [ -s "${chunks_file}.tmp" ] || { echo "[X]  CANCELLED: $step_name failed for ${CHR}"; exit 1; }
    mv "${chunks_file}.tmp" "$chunks_file"
    echo "[>]  $chunks_file"
    n_done=$(( n_done + 1 ))
  done

  echo "[i]  Chromosomes chunked: $n_done (cached: $n_skip)"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 2: GLIMPSE2_split_reference — build binary reference panels per chunk (cached)
step2_split_reference() {
  local step_name="Step 2: Split Reference Panels"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  local n_done=0 n_skip=0

  for CHR in $CANON_CHROMS; do
    local ref_panel_chr="${ref_panel}/reference_panel.${CHR}.bcf"
    local gmap_chr="${ref_gmap}/${CHR}.b38.gmap.gz"
    local chunks_file="${CACHE_DIR}/chunks.${CHR}.txt"

    [ -f "$chunks_file" ] || continue

    while IFS=$'\t' read -r idx chr input_region output_region; do
      local bin_prefix="${CACHE_DIR}/ref_panel.${CHR}.chunk${idx}"
      local bin_file="${bin_prefix}.bin"

      if [ -s "$bin_file" ]; then
        n_skip=$(( n_skip + 1 ))
        continue
      fi

      echo "[&]  Splitting reference: ${CHR} chunk ${idx}"
      set -o xtrace
      GLIMPSE2_split_reference \
        --reference "$ref_panel_chr" \
        --map "$gmap_chr" \
        --input-region "$input_region" \
        --output-region "$output_region" \
        --output "$bin_prefix" \
        --threads "$njobs"
      set +o xtrace

      [ -s "$bin_file" ] || { echo "[X]  CANCELLED: $step_name failed for ${CHR} chunk ${idx}"; exit 1; }
      echo "[>]  $bin_file"
      n_done=$(( n_done + 1 ))
    done < "$chunks_file"
  done

  echo "[i]  Reference chunks built: $n_done (cached: $n_skip)"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 3: GLIMPSE2_phase — impute and phase per chunk per sample
step3_phase_impute() {
  local step_name="Step 3: Phase and Impute"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  local n_done=0 n_skip=0

  for CHR in $CANON_CHROMS; do
    local chunks_file="${CACHE_DIR}/chunks.${CHR}.txt"
    [ -f "$chunks_file" ] || continue

    # Skip whole chromosome if step 4 already completed
    local ligated="${CHROM_DIR}/${SAMPLE_NAME}.imputed.${CHR}.vcf.gz"
    if [ -s "$ligated" ]; then
      echo "[i]  Already ligated: ${CHR} — skipping phase"
      continue
    fi

    echo "[&]  Phasing: ${CHR}"

    while IFS=$'\t' read -r idx chr input_region output_region; do
      local bin_file="${CACHE_DIR}/ref_panel.${CHR}.chunk${idx}.bin"
      local outfile="${CHUNKS_WORKDIR}/${SAMPLE_NAME}.imputed.${CHR}.chunk${idx}.bcf"

      [ -f "$bin_file" ] || { echo "[X]  Missing binary reference: $bin_file"; exit 1; }

      if [ -s "$outfile" ]; then
        n_skip=$(( n_skip + 1 ))
        continue
      fi

      local outfile_tmp="${CHUNKS_WORKDIR}/${SAMPLE_NAME}.imputed.${CHR}.chunk${idx}.tmp.bcf"
      set -o xtrace
      GLIMPSE2_phase \
        --bam-file "$BAM_FILE" \
        --reference "$bin_file" \
        --output "$outfile_tmp" \
        --threads "$njobs"
      set +o xtrace

      [ -s "$outfile_tmp" ] || { echo "[X]  CANCELLED: $step_name failed for ${CHR} chunk ${idx}"; exit 1; }
      mv "$outfile_tmp" "$outfile"
      echo "[>]  $outfile"
      n_done=$(( n_done + 1 ))
    done < "$chunks_file"
  done

  echo "[i]  Chunks phased: $n_done (skipped: $n_skip)"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 4: GLIMPSE2_ligate — merge chunks into per-chromosome VCF.gz per sample
step4_ligate_chromosomes() {
  local step_name="Step 4: Ligate Chromosomes"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  local n_done=0 n_skip=0

  for CHR in $CANON_CHROMS; do
    local chunks_file="${CACHE_DIR}/chunks.${CHR}.txt"
    [ -f "$chunks_file" ] || continue

    local outfile="${CHROM_DIR}/${SAMPLE_NAME}.imputed.${CHR}.vcf.gz"

    if [ -s "$outfile" ]; then
      echo "[i]  Already completed: ${CHR}"
      n_skip=$(( n_skip + 1 ))
      continue
    fi

    # Build ordered list of chunk BCFs
    local list_file="${CHUNKS_WORKDIR}/${SAMPLE_NAME}.chunks.${CHR}.txt"
    : > "$list_file"
    while IFS=$'\t' read -r idx chr input_region output_region; do
      local chunk_bcf="${CHUNKS_WORKDIR}/${SAMPLE_NAME}.imputed.${CHR}.chunk${idx}.bcf"
      [ -s "$chunk_bcf" ] || { echo "[X]  Missing chunk BCF for ligate: $chunk_bcf"; exit 1; }
      echo "$chunk_bcf" >> "$list_file"
    done < "$chunks_file"

    echo "[&]  Ligating: ${CHR}"
    local ligated_bcf="${CHUNKS_WORKDIR}/${SAMPLE_NAME}.ligated.${CHR}.tmp.bcf"

    set -o xtrace
    GLIMPSE2_ligate \
      --input "$list_file" \
      --output "$ligated_bcf" \
      --threads "$njobs"

    bcftools view "$ligated_bcf" \
      --threads "$njobs" \
      -Oz -o "${outfile}.tmp"
    set +o xtrace

    rm -f "$ligated_bcf"
    mv "${outfile}.tmp" "$outfile"
    bcftools index -t "$outfile"

    [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed for ${CHR}, output missing: $outfile"; exit 1; }
    echo "[>]  $outfile"
    n_done=$(( n_done + 1 ))
  done

  echo "[i]  Chromosomes ligated: $n_done (skipped: $n_skip)"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 5: Collect per-chromosome VCFs into single whole-genome VCF
step5_collect_chromosomes() {
  local step_name="Step 5: Collect Chromosomes"
  local step_timestamp=$(date +%s)
  local outfile="$FINAL_FILE"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -s "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; return 0; }

  local chr_vcfs=()
  for CHR in $CANON_CHROMS; do
    local chr_vcf="${CHROM_DIR}/${SAMPLE_NAME}.imputed.${CHR}.vcf.gz"
    [ -s "$chr_vcf" ] && chr_vcfs+=( "$chr_vcf" )
  done

  if [ ${#chr_vcfs[@]} -eq 0 ]; then
    echo "[X]  No per-chromosome VCFs found in: $CHROM_DIR"
    exit 1
  fi

  echo "[i]  Merging ${#chr_vcfs[@]} chromosomes..."

  set -o xtrace
  bcftools concat \
    "${chr_vcfs[@]}" \
    --threads "$njobs" \
    -Oz -o "${outfile}.tmp"
  set +o xtrace

  mv "${outfile}.tmp" "$outfile"
  bcftools index -t "$outfile"

  [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
  echo "[>]  $outfile"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Housekeeping: Remove intermediate chunk BCFs (optional)
housekeeping() {
  echo
  echo "[*]  Housekeeping: Remove intermediate files"

  [ "$HOUSEKEEP" = true ] || { echo "[i]  Skipped (HOUSEKEEP=false)"; return 0; }

  if [ -s "$FINAL_FILE" ]; then
    find "$CHUNKS_WORKDIR" \
      -name "${SAMPLE_NAME}.imputed.*.chunk*.bcf" \
      -delete -print \
      | sed 's/^/[W]  Removed: /'
    echo "[W]  Intermediate chunk BCFs removed."
  else
    echo "[X]  Final file not found. Intermediate files not removed."
  fi
}

## Finisher: check final output and report
finisher() {
  if [ -s "$FINAL_FILE" ]; then
    echo
    echo "[>]  $FINAL_FILE"
    echo "[$] GLIMPSE2 Imputation [FENIX] completed successfully!"
    local elapsed=$(( EPOCHSECONDS - script_timestamp ))
    printf "[&]  Total time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
    exit 0
  else
    echo
    echo "[X]  GLIMPSE2 Imputation [FENIX] INCOMPLETE. Final output not found: ${FINAL_FILE}"
    local elapsed=$(( EPOCHSECONDS - script_timestamp ))
    printf "[&]  Total time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
    exit 1
  fi
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

main() {
  step1_chunk_chromosomes
  step2_split_reference
  step3_phase_impute
  step4_ligate_chromosomes
  step5_collect_chromosomes
  housekeeping
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
