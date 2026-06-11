#!/usr/bin/env bash
#SBATCH --job-name=GLIMPSE2_PREP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=END,FAIL

# =============================================================================
# Title       : GLIMPSE2 Reference Panel Prep [FENIX]
# Description : QC, chunk, and split a phased reference panel VCF into the
#               binary format required by GLIMPSE2_phase. One VCF per
#               chromosome; submit one job per chromosome in parallel.
#               Designed for LAVIS-FENIX with /scratch I/O isolation.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-06-10
# Version     : 1.2
# Usage       : sbatch bin/00_glimpse2_ref_panel_prep.sh <panel.vcf.gz> [output_panel_dir] [run_qc]
#             : bash   bin/00_glimpse2_ref_panel_prep.sh <panel.vcf.gz> [output_panel_dir] [run_qc]
#             : panel.vcf.gz       — phased reference panel VCF/BCF for one chromosome
#             : output_panel_dir   — destination for BCF + binary chunks (default: config ref_panel)
#             : run_qc             — true (default) | false  skip norm+SNP-filter if panel is pre-QC'd
# Output      : <output_panel_dir>/reference_panel.<CHR>.bcf[.csi]
#             : <output_panel_dir>/../glimpse2_cache/chunks.<CHR>.txt
#             : <output_panel_dir>/../glimpse2_cache/ref_panel.<CHR>.chunk<N>.bin
# Source      : GLIMPSE2 — https://odelaneau.github.io/GLIMPSE/
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

VCF_INPUT="${1:?Usage: $(basename "$0") <panel.vcf.gz> [output_panel_dir] [run_qc]}"
OUTPUT_PANEL_ARG="${2:-}"
RUN_QC="${3:-true}"   # false = skip norm+SNP-filter (use when panel is already pre-QC'd)

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GLIMPSE2 Reference Panel Prep [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

# Validate input
if [ ! -f "$VCF_INPUT" ]; then
  echo "[X]  CANCELLED. Input file not found: ${VCF_INPUT}"
  exit 1
fi
echo "[<]  ${VCF_INPUT}"

# Detect chromosome from filename (e.g. H1K2-v1.0.chr1.ARPnorm.vcf.gz → chr1)
VCF_BASE=$(basename "$VCF_INPUT")
CHR=$(echo "$VCF_BASE" | grep -oP 'chr[0-9XYM]+' | head -1)
if [ -z "$CHR" ]; then
  echo "[X]  CANCELLED. Could not detect chromosome from filename: ${VCF_BASE}"
  echo "[X]  Filename must contain a token matching chrN, chrX, chrY, or chrM."
  exit 1
fi
echo "[i]  Chromosome: ${CHR}"

# Config file resolution (in priority order):
# 1. PIPELINE_CONFIG env var — explicit override, required when submitting from outside the repo
# 2. SLURM_SUBMIT_DIR — works when sbatch is called from bin/supp/
# 3. $0-relative — interactive runs from the script's own directory
if [[ -n "${PIPELINE_CONFIG:-}" ]]; then
  CONFIG_FILE="$PIPELINE_CONFIG"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/config.yaml" ]]; then
  CONFIG_FILE="${SLURM_SUBMIT_DIR}/config.yaml"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  CONFIG_FILE="${SLURM_SUBMIT_DIR}/../../config/config.yaml"
else
  CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../../config/config.yaml"
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "[X]  Config not found: ${CONFIG_FILE}"
  echo "[X]  Options:"
  echo "[X]    1. Copy or symlink config.yaml into your working directory"
  echo "[X]    2. Set PIPELINE_CONFIG=/absolute/path/to/config/config.yaml"
  exit 1
fi

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

# Resolve output directories
OUTPUT_PANEL="${OUTPUT_PANEL_ARG:-${ref_panel}}"
CACHE_DIR="${OUTPUT_PANEL}/../glimpse2_cache"

echo "[i]  Panel dir  : ${OUTPUT_PANEL}"
echo "[i]  Cache dir  : ${CACHE_DIR}"

# Resolve genetic map for this chromosome
GMAP_FILE="${ref_gmap}/${CHR}.b38.gmap.gz"
if [ ! -f "$GMAP_FILE" ]; then
  echo "[X]  CANCELLED. Genetic map not found: ${GMAP_FILE}"
  exit 1
fi
echo "[i]  Genetic map: ${GMAP_FILE}"

# Scratch setup — fall back to /tmp when scratch_root is unset (local dev)
JOB_ID="${SLURM_JOB_ID:-$(date +%s)}"
if [ -n "${scratch_root:-}" ]; then
  SCRATCH_JOB="${scratch_root}/${USER}/job_${JOB_ID}"
else
  SCRATCH_JOB="${TMPDIR:-/tmp}/glimpse2_prep_${JOB_ID}"
fi
export TMPDIR="${SCRATCH_JOB}/tmp"

echo "[i]  Scratch dir: ${SCRATCH_JOB}"
mkdir -p "$SCRATCH_JOB" "$TMPDIR"

# Load modules on remote
if [[ "$env_type" == "remote" ]]; then
  echo "[i]  Loading modules..."
  module load glimpse 2>/dev/null \
    || echo "[W]  glimpse module unavailable — checking PATH"
  module load bcftools
fi

# Guard: required binaries
for bin in bcftools GLIMPSE2_chunk_static GLIMPSE2_split_reference_static; do
  command -v "$bin" > /dev/null 2>&1 \
    || { echo "[X]  ${bin} not found in PATH"; exit 1; }
done

# Final output sentinel — skip entire job if already complete
SENTINEL="${OUTPUT_PANEL}/reference_panel.${CHR}.bcf"
if [ -s "$SENTINEL" ] && \
   [ -s "${CACHE_DIR}/chunks.${CHR}.txt" ]; then
  # Count expected vs present binary files
  n_chunks=$(wc -l < "${CACHE_DIR}/chunks.${CHR}.txt")
  n_bins=$(find "$CACHE_DIR" -name "ref_panel.${CHR}.chunk*.bin" 2>/dev/null | wc -l)
  if [ "$n_bins" -ge "$n_chunks" ]; then
    echo "[i]  All outputs already present for ${CHR} — nothing to do."
    exit 0
  fi
fi

# Destination BCF and sites paths (persistent storage)
DEST_BCF="${OUTPUT_PANEL}/reference_panel.${CHR}.bcf"
DEST_SITES="${OUTPUT_PANEL}/reference_panel.${CHR}.sites.vcf.gz"
DEST_CHUNKS="${CACHE_DIR}/chunks.${CHR}.txt"

mkdir -p "$OUTPUT_PANEL" "$CACHE_DIR"

# Scratch-local working paths
SCR_VCF="${SCRATCH_JOB}/$(basename "$VCF_INPUT")"
SCR_BCF="${SCRATCH_JOB}/reference_panel.${CHR}.bcf"
SCR_SITES="${SCRATCH_JOB}/reference_panel.${CHR}.sites.vcf.gz"
SCR_CHUNKS="${SCRATCH_JOB}/chunks.${CHR}.txt"
SCR_BIN_PREFIX="${SCRATCH_JOB}/ref_panel.${CHR}"

njobs="${SLURM_CPUS_PER_TASK:-4}"

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## Step 0: Stage input VCF to scratch
step0_stage_input() {
  local step_name="Step 0: Stage Input to Scratch"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  if [ -f "$SCR_VCF" ]; then
    echo "[i]  Already staged ($(basename "$SCR_VCF"))"
  else
    echo "[i]  Copying $(basename "$VCF_INPUT") → scratch..."
    set -o xtrace
    cp "$VCF_INPUT" "$SCR_VCF"
    # Copy index if present
    for idx_ext in .tbi .csi; do
      [ -f "${VCF_INPUT}${idx_ext}" ] && cp "${VCF_INPUT}${idx_ext}" "${SCR_VCF}${idx_ext}"
    done
    set +o xtrace
  fi

  # Index if no index present
  if [ ! -f "${SCR_VCF}.tbi" ] && [ ! -f "${SCR_VCF}.csi" ]; then
    echo "[i]  Indexing input VCF in scratch..."
    set -o xtrace
    bcftools index -t --threads "$njobs" "$SCR_VCF"
    set +o xtrace
  fi

  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 1: QC — normalize, biallelic SNPs only → BCF
## Set RUN_QC=false to skip when the panel is already pre-QC'd; the input VCF
## is converted to BCF in-place so downstream steps have a consistent source.
step1_qc_panel() {
  local step_name="Step 1: QC Reference Panel"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  if [ "$RUN_QC" != "true" ]; then
    echo "[i]  RUN_QC=false — converting input to BCF without norm/filter..."

    if [ -s "$DEST_BCF" ]; then
      echo "[SKIP]  ${DEST_BCF} already exists"
      [ -s "$SCR_BCF" ] || cp "$DEST_BCF" "$SCR_BCF"
      [ -f "${DEST_BCF}.csi" ] && { [ -f "${SCR_BCF}.csi" ] || cp "${DEST_BCF}.csi" "${SCR_BCF}.csi"; }
    else
      set -o xtrace
      bcftools view \
        --threads "$njobs" \
        --output-type b \
        --output "${SCR_BCF}.tmp" \
        "$SCR_VCF"
      set +o xtrace

      [ -s "${SCR_BCF}.tmp" ] || { echo "[X]  CANCELLED: $step_name produced empty output"; exit 1; }
      mv "${SCR_BCF}.tmp" "$SCR_BCF"
      bcftools index --threads "$njobs" "$SCR_BCF"

      cp "$SCR_BCF" "$DEST_BCF"
      cp "${SCR_BCF}.csi" "${DEST_BCF}.csi"
      echo "[>]  $DEST_BCF"
    fi

    echo "[!]  $step_name (QC skipped)"
    local elapsed=$(( EPOCHSECONDS - step_timestamp ))
    printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
    return 0
  fi

  if [ -s "$DEST_BCF" ]; then
    echo "[SKIP]  ${DEST_BCF} already exists"
    [ -s "$SCR_BCF" ] || cp "$DEST_BCF" "$SCR_BCF"
    [ -f "${DEST_BCF}.csi" ] && { [ -f "${SCR_BCF}.csi" ] || cp "${DEST_BCF}.csi" "${SCR_BCF}.csi"; }
    echo "[!]  $step_name"
    return 0
  fi

  echo "[i]  Normalizing and filtering → BCF..."
  set -o xtrace
  bcftools norm \
    --multiallelics -any \
    --threads "$njobs" \
    "$SCR_VCF" \
    --output-type u |
  bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    --types snps \
    --threads "$njobs" \
    --output-type b \
    --output "${SCR_BCF}.tmp"
  set +o xtrace

  [ -s "${SCR_BCF}.tmp" ] || { echo "[X]  CANCELLED: $step_name produced empty output"; exit 1; }
  mv "${SCR_BCF}.tmp" "$SCR_BCF"
  bcftools index --threads "$njobs" "$SCR_BCF"

  echo "[i]  Copying BCF to persistent storage..."
  cp "$SCR_BCF" "$DEST_BCF"
  cp "${SCR_BCF}.csi" "${DEST_BCF}.csi"

  echo "[>]  $DEST_BCF"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 2: Extract sites-only VCF (used by GLIMPSE2_chunk)
step2_extract_sites() {
  local step_name="Step 2: Extract Sites"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  if [ -s "$DEST_SITES" ]; then
    echo "[SKIP]  ${DEST_SITES} already exists"
    [ -s "$SCR_SITES" ] || cp "$DEST_SITES" "$SCR_SITES"
    [ -f "${DEST_SITES}.tbi" ] && { [ -f "${SCR_SITES}.tbi" ] || cp "${DEST_SITES}.tbi" "${SCR_SITES}.tbi"; }
    echo "[!]  $step_name"
    return 0
  fi

  [ -s "$SCR_BCF" ] || { echo "[X]  CANCELLED: QC BCF not found in scratch: $SCR_BCF"; exit 1; }

  echo "[i]  Extracting sites-only VCF..."
  set -o xtrace
  bcftools view \
    --drop-genotypes \
    --threads "$njobs" \
    --output-type z \
    --output "${SCR_SITES}.tmp" \
    "$SCR_BCF"
  set +o xtrace

  [ -s "${SCR_SITES}.tmp" ] || { echo "[X]  CANCELLED: $step_name produced empty output"; exit 1; }
  mv "${SCR_SITES}.tmp" "$SCR_SITES"
  bcftools index -t --threads "$njobs" "$SCR_SITES"

  echo "[i]  Copying sites VCF to persistent storage..."
  cp "$SCR_SITES" "$DEST_SITES"
  cp "${SCR_SITES}.tbi" "${DEST_SITES}.tbi"

  echo "[>]  $DEST_SITES"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 3: GLIMPSE2_chunk — define imputation windows
step3_chunk() {
  local step_name="Step 3: Chunk Chromosome"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  if [ -s "$DEST_CHUNKS" ]; then
    echo "[SKIP]  ${DEST_CHUNKS} already exists ($(wc -l < "$DEST_CHUNKS") chunks)"
    [ -s "$SCR_CHUNKS" ] || cp "$DEST_CHUNKS" "$SCR_CHUNKS"
    echo "[!]  $step_name"
    return 0
  fi

  [ -s "$SCR_SITES" ] || { echo "[X]  CANCELLED: sites VCF not found in scratch: $SCR_SITES"; exit 1; }

  echo "[i]  Running GLIMPSE2_chunk for ${CHR}..."
  set -o xtrace
  GLIMPSE2_chunk_static \
    --input "$SCR_SITES" \
    --region "$CHR" \
    --sequential \
    --map "$GMAP_FILE" \
    --output "${SCR_CHUNKS}.tmp"
  set +o xtrace

  [ -s "${SCR_CHUNKS}.tmp" ] || { echo "[X]  CANCELLED: $step_name produced empty chunk file"; exit 1; }
  mv "${SCR_CHUNKS}.tmp" "$SCR_CHUNKS"

  echo "[i]  Copying chunk file to persistent storage..."
  cp "$SCR_CHUNKS" "$DEST_CHUNKS"

  local n_chunks
  n_chunks=$(wc -l < "$SCR_CHUNKS")
  echo "[>]  $DEST_CHUNKS  (${n_chunks} chunks)"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Step 4: GLIMPSE2_split_reference — build binary reference panel per chunk
step4_split_reference() {
  local step_name="Step 4: Split Reference"
  local step_timestamp=$(date +%s)

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ -s "$SCR_CHUNKS" ] || { echo "[X]  CANCELLED: chunk file not found in scratch: $SCR_CHUNKS"; exit 1; }
  [ -s "$SCR_BCF" ]    || { echo "[X]  CANCELLED: QC BCF not found in scratch: $SCR_BCF"; exit 1; }

  local n_done=0 n_skip=0

  while IFS=$'\t' read -r idx _chr input_region output_region; do
    local dest_bin="${CACHE_DIR}/ref_panel.${CHR}.chunk${idx}.bin"

    if [ -s "$dest_bin" ]; then
      n_skip=$(( n_skip + 1 ))
      continue
    fi

    echo "[&]  Splitting: ${CHR} chunk ${idx}  (${output_region})"
    set -o xtrace
    GLIMPSE2_split_reference_static \
      --reference "$SCR_BCF" \
      --map "$GMAP_FILE" \
      --input-region "$input_region" \
      --output-region "$output_region" \
      --output "${SCR_BIN_PREFIX}.chunk${idx}" \
      --threads "$njobs"
    set +o xtrace

    # GLIMPSE2 appends region coordinates to the output filename (e.g. chunk0_chr22_1_17316684.bin)
    local scr_bin_actual
    scr_bin_actual=$(find "$SCRATCH_JOB" -maxdepth 1 -name "ref_panel.${CHR}.chunk${idx}*.bin" 2>/dev/null | head -1)
    [ -n "$scr_bin_actual" ] || { echo "[X]  CANCELLED: $step_name — no .bin output found for chunk ${idx}"; exit 1; }

    echo "[i]  Copying binary chunk to persistent storage..."
    cp "$scr_bin_actual" "$dest_bin"

    echo "[>]  $dest_bin"
    n_done=$(( n_done + 1 ))
  done < "$SCR_CHUNKS"

  echo "[i]  Binary chunks built: $n_done (cached: $n_skip)"
  echo "[!]  $step_name"
  local elapsed=$(( EPOCHSECONDS - step_timestamp ))
  printf "[&]  Step time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

## Cleanup: remove scratch directory
cleanup_scratch() {
  echo
  echo "[*]  Scratch Cleanup"
  echo "[i]  Removing: ${SCRATCH_JOB}"
  cd /
  rm -rf "$SCRATCH_JOB"
  echo "[W]  Scratch removed."
}

## Finisher
finisher() {
  local n_bins
  n_bins=$(find "$CACHE_DIR" -name "ref_panel.${CHR}.chunk*.bin" 2>/dev/null | wc -l)

  echo
  if [ -s "$DEST_BCF" ] && [ -s "$DEST_CHUNKS" ] && [ "$n_bins" -gt 0 ]; then
    echo "[>]  BCF      : ${DEST_BCF}"
    echo "[>]  Chunks   : ${DEST_CHUNKS}"
    echo "[>]  Bin files: ${n_bins} in ${CACHE_DIR}"
    echo "[$] GLIMPSE2 Reference Panel Prep [FENIX] completed successfully!"
  else
    echo "[X]  GLIMPSE2 Reference Panel Prep [FENIX] INCOMPLETE."
    echo "[X]  BCF    : ${DEST_BCF}"
    echo "[X]  Chunks : ${DEST_CHUNKS}"
    echo "[X]  Bins   : ${n_bins} found in ${CACHE_DIR}"
  fi
  local elapsed=$(( EPOCHSECONDS - script_timestamp ))
  printf "[&]  Total time: %02d:%02d:%02d\n" $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

main() {
  step0_stage_input
  step1_qc_panel
  step2_extract_sites
  step3_chunk
  step4_split_reference
  cleanup_scratch
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
