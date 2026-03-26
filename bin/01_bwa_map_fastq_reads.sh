#!/usr/bin/env bash

# =============================================================================
# Title       : BWA FASTQ Reads Mapper [FENIX]
# Description : Maps FASTQ reads to a reference genome to produce mapped BAMs. Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-03-24
# Version     : 1.0
# Usage       : 01_bwa_map_fastq_reads_WIP.sh [input_dir] [output_path]
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
#             : https://gatk.broadinstitute.org/hc/en-us/articles/360035889471
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

INPUT_DIR="${1:-$PWD}"
OUTPUT_PATH="${2:-$PWD}"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] BWA FASTQ Reads Mapper [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

SAMPLE_NAME=$(basename "$INPUT_DIR")

echo
echo "[i]  Checking input files..."

if find . -maxdepth 1 -name "${SAMPLE_NAME}*.f*q.gz" | grep -q .; then
  echo "[<]  ${INPUT_DIR}"
  echo "[i]  Sample : ${SAMPLE_NAME}"
  echo "[i]  Output : ${OUTPUT_PATH}"
else
  echo "[X]  CANCELLED. FASTQ files starting with ${SAMPLE_NAME}* not found."
  exit 1
fi

# Options
njobs=4
RUN_FASTQC=true
RUN_STATS=true

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
else
  env_type="local"
  njobs=8
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
  module load fastqc
  module load bwa
  module load samtools
fi

# Guard: reference file paths
if [ -z "$ref_gnm" ]; then
  echo "[X]  Missing required reference paths. Check config: $CONFIG_FILE"
  exit 1
fi

# Guard: reference file availability
echo "[i]  References:"
[ -f "${ref_gnm}" ] || { echo "[X]  Reference genome not found: ${ref_gnm}"; exit 1; }
echo "[i]    Genome  : ${ref_gnm}"
[ -f "${ref_gnm}.bwt" ] || { echo "[X]  BWA index not found: ${ref_gnm}.bwt — run: bwa index ${ref_gnm}"; exit 1; }

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## Step 0: Find and map input files per read group
step0_map_reads_per_sample() {
  local step_name="Step 0: Map Input Files"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  # Find R1 files: any sample file ending in _[R]1.f*q.gz (R0/index reads are naturally excluded)
  # R2 is derived by replacing the trailing 1 with 2 (preserving the optional R prefix)
  local r1_files r1_count
  r1_files=$(find . -maxdepth 1 -name "${SAMPLE_NAME}*.f*q.gz" \
    | grep -E '_R?1\.f[^.]*q\.gz$' | sort)
  r1_count=$(echo "$r1_files" | wc -l | tr -d ' ')

  [ -z "$r1_files" ] && { echo "[X]  No R1 files found for: ${SAMPLE_NAME}"; exit 1; }

  echo "[i]  Pairs found: $r1_count"
  read_groups=""

  if [ "$r1_count" -eq 1 ]; then
    # Single pair — no read group key; BAM will be named ${SAMPLE_NAME}.bam
    local R1 R2
    R1="$r1_files"
    R2=$(echo "$R1" | sed -E 's/1(\.f[^.]*q\.gz)$/2\1/')
    [ -f "$R2" ] || { echo "[X]  R2 not found: $R2"; exit 1; }
    echo "[i]    R1: $R1"
    echo "[i]    R2: $R2"
    printf "%s\n%s\n" "$R1" "$R2" > "${OUTPUT_PATH}/${SAMPLE_NAME}_bwa_inputs.txt"
    echo "[>]  ${OUTPUT_PATH}/${SAMPLE_NAME}_bwa_inputs.txt"
  else
    # Multiple pairs — extract PLATE_LANE as read group key; BAM named ${SAMPLE_NAME}_${rg}.bam
    while IFS= read -r R1; do
      local R2 base rg
      R2=$(echo "$R1" | sed -E 's/1(\.f[^.]*q\.gz)$/2\1/')
      [ -f "$R2" ] || { echo "[X]  R2 not found: $R2"; exit 1; }

      base=$(basename "$R1")
      base="${base#"${SAMPLE_NAME}_"}"
      rg=$(echo "$base" | sed -E 's/[A-Za-z0-9]*-1A_//; s/_?R?1\.f[^.]*q\.gz$//')

      read_groups="${read_groups:+$read_groups }$rg"
      echo "[i]    ${rg}:"
      echo "[i]      R1: $R1"
      echo "[i]      R2: $R2"

      printf "%s\n%s\n" "$R1" "$R2" > "${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}_bwa_inputs.txt"
      echo "[>]  ${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}_bwa_inputs.txt"
    done <<< "$r1_files"
  fi

  echo "[!]  $step_name"
}

## Step 0b: Annotate Read Groups (optional, BUILD_RG=false)
step0b_annotate_read_groups() {
  local step_name="Step 0b: Annotate Read Groups"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  echo "[i]  Prepending @RG string to inputs files..."

  if [ -z "$read_groups" ]; then
    # Single pair — construct default @RG from sample name and R1 filename
    local infile="${OUTPUT_PATH}/${SAMPLE_NAME}_bwa_inputs.txt"
    [ -f "$infile" ] || { echo "[X]  Missing inputs file: $infile"; exit 1; }
    local r1 pu rg_line
    r1=$(head -1 "$infile")
    pu=$(basename "$r1" | sed -E 's/\.f[^.]*q\.gz$//')
    rg_line="@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${pu}\tLB:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}"
    { echo "$rg_line"; cat "$infile"; } > "${infile}.tmp"
    mv "${infile}.tmp" "$infile"
    echo "[>]  $infile"
  else
    for rg in $read_groups; do
      local infile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}_bwa_inputs.txt"
      [ -f "$infile" ] || { echo "[X]  Missing inputs file: $infile"; exit 1; }

      # Read R1 path from inputs file and extract barcode from filename
      local r1 r2 barcode rgID rg_line
      { read -r r1; read -r r2; } < "$infile"
      barcode=$(basename "$r1" | sed -nE 's/.*_([A-Za-z0-9]+)-1A_.*/\1/p')
      barcode="${barcode:-$SAMPLE_NAME}"
      rgID="${rg/_/.}"

      rg_line="@RG\tID:${rgID}\tPL:ILLUMINA\tPU:${rgID}.${barcode}\tLB:${barcode}\tSM:${SAMPLE_NAME}"

      { echo "$rg_line"; echo "$r1"; echo "$r2"; } > "${infile}.tmp"
      mv "${infile}.tmp" "$infile"
      echo "[>]  $infile"
    done
  fi

  echo "[!]  $step_name"
}

## Step 1: Check Reads Quality (FastQC)
step1_reads_quality() {
  local step_name="Step 1: Check Reads Quality"
  local outdir="${OUTPUT_PATH}/${SAMPLE_NAME}-fastqc/"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ "$RUN_FASTQC" = true ] || { echo "[i]  Skipped (RUN_FASTQC=false)"; return 0; }

  mkdir -p "$outdir"
  echo "[i]  Output folder: ${outdir}"

  echo "[&]  Running FastQC for all $SAMPLE_NAME reads..."
  if JAVA_TOOL_OPTIONS="-Djava.awt.headless=true" fastqc \
      --outdir "$outdir" \
      --threads "$njobs" \
      --noextract \
      --quiet \
      $(find . -maxdepth 1 -name "${SAMPLE_NAME}*.f*q.gz" | grep -E '_R?[12]\.f[^.]*q\.gz$' | sort); then
    echo "[>]  $outdir"
  else
    echo "[W]  FastQC failed (non-fatal) — skipping quality check, continuing pipeline"
  fi

  echo "[!]  $step_name"
}

## Step 2: BWA Mapping to Reference (per read group)
step2_bwa_mapping_per_readgroup() {
  local step_name="Step 2: BWA Mapping to Reference [Iterative]"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  local items
  if [ -z "$read_groups" ]; then items=(""); else items=($read_groups); fi

  for rg in "${items[@]}"; do
    echo
    echo "[&]  Processing: $SAMPLE_NAME${rg:+ [$rg]}..."
    local infile outfile
    if [ -z "$rg" ]; then
      infile="${OUTPUT_PATH}/${SAMPLE_NAME}_bwa_inputs.txt"
      outfile="${OUTPUT_PATH}/${SAMPLE_NAME}.bam"
    else
      infile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}_bwa_inputs.txt"
      outfile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}.bam"
    fi

    [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }

    # Inputs array: optional @RG line, then R1, R2
    local inputs_array=()
    while IFS= read -r line; do
      inputs_array+=( "$line" )
    done < "$infile"

    local bwa_rg_args=() r1 r2
    if [[ "${inputs_array[0]}" == @RG* ]]; then
      bwa_rg_args=(-R "${inputs_array[0]}")
      r1="${inputs_array[1]}"
      r2="${inputs_array[2]}"
    else
      r1="${inputs_array[0]}"
      r2="${inputs_array[1]}"
    fi

    [ -f "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; continue; }

    echo "[&]  Mapping reads to: ${outfile}"

    set -o xtrace
    bwa mem \
      -t "$njobs" \
      "${bwa_rg_args[@]}" \
      -v 1 \
      "$ref_gnm" \
      "$r1" "$r2" |
    samtools view \
      --threads "$njobs" \
      --bam \
      --output "${outfile}.tmp" -
    set +o xtrace

    mv "${outfile}.tmp" "$outfile"

    [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo "[>]  $outfile"
  done

  echo "[!]  $step_name"
}

## Step 4: Sort and Index BAM files
step3_sort_mapped_bams() {
  local step_name="Step 3: Sort and Index BAM Files"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  local items
  if [ -z "$read_groups" ]; then items=(""); else items=($read_groups); fi

  for rg in "${items[@]}"; do
    local infile outfile
    if [ -z "$rg" ]; then
      infile="${OUTPUT_PATH}/${SAMPLE_NAME}.bam"
      outfile="${OUTPUT_PATH}/${SAMPLE_NAME}.sort.bam"
    else
      infile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}.bam"
      outfile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}.sort.bam"
    fi

    [ -f "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; continue; }

    echo "[&]  Sorting and indexing: $(basename "$infile")"
    set -o xtrace
    samtools sort "$infile" -O bam -o "${outfile}.tmp" -@ "$njobs"

    mv "${outfile}.tmp" "$outfile"

    samtools index "$outfile" -@ "$njobs"
    set +o xtrace

    echo "[>]  $outfile"
  done

  echo "[!]  $step_name"
}

## Step 4: BAM Statistics (samtools stats)
step4_bam_stats() {
  local step_name="Step 4: BAM Statistics"

  echo
  echo "[*]  $step_name"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  [ "$RUN_STATS" = true ] || { echo "[i]  Skipped (RUN_STATS=false)"; return 0; }

  local items
  if [ -z "$read_groups" ]; then items=(""); else items=($read_groups); fi

  for rg in "${items[@]}"; do
    local infile outfile
    if [ -z "$rg" ]; then
      infile="${OUTPUT_PATH}/${SAMPLE_NAME}.sort.bam"
      outfile="${OUTPUT_PATH}/${SAMPLE_NAME}.sort.stats"
    else
      infile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}.sort.bam"
      outfile="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}.sort.stats"
    fi

    [ -f "$infile" ] || { echo "[X]  Missing input: $infile"; exit 1; }
    [ -f "$outfile" ] && { echo "[i]  Already completed ($outfile exists)"; continue; }

    echo "[&]  Running samtools stats: $(basename "$infile")"
    set -o xtrace
    samtools stats \
      --reference "$ref_gnm" \
      --threads "$njobs" \
      "$infile" \
      > "${outfile}.tmp"
    set +o xtrace

    mv "${outfile}.tmp" "$outfile"

    [ -s "$outfile" ] || { echo "[X]  CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo "[>]  $outfile"
  done

  echo "[!]  $step_name"
}

## Finisher: check final output and report
finisher() {
  local items all_ok=true
  if [ -z "$read_groups" ]; then items=(""); else items=($read_groups); fi

  echo
  for rg in "${items[@]}"; do
    local bam
    [ -z "$rg" ] && bam="${OUTPUT_PATH}/${SAMPLE_NAME}.sort.bam" \
                 || bam="${OUTPUT_PATH}/${SAMPLE_NAME}_${rg}.sort.bam"
    if [ -f "$bam" ]; then
      echo "[>]  $bam"
    else
      echo "[X]  INCOMPLETE. Output not found: $bam"
      all_ok=false
    fi
  done

  local elapsed=$(( EPOCHSECONDS - script_timestamp ))
  local h=$(( elapsed / 3600 )) m=$(( (elapsed % 3600) / 60 )) s=$(( elapsed % 60 ))
  printf "[&]  Total time: %02d:%02d:%02d\n" "$h" "$m" "$s"

  if [ "$all_ok" = true ]; then
    echo "[$] BWA FASTQ Reads Mapper [FENIX] completed successfully!"
    exit 0
  else
    exit 1
  fi
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

main() {
  step0_map_reads_per_sample
  step0b_annotate_read_groups
  step1_reads_quality
  step2_bwa_mapping_per_readgroup
  step3_sort_mapped_bams
  step4_bam_stats
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
