#!/usr/bin/env bash

# =============================================================================
# Title       : GATK GenotypeGVCFs [FENIX]
# Description : Joint genotyping. Reads the per-chromosome GenomicsDB
#               workspaces built by Step 04, runs GenotypeGVCFs against each
#               (gendb://...), then gathers the per-chromosome joint VCFs into
#               one cohort-level callset. Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-09-12
# Version     : 1.0
# Usage       : 05_gatk_GenotypeGVCFs.sh <genomicsdb_path> <output_path> <chrom> [cohort_name]
#             : genomicsdb_path — the <output_path> given to Step 04; workspaces
#             :                   are read from <genomicsdb_path>/genomicsdb/<chrom>
#             : output_path     — per-chrom VCFs go to <output_path>/chrom_vcf/,
#             :                   the gathered callset to <output_path>/<cohort>.joint.vcf.gz
#             : chrom           — chr1..chr22 | chrX | chrY | chrM | autosomes | all
#             : cohort_name     — label used in output filenames (default: "cohort")
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# =============================================================================
#
# ---------------------------------------------------------------------------
#  STATUS
# ---------------------------------------------------------------------------
#
#   DRAFT — written against synthetic-scale reasoning only, not yet run on
#   FENIX or against a real GenomicsDB workspace. Before production use:
#     - confirm Java heap sizing (GENO_JAVA_MEM below is a placeholder copied
#       from Step 04's measured GenomicsDBImport footprint; GenotypeGVCFs'
#       own memory profile against a multi-fragment, whole-chr1-size DB is
#       unmeasured)
#     - confirm whether reading a GenomicsDB during GenotypeGVCFs benefits
#       from /scratch the way Step 04's TileDB writes did (unlike Step 04,
#       this script reads the workspace in place over NFS — no scratch
#       staging is implemented here pending that measurement)
#     - GenotypeGVCFs is effectively single-threaded per contig; parallelism
#       comes only from scattering chromosomes across separate job submissions
#       (as Step 04 does) — no per-chromosome SLURM launcher exists yet
#
# ---------------------------------------------------------------------------
#  ALLELE-SPECIFIC ANNOTATIONS — deliberately NOT enabled
# ---------------------------------------------------------------------------
#
#   This runs GenotypeGVCFs with GATK's default (non-AS) annotation set.
#   Allele-specific VQSR needs -G AS_StandardAnnotation at BOTH HaplotypeCaller
#   (Step 03a) and GenotypeGVCFs — Step 03a does not currently request it, so
#   turning it on only here would produce incomplete AS annotations. Whether
#   Step 06 uses allele-specific VQSR, classic VQSR, or hard-filtering is still
#   undecided (see docs/PIPELINE_STATUS.md); revisit Step 03a first if
#   allele-specific VQSR is chosen.
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

GDB_BASE="${1:?Usage: $(basename "$0") <genomicsdb_path> <output_path> <chrom> [cohort_name]}"
OUTPUT_PATH="${2:?Usage: $(basename "$0") <genomicsdb_path> <output_path> <chrom> [cohort_name]}"
CHROM_ARG="${3:?Usage: $(basename "$0") <genomicsdb_path> <output_path> <chrom> [cohort_name]}"
COHORT="${4:-cohort}"

[[ "$COHORT" == "cohort" ]] && echo "[!]  No cohort_name given — using generic 'cohort' for output filenames." >&2

# <\ARGS> -----------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GATK GenotypeGVCFs [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking inputs..."

GDB_BASE="$(readlink -f "$GDB_BASE")"
[ -d "${GDB_BASE}/genomicsdb" ] || { echo "[X]  CANCELLED. No genomicsdb/ under: $GDB_BASE (run Step 04 first)"; exit 1; }
mkdir -p "$OUTPUT_PATH"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"

echo "[<]  GenomicsDB : ${GDB_BASE}/genomicsdb/<chrom>"
echo "[i]  Output     : ${OUTPUT_PATH}"
echo "[i]  Cohort     : ${COHORT}"

# Options (env-overridable so a launcher can vary resources per run)
#   GENO_JAVA_MEM   — Java -Xms/-Xmx                       (default 6G remote / 4G local; PLACEHOLDER, see STATUS above)
#   GENO_VERBOSITY  — GenotypeGVCFs --verbosity             (default ERROR)
#   GENO_THREADS    — bcftools threads for the gather step  (default: $SLURM_CPUS_PER_TASK, else 2)
VERBOSITY="${GENO_VERBOSITY:-ERROR}"
njobs="${GENO_THREADS:-${SLURM_CPUS_PER_TASK:-2}}"

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
  MEM="${GENO_JAVA_MEM:-6G}"
else
  env_type="local"
  MEM="${GENO_JAVA_MEM:-4G}"
fi
echo "[i]  Environment: $env_type"
echo "[i]  Knobs      : java-mem=${MEM}  verbosity=${VERBOSITY}  gather-threads=${njobs}"

# Parse YAML config into Bash variables (embedded parser — no external tool)
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
  module unload oracle-java 2>/dev/null || true
  module load oracle-java/25.0.2
  module load gatk
  module load bcftools
fi

# Guard: reference file paths
if [ -z "${ref_gnm:-}" ] || [ -z "${ref_vars:-}" ]; then
  echo "[X]  Missing required reference paths (ref_gnm / ref_vars). Check config."
  exit 1
fi
[ -f "${ref_gnm}" ]  || { echo "[X]  Reference genome not found: ${ref_gnm}"; exit 1; }
[ -f "${ref_vars}" ] || { echo "[X]  Reference variants (dbSNP) not found: ${ref_vars}"; exit 1; }
echo "[i]  References : genome=${ref_gnm}  dbsnp=${ref_vars}"

TMPDIR_LOCAL="${OUTPUT_PATH}/.geno_tmp_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR_LOCAL"
export TMPDIR="$TMPDIR_LOCAL"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

# <\ENV> --------------------------------------------------------------------

# <FUNCTIONS> ---------------------------------------------------------------

AUTOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
ALLCHROMS=("${AUTOSOMES[@]}" chrX chrY chrM)
CHROM_VCFS=()   # per-chromosome outputs, in canonical order, fed to the gather step

## Resolve the requested selector into a list of chromosomes.
resolve_chrom_list() {
  case "$CHROM_ARG" in
    autosomes) CHROMS=("${AUTOSOMES[@]}") ;;
    all)       CHROMS=("${ALLCHROMS[@]}") ;;
    chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY|chrM|chrMT)
               CHROMS=("$CHROM_ARG") ;;
    *) echo "[X]  Invalid chrom selector: '$CHROM_ARG'"
       echo "[i]  Expected: chr1..chr22 | chrX | chrY | chrM | autosomes | all"
       exit 1 ;;
  esac
}

## Match a requested chrom to the workspace directory Step 04 actually built
## (mirrors the chrM/chrMT fallback in 04_gatk_GenomicsDB_import.sh's finisher).
resolve_workspace_chrom() {
  local want="$1"
  [[ -d "${GDB_BASE}/genomicsdb/${want}" ]] && { echo "$want"; return 0; }
  if [[ "$want" == "chrM" || "$want" == "chrMT" ]]; then
    local alt="chrMT"; [[ "$want" == "chrMT" ]] && alt="chrM"
    if [[ -d "${GDB_BASE}/genomicsdb/${alt}" ]]; then
      echo "[!]  WARNING: requested '$want' but workspace uses '$alt' notation — using '$alt'" >&2
      echo "$alt"; return 0
    fi
  fi
  echo "$want"   # let the caller's existence guard produce the missing-workspace error
}

## Joint-genotype one chromosome's GenomicsDB workspace.
genotype_one_chrom() {
  local req="$1"
  local step_timestamp=$EPOCHSECONDS
  local chr; chr="$(resolve_workspace_chrom "$req")"
  local ws="${GDB_BASE}/genomicsdb/${chr}"
  local outdir="${OUTPUT_PATH}/chrom_vcf"
  local outfile="${outdir}/${COHORT}.joint.${chr}.vcf.gz"
  local tmp="${outdir}/.${COHORT}.joint.${chr}.tmp.vcf.gz"

  echo
  echo "[*]  Chromosome ${req}$( [[ "$chr" != "$req" ]] && echo "  (effective: ${chr})" )"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ -s "$outfile" && -s "${outfile}.tbi" ]]; then
    echo "[i]  Already completed (${outfile} exists)"
    CHROM_VCFS+=("$outfile")
    return 0
  fi

  [ -s "${ws}/callset.json" ] || { echo "[X]  No GenomicsDB workspace for ${req}: ${ws} (run Step 04 first)"; exit 1; }

  mkdir -p "$outdir"
  rm -f "$tmp" "${tmp}.tbi"

  set -o xtrace
  gatk --java-options "-Xms${MEM} -Xmx${MEM}" GenotypeGVCFs \
    --variant "gendb://${ws}" \
    --reference "$ref_gnm" \
    --dbsnp "$ref_vars" \
    --intervals "$chr" \
    --tmp-dir "$TMPDIR" \
    --verbosity "$VERBOSITY" \
    --create-output-variant-index \
    --output "$tmp"
  set +o xtrace

  [ -s "$tmp" ] || { echo "[X]  CANCELLED: GenotypeGVCFs produced no output for ${chr}"; exit 1; }
  mv "$tmp" "$outfile"
  mv "${tmp}.tbi" "${outfile}.tbi"
  CHROM_VCFS+=("$outfile")

  echo "[>]  $outfile"
  echo "[!]  Chromosome ${chr} — genotyped"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Gather all per-chromosome joint VCFs into one cohort-level callset.
gather_cohort_vcf() {
  local step_timestamp=$EPOCHSECONDS
  local final="${OUTPUT_PATH}/${COHORT}.joint.vcf.gz"

  echo
  echo "[*]  Gather: cohort callset"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  (( ${#CHROM_VCFS[@]} >= 1 )) || { echo "[X]  CANCELLED: no per-chromosome VCFs to gather"; exit 1; }

  if [[ -s "$final" && -s "${final}.tbi" ]]; then
    echo "[i]  Already completed (${final} exists)"
    return 0
  fi

  if (( ${#CHROM_VCFS[@]} == 1 )); then
    cp "${CHROM_VCFS[0]}" "$final"
    cp "${CHROM_VCFS[0]}.tbi" "${final}.tbi"
  else
    # Guard: every per-chromosome VCF must carry the same sample set. An uneven
    # Step 04 wave import (a sample added to some chromosomes' workspaces but
    # not others) would otherwise silently corrupt the concat.
    local ref_samples; ref_samples="$(bcftools query -l "${CHROM_VCFS[0]}")"
    local f samples
    for f in "${CHROM_VCFS[@]:1}"; do
      samples="$(bcftools query -l "$f")"
      if [[ "$samples" != "$ref_samples" ]]; then
        echo "[X]  CANCELLED: sample set in $(basename "$f") differs from $(basename "${CHROM_VCFS[0]}")"
        echo "[i]  Likely an uneven Step 04 wave import across chromosomes — reconcile before gathering."
        exit 1
      fi
    done

    local tmp="${OUTPUT_PATH}/.${COHORT}.joint.tmp.vcf.gz"
    rm -f "$tmp" "${tmp}.tbi"
    bcftools concat "${CHROM_VCFS[@]}" \
      --output-type z \
      --threads "$njobs" \
      --write-index=tbi \
      --output "$tmp"
    [ -s "$tmp" ] || { echo "[X]  CANCELLED: bcftools concat produced no output"; exit 1; }
    mv "$tmp" "$final"
    mv "${tmp}.tbi" "${final}.tbi"
  fi

  [ -s "$final" ] || { echo "[X]  CANCELLED: gather failed, output missing: $final"; exit 1; }
  echo "[>]  $final"
  echo "[!]  Gather done (${#CHROM_VCFS[@]} chromosome$( [[ ${#CHROM_VCFS[@]} -gt 1 ]] && echo "s" ))"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Finisher: verify outputs and report.
finisher() {
  echo
  local missing=0
  for req in "${CHROMS[@]}"; do
    local chr; chr="$(resolve_workspace_chrom "$req")"
    local f="${OUTPUT_PATH}/chrom_vcf/${COHORT}.joint.${chr}.vcf.gz"
    if [[ -s "$f" ]]; then
      echo "[>]  $f"
    else
      echo "[X]  MISSING per-chromosome VCF for ${req}"
      missing=1
    fi
  done

  local final="${OUTPUT_PATH}/${COHORT}.joint.vcf.gz"
  if [[ -s "$final" ]]; then
    echo "[>]  $final"
  else
    echo "[X]  MISSING gathered cohort VCF: $final"
    missing=1
  fi

  echo
  if (( missing == 0 )); then
    echo "[$] GATK GenotypeGVCFs [FENIX] completed successfully!"
    echo "[i]  Sanity check before Step 06, e.g.:"
    echo "[i]    bcftools stats '${final}' | grep 'number of samples:'"
    echo "[i]  Next: Step 06 (VQSR / hard-filtering) reads '${final}'."
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 0
  else
    echo "[X]  GATK GenotypeGVCFs [FENIX] INCOMPLETE — see MISSING entries above."
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 1
  fi
}

# <\FUNCTIONS> --------------------------------------------------------------

# <MAIN> ------------------------------------------------------------------------

main() {
  resolve_chrom_list

  echo
  echo "[i]  Chromosomes to process (${#CHROMS[@]}): ${CHROMS[*]}"
  [[ ${#CHROMS[@]} -gt 1 ]] && echo "[i]  Note: chromosomes run serially in this process — no per-chromosome" \
                            && echo "[i]        SLURM launcher exists yet for this step (see STATUS in header)."

  for req in "${CHROMS[@]}"; do
    genotype_one_chrom "$req"
  done

  gather_cohort_vcf
  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
