#!/usr/bin/env bash

# =============================================================================
# Title       : GATK GenotypeGVCFs [FENIX]
# Description : Joint genotyping of one chromosome: reads the Step 04 GenomicsDB
#               workspace for that chromosome and writes one multi-sample VCF
#               (all cohort samples) for Step 06 (variant filtering).
#               Optional sub-region for quick tests (e.g. 500 kb of chr22).
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-11-24
# Version     : 1.0
# Usage       : 05_gatk_GenotypeGVCFs.sh <genomicsdb_root> <output_path> <chrom> [region]
#             : genomicsdb_root — Step 04's <output_path>/genomicsdb (one workspace
#             :                   per chromosome: <genomicsdb_root>/<chrom>)
#             : output_path     — VCF goes to <output_path>/<prefix>.<chrom>.vcf.gz
#             :                   (+ .tbi). Group storage only, never $HOME.
#             : chrom           — chr1..chr22 | chrX | chrY | chrM
#             : region          — optional, inside <chrom>: chr22:20000001-20500000
#             :                   Output then ends in .<chrom>_<start>-<end>.vcf.gz
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# =============================================================================
#
# Options (environment variables):
#   GTGVCF_PREFIX         output file prefix                 (default: cohort)
#   GTGVCF_JAVA_MEM       Java -Xmx                          (default 12G remote / 4G local)
#                         GenomicsDB reads through native (non-heap) memory too,
#                         so give the SLURM job ~2x this.
#   GTGVCF_SCRATCH_ROOT   overrides config `scratch_root` for --tmp-dir; use the
#                         PRIMARY-group volume /scratch/groups/$(id -gn)
#   GTGVCF_DBSNP          true (default): fill the ID column from config ref_vars
#   GTGVCF_PROGRESS_SECS  seconds between GATK progress lines (default 600) —
#                         the log then shows how far along the chromosome it is
#
# Skip logic: an existing non-empty <output>.vcf.gz WITH its .tbi is kept. The
# VCF is written as <output>.tmp.vcf.gz and renamed only after it passes the
# checks (BGZF end-of-file marker, sample count = workspace sample count).
# =============================================================================

set -euo pipefail

# <ARGUMENTS> ---------------------------------------------------------------

USAGE="Usage: $(basename "$0") <genomicsdb_root> <output_path> <chrom> [region]"
GDB_ROOT="${1:?$USAGE}"
OUTPUT_PATH="${2:?$USAGE}"
CHROM="${3:?$USAGE}"
REGION="${4:-}"

# <\ARGS> -------------------------------------------------------------------

# <ENVIRONMENT> -------------------------------------------------------------

echo
echo "[$] GATK GenotypeGVCFs [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

PREFIX="${GTGVCF_PREFIX:-cohort}"
USE_DBSNP="${GTGVCF_DBSNP:-true}"
PROGRESS_SECS="${GTGVCF_PROGRESS_SECS:-600}"

case " chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM " in
  *" $CHROM "*) ;;
  *) echo "[X]  CANCELLED. Invalid chromosome: '$CHROM' (chr1..chr22 | chrX | chrY | chrM)"; exit 1 ;;
esac

# Interval: the whole chromosome, or a sub-region of it (quick tests).
if [[ -n "$REGION" ]]; then
  if ! [[ "$REGION" =~ ^${CHROM}:([0-9]+)-([0-9]+)$ ]]; then
    echo "[X]  CANCELLED. Region must look like ${CHROM}:<start>-<end> (got '$REGION')"; exit 1
  fi
  r_start="${BASH_REMATCH[1]}"; r_end="${BASH_REMATCH[2]}"
  (( r_start >= 1 && r_end >= r_start )) || { echo "[X]  CANCELLED. Bad region bounds: $REGION"; exit 1; }
  INTERVAL="$REGION"
  OUT_STEM="${PREFIX}.${CHROM}_${r_start}-${r_end}"
else
  INTERVAL="$CHROM"
  OUT_STEM="${PREFIX}.${CHROM}"
fi

[[ -d "$GDB_ROOT" ]] || { echo "[X]  CANCELLED. genomicsdb_root not found: $GDB_ROOT"; exit 1; }
GDB_ROOT="$(readlink -f "$GDB_ROOT")"
WS="${GDB_ROOT}/${CHROM}"
[[ -s "$WS/callset.json" && -s "$WS/vidmap.json" ]] \
  || { echo "[X]  CANCELLED. No complete Step 04 workspace for ${CHROM}: $WS (callset.json + vidmap.json)"; exit 1; }

# Samples in the workspace (callset.json holds one "row_idx" per sample).
N_WS=$(tr -d '\n' < "$WS/callset.json" | grep -o '"row_idx"' | wc -l | tr -d '[:space:]')
(( N_WS >= 1 )) || { echo "[X]  CANCELLED. No samples listed in $WS/callset.json"; exit 1; }

home_guard() {
  case "$1/" in
    "${HOME}"/*)
      echo "[X]  CANCELLED. Output is under \$HOME (${HOME}) — home has a hard quota on FENIX."
      exit 1 ;;
  esac
}
home_guard "$OUTPUT_PATH"                 # before mkdir: leave nothing behind
mkdir -p "$OUTPUT_PATH"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"
home_guard "$OUTPUT_PATH"                 # again, after resolving symlinks
outfile="${OUTPUT_PATH}/${OUT_STEM}.vcf.gz"
tmpfile="${OUTPUT_PATH}/${OUT_STEM}.tmp.vcf.gz"

echo "[<]  Workspace  : gendb://${WS}  (${N_WS} samples)"
echo "[i]  Interval   : ${INTERVAL}"
echo "[>]  Output     : ${outfile}"

# Skip if already done
if [[ -s "$outfile" && -s "${outfile}.tbi" ]]; then
  echo "[SKIP] Output already exists (with index): $outfile"
  echo "[i]  Delete it (and its .tbi) to genotype ${INTERVAL} again."
  exit 0
fi

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
  MEM="${GTGVCF_JAVA_MEM:-12G}"
else
  env_type="local"
  MEM="${GTGVCF_JAVA_MEM:-4G}"
fi
echo "[i]  Environment: $env_type   java-mem=${MEM}   dbsnp=${USE_DBSNP}"

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

if [[ "$env_type" == "remote" ]]; then
  echo "[i]  Loading modules..."
  module unload oracle-java 2>/dev/null || true
  module load oracle-java/25.0.2
  module load gatk
  module load bcftools
fi

[[ -n "${ref_gnm:-}" && -f "${ref_gnm}" ]] || { echo "[X]  Reference genome not found: ${ref_gnm:-<unset>}"; exit 1; }
echo "[i]  Reference  : ${ref_gnm}"
DBSNP_ARGS=()
if [[ "$USE_DBSNP" == true ]]; then
  [[ -n "${ref_vars:-}" && -f "${ref_vars}" ]] || { echo "[X]  dbSNP (config ref_vars) not found: ${ref_vars:-<unset>}  (or set GTGVCF_DBSNP=false)"; exit 1; }
  DBSNP_ARGS=(--dbsnp "$ref_vars")
  echo "[i]  dbSNP      : ${ref_vars}"
fi

# Temp space: scratch when available (GATK spills sort/temp files there).
scratch_root="${GTGVCF_SCRATCH_ROOT:-${scratch_root:-}}"
if [[ -n "$scratch_root" && "$env_type" == "remote" ]]; then
  SCRATCH_BASE="${scratch_root}/${USER}/gtgvcf_${SLURM_JOB_ID:-$$}_$(date +%s)"
else
  SCRATCH_BASE="${OUTPUT_PATH}/.gtgvcf_work_${SLURM_JOB_ID:-$$}_$(date +%s)"
fi
mkdir -p "$SCRATCH_BASE/tmp" || { echo "[X]  Cannot create temp dir: $SCRATCH_BASE/tmp"; exit 1; }
export TMPDIR="$SCRATCH_BASE/tmp"
echo "[i]  Temp       : $TMPDIR"

cleanup() { rm -rf "$SCRATCH_BASE"; }
trap cleanup EXIT

# <\ENV> --------------------------------------------------------------------

# <FUNCTIONS> ---------------------------------------------------------------

## A well-formed BGZF file ends with a fixed 28-byte empty block (htslib's test).
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
bgzf_complete() {
  [[ "$(tail -c 28 "$1" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')" == "$BGZF_EOF" ]]
}

step1_genotype() {
  local step_timestamp; step_timestamp=$(date +%s)
  local only_starting=()
  [[ -n "$REGION" ]] && only_starting=(--only-output-calls-starting-in-intervals)

  echo
  echo "<START> GenotypeGVCFs ${INTERVAL}  ($(date))"
  rm -f "$tmpfile" "${tmpfile}.tbi"

  set -o xtrace
  gatk --java-options "-Xms2G -Xmx${MEM}" GenotypeGVCFs \
    --reference "$ref_gnm" \
    --variant "gendb://${WS}" \
    --intervals "$INTERVAL" \
    "${only_starting[@]+"${only_starting[@]}"}" \
    "${DBSNP_ARGS[@]+"${DBSNP_ARGS[@]}"}" \
    --tmp-dir "$TMPDIR" \
    --genomicsdb-shared-posixfs-optimizations true \
    --seconds-between-progress-updates "$PROGRESS_SECS" \
    --create-output-variant-index true \
    --output "$tmpfile"
  { set +o xtrace; } 2>/dev/null

  echo "<SUCCESS> GenotypeGVCFs  [$(echo $(( $(date +%s) - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')]"
}

step2_verify_and_rename() {
  echo
  echo "[i]  Verifying output..."
  [[ -s "$tmpfile" && -s "${tmpfile}.tbi" ]] || { echo "<ERROR> GenotypeGVCFs left no VCF + index: $tmpfile"; exit 1; }
  bgzf_complete "$tmpfile" || { echo "<ERROR> VCF looks truncated (no BGZF EOF marker): $tmpfile"; exit 1; }

  local n_vcf n_rec
  n_vcf=$(bcftools query -l "$tmpfile" | wc -l | tr -d '[:space:]')
  if (( n_vcf != N_WS )); then
    echo "<ERROR> VCF has ${n_vcf} samples, workspace has ${N_WS}: $tmpfile"; exit 1
  fi
  # Full read: counts records (GATK's .tbi carries no counts) AND decompresses
  # every block, so a corrupt block fails here instead of in Step 06.
  # Minutes on a whole chromosome.
  if ! n_rec=$(bcftools view -H "$tmpfile" | wc -l | tr -d '[:space:]'); then
    echo "<ERROR> VCF could not be read to the end (corrupt block?): $tmpfile"; exit 1
  fi

  mv -f "${tmpfile}.tbi" "${outfile}.tbi"
  mv -f "$tmpfile" "$outfile"
  echo "[DONE] ${outfile}"
  echo "[i]  Samples: ${n_vcf}   Records: ${n_rec}   Size: $(du -h "$outfile" | cut -f1)"
}

# <\FUNCTIONS> --------------------------------------------------------------

# <MAIN> --------------------------------------------------------------------

step1_genotype
step2_verify_and_rename

echo
echo "[$] GATK GenotypeGVCFs [FENIX] completed successfully!"
echo "[i]  Quick look:  bcftools stats '${outfile}' | grep '^SN'"
echo "[i]  Next: Step 06 (variant filtering) on the per-chromosome VCFs."
echo "[&]  Total time: $(echo $(( $(date +%s) - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"

# <\MAIN> -------------------------------------------------------------------

#EOF
