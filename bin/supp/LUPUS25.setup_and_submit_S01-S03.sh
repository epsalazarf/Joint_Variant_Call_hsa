#!/usr/bin/env bash

# =============================================================================
# Title       : LUPUS 2025 — Batch setup + S01→S03 launcher [FENIX]
# Description : Stages raw FASTQ pairs for the pending Lupus2025 samples into the
#               BAMQC working tree and submits Steps 01→02→03 as chained SLURM
#               jobs (one self-sequencing chain per sample) via the repo's
#               bin/PIPELINE.single_sample_01-03.sh launcher.
# Author      : generated for epsalazarf
# Date        : 2026-09-07
# Usage       : LUPUS25.setup_and_submit_S01-S03.sh [batchlist.tsv] [action]
#                 action : setup   — create dirs + symlinks + status report (default)
#                          submit  — setup, then submit the SLURM chains
#                          list    — parse + report only, touch nothing
#               Env overrides: REPO_DIR, RAW_ROOT, WORK_ROOT, ASSUME_YES=1
# =============================================================================

set -euo pipefail

# <CONFIG> ------------------------------------------------------------------------

BATCHLIST="${1:-Lupus25_batchlist.tsv}"
ACTION="${2:-setup}"

# Repo root (auto-detected from this script's location if it sits inside the repo,
# else falls back to this guess — override with REPO_DIR=... if wrong).
_self="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO_DIR="${REPO_DIR:-}"
if [[ -z "$REPO_DIR" ]]; then
  for cand in "$_self" "$_self/.." "$_self/../.." "$HOME/GitHub/lambda/pipelines/Joint_Variant_Call_hsa"; do
    if [[ -f "$cand/bin/PIPELINE.single_sample_01-03.sh" ]]; then
      REPO_DIR="$(cd "$cand" && pwd)"; break
    fi
  done
fi

RAW_ROOT="${RAW_ROOT:-/mnt/data/amedina/amedina/Lupus2025/01.RawData}"
WORK_ROOT="${WORK_ROOT:-/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC}"
LAUNCHER="${REPO_DIR:-/nonexistent}/bin/PIPELINE.single_sample_01-03.sh"

# Final per-sample output that means "this sample is already done" (Step 03).
DONE_GLOB='*.raw_variants.canon_chr.g.vcf.gz'

# <\CONFIG> ---------------------------------------------------------------------

# <CHECKS> --------------------------------------------------------------------

die() { echo "<ERROR> $*" >&2; exit 1; }

[[ -f "$BATCHLIST" ]]            || die "Batch list not found: $BATCHLIST"
[[ -n "${REPO_DIR:-}" && -d "$REPO_DIR" ]] || die "Repo dir not found — set REPO_DIR=/path/to/Joint_Variant_Call_hsa"
[[ -f "$LAUNCHER" ]]            || die "Launcher not found: $LAUNCHER"

if [[ "$ACTION" != "list" ]]; then
  [[ -d "$RAW_ROOT" ]]  || die "Raw FASTQ root not found: $RAW_ROOT"
  mkdir -p "$WORK_ROOT" || die "Cannot create work root: $WORK_ROOT"
fi
if [[ "$ACTION" == "submit" ]]; then
  command -v sbatch &>/dev/null || die "sbatch not found — run 'submit' on a SLURM login node."
fi

echo "[i] repo      : $REPO_DIR"
echo "[i] launcher  : $LAUNCHER"
echo "[i] raw root  : $RAW_ROOT"
echo "[i] work root : $WORK_ROOT"
echo "[i] batchlist : $BATCHLIST"
echo "[i] action    : $ACTION"
echo

# <\CHECKS> -----------------------------------------------------------------------

# <PARSE> -------------------------------------------------------------------------
# Columns (TAB): 1 SAMPLE_ID  2 LIBRARY  3 SUBSAMPLE_ID  4 RG_identifier
#                5 MATCHED_FASTQ  6 MAPPED_BAM  7 RAW_FILE_PREFIX  8 OUTPUT_NAME  9 BAMQC_BATCH
#
# RAW_FILE_PREFIX is a shell glob for the pair, e.g. "L065_[12].fq.gz" or
# "L061_CKDN250015297-1A_22W5WYLT4_L3_[12].fq.gz".

ROWS=()
while IFS= read -r line; do
  ROWS+=("$line")
done < <(
  awk -F'\t' '
    { sub(/\r$/, "") }
    NR == 1 { next }                              # header
    $1 == "" || $1 ~ /^#/ || $1 == "SAMPLE_ID" { next }
    {
      batch = $9; sub(/^B/, "", batch)            # B16 -> 16
      printf "%s\t%s\t%s\n", $1, $7, batch
    }
  ' "$BATCHLIST"
)

[[ ${#ROWS[@]} -gt 0 ]] || die "No sample rows parsed from $BATCHLIST"

printf '%-10s  %-8s  %-46s  %-9s  %s\n' SAMPLE BATCH FASTQ_GLOB STATUS DEST
printf '%.0s-' {1..110}; echo

SETUP_OK=(); SETUP_DONE=(); SETUP_FAIL=()
SAMPLE_DIRS=()

for row in "${ROWS[@]}"; do
  IFS=$'\t' read -r sample glob bnum <<< "$row"
  batchdir="Batch$(printf '%02d' "$bnum")"
  dest="$WORK_ROOT/$batchdir/$sample"
  src_dir="$RAW_ROOT/$sample"

  # already-done?
  status="PENDING"
  if compgen -G "$dest/$DONE_GLOB" >/dev/null 2>&1 || [[ -d "$dest/chrom_gvcf" ]]; then
    status="DONE"
  fi

  printf '%-10s  %-8s  %-46s  %-9s  %s\n' "$sample" "$batchdir" "$glob" "$status" "$dest"

  SAMPLE_DIRS+=("$dest")

  [[ "$ACTION" == "list" ]] && continue
  if [[ "$status" == "DONE" ]]; then SETUP_DONE+=("$sample"); continue; fi

  # resolve the pair from the raw tree
  shopt -s nullglob
  pair=( "$src_dir"/$glob )
  shopt -u nullglob
  if [[ ${#pair[@]} -lt 2 ]]; then
    # fall back to the generic <sample>_[12].fq.gz convention
    shopt -s nullglob
    pair=( "$src_dir"/${sample}_*[12].f*q.gz )
    shopt -u nullglob
  fi
  if [[ ${#pair[@]} -lt 2 ]]; then
    echo "    <WARN> could not resolve a FASTQ pair in $src_dir (glob: $glob) — skipped"
    SETUP_FAIL+=("$sample")
    continue
  fi

  mkdir -p "$dest"
  for f in "${pair[@]}"; do ln -sf "$f" "$dest/"; done
  SETUP_OK+=("$sample")
done

echo
echo "[i] pending & staged : ${#SETUP_OK[@]}   ${SETUP_OK[*]:-}"
echo "[i] already done     : ${#SETUP_DONE[@]}   ${SETUP_DONE[*]:-}"
echo "[i] unresolved       : ${#SETUP_FAIL[@]}   ${SETUP_FAIL[*]:-}"
echo

[[ "$ACTION" == "list"  ]] && exit 0
[[ "$ACTION" == "setup" ]] && {
  echo "[i] Setup complete. Review the table above, then submit with:"
  echo "      $0 $BATCHLIST submit"
  exit 0
}

# <\PARSE> ----------------------------------------------------------------------

# <SUBMIT> ------------------------------------------------------------------------

if [[ "${ASSUME_YES:-0}" != "1" ]]; then
  n=$(( ${#SETUP_OK[@]} ))
  read -r -p "[?] Submit S01->S02->S03 chains for $n sample(s) (~$((n*3)) SLURM jobs)? [y/N] " ans
  [[ "$ans" == [yY] ]] || { echo "[i] Aborted — nothing submitted."; exit 0; }
fi

submitted=0
for dest in "${SAMPLE_DIRS[@]}"; do
  sample="$(basename "$dest")"
  # skip done / unstaged
  printf '%s\n' "${SETUP_OK[@]}" | grep -qxF "$sample" || continue
  echo
  echo "=== $sample :: $dest ==="
  bash "$LAUNCHER" "$dest"
  submitted=$(( submitted + 1 ))
done

echo
echo "[i] launcher invoked for $submitted sample(s)."
echo "[i] monitor:  squeue -u \$USER"
echo "[i] per-sample logs:  $WORK_ROOT/Batch*/<SAMPLE>/log/"

# <\SUBMIT> ---------------------------------------------------------------------
#EOF
