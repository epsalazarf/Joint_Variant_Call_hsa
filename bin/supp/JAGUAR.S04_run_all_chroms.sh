#!/usr/bin/env bash

# =============================================================================
# Title       : JAGUAR — Step 04 full-cohort import, all chromosomes [FENIX]
# Description : One-shot launcher for JAGUAR's first full cohort merge into
#               GenomicsDB. Builds the sample map from the fixed cohort folder
#               below, then submits ONE SLURM array PER CHROMOSOME-SIZE CLASS
#               (chr1/chr2, chr3-8+chrX, the rest) — each chromosome writes an
#               independent workspace, so all 25 run in parallel instead of
#               serially. `action=create`, single wave: the cohort is fixed at
#               93 samples for this run, no update wave needed.
#
#               Meant to be started ONCE and left unattended (e.g. over a
#               weekend). Bounded by the slowest class (chr1/chr2, ~20h from
#               the chr1 test — see docs/S04_GenomicsDBImport_design.md), not
#               by the sum of all 25 chromosomes (~8 days if run serially).
#
#               Safe to re-run: Step 04 skips a chromosome whose workspace
#               already has callset.json + vidmap.json, so if some tasks fail
#               (bad node, hit the wall) just run this script again — only the
#               missing chromosomes get resubmitted.
#
# Usage       : bash JAGUAR.S04_run_all_chroms.sh [--dry-run] [--status]
#               No positional arguments — cohort + output are fixed below for
#               this run. Edit the <CONFIG> block if either path changes.
#                 --dry-run   print the sbatch commands, submit nothing
#                 --status    report which chromosome workspaces exist vs are
#                             still missing, submit nothing
#
# Coordinate with the maintainer before launching — this is a cohort-level
# step, run once after ALL samples have completed Step 03.
# =============================================================================

set -euo pipefail

# <CONFIG> ----------------------------------------------------------------
# Fixed for this run — JAGUAR cohort, 93 samples, all already through Step 03.
# OUTPUT_ROOT is the output_path Step 04 itself takes — it creates workspaces
# at <OUTPUT_ROOT>/genomicsdb/<chrom> on its own, so do NOT add /genomicsdb here.
COHORT_ROOT="/mnt/data/amedina/mramirezc/JAGUAR_JVC/bam"
OUTPUT_ROOT="/mnt/data/amedina/${USER:-mramirezc}/JAGUAR_JVC"
TAG="jaguar"
ACTION="create"     # single wave: all 93 samples go in at once

CANON_CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
              chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
              chrX chrY chrM)

# Resource classes (see docs/S04_GenomicsDBImport_design.md — memory/time scale
# with CONTIG SIZE, not sample count. --cpus stays 2: GenomicsDBImport is
# serial + I/O-bound, more cores don't help.) ARRAY_CONC caps how many tasks
# in that class's array run at once, to keep concurrent NFS reads + BeeGFS
# scratch writes from stepping on each other. Lower it if the cluster looks
# saturated; the array will just take longer, not fail.
declare -A CLASS_CHROMS=(
  [big]="chr1 chr2"
  [mid]="chr3 chr4 chr5 chr6 chr7 chr8 chrX"
  [small]="chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrY chrM"
)
declare -A CLASS_MEM=([big]="64G" [mid]="32G" [small]="16G")
declare -A CLASS_HOURS=([big]=20 [mid]=12 [small]=6)
declare -A CLASS_CONC=([big]=2 [mid]=4 [small]=6)
# <\CONFIG> -----------------------------------------------------------------

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S04="$REPO/bin/04_gatk_GenomicsDB_import.sh"
SAMPLE_MAP="${OUTPUT_ROOT}/${TAG}.sample_map.tsv"
LOG_DIR="${OUTPUT_ROOT}/logs"

DRY=false
STATUS=false
for a in "$@"; do
  case "$a" in
    --dry-run) DRY=true ;;
    --status)  STATUS=true ;;
    -h|--help) sed -n '3,32p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
    *)         echo "[X]  unknown argument: $a  (--help)"; exit 2 ;;
  esac
done

# --- status mode: report and exit, no submission ---------------------------
if $STATUS; then
  echo "[&]  ${TAG} Step 04 status  — ${OUTPUT_ROOT}/genomicsdb/<chrom>"
  echo
  done_n=0
  for c in "${CANON_CHROMS[@]}"; do
    ws="${OUTPUT_ROOT}/genomicsdb/${c}"
    if [[ -s "${ws}/callset.json" && -s "${ws}/vidmap.json" ]]; then
      printf '  [DONE]    %s\n' "$c"; (( ++done_n ))
    else
      printf '  [MISSING] %s\n' "$c"
    fi
  done
  echo
  echo "[i]  ${done_n}/${#CANON_CHROMS[@]} chromosome workspaces complete."
  (( done_n < ${#CANON_CHROMS[@]} )) && echo "[i]  Re-run without --status to submit the missing ones."
  exit 0
fi

# <CHECKS> --------------------------------------------------------------------

[[ -f "$S04" ]] || { echo "[X]  Step 04 script not found: $S04"; exit 1; }
[[ -d "$COHORT_ROOT" ]] || { echo "[X]  cohort_root not a directory: $COHORT_ROOT"; exit 1; }
$DRY || command -v sbatch >/dev/null || { echo "[X]  sbatch not found — SLURM environment required."; exit 1; }

case "${OUTPUT_ROOT}/" in
  "${HOME}"/*)
    echo "[X]  output_root is under \$HOME ($HOME) — FENIX home has a hard quota."
    echo "[i]  Edit OUTPUT_ROOT in this script's <CONFIG> block to group storage."
    exit 1 ;;
esac

mkdir -p "$(dirname "$SAMPLE_MAP")" "$LOG_DIR"

# <\CHECKS> ---------------------------------------------------------------------

# --- Step 1: sample map (built once, reused for every chromosome) ----------
# One map covers ALL chromosomes: Step 04 resolves the per-chromosome GVCF by
# globbing *.<chrom>.g.vcf.gz inside each sample's chrom_gvcf/ dir, so the map
# only needs to point at that directory (see 04_gatk_GenomicsDB_import.sh).
if [[ -s "$SAMPLE_MAP" ]]; then
  n_samples=$(grep -cvE '^[[:space:]]*(#|$)' "$SAMPLE_MAP")
  echo "[i]  Using existing sample map: $SAMPLE_MAP  (${n_samples} samples)"
  echo "[i]  Delete it first if the cohort folder has changed and you want it rebuilt."
else
  tmp="${SAMPLE_MAP}.tmp"
  : > "$tmp"
  for d in "$COHORT_ROOT"/*/chrom_gvcf; do
    [[ -d "$d" ]] || continue
    s="$(basename "$(dirname "$d")")"
    printf '%s\t%s\n' "$s" "$(readlink -f "$d")" >> "$tmp"
  done
  n_samples=$(wc -l < "$tmp" | tr -d '[:space:]')   # BSD wc pads its count with spaces
  if (( n_samples == 0 )); then
    rm -f "$tmp"
    echo "[X]  No */chrom_gvcf directories found under $COHORT_ROOT"
    exit 1
  fi
  mv "$tmp" "$SAMPLE_MAP"
  echo "[i]  Built sample map: $SAMPLE_MAP  (${n_samples} samples)"
fi
echo "[i]  First 3 lines:"
head -3 "$SAMPLE_MAP" | sed 's/^/       /'
echo

# --- Step 2: submit one array per resource class ----------------------------

echo "[&]  ${TAG} Step 04 — full cohort, all chromosomes   ($(date))"
echo "[i]  cohort_root : $COHORT_ROOT"
echo "[i]  output_root : $OUTPUT_ROOT   (workspaces -> \${output_root}/genomicsdb/<chrom>)"
echo "[i]  sample_map  : $SAMPLE_MAP  (${n_samples} samples)"
echo "[i]  action      : $ACTION"
echo "[i]  dry-run     : $DRY"
echo

submit_class() {
  local class="$1"
  local all_chrs=(${CLASS_CHROMS[$class]})
  local mem="${CLASS_MEM[$class]}" hours="${CLASS_HOURS[$class]}" conc="${CLASS_CONC[$class]}"

  # Only submit chromosomes that don't already have a complete workspace, so a
  # re-run after a partial failure resubmits just the stragglers instead of
  # re-queueing (and re-skipping, inside Step 04) everything in the class.
  local chrs=()
  local c
  for c in "${all_chrs[@]}"; do
    local ws="${OUTPUT_ROOT}/genomicsdb/${c}"
    [[ -s "${ws}/callset.json" && -s "${ws}/vidmap.json" ]] && continue
    chrs+=("$c")
  done
  if (( ${#chrs[@]} == 0 )); then
    echo "[SKIP] class=${class} — all ${#all_chrs[@]} chromosome(s) already built"
    return 0
  fi

  local n=${#chrs[@]}
  local jobname="JVC-GDBI-${TAG}-${class}"
  local logfile="${LOG_DIR}/${TAG}-${class}-%A_%a.out"
  local wrap="
set -euo pipefail
_chrs=(${chrs[@]})
_chrom=\"\${_chrs[\$((SLURM_ARRAY_TASK_ID - 1))]:-}\"
[[ -z \"\$_chrom\" ]] && { echo \"<ERROR> no chromosome for array index \$SLURM_ARRAY_TASK_ID\" >&2; exit 1; }
echo \"[i] array task \$SLURM_ARRAY_TASK_ID -> \$_chrom  (class ${class})\"
env GENDBI_READER_THREADS=2 GENDBI_CONSOLIDATE=true \
  bash '$S04' '$SAMPLE_MAP' '$OUTPUT_ROOT' \"\$_chrom\" '$ACTION'
"

  echo "[*]  class=${class}  chroms=(${chrs[*]})  ${mem} / ${hours}h  %${conc} concurrent"

  if $DRY; then
    echo "     sbatch --job-name=$jobname --array=1-${n}%${conc} --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=$mem --time=${hours}:00:00 --output=$logfile --wrap \"...\""
  else
    local jid
    jid=$(sbatch --parsable \
      --job-name="$jobname" \
      --array="1-${n}%${conc}" \
      --nodes=1 --ntasks=1 --cpus-per-task=2 \
      --mem="$mem" --time="${hours}:00:00" \
      --output="$logfile" \
      --wrap "$wrap")
    echo "     submitted: job $jid   log: ${logfile/\%A_\%a/${jid}_*}"
  fi
}

for class in big mid small; do
  submit_class "$class"
done

echo
if $DRY; then
  echo "[i]  dry run — nothing submitted."
else
  echo "[i]  Monitor:   squeue -u \$USER | grep JVC-GDBI-${TAG}"
  echo "[i]  Status:    bash $(basename "$0") --status"
  echo "[i]  Logs:      ${LOG_DIR}/${TAG}-<class>-<jobid>_<task>.out"
  echo "[i]  A failed task leaves that chromosome's workspace missing; re-run"
  echo "[i]  this script (no flags) once the cause is fixed — it only submits"
  echo "[i]  what --status still shows as MISSING."
  echo
  echo "[i]  Once all 25 report DONE, verify one workspace is readable before"
  echo "[i]  moving to Step 05, e.g.:"
  echo "[i]    gatk SelectVariants -R <ref_gnm> -V gendb://${OUTPUT_ROOT}/genomicsdb/chr22 \\"
  echo "[i]      -L chr22:1-200000 -O /tmp/check_chr22.vcf.gz"
  echo "[i]    bcftools query -l /tmp/check_chr22.vcf.gz | wc -l   # expect ${n_samples}"
fi

#EOF
