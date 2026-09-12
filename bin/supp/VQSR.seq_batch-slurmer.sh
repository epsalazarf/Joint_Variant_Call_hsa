#!/usr/bin/env bash
# =============================================================================
# Title : Step 06 (VQSR / Hard-Filtering) SLURM launcher
# About : Submits ONE sbatch job running 06_gatk_vqsr.sh against the cohort's
#         joint VCF. Unlike Step 04/05, this step is per-cohort, not
#         per-chromosome (docs/PIPELINE_STATUS.md's step table already scopes
#         S06 that way) — there is nothing to scatter across separate jobs at
#         the chromosome level, so one job is the natural granularity, not a
#         simplification made for this launcher's sake.
#
#         NOT DONE HERE: vqsr mode's own two independent VariantRecalibrator
#         calls (SNP, INDEL) — and hard-filter mode's two independent
#         SelectVariants+VariantFiltration branches — could each run as two
#         parallel jobs before a final ApplyVQSR/merge job, the same
#         scatter-then-gather shape as Step 05's launcher. 06_gatk_vqsr.sh
#         does not expose a phase switch for that (no GENO_PHASE-style
#         mechanism), so it isn't wired up here; revisit once Step 06 has run
#         for real and it's clear that split is worth the added complexity.
#
# Usage :
#   bash VQSR.seq_batch-slurmer.sh [options] <joint_vcf> <output_path> [mode] [cohort_name]
#   bash VQSR.seq_batch-slurmer.sh --dry-run [options] <joint_vcf> <output_path> [mode] [cohort_name]
#   bash VQSR.seq_batch-slurmer.sh report <manifest_file>
#
#   --cpus  N     --cpus-per-task (default: 2 — GATK's VQSR/filtering tools used
#                 here are not internally multi-threaded)
#   --mem   SIZE  --mem (default: 16G — PLACEHOLDER, unmeasured; see Step 06's
#                 STATUS block. A whole-genome joint VCF is a single file (no
#                 per-chromosome split the way Step 04/05's GenomicsDB/VCF
#                 outputs are), so this may need to scale with cohort size and
#                 total variant count in ways there is no FENIX data for yet)
#   --hours H     --time in hours (default: 24 — also a placeholder)
#
#   mode          vqsr | hard-filter | auto (forwarded to 06_gatk_vqsr.sh; default: auto)
#   joint_vcf     the gathered <cohort>.joint.vcf.gz from Step 05
#   output_path   where Step 06 writes vqsr_work/ or hardfilter_work/, logs/, and
#                 the final <cohort>.filtered.vcf.gz
#   cohort_name   optional; forwarded as-is (06_gatk_vqsr.sh derives its own
#                 default from joint_vcf's filename when omitted)
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S06="$REPO/bin/06_gatk_vqsr.sh"

# ---------- report mode (before option parsing) ------------------------------
if [[ "${1:-}" == "report" ]]; then
  manifest="${2:?Usage: $(basename "$0") report <manifest_file>}"
  [[ -s "$manifest" ]] || { echo "[X]  manifest not found: $manifest"; exit 1; }
  echo "[i]  manifest: $manifest"; echo
  jids=$(awk -F'\t' '/^vqsr/{print $2}' "$manifest" | paste -sd, -)
  [[ -n "$jids" ]] || { echo "[X]  no job ids in manifest"; exit 1; }
  if command -v sacct >/dev/null; then
    sacct -j "$jids" --units=G \
      --format=JobID%14,JobName%24,State%12,Elapsed%12,TotalCPU%12,MaxRSS%10,AllocCPUS%9,ReqMem%8,Start%20,End%20
  else
    echo "[!]  sacct not available; raw manifest:"; cat "$manifest"
  fi
  exit 0
fi

fmt_time() { printf '%d:00:00' "$1"; }   # H -> H:00:00

# ---------- options ------------------------------------------------------------
CPUS=2 ; MEM="16G" ; HOURS=24 ; DRY=false
pos=()
while (( $# )); do
  case "$1" in
    --dry-run)  DRY=true; shift ;;
    --cpus)     CPUS="${2:?}"; shift 2 ;;
    --mem)      MEM="${2:?}"; shift 2 ;;
    --hours)    HOURS="${2:?}"; shift 2 ;;
    -h|--help)  sed -n '3,39p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
    -*)         echo "[X]  unknown option: $1  (--help)"; exit 2 ;;
    *)          pos+=("$1"); shift ;;
  esac
done
JOINT_VCF="${pos[0]:?Usage: $(basename "$0") [options] <joint_vcf> <output_path> [mode] [cohort_name]}"
OUTPUT_PATH="${pos[1]:?Usage: $(basename "$0") [options] <joint_vcf> <output_path> [mode] [cohort_name]}"
MODE="${pos[2]:-auto}"
COHORT="${pos[3]:-}"

case "$MODE" in
  vqsr|hard-filter|auto) ;;
  *) echo "[X]  Invalid mode: '$MODE' (expected: vqsr | hard-filter | auto)"; exit 1 ;;
esac

# ---------- checks --------------------------------------------------------
[[ -f "$S06" ]] || { echo "[X]  Step 06 script not found: $S06"; exit 1; }
[[ -s "$JOINT_VCF" ]] || { echo "[X]  Joint VCF not found or empty: $JOINT_VCF"; exit 1; }
JOINT_VCF="$(readlink -f "$JOINT_VCF")"
$DRY || command -v sbatch >/dev/null || { echo "[X]  sbatch not found"; exit 1; }
mkdir -p "$OUTPUT_PATH"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"
LOG_DIR="${OUTPUT_PATH}/logs"
mkdir -p "$LOG_DIR"

COHORT_LABEL="${COHORT:-$(basename "$JOINT_VCF" | sed -E 's/\.vcf\.gz$//; s/\.joint$//')}"
[[ -n "$COHORT_LABEL" ]] || COHORT_LABEL="cohort"

# ---------- submit -------------------------------------------------------------
TS="$(date +%Y%m%d-%H%M%S)"
MANIFEST="$LOG_DIR/vqsr_${COHORT_LABEL}_run_${TS}.manifest"
{
  echo -e "# Step 06 (VQSR/hard-filter) run\t$TS"
  echo -e "# joint_vcf\t$JOINT_VCF"
  echo -e "# output_path\t$OUTPUT_PATH"
  echo -e "# mode\t$MODE"
  echo -e "# cohort\t$COHORT_LABEL"
  echo -e "# s06\t$S06"
  echo -e "# columns: step\tjobid\tcpus\tmem\ttime\tlogfile"
} > "$MANIFEST"

echo "[&]  ${COHORT_LABEL} — Step 06 launcher   ($TS)"
echo "[i]  joint_vcf   : $JOINT_VCF"
echo "[i]  output_path : $OUTPUT_PATH"
echo "[i]  mode        : $MODE"
echo "[i]  cohort      : $COHORT_LABEL"
echo "[i]  job         : ${CPUS} CPU / ${MEM} / ${HOURS}h"
echo "[i]  logs        : $LOG_DIR"
echo "[i]  manifest    : $MANIFEST"
echo "[i]  dry-run     : $DRY"
echo

jobname="JVC-VQSR-${COHORT_LABEL}"
logfile="$LOG_DIR/vqsr-${COHORT_LABEL}-%j.out"
cohort_arg=""
[[ -n "$COHORT" ]] && cohort_arg=" '$COHORT'"
wrap="bash '$S06' '$JOINT_VCF' '$OUTPUT_PATH' '$MODE'${cohort_arg}"

set -- --parsable --job-name="$jobname" \
       --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --mem="$MEM" --time="$(fmt_time "$HOURS")" \
       --output="$logfile"

echo "[*]  vqsr  (mode: $MODE)"

if $DRY; then
  echo "       sbatch $* --wrap \"$wrap\""
  jid="DRYRUN-vqsr"
else
  jid=$(sbatch "$@" --wrap "$wrap")
  echo "       submitted: job $jid   log: ${logfile/\%j/$jid}"
fi
echo -e "vqsr\t$jid\t${CPUS}\t${MEM}\t$(fmt_time "$HOURS")\t${logfile/\%j/$jid}" >> "$MANIFEST"

echo
echo "[i]  manifest : $MANIFEST"
echo "[i]  report   : bash '$0' report '$MANIFEST'"

#EOF
