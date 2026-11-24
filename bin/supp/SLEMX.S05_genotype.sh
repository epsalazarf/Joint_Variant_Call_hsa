#!/usr/bin/env bash

# =============================================================================
# Title       : SLEmx (lupus) — Step 05 joint genotyping launcher [FENIX]
# Description : Checks the Step 04 GenomicsDB workspaces and submits Step 05
#               (bin/05_gatk_GenotypeGVCFs.sh), one SLURM job per chromosome.
#
#               --check (no submission, login node, ~1 min): for every
#               autosome workspace, reports complete / sample count vs the
#               frozen S04 map / same sample names on every chromosome /
#               fragment count / group-readable / size, and lists leftover
#               Step 04 scratch dirs. Fails loudly on anything wrong.
#
#               Submit: --chroms chr22,chr1 (required — no "all" by accident).
#               The same check runs first for those chromosomes. A chromosome
#               whose VCF (+ .tbi) exists, or whose job is already queued, is
#               skipped, so re-running after a failure is safe.
#
#               --region chr22:20000001-20500000 runs a small piece of ONE
#               chromosome (minutes): a cheap read test before the long jobs.
#
# Usage       : bash SLEMX.S05_genotype.sh --check
#               bash SLEMX.S05_genotype.sh --chroms chr22 --region chr22:20000001-20500000
#               bash SLEMX.S05_genotype.sh --chroms chr22,chr1 [--dry-run]
#               bash SLEMX.S05_genotype.sh --status
#                 --out DIR   cohort root (default /mnt/data/amedina/esalazarf/SLEmx_JVC)
#                             reads <out>/genomicsdb/<chrom>, writes <out>/genotyped/
# =============================================================================

set -euo pipefail

# <CONFIG> ----------------------------------------------------------------
OUTPUT_ROOT="/mnt/data/amedina/esalazarf/SLEmx_JVC"
TAG="slemx"
DATA_GROUP="amedina"
AUTOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
# Resources. GenotypeGVCFs is single-threaded; time scales with contig length x
# samples. No measurement yet for 164 samples — these are generous first
# guesses (defq has no time limit; an early finish frees the slot). Check
# sacct MaxRSS/Elapsed after chr22 and chr1 and tighten. Job memory = 2x Java
# heap: GenomicsDB reads through native memory too.
JOB_MEM="32G";  JAVA_MEM="12G"
REGION_MEM="16G"; REGION_JAVA_MEM="8G"; REGION_HOURS=12
hours_for() {   # $1 chrom -> walltime hours
  case "$1" in
    chr1|chr2) echo 168 ;;
    chr3|chr4|chr5|chr6|chr7|chr8) echo 120 ;;
    *) echo 72 ;;
  esac
}
# <\CONFIG> -----------------------------------------------------------------

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S05="$REPO/bin/05_gatk_GenotypeGVCFs.sh"
CONFIG_FILE="$REPO/config/config.yaml"

MODE="submit"; CHROMS=""; REGION=""
while [ $# -gt 0 ]; do
  case "$1" in
    --check)   MODE="check";  shift ;;
    --status)  MODE="status"; shift ;;
    --dry-run) MODE="dry";    shift ;;
    --chroms)  CHROMS="$(echo "${2:?--chroms needs a list, e.g. chr22,chr1}" | tr ',' ' ')"; shift 2 ;;
    --region)  REGION="${2:?--region needs chrN:start-end}"; shift 2 ;;
    --out)     OUTPUT_ROOT="${2:?--out needs a directory}"; shift 2 ;;
    -h|--help) sed -n '4,27p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)         echo "[X]  unknown argument: $1  (--help)"; exit 2 ;;
  esac
done
OUTPUT_ROOT="${OUTPUT_ROOT%/}"
GDB_ROOT="$OUTPUT_ROOT/genomicsdb"
VCF_DIR="$OUTPUT_ROOT/genotyped"
LOG_DIR="$OUTPUT_ROOT/logs"
FROZEN_MAP="$OUTPUT_ROOT/${TAG}.sample_map.S04.tsv"

# Output name, same rule as 05_gatk_GenotypeGVCFs.sh (prefix TAG).
vcf_for() {   # $1 chrom [$2 region]
  if [[ -n "${2:-}" ]]; then local se="${2#*:}"; echo "$VCF_DIR/${TAG}.$1_${se}.vcf.gz"
  else echo "$VCF_DIR/${TAG}.$1.vcf.gz"; fi
}
vcf_done() { [[ -s "$1" && -s "$1.tbi" ]]; }
jobname_for() { if [[ -n "${2:-}" ]]; then echo "JVC-GTGVCF-${TAG}-$1-region"; else echo "JVC-GTGVCF-${TAG}-$1"; fi; }
queued() { command -v squeue >/dev/null && squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -qx "$1"; }

# --- status -----------------------------------------------------------------
if [[ "$MODE" == status ]]; then
  echo "[&]  ${TAG} Step 05 status — $VCF_DIR"
  for c in $AUTOSOMES; do
    f="$(vcf_for "$c")"
    if vcf_done "$f"; then printf '  [DONE]    %-6s %s\n' "$c" "$(du -h "$f" | cut -f1)"
    elif queued "$(jobname_for "$c")"; then printf '  [QUEUED]  %s\n' "$c"
    else printf '  [ ---- ]  %s\n' "$c"; fi
  done
  shopt -s nullglob; for f in "$VCF_DIR"/"${TAG}".chr*_*-*.vcf.gz; do echo "  [region]  $(basename "$f")"; done; shopt -u nullglob
  command -v squeue >/dev/null && squeue -u "$USER" -h -o '%.18i %.34j %.8T %.12M %R' | grep "JVC-GTGVCF-${TAG}" || true
  exit 0
fi

# --- what to check / submit ---------------------------------------------------
if [[ "$MODE" == check ]]; then
  CHECK_CHROMS="$AUTOSOMES"
else
  [[ -n "$CHROMS" ]] || { echo "[X]  --chroms is required to submit (e.g. --chroms chr22,chr1). Use --check to only check."; exit 2; }
  for c in $CHROMS; do
    case " $AUTOSOMES " in *" $c "*) ;; *) echo "[X]  not an autosome: $c"; exit 2 ;; esac
  done
  if [[ -n "$REGION" ]]; then
    set -- $CHROMS
    (( $# == 1 )) || { echo "[X]  --region works with exactly one chromosome in --chroms"; exit 2; }
    [[ "$REGION" =~ ^$1:[0-9]+-[0-9]+$ ]] || { echo "[X]  --region must be $1:<start>-<end> (got '$REGION')"; exit 2; }
  fi
  CHECK_CHROMS="$CHROMS"
fi

# <CHECKS> ------------------------------------------------------------------
ERR=0
fail() { echo "[X]  $*"; ERR=$((ERR + 1)); }
warn() { echo "[!]  $*"; }

[[ -f "$S05" ]] || { echo "[X]  Step 05 script not found: $S05"; exit 1; }
[[ -d "$GDB_ROOT" ]] || { echo "[X]  No Step 04 output: $GDB_ROOT"; exit 1; }
N_MAP=""
if [[ -s "$FROZEN_MAP" ]]; then N_MAP=$(grep -cvE '^[[:space:]]*(#|$)' "$FROZEN_MAP" || true)
else warn "frozen S04 map not found ($FROZEN_MAP) — sample counts only compared across chromosomes"; fi

echo "[&]  ${TAG} Step 04 workspace check   ($(date))"
echo "[i]  workspaces : $GDB_ROOT/<chrom>"
echo "[i]  S04 map    : ${N_MAP:-?} samples"
echo
printf '  %-6s %-8s %7s %5s %9s  %s\n' chrom state samples frags size note
ref_names_ck=""; ref_chrom=""; total_kb=0
for c in $CHECK_CHROMS; do
  ws="$GDB_ROOT/$c"; note=""; err0=$ERR
  if [[ ! -s "$ws/callset.json" || ! -s "$ws/vidmap.json" ]]; then
    printf '  %-6s %-8s\n' "$c" MISSING; fail "$c: no complete workspace (callset.json + vidmap.json)"; continue
  fi
  names="$(tr -d '\n' < "$ws/callset.json" | grep -oE '"sample_name"[[:space:]]*:[[:space:]]*"[^"]*"' | sed -E 's/.*"([^"]*)"$/\1/' | sort)"
  n=$(printf '%s\n' "$names" | grep -c . || true)
  dups="$(printf '%s\n' "$names" | uniq -d | paste -sd' ' -)"
  ck="$(printf '%s\n' "$names" | cksum | cut -d' ' -f1)"
  shopt -s nullglob; arr=( "$ws/$c"\$* ); shopt -u nullglob
  frags="?"
  if (( ${#arr[@]} == 1 )); then
    frags=$(find "${arr[0]}" -mindepth 1 -maxdepth 1 -type d -name '__*' | wc -l | tr -d '[:space:]')
  else
    fail "$c: expected one array dir '$c\$1\$<length>' in $ws, found ${#arr[@]}"
  fi
  kb=$(du -sk "$ws" 2>/dev/null | cut -f1); total_kb=$((total_kb + ${kb:-0}))
  [[ -n "$N_MAP" && "$n" != "$N_MAP" ]] && fail "$c: $n samples in workspace, S04 map has $N_MAP"
  [[ -n "$dups" ]] && fail "$c: duplicate sample names: $dups"
  if [[ -z "$ref_names_ck" ]]; then ref_names_ck="$ck"; ref_chrom="$c"
  elif [[ "$ck" != "$ref_names_ck" ]]; then fail "$c: sample names differ from $ref_chrom (chromosomes must hold the same cohort)"; note="names differ"; fi
  [[ "$frags" == 0 ]] && fail "$c: array has no fragments (empty workspace?)"
  if [[ -n "$(find "$ws" ! -perm -g=r -print -quit 2>/dev/null)" ]]; then
    warn "$c: some files are not group-readable — fix: chgrp -R $DATA_GROUP '$ws' && chmod -R g+rX '$ws'"; note="${note:+$note, }not group-readable"
  fi
  state=OK; (( ERR > err0 )) && state=BAD
  printf '  %-6s %-8s %7s %5s %8sG  %s\n' "$c" "$state" "$n" "$frags" "$(( ${kb:-0} / 1024 / 1024 ))" "$note"
done
echo
echo "[i]  total size : $(( total_kb / 1024 / 1024 )) GB   (fragments: 1 = consolidated; >1 is fine, Step 05 reads both)"
odd="$(printf '%s\n' "${names:-}" | grep '_' | head -3 | paste -sd' ' - || true)"
[[ -n "$odd" ]] && echo "[i]  names with '_' (header name, not map label): $odd — expected for Q_040"

# Leftover Step 04 scratch (an OOM kill skips the cleanup trap).
scr="/scratch/groups/$(id -gn)/$USER"
shopt -s nullglob; left=( "$scr"/gendbi_* ); shopt -u nullglob
if (( ${#left[@]} )); then
  warn "${#left[@]} leftover Step 04 scratch dir(s) in $scr (safe to delete once no S04 job runs):"
  for d in "${left[@]}"; do echo "       $(du -sh "$d" 2>/dev/null | cut -f1)  $d"; done
  echo "       rm -rf $scr/gendbi_*"
fi

echo
if (( ERR )); then echo "[X]  Check FAILED: $ERR problem(s) above. Nothing submitted."; exit 1; fi
echo "[OK] ${CHECK_CHROMS// /, }: workspaces look complete and consistent."
[[ "$MODE" == check ]] && { echo "[i]  --check: nothing submitted."; exit 0; }
# <\CHECKS> -----------------------------------------------------------------

# --- submit -------------------------------------------------------------------
if [[ "$MODE" != dry ]] && ! command -v sbatch >/dev/null; then
  echo "[X]  sbatch not found — run this on a FENIX login node."; exit 1
fi
mkdir -p "$VCF_DIR" "$LOG_DIR"
exclude="$(awk '/^remote:/{r=1;next} /^[^[:space:]]/{r=0} r && $1=="sbatch_exclude:"{gsub(/"/,"",$2); print $2}' "$CONFIG_FILE")"

echo
echo "[&]  ${TAG} Step 05 — GenotypeGVCFs   (output: $VCF_DIR)"
for c in $CHROMS; do
  out="$(vcf_for "$c" "$REGION")"; jn="$(jobname_for "$c" "$REGION")"
  if vcf_done "$out"; then echo "[SKIP] $c — $(basename "$out") already exists"; continue; fi
  if queued "$jn"; then echo "[SKIP] $c — job $jn already queued/running"; continue; fi
  if [[ -n "$REGION" ]]; then
    mem="$REGION_MEM"; jmem="$REGION_JAVA_MEM"; hours="$REGION_HOURS"; label="$REGION"
    log="$LOG_DIR/${TAG}-S05-${c}-region-%j.out"
  else
    mem="$JOB_MEM"; jmem="$JAVA_MEM"; hours="$(hours_for "$c")"; label="$c"
    log="$LOG_DIR/${TAG}-S05-${c}-%j.out"
  fi
  wrap="
set -euo pipefail
umask 0027
_scr=\"/scratch/groups/\$(id -gn)\"
mkdir -p \"\$_scr/\$USER\" || { echo \"<ERROR> scratch not writable on \$(hostname): \$_scr/\$USER\" >&2; exit 1; }
echo \"[i] ${label} on node \$(hostname), temp in \$_scr\"
env GTGVCF_SCRATCH_ROOT=\"\$_scr\" GTGVCF_PREFIX='${TAG}' GTGVCF_JAVA_MEM='${jmem}' \\
  bash '$S05' '$GDB_ROOT' '$VCF_DIR' '$c' '$REGION'
chgrp '$DATA_GROUP' '$out' '$out.tbi' 2>/dev/null && chmod g+r '$out' '$out.tbi' || true
"
  args=(--parsable --job-name="$jn" --nodes=1 --ntasks=1 --cpus-per-task=2
        --mem="$mem" --time="${hours}:00:00" --output="$log")
  [[ -n "$exclude" ]] && args+=(--exclude="$exclude")
  echo "[*]  $label   $mem / ${hours}h   -> $(basename "$out")"
  if [[ "$MODE" == dry ]]; then
    echo "     sbatch ${args[*]} --wrap \"...\""
  else
    jid=$(sbatch "${args[@]}" --wrap "$wrap")
    echo "     submitted: job $jid   log: ${log/\%j/$jid}"
  fi
done

echo
if [[ "$MODE" == dry ]]; then echo "[i]  --dry-run: nothing submitted."
else
  echo "[i]  Monitor : squeue -u \$USER | grep JVC-GTGVCF-${TAG}"
  echo "[i]  Status  : bash $0 --status"
  echo "[i]  Progress: tail -f the log — GATK prints its position every 10 min."
  echo "[i]  After it ends: sacct -j <jobid> -o JobID,State,Elapsed,MaxRSS  (to size the other chromosomes)"
fi

#EOF
