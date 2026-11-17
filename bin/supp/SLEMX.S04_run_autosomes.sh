#!/usr/bin/env bash

# =============================================================================
# Title       : SLEmx (lupus) — Step 04 cohort import, autosomes [FENIX]
# Description : Launcher for the lupus cohort's first GenomicsDB build. Takes
#               the curated cohort sample map (SLEmx_JVC.sample_map.2026-10-23,
#               164 samples), checks it hard, freezes a copy, then submits ONE
#               SLURM array PER CHROMOSOME-SIZE CLASS over chr1..chr22. Each
#               chromosome writes its own workspace, so they run in parallel.
#               action=create, single wave (all samples at once).
#
#               MERGED SAMPLES — the 11 duplicate IDs (two sequencing batches
#               each) were re-merged at Step 02 into BAMQC/Batch25/<ID>. Each
#               must enter the database ONCE, from Batch25. Pre-flight fails if
#               a merged ID is missing, does not point at Batch25, or if one of
#               its original per-batch copies (listed in the IDCHECK pairs
#               file) is still in the map. Excluded samples (Qro113, Qro115,
#               Q038, Q002) also fail the check.
#
#               Pre-flight also checks: unique labels and dirs, 22 autosomal
#               GVCFs + index per sample, one unique sample name per GVCF
#               header (chr22, needs bcftools), free space, no job already
#               queued.
#
#               Safe to re-run: a chromosome whose workspace has callset.json +
#               vidmap.json is not resubmitted, so after a failure (bad node,
#               walltime) run the script again and only the missing ones go.
#
# Usage       : bash SLEMX.S04_run_autosomes.sh [--map FILE] [--out DIR]
#                                                [--build-map | --check | --dry-run | --status]
#                 --build-map regenerate the cohort map from BAMQC/Batch*/*/chrom_gvcf
#                             with the 2026-10-23 selection rule (all dirs, minus
#                             the original copies of the merged IDs, minus
#                             EXCLUDED), write it to the --map path, stop.
#                             Review it, then run --check.
#                 --map FILE  cohort sample map (label <TAB> sample or chrom_gvcf
#                             dir). Default: <out>/SLEmx_JVC.sample_map.tsv
#                 --out DIR   output root (default /mnt/data/amedina/esalazarf/SLEmx_JVC)
#                             -> <out>/genomicsdb/chr1..chr22, <out>/logs/
#                 --check     pre-flight only, submit nothing
#                 --dry-run   pre-flight + print sbatch commands, submit nothing
#                 --status    which chromosome workspaces exist / are missing
#
# Run it from a FENIX login node, after `git pull` in the production repo:
#   bash /mnt/data/amedina/esalazarf/Joint_Variant_Call_hsa/bin/supp/SLEMX.S04_run_autosomes.sh --check
# =============================================================================

set -euo pipefail

# <CONFIG> ----------------------------------------------------------------
OUTPUT_ROOT="/mnt/data/amedina/esalazarf/SLEmx_JVC"
MAP_NAME="SLEmx_JVC.sample_map.tsv"   # same selection as SLEmx_JVC.sample_map.2026-10-23.tsv
BAMQC_ROOT="${SLEMX_BAMQC_ROOT:-/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC}"   # env override: local tests only
MERGED_BATCH="Batch25"           # Step 02 re-merge of the duplicate IDs
EXPECTED_N=164                   # warn (not fail) if the map has a different count
EXCLUDED=(Qro113 Qro115 Q038 Q002)
TAG="slemx"
ACTION="create"
MIN_FREE_GB="${SLEMX_MIN_FREE_GB:-1500}"                # 22 autosomes ~0.7-0.85 TiB, plus headroom
DATA_GROUP="amedina"             # group that must be able to read the database

AUTOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

# Resource classes: memory/time scale with CONTIG SIZE (docs/S04_GenomicsDBImport_design.md).
# 2 CPUs: GenomicsDBImport is serial + I/O-bound. Walltimes are 2x the JAGUAR
# launcher's (48/24/12 h), because this cohort has 164 samples (JAGUAR: 93) and
# many are >20x. defq has no time limit; a job that ends early frees its slot.
declare -A CLASS_CHROMS=(
  [big]="chr1 chr2"
  [mid]="chr3 chr4 chr5 chr6 chr7 chr8"
  [small]="chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
)
#
# 2026-11-13 run (jobs 457714-16): big (64G) and chr21/22 passed. ALL 18 other
# tasks imported all 164 samples, then were OOM-killed in the final
# `--consolidate` pass (mid at 32G, small at 16G). Consolidation rewrites the
# whole array, and its memory exceeds the import. It is not needed:
# GenotypeGVCFs reads a multi-fragment workspace (164 samples / batch 50 = 4
# fragments). So mid/small now skip it (CLASS_CONSOL=false). chr1/2/21/22 stay
# as built (consolidated). Memory/time get headroom over the import alone
# (small import took up to ~20.5 h of the old 24 h wall).
declare -A CLASS_MEM=([big]="64G" [mid]="48G" [small]="24G")
declare -A CLASS_HOURS=([big]=96 [mid]=72 [small]=48)
declare -A CLASS_CONC=([big]=2 [mid]=4 [small]=6)
declare -A CLASS_CONSOL=([big]=true [mid]=false [small]=false)
# <\CONFIG> -----------------------------------------------------------------

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S04="$REPO/bin/04_gatk_GenomicsDB_import.sh"
PAIRS="${SLEMX_PAIRS:-$REPO/test/idcheck/SLEmx_duplicates.pairs.tsv}"
CONFIG_FILE="$REPO/config/config.yaml"

MAP_IN=""
MODE="submit"
while [ $# -gt 0 ]; do
  case "$1" in
    --map)     MAP_IN="${2:?--map needs a file}"; shift 2 ;;
    --out)     OUTPUT_ROOT="${2:?--out needs a directory}"; shift 2 ;;
    --build-map) MODE="build"; shift ;;
    --check)   MODE="check";  shift ;;
    --dry-run) MODE="dry";    shift ;;
    --status)  MODE="status"; shift ;;
    -h|--help) sed -n '4,45p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)         echo "[X]  unknown argument: $1  (--help)"; exit 2 ;;
  esac
done
OUTPUT_ROOT="${OUTPUT_ROOT%/}"
MAP_IN="${MAP_IN:-$OUTPUT_ROOT/$MAP_NAME}"
FROZEN_MAP="$OUTPUT_ROOT/${TAG}.sample_map.S04.tsv"   # what the jobs actually read
LOG_DIR="$OUTPUT_ROOT/logs"

ws_done() { [[ -s "$OUTPUT_ROOT/genomicsdb/$1/callset.json" && -s "$OUTPUT_ROOT/genomicsdb/$1/vidmap.json" ]]; }

# --- status mode ------------------------------------------------------------
if [[ "$MODE" == status ]]; then
  echo "[&]  ${TAG} Step 04 status — $OUTPUT_ROOT/genomicsdb/<chrom>"
  n=0
  for c in "${AUTOSOMES[@]}"; do
    if ws_done "$c"; then printf '  [DONE]    %s\n' "$c"; n=$((n + 1)); else printf '  [MISSING] %s\n' "$c"; fi
  done
  echo "[i]  ${n}/${#AUTOSOMES[@]} autosome workspaces complete."
  command -v squeue >/dev/null && squeue -u "$USER" -h -o '%.18i %.30j %.8T %.10M' | grep "JVC-GDBI-${TAG}" || true
  exit 0
fi

# --- build-map mode: regenerate the cohort map, then stop --------------------
# Rule (reproduces the 2026-10-23 selection from the 2026-10-07 full scan:
# 178 dirs - 22 originals - 3 excluded + 11 Batch25 merges = 164):
#   every BAMQC/Batch*/<sample>/chrom_gvcf, label = <sample> dir name,
#   minus both original per-batch dirs of each merged ID in the pairs file,
#   minus EXCLUDED labels. Batch25/<ID> is picked up by the scan itself.
if [[ "$MODE" == build ]]; then
  [[ -f "$PAIRS" ]] || { echo "[X]  Duplicate-ID pairs file not found: $PAIRS"; exit 1; }
  [[ -d "$BAMQC_ROOT" ]] || { echo "[X]  BAMQC root not found: $BAMQC_ROOT"; exit 1; }
  [[ -e "$MAP_IN" ]] && { echo "[X]  $MAP_IN already exists — move it away first, or pass --map NEW_FILE"; exit 1; }
  declare -A DROP=()
  while IFS=$'\t' read -r label dirA dirB _; do
    case "$label" in ''|'#'*|*_vs_*) continue ;; esac
    for o in "$dirA" "$dirB"; do DROP[$(readlink -f "$o" 2>/dev/null || echo "$o")]="$label"; done
  done < "$PAIRS"
  mkdir -p "$(dirname "$MAP_IN")"
  tmp="$MAP_IN.tmp"
  printf '# SLEmx S04 cohort map, rebuilt %s by %s --build-map\n# rule: BAMQC/Batch*/*/chrom_gvcf - merged-ID originals (%s) - EXCLUDED (%s)\n' \
    "$(date +%F)" "$(basename "$0")" "$(basename "$PAIRS")" "${EXCLUDED[*]}" > "$tmp"
  n=0; n_drop=0; n_excl=0
  declare -A PER_BATCH=()
  for d in "$BAMQC_ROOT"/Batch*/*/chrom_gvcf; do
    [[ -d "$d" ]] || continue
    sdir="$(dirname "$d")"; label="$(basename "$sdir")"; batch="$(basename "$(dirname "$sdir")")"
    rs="$(readlink -f "$sdir")"
    if [[ -n "${DROP[$rs]:-}" || -n "${DROP[$sdir]:-}" ]]; then
      echo "[-]  original copy of merged ${DROP[$rs]:-${DROP[$sdir]}}: $batch/$label"; n_drop=$((n_drop + 1)); continue
    fi
    if [[ " ${EXCLUDED[*]} " == *" $label "* ]]; then
      echo "[-]  excluded: $batch/$label"; n_excl=$((n_excl + 1)); continue
    fi
    printf '%s\t%s\n' "$label" "$(readlink -f "$d")" >> "$tmp"
    PER_BATCH[$batch]=$(( ${PER_BATCH[$batch]:-0} + 1 )); n=$((n + 1))
  done
  mv "$tmp" "$MAP_IN"
  echo
  echo "[i]  per batch: $(for b in $(printf '%s\n' "${!PER_BATCH[@]}" | sort -V); do printf '%s=%s ' "$b" "${PER_BATCH[$b]}"; done)"
  echo "[i]  dropped $n_drop merged-ID originals (expect 22), $n_excl excluded (expect 3)"
  dups="$(grep -v '^#' "$MAP_IN" | cut -f1 | sort | uniq -d | paste -sd' ' -)"
  [[ -n "$dups" ]] && echo "[!]  labels listed more than once (fix before --check): $dups"
  if (( n == EXPECTED_N )); then echo "[OK] wrote $MAP_IN — $n samples (matches the 2026-10-23 selection count)"
  else echo "[!]  wrote $MAP_IN — $n samples, EXPECTED $EXPECTED_N: new/removed dirs since 2026-10-23? Review before --check."; fi
  echo "[i]  Next: bash $0 --check"
  exit 0
fi

# <CHECKS> ------------------------------------------------------------------
ERR=0
fail() { echo "[X]  $*"; ERR=$((ERR + 1)); }
warn() { echo "[!]  $*"; }

echo "[&]  ${TAG} Step 04 pre-flight   ($(date))"
echo "[i]  sample map  : $MAP_IN"
echo "[i]  output root : $OUTPUT_ROOT"
echo

[[ -f "$S04" ]]   || { echo "[X]  Step 04 script not found: $S04"; exit 1; }
[[ -f "$PAIRS" ]] || { echo "[X]  Duplicate-ID pairs file not found: $PAIRS"; exit 1; }
[[ -s "$MAP_IN" ]] || { echo "[X]  Sample map not found or empty: $MAP_IN"; echo "[i]  Regenerate it: bash $0 --build-map   (or pass --map FILE)"; exit 1; }
case "$OUTPUT_ROOT/" in
  "$HOME"/*) echo "[X]  Output root is under \$HOME — FENIX home has a hard quota. Use group storage."; exit 1 ;;
esac
[[ "$MODE" == check ]] || [[ "$MODE" == dry ]] || command -v sbatch >/dev/null \
  || { echo "[X]  sbatch not found — run this on a FENIX login node."; exit 1; }

# 1) Normalise the map: label <TAB> absolute chrom_gvcf dir (column 2 may be the
#    sample dir or its chrom_gvcf/). Labels and dirs must be unique.
declare -A DIR_OF=() LABEL_OF_DIR=()
LABELS=()
lineno=0
while IFS=$'\t' read -r label dir _rest || [[ -n "${label:-}" ]]; do
  lineno=$((lineno + 1))
  label="${label%$'\r'}"; dir="${dir%$'\r'}"
  [[ -z "${label// }" || "$label" == \#* ]] && continue
  if [[ -z "${dir:-}" ]]; then fail "map line $lineno: need 2 TAB-separated columns: '$label'"; continue; fi
  dir="${dir%/}"
  [[ "$(basename "$dir")" != chrom_gvcf && -d "$dir/chrom_gvcf" ]] && dir="$dir/chrom_gvcf"
  if [[ ! -d "$dir" ]]; then fail "$label: chrom_gvcf dir not found: $dir"; continue; fi
  dir="$(readlink -f "$dir")"
  if [[ -n "${DIR_OF[$label]:-}" ]]; then fail "$label: label listed twice in the map"; continue; fi
  if [[ -n "${LABEL_OF_DIR[$dir]:-}" ]]; then fail "$label: same dir as '${LABEL_OF_DIR[$dir]}': $dir"; continue; fi
  DIR_OF[$label]="$dir"; LABEL_OF_DIR[$dir]="$label"; LABELS+=("$label")
done < "$MAP_IN"
N=${#LABELS[@]}
echo "[i]  map rows    : $N samples"
(( N == EXPECTED_N )) || warn "expected $EXPECTED_N samples (2026-10-23 selection), map has $N"

# 2) Excluded samples must not be in the map.
for x in "${EXCLUDED[@]}"; do
  [[ -n "${DIR_OF[$x]:-}" ]] && fail "$x is EXCLUDED from S04 (QC decision) but is in the map"
done

# 3) Merged duplicate IDs: present once, from Batch25; neither original copy in the map.
n_merged=0
while IFS=$'\t' read -r label dirA dirB _; do
  case "$label" in ''|'#'*|*_vs_*) continue ;; esac
  n_merged=$((n_merged + 1))
  want="$(readlink -f "$BAMQC_ROOT/$MERGED_BATCH/$label/chrom_gvcf" 2>/dev/null || echo "$BAMQC_ROOT/$MERGED_BATCH/$label/chrom_gvcf")"
  if [[ -z "${DIR_OF[$label]:-}" ]]; then
    fail "merged sample $label missing from the map (expected $want)"
  elif [[ "${DIR_OF[$label]}" != "$want" ]]; then
    fail "merged sample $label points at ${DIR_OF[$label]} — must be the $MERGED_BATCH merge: $want"
  fi
  for orig in "$dirA" "$dirB"; do
    o="$(readlink -f "$orig/chrom_gvcf" 2>/dev/null || true)"
    [[ -n "$o" && -n "${LABEL_OF_DIR[$o]:-}" ]] \
      && fail "original copy of merged $label is in the map as '${LABEL_OF_DIR[$o]}': $o (would import $label twice)"
  done
done < "$PAIRS"
echo "[i]  merged IDs  : $n_merged checked against $MERGED_BATCH (from $(basename "$PAIRS"))"

# 4) Every sample has exactly one indexed GVCF per autosome.
n_bad=0
for label in "${LABELS[@]}"; do
  d="${DIR_OF[$label]}"; miss=()
  for c in "${AUTOSOMES[@]}"; do
    shopt -s nullglob; hits=( "$d"/*."$c".g.vcf.gz ); shopt -u nullglob
    if (( ${#hits[@]} != 1 )); then miss+=("$c(${#hits[@]})")
    elif [[ ! -s "${hits[0]}.tbi" && ! -s "${hits[0]}.csi" ]]; then miss+=("$c(no index)"); fi
  done
  if (( ${#miss[@]} )); then fail "$label: autosome GVCF problems: ${miss[*]}"; n_bad=$((n_bad + 1)); fi
done
echo "[i]  GVCF files  : $((N - n_bad))/$N samples have all 22 autosomes, indexed"

# 5) GVCF header sample names: one per file, unique across the cohort (chr22).
if ! command -v bcftools >/dev/null && command -v module >/dev/null; then
  module load bcftools >/dev/null 2>&1 || true
fi
if command -v bcftools >/dev/null; then
  declare -A SEEN_SM=()
  for label in "${LABELS[@]}"; do
    shopt -s nullglob; g=( "${DIR_OF[$label]}"/*.chr22.g.vcf.gz ); shopt -u nullglob
    (( ${#g[@]} == 1 )) || continue
    sms="$(bcftools query -l "${g[0]}" 2>/dev/null | paste -sd, -)"
    case "$sms" in
      "")  fail "$label: cannot read sample name from ${g[0]}" ;;
      *,*) fail "$label: GVCF holds more than one sample ($sms)" ;;
      *)   if [[ -n "${SEEN_SM[$sms]:-}" ]]; then fail "$label: header sample '$sms' already used by '${SEEN_SM[$sms]}'"
           else SEEN_SM[$sms]="$label"; [[ "$sms" == "$label" ]] || warn "$label: GVCF header says '$sms' (S04 keeps the header name)"; fi ;;
    esac
  done
  echo "[i]  GVCF header : ${#SEEN_SM[@]} unique sample names"
else
  warn "bcftools not on PATH — header sample-name check skipped (S04 repeats it in every job)"
fi

# 6) Frozen map: all 22 chromosomes MUST be built from the same sample set.
mkdir -p "$OUTPUT_ROOT" "$LOG_DIR"
tmp_map="$(mktemp "$OUTPUT_ROOT/.map.XXXXXX")"
for label in "${LABELS[@]}"; do printf '%s\t%s\n' "$label" "${DIR_OF[$label]}"; done > "$tmp_map"
if [[ -s "$FROZEN_MAP" ]] && ! cmp -s "$tmp_map" "$FROZEN_MAP"; then
  done_any=false; for c in "${AUTOSOMES[@]}"; do ws_done "$c" && done_any=true; done
  if $done_any; then
    fail "map differs from $FROZEN_MAP, which built existing workspaces — chromosomes would hold different samples. Restore the old map, or delete genomicsdb/ + the frozen map to rebuild."
  else
    warn "replacing previous frozen map (no workspace built from it yet)"
  fi
fi

# 7) Free space and queue.
avail_gb=$(( $(df -Pk "$OUTPUT_ROOT" | awk 'NR==2{print $4}') / 1024 / 1024 ))
echo "[i]  free space  : ${avail_gb} GB on $(df -P "$OUTPUT_ROOT" | awk 'NR==2{print $6}')"
(( avail_gb >= MIN_FREE_GB )) || fail "less than ${MIN_FREE_GB} GB free on the output filesystem"
if command -v squeue >/dev/null && squeue -u "$USER" -h -o '%j' | grep -q "^JVC-GDBI-${TAG}-"; then
  fail "JVC-GDBI-${TAG} jobs are already queued/running — wait or scancel them first (no double submission)"
fi

echo
if (( ERR )); then
  rm -f "$tmp_map"
  echo "[X]  Pre-flight FAILED: $ERR problem(s) above. Nothing submitted."
  exit 1
fi
mv "$tmp_map" "$FROZEN_MAP"
echo "[OK] Pre-flight passed. Frozen map: $FROZEN_MAP ($N samples)"

# Shared group read access to the database (workspace files default to 0700).
if chgrp "$DATA_GROUP" "$OUTPUT_ROOT" 2>/dev/null && chmod g+rxs "$OUTPUT_ROOT" 2>/dev/null; then
  echo "[i]  $OUTPUT_ROOT is setgid $DATA_GROUP (new files inherit the group)"
else
  warn "could not set group $DATA_GROUP + setgid on $OUTPUT_ROOT — labmates may not read the database"
fi
[[ "$MODE" == check ]] && { echo "[i]  --check: nothing submitted."; exit 0; }
# <\CHECKS> -----------------------------------------------------------------

# --- submit one array per resource class ------------------------------------
exclude="$(awk '/^remote:/{r=1;next} /^[^[:space:]]/{r=0} r && $1=="sbatch_exclude:"{gsub(/"/,"",$2); print $2}' "$CONFIG_FILE")"
EXCL_ARG=(); [[ -n "$exclude" ]] && EXCL_ARG=(--exclude="$exclude")

echo
echo "[&]  ${TAG} Step 04 — $N samples, autosomes, action=$ACTION"
submit_class() {
  local class="$1" c chrs=()
  for c in ${CLASS_CHROMS[$class]}; do ws_done "$c" || chrs+=("$c"); done
  if (( ${#chrs[@]} == 0 )); then echo "[SKIP] class=$class — already built"; return 0; fi

  local mem="${CLASS_MEM[$class]}" hours="${CLASS_HOURS[$class]}" conc="${CLASS_CONC[$class]}"
  local consol="${CLASS_CONSOL[$class]}"
  local jobname="JVC-GDBI-${TAG}-${class}"
  local logfile="${LOG_DIR}/${TAG}-${class}-%A_%a.out"
  # Scratch = PRIMARY-group volume, resolved on the compute node (reachable on
  # every node; the guest amedina tree is not). umask 0027: group can read.
  local wrap="
set -euo pipefail
umask 0027
_chrs=(${chrs[*]})
_chrom=\"\${_chrs[\$((SLURM_ARRAY_TASK_ID - 1))]:-}\"
[[ -n \"\$_chrom\" ]] || { echo \"<ERROR> no chromosome for array index \$SLURM_ARRAY_TASK_ID\" >&2; exit 1; }
_scr=\"/scratch/groups/\$(id -gn)\"
mkdir -p \"\$_scr/\$USER\" || { echo \"<ERROR> scratch not writable on \$(hostname): \$_scr/\$USER\" >&2; exit 1; }
echo \"[i] array task \$SLURM_ARRAY_TASK_ID -> \$_chrom  (class ${class}, node \$(hostname), scratch \$_scr)\"
env GENDBI_SCRATCH_ROOT=\"\$_scr\" GENDBI_READER_THREADS=2 GENDBI_CONSOLIDATE=${consol} \\
  bash '$S04' '$FROZEN_MAP' '$OUTPUT_ROOT' \"\$_chrom\" '$ACTION'
# TileDB writes 0700 files; open the finished workspace to the data group.
chgrp -R '$DATA_GROUP' '$OUTPUT_ROOT/genomicsdb/'\"\$_chrom\" 2>/dev/null \\
  && chmod -R g+rX '$OUTPUT_ROOT/genomicsdb/'\"\$_chrom\" \\
  || echo \"[!] could not open \$_chrom workspace to group $DATA_GROUP — run: chgrp -R $DATA_GROUP ...; chmod -R g+rX ...\"
"
  echo "[*]  class=$class  chroms=(${chrs[*]})  $mem / ${hours}h  %$conc concurrent  consolidate=$consol"
  local args=(--parsable --job-name="$jobname" --array="1-${#chrs[@]}%${conc}"
              --nodes=1 --ntasks=1 --cpus-per-task=2 --mem="$mem" --time="${hours}:00:00"
              --output="$logfile" ${EXCL_ARG[@]+"${EXCL_ARG[@]}"})
  if [[ "$MODE" == dry ]]; then
    echo "     sbatch ${args[*]} --wrap \"...\""
  else
    local jid; jid=$(sbatch "${args[@]}" --wrap "$wrap")
    echo "     submitted: job $jid   logs: ${LOG_DIR}/${TAG}-${class}-${jid}_*.out"
  fi
}
for class in big mid small; do submit_class "$class"; done

echo
if [[ "$MODE" == dry ]]; then
  echo "[i]  --dry-run: nothing submitted."
else
  echo "[i]  Monitor : squeue -u \$USER | grep JVC-GDBI-${TAG}"
  echo "[i]  Status  : bash $0 --status"
  echo "[i]  A failed task leaves its chromosome MISSING; fix the cause, then re-run"
  echo "[i]  this script (no flags) — it only submits what is still missing."
  echo "[i]  When all 22 are DONE, check one workspace before Step 05:"
  echo "[i]    gatk SelectVariants -R <ref_gnm> -V gendb://${OUTPUT_ROOT}/genomicsdb/chr22 -L chr22:1-200000 -O ${OUTPUT_ROOT}/check_chr22.vcf.gz"
  echo "[i]    bcftools query -l ${OUTPUT_ROOT}/check_chr22.vcf.gz | wc -l   # expect $N"
fi

#EOF
