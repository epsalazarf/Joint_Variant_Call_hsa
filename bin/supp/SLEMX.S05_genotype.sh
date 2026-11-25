#!/usr/bin/env bash

# =============================================================================
# Title       : SLEmx (lupus) — Step 05 joint genotyping launcher [FENIX]
# Description : Checks the Step 04 GenomicsDB workspaces and submits Step 05
#               (bin/05_gatk_GenotypeGVCFs.sh) SCATTERED: each chromosome is cut
#               into fixed windows (default 5 Mb), one SLURM array task per
#               window, then a gather job joins the windows into
#               genotyped/slemx.<chrom>.vcf.gz (+ .tbi).
#
#               Why scatter: on 164 samples a whole-chromosome GenotypeGVCFs
#               grows its native (GenomicsDB, non-Java) memory as it walks the
#               chromosome — chr22 and chr1 were both OOM-killed at 32G after
#               ~2.5-3 h. A window job ends before it grows; the Broad's own
#               joint-calling workflow scatters the same way. Windows use
#               --only-output-calls-starting-in-intervals, so no site is
#               written twice at a boundary.
#
#               --check (no submission, login node, ~1 min): for every
#               autosome workspace, reports complete / sample count vs the
#               frozen S04 map / same sample names on every chromosome /
#               fragment count / group-readable / size, and lists leftovers
#               (Step 04/05 scratch, unfinished .tmp VCFs). Fails loudly.
#
#               Submit: --chroms chr22,chr1 (required — no "all" by accident).
#               The same check runs first. Finished windows are kept: a re-run
#               submits only the missing windows (+ the gather), so re-running
#               after a failure is safe. A chromosome whose final VCF exists,
#               or whose jobs are queued, is skipped.
#
#               The gather (GATK GatherVcfs, block copy) checks every window is
#               there and complete, drops header-only windows (all-N, e.g.
#               centromeres), re-indexes with bcftools, and compares the record
#               count with the sum of the windows. Only then are the windows
#               deleted (keep them: SLEMX_S05_KEEP_CHUNKS=true).
#
#               --region chr22:20000001-20500000 runs a small piece of ONE
#               chromosome (minutes): a cheap read test, not scattered.
#
# Usage       : bash SLEMX.S05_genotype.sh --check
#               bash SLEMX.S05_genotype.sh --chroms chr22 --region chr22:20000001-20500000
#               bash SLEMX.S05_genotype.sh --chroms chr22,chr1 [--dry-run]
#               bash SLEMX.S05_genotype.sh --status
#                 --out DIR      cohort root (default /mnt/data/amedina/esalazarf/SLEmx_JVC)
#                                reads <out>/genomicsdb/<chrom>, writes <out>/genotyped/
#                 --chunk-mb N   window size in Mb (default 5). Keep it the same
#                                until a chromosome is gathered.
#                 --par N        max running windows per chromosome (default 10)
# =============================================================================

set -euo pipefail

# <CONFIG> ----------------------------------------------------------------
OUTPUT_ROOT="/mnt/data/amedina/esalazarf/SLEmx_JVC"
TAG="slemx"
DATA_GROUP="amedina"
AUTOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
# Resources. Measured (164 samples, 5 Mb windows, 8G Java heap; 2026-11-25,
# chr1/6/19/20/21/22 = 131 windows, 0 failed): peak MaxRSS 15.7G (chr1
# 35-40 Mb), 12-15G on dense windows, ~1.5G on all-N ones; longest window
# 56 min. Whole chromosomes (no scatter) grew past 32G in ~2.5-3 h — native
# GenomicsDB memory, not the Java heap. 24G = ~1.5x the worst window seen;
# 4 h = ~4x the longest (a short limit helps SLURM backfill).
CHUNK_MB=5
CHUNK_MEM="24G";  CHUNK_JAVA_MEM="8G";  CHUNK_HOURS=4
PAR=10                                   # running windows per chromosome
GATHER_MEM="8G";  GATHER_HOURS=24
REGION_MEM="16G"; REGION_JAVA_MEM="8G"; REGION_HOURS=12
KEEP_CHUNKS="${SLEMX_S05_KEEP_CHUNKS:-false}"
# <\CONFIG> -----------------------------------------------------------------

SELF="$(readlink -f "$0")"
HERE="$(cd "$(dirname "$SELF")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S05="$REPO/bin/05_gatk_GenotypeGVCFs.sh"
CONFIG_FILE="$REPO/config/config.yaml"

MODE="submit"; CHROMS=""; REGION=""; GATHER_CHROM=""
while [ $# -gt 0 ]; do
  case "$1" in
    --check)    MODE="check";  shift ;;
    --status)   MODE="status"; shift ;;
    --dry-run)  MODE="dry";    shift ;;
    --chroms)   CHROMS="$(echo "${2:?--chroms needs a list, e.g. chr22,chr1}" | tr ',' ' ')"; shift 2 ;;
    --region)   REGION="${2:?--region needs chrN:start-end}"; shift 2 ;;
    --out)      OUTPUT_ROOT="${2:?--out needs a directory}"; shift 2 ;;
    --chunk-mb) CHUNK_MB="${2:?--chunk-mb needs a number}"; shift 2 ;;
    --par)      PAR="${2:?--par needs a number}"; shift 2 ;;
    --gather)   MODE="gather"; GATHER_CHROM="${2:?--gather needs a chromosome}"; shift 2 ;;  # used by the gather job
    -h|--help)  sed -n '4,49p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)          echo "[X]  unknown argument: $1  (--help)"; exit 2 ;;
  esac
done
[[ "$CHUNK_MB" =~ ^[1-9][0-9]*$ ]] || { echo "[X]  --chunk-mb must be a whole number of Mb (got '$CHUNK_MB')"; exit 2; }
[[ "$PAR" =~ ^[1-9][0-9]*$ ]]      || { echo "[X]  --par must be a positive number (got '$PAR')"; exit 2; }
CHUNK_BP="${SLEMX_S05_CHUNK_BP:-$(( CHUNK_MB * 1000000 ))}"   # env: tests only
OUTPUT_ROOT="${OUTPUT_ROOT%/}"
GDB_ROOT="$OUTPUT_ROOT/genomicsdb"
VCF_DIR="$OUTPUT_ROOT/genotyped"
LOG_DIR="$OUTPUT_ROOT/logs"
FROZEN_MAP="$OUTPUT_ROOT/${TAG}.sample_map.S04.tsv"

# Output names, same rule as 05_gatk_GenotypeGVCFs.sh (prefix TAG).
vcf_for() {   # $1 chrom [$2 region]
  if [[ -n "${2:-}" ]]; then local se="${2#*:}"; echo "$VCF_DIR/${TAG}.$1_${se}.vcf.gz"
  else echo "$VCF_DIR/${TAG}.$1.vcf.gz"; fi
}
chunk_dir() { echo "$VCF_DIR/chunks/$1"; }
chrom_len() {   # $1 chrom -> length, from the workspace array dir '<chrom>$1$<length>'
  local d; for d in "$GDB_ROOT/$1/$1"\$1\$*; do [[ -d "$d" ]] && { echo "${d##*\$}"; return 0; }; done
  return 1
}
n_chunks() { echo $(( ($1 + CHUNK_BP - 1) / CHUNK_BP )); }        # $1 length
chunk_region() {   # $1 chrom $2 length $3 index(1-based) -> chrN:start-end
  local s=$(( ($3 - 1) * CHUNK_BP + 1 )) e=$(( $3 * CHUNK_BP ))
  if (( e > $2 )); then e=$2; fi
  echo "$1:${s}-${e}"
}
chunk_vcf() { local r; r="$(chunk_region "$1" "$2" "$3")"; echo "$(chunk_dir "$1")/${TAG}.$1_${r#*:}.vcf.gz"; }
vcf_done() { [[ -s "$1" && -s "$1.tbi" ]]; }
jobname_for() { if [[ -n "${2:-}" ]]; then echo "JVC-GTGVCF-${TAG}-$1-$2"; else echo "JVC-GTGVCF-${TAG}-$1"; fi; }
queued() { command -v squeue >/dev/null && squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -qx "$1"; }
chunks_done() {   # $1 chrom $2 length -> number of finished windows
  local i k=0 n; n=$(n_chunks "$2")
  for (( i = 1; i <= n; i++ )); do vcf_done "$(chunk_vcf "$1" "$2" "$i")" && k=$((k + 1)); done
  echo "$k"
}
hms() { echo $(( $(date +%s) - $1 )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp'; }
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
bgzf_complete() { [[ "$(tail -c 28 "$1" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')" == "$BGZF_EOF" ]]; }

# --- status -----------------------------------------------------------------
if [[ "$MODE" == status ]]; then
  echo "[&]  ${TAG} Step 05 status — $VCF_DIR   (windows: ${CHUNK_MB} Mb)"
  for c in $AUTOSOMES; do
    f="$(vcf_for "$c")"; q=""
    queued "$(jobname_for "$c")" && q="  windows queued/running"
    queued "$(jobname_for "$c" gather)" && q="${q}  gather queued/running"
    if vcf_done "$f"; then printf '  [DONE]    %-6s %s\n' "$c" "$(du -h "$f" | cut -f1)"; continue; fi
    len="$(chrom_len "$c" || true)"
    if [[ -n "$len" && -d "$(chunk_dir "$c")" ]]; then
      printf '  [%3s/%-3s] %-6s windows done%s\n' "$(chunks_done "$c" "$len")" "$(n_chunks "$len")" "$c" "$q"
    elif [[ -n "$q" ]]; then printf '  [QUEUED]  %-6s%s\n' "$c" "$q"
    else printf '  [ ---- ]  %s\n' "$c"; fi
  done
  shopt -s nullglob; for f in "$VCF_DIR"/"${TAG}".chr*_*-*.vcf.gz; do echo "  [region]  $(basename "$f")"; done; shopt -u nullglob
  command -v squeue >/dev/null && squeue -u "$USER" -h -o '%.18i %.40j %.8T %.12M %R' | grep "JVC-GTGVCF-${TAG}" || true
  exit 0
fi

# --- gather (runs inside the SLURM gather job) --------------------------------
if [[ "$MODE" == gather ]]; then
  t0=$(date +%s); c="$GATHER_CHROM"
  echo "[$] ${TAG} Step 05 gather — $c   ($(date), node $(hostname))"
  final="$(vcf_for "$c")"; tmp="${final%.vcf.gz}.tmp.vcf.gz"
  if vcf_done "$final"; then echo "[SKIP] $final already exists"; exit 0; fi
  len="$(chrom_len "$c")" || { echo "<ERROR> no workspace array for $c in $GDB_ROOT/$c"; exit 1; }
  n=$(n_chunks "$len")
  if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]] && command -v module >/dev/null 2>&1; then
    module unload oracle-java 2>/dev/null || true
    module load oracle-java/25.0.2; module load gatk; module load bcftools
  fi
  echo "[i]  $c: $len bp, $n windows of ${CHUNK_MB} Mb -> $final"
  inputs=(); n_empty=0; n_sum=0; n_samp=""
  for (( i = 1; i <= n; i++ )); do
    f="$(chunk_vcf "$c" "$len" "$i")"
    vcf_done "$f" || { echo "<ERROR> window $i missing (no VCF + .tbi): $f"; exit 1; }
    bgzf_complete "$f" || { echo "<ERROR> window $i truncated (no BGZF EOF): $f"; exit 1; }
    s=$(bcftools query -l "$f" | wc -l | tr -d '[:space:]')
    if [[ -z "$n_samp" ]]; then n_samp=$s
    elif (( s != n_samp )); then echo "<ERROR> window $i has $s samples, window 1 has $n_samp: $f"; exit 1; fi
    # Full read: record count, and every block decompresses.
    r=$(bcftools view -H "$f" | wc -l | tr -d '[:space:]')
    if (( r == 0 )); then n_empty=$((n_empty + 1)); echo "[i]  window $i: no records (all-N?) — left out"; continue; fi
    n_sum=$((n_sum + r)); inputs+=(-I "$f")
  done
  (( ${#inputs[@]} )) || { echo "<ERROR> all $n windows are empty for $c"; exit 1; }
  echo "[i]  windows OK: $n ($n_empty empty), $n_samp samples, $n_sum records   [$(hms "$t0")]"

  _scr="/scratch/groups/$(id -gn)/$USER"; mkdir -p "$_scr" 2>/dev/null || _scr="$VCF_DIR"
  rm -f "$tmp" "$tmp.tbi"
  set -o xtrace
  gatk --java-options "-Xmx4G" GatherVcfs --TMP_DIR "$_scr" "${inputs[@]}" -O "$tmp"
  bcftools index -t --threads 2 -f "$tmp"
  { set +o xtrace; } 2>/dev/null
  bgzf_complete "$tmp" || { echo "<ERROR> gathered VCF truncated: $tmp"; exit 1; }
  s=$(bcftools query -l "$tmp" | wc -l | tr -d '[:space:]')
  (( s == n_samp )) || { echo "<ERROR> gathered VCF has $s samples, windows have $n_samp"; exit 1; }
  r=$(bcftools index -n "$tmp")
  (( r == n_sum )) || { echo "<ERROR> gathered VCF has $r records, windows sum to $n_sum — windows kept"; exit 1; }
  mv -f "$tmp.tbi" "$final.tbi"; mv -f "$tmp" "$final"
  chgrp "$DATA_GROUP" "$final" "$final.tbi" 2>/dev/null && chmod g+r "$final" "$final.tbi" || true
  echo "[DONE] $final   ($s samples, $r records, $(du -h "$final" | cut -f1))"
  if [[ "$KEEP_CHUNKS" == true ]]; then echo "[i]  windows kept: $(chunk_dir "$c")"
  else rm -rf "$(chunk_dir "$c")"; echo "[i]  windows deleted: $(chunk_dir "$c")  (SLEMX_S05_KEEP_CHUNKS=true keeps them)"; fi
  echo "[&]  Total time: $(hms "$t0")"
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

# Leftovers: an OOM kill skips the cleanup traps of Steps 04 and 05.
scr="/scratch/groups/$(id -gn)/$USER"
shopt -s nullglob
left=( "$scr"/gendbi_* "$scr"/gtgvcf_* )
tmps=( "$VCF_DIR"/*.tmp.vcf.gz "$VCF_DIR"/*.tmp.vcf.gz.tbi "$VCF_DIR"/chunks/*/*.tmp.vcf.gz "$VCF_DIR"/chunks/*/*.tmp.vcf.gz.tbi )
shopt -u nullglob
if (( ${#left[@]} )); then
  warn "${#left[@]} leftover Step 04/05 scratch dir(s) in $scr (safe to delete once none of your S04/S05 jobs run):"
  for d in "${left[@]}"; do echo "       $(du -sh "$d" 2>/dev/null | cut -f1)  $d"; done
  echo "       rm -rf $scr/gendbi_* $scr/gtgvcf_*"
fi
if (( ${#tmps[@]} )); then
  warn "${#tmps[@]} unfinished .tmp VCF file(s) (from killed jobs; rewritten on the next run, safe to delete once no S05 job runs):"
  for f in "${tmps[@]}"; do echo "       $(du -h "$f" 2>/dev/null | cut -f1)  $f"; done
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
common=(--parsable --nodes=1 --ntasks=1)
[[ -n "$exclude" ]] && common+=(--exclude="$exclude")
submit() {   # sbatch args... (last: --wrap text) -> prints job id (fake id on --dry-run)
  if [[ "$MODE" == dry ]]; then echo "DRY"; else sbatch "$@"; fi
}

# The scratch line shared by window and region jobs.
scratch_wrap='
_scr="/scratch/groups/$(id -gn)"
mkdir -p "$_scr/$USER" || { echo "<ERROR> scratch not writable on $(hostname): $_scr/$USER" >&2; exit 1; }'

echo
echo "[&]  ${TAG} Step 05 — GenotypeGVCFs   (output: $VCF_DIR)"

# Region test: one job, not scattered.
if [[ -n "$REGION" ]]; then
  c="$CHROMS"; out="$(vcf_for "$c" "$REGION")"; jn="$(jobname_for "$c" region)"
  if vcf_done "$out"; then echo "[SKIP] $c — $(basename "$out") already exists"
  elif queued "$jn"; then echo "[SKIP] $c — job $jn already queued/running"
  else
    log="$LOG_DIR/${TAG}-S05-${c}-region-%j.out"
    wrap="set -euo pipefail; umask 0027
$scratch_wrap
echo \"[i] ${REGION} on node \$(hostname), temp in \$_scr\"
env GTGVCF_SCRATCH_ROOT=\"\$_scr\" GTGVCF_PREFIX='${TAG}' GTGVCF_JAVA_MEM='${REGION_JAVA_MEM}' \\
  bash '$S05' '$GDB_ROOT' '$VCF_DIR' '$c' '$REGION'
chgrp '$DATA_GROUP' '$out' '$out.tbi' 2>/dev/null && chmod g+r '$out' '$out.tbi' || true"
    echo "[*]  $REGION   $REGION_MEM / ${REGION_HOURS}h   -> $(basename "$out")"
    jid=$(submit "${common[@]}" --job-name="$jn" --cpus-per-task=2 --mem="$REGION_MEM" \
            --time="${REGION_HOURS}:00:00" --output="$log" --wrap "$wrap")
    echo "     job $jid   log: ${log/\%j/$jid}"
  fi
  echo; [[ "$MODE" == dry ]] && echo "[i]  --dry-run: nothing submitted."
  exit 0
fi

for c in $CHROMS; do
  final="$(vcf_for "$c")"; jn="$(jobname_for "$c")"; jg="$(jobname_for "$c" gather)"
  if vcf_done "$final"; then echo "[SKIP] $c — $(basename "$final") already exists"; continue; fi
  if queued "$jn" || queued "$jg"; then echo "[SKIP] $c — $jn / gather already queued/running"; continue; fi
  len="$(chrom_len "$c")"; n=$(n_chunks "$len"); cdir="$(chunk_dir "$c")"

  # Windows of another size in the chunk dir would never be gathered: refuse.
  if [[ -d "$cdir" ]]; then
    stray=0
    shopt -s nullglob
    for f in "$cdir"/"${TAG}".*.vcf.gz; do
      case "$f" in *.tmp.vcf.gz) continue ;; esac
      ok=0; for (( i = 1; i <= n; i++ )); do [[ "$f" == "$(chunk_vcf "$c" "$len" "$i")" ]] && { ok=1; break; }; done
      (( ok )) || stray=$((stray + 1))
    done
    shopt -u nullglob
    (( stray == 0 )) || { echo "[X]  $c: $stray window(s) in $cdir are not ${CHUNK_MB} Mb windows — use the --chunk-mb they were made with, or delete them."; exit 1; }
  fi

  todo=""; k=0
  for (( i = 1; i <= n; i++ )); do
    if vcf_done "$(chunk_vcf "$c" "$len" "$i")"; then k=$((k + 1)); else todo="${todo:+$todo,}$i"; fi
  done
  echo "[*]  $c   $len bp = $n windows of ${CHUNK_MB} Mb   ($k done, $((n - k)) to run)"

  dep=()
  if [[ -n "$todo" ]]; then
    if [[ "$MODE" != dry ]]; then mkdir -p "$cdir"; fi
    log="$LOG_DIR/${TAG}-S05-${c}-w%a-%A.out"
    wrap="set -euo pipefail; umask 0027
$scratch_wrap
i=\$SLURM_ARRAY_TASK_ID
s=\$(( (i - 1) * $CHUNK_BP + 1 )); e=\$(( i * $CHUNK_BP ))
if (( e > $len )); then e=$len; fi
echo \"[i] ${c} window \$i/${n} = ${c}:\$s-\$e on node \$(hostname), temp in \$_scr\"
env GTGVCF_SCRATCH_ROOT=\"\$_scr\" GTGVCF_PREFIX='${TAG}' GTGVCF_JAVA_MEM='${CHUNK_JAVA_MEM}' \\
  bash '$S05' '$GDB_ROOT' '$cdir' '$c' \"${c}:\$s-\$e\"
out='$cdir/${TAG}.${c}_'\"\$s-\$e\".vcf.gz
chgrp '$DATA_GROUP' \"\$out\" \"\$out.tbi\" 2>/dev/null && chmod g+r \"\$out\" \"\$out.tbi\" || true"
    jid=$(submit "${common[@]}" --job-name="$jn" --array="${todo}%${PAR}" --cpus-per-task=2 \
            --mem="$CHUNK_MEM" --time="${CHUNK_HOURS}:00:00" --output="$log" --wrap "$wrap")
    echo "     windows : array job $jid  [$todo]%${PAR}   $CHUNK_MEM / ${CHUNK_HOURS}h each"
    echo "               logs: ${log/\%A/$jid}"
    dep=(--dependency="afterok:$jid" --kill-on-invalid-dep=yes)
  fi

  log="$LOG_DIR/${TAG}-S05-${c}-gather-%j.out"
  wrap="set -euo pipefail; umask 0027
env SLEMX_S05_KEEP_CHUNKS='${KEEP_CHUNKS}' bash '$SELF' --gather '$c' --out '$OUTPUT_ROOT' --chunk-mb '$CHUNK_MB'"
  gid=$(submit "${common[@]}" --job-name="$jg" "${dep[@]+"${dep[@]}"}" --cpus-per-task=2 \
          --mem="$GATHER_MEM" --time="${GATHER_HOURS}:00:00" --output="$log" --wrap "$wrap")
  echo "     gather  : job $gid  (after all windows succeed)  -> $(basename "$final")"
  echo "               log: ${log/\%j/$gid}"
done

echo
if [[ "$MODE" == dry ]]; then echo "[i]  --dry-run: nothing submitted."
else
  echo "[i]  Monitor : squeue -u \$USER | grep JVC-GTGVCF-${TAG}"
  echo "[i]  Status  : bash $0 --status        (windows done per chromosome)"
  echo "[i]  A window fails -> its gather is cancelled; just run the same command again (only missing windows rerun)."
  echo "[i]  Sizing  : sacct -j <array jobid> -o JobID,State,Elapsed,MaxRSS"
fi

#EOF
