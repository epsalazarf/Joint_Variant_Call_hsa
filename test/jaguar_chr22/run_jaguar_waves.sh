#!/usr/bin/env bash
# =============================================================================
# Title : JAGUAR chr22 — sequential 3-wave launcher for Step 04 (GenomicsDBImport)
# About : Submits three chained SLURM jobs:
#           wave 1  create   (fresh GenomicsDB for chr22)
#           wave 2  update   (afterok wave 1)   — first incremental import
#           wave 3  update   (afterok wave 2)   — second incremental import
#         Each wave requests a DIFFERENT resource envelope so we can see whether
#         more CPU / RAM actually shortens GenomicsDBImport (it is mostly native,
#         mildly parallel — this is the experiment). Wave 3 doubles wave 2.
#         All SLURM logs go to LOG_DIR and a run manifest is written for
#         debugging + the resource-vs-runtime comparison afterwards.
#
# Usage :
#   bash run_jaguar_waves.sh [cohort_root] [output_dir]         # submit
#   bash run_jaguar_waves.sh --dry-run [cohort_root] [output_dir]
#   bash run_jaguar_waves.sh report [manifest_file]             # timings, post-run
#
#     cohort_root — dir containing bam/<SAMPLE>/chrom_gvcf/...
#                   default: /mnt/data/amedina/mramirezc/JAGUAR_JVC  (FENIX)
#     output_dir  — where genomicsdb/chr22 is built (default: ./out next to script)
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S04="$REPO/bin/04_gatk_GenomicsDB_import.sh"
CHROM="chr22"
MAPS_DIR="$HERE"
LOG_DIR="$HERE/logs"

# --- per-wave resource envelope (edit here) --------------------------------------
#            wave:      1        2         3
WAVE_ACTION=( -    create   update    update  )
WAVE_CPUS=(   -    8        8         16      )   # --cpus-per-task
WAVE_MEM=(    -    32G      64G       128G    )   # --mem
WAVE_TIME=(   -    12:00:00 24:00:00  48:00:00)   # --time
WAVE_RTHREADS=(-   8        8         16      )   # GENDBI_READER_THREADS
WAVE_BATCH=(  -    50       50        100     )   # GENDBI_BATCH_SIZE
N_WAVES=3
# NOTE: GenomicsDBImport keeps a small Java heap; big --mem is intentional
# head-room for the experiment, not a real requirement. The group tolerates
# large asks (ref: ADMIXTURE at 8 CPU / 192G / 96h scheduled without delay).

# =============================================================================

# ---------- report mode -------------------------------------------------------
if [[ "${1:-}" == "report" ]]; then
  manifest="${2:-$(ls -t "$LOG_DIR"/jaguar_run_*.manifest 2>/dev/null | head -1 || true)}"
  [[ -s "${manifest:-}" ]] || { echo "[X]  no manifest found (looked in $LOG_DIR)"; exit 1; }
  echo "[i]  manifest: $manifest"; echo
  jids=$(awk -F'\t' '/^wave/{print $3}' "$manifest" | paste -sd, -)
  [[ -n "$jids" ]] || { echo "[X]  no job ids in manifest"; exit 1; }
  if command -v sacct >/dev/null; then
    sacct -j "$jids" --units=G \
      --format=JobID%14,JobName%16,State%12,Elapsed%12,TotalCPU%12,MaxRSS%10,AllocCPUS%9,ReqMem%8,Start%20,End%20
  else
    echo "[!]  sacct not available; raw manifest:"; cat "$manifest"
  fi
  exit 0
fi

# ---------- submit mode ------------------------------------------------------
DRY=false
if [[ "${1:-}" == "--dry-run" ]]; then DRY=true; shift; fi

COHORT_ROOT="${1:-/mnt/data/amedina/mramirezc/JAGUAR_JVC}"   # JAGUAR test default
OUTPUT_DIR="${2:-$HERE/out}"

[[ -f "$S04" ]] || { echo "[X]  Step 04 script not found: $S04"; exit 1; }
[[ -d "$COHORT_ROOT" ]] || { echo "[X]  cohort_root not a directory: $COHORT_ROOT"; exit 1; }
COHORT_ROOT="$(readlink -f "$COHORT_ROOT")"
$DRY || command -v sbatch >/dev/null || { echo "[X]  sbatch not found"; exit 1; }
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
OUTPUT_DIR="$(readlink -f "$OUTPUT_DIR")"

# maps: build if absent
need_maps=false
for (( w=1; w<=N_WAVES; w++ )); do
  [[ -s "$MAPS_DIR/jaguar_${CHROM}_wave${w}.sample_map.tsv" ]] || need_maps=true
done
if $need_maps; then
  echo "[i]  wave maps missing — generating with make_wave_maps.sh"
  bash "$HERE/make_wave_maps.sh" "$COHORT_ROOT" "$MAPS_DIR/JAGUAR_JVC_${CHROM}_list.txt" "$N_WAVES" "$MAPS_DIR"
  echo
fi

TS="$(date +%Y%m%d-%H%M%S)"
MANIFEST="$LOG_DIR/jaguar_run_${TS}.manifest"
{
  echo -e "# JAGUAR chr22 Step-04 wave run\t$TS"
  echo -e "# cohort_root\t$COHORT_ROOT"
  echo -e "# output_dir\t$OUTPUT_DIR"
  echo -e "# s04\t$S04"
  echo -e "# columns: wave<TAB>action<TAB>jobid<TAB>cpus<TAB>mem<TAB>time<TAB>reader_threads<TAB>batch<TAB>map<TAB>logfile"
} > "$MANIFEST"

echo "[&]  JAGUAR chr22 — Step 04 wave launcher   ($TS)"
echo "[i]  cohort_root : $COHORT_ROOT"
echo "[i]  output_dir  : $OUTPUT_DIR"
echo "[i]  logs        : $LOG_DIR"
echo "[i]  manifest    : $MANIFEST"
echo "[i]  dry-run     : $DRY"
echo

dep=""
declare -a SUBMITTED
for (( w=1; w<=N_WAVES; w++ )); do
  action="${WAVE_ACTION[$w]}"
  map="$MAPS_DIR/jaguar_${CHROM}_wave${w}.sample_map.tsv"
  [[ -s "$map" ]] || { echo "[X]  missing map: $map"; exit 1; }
  nsamp=$(grep -cvE '^[[:space:]]*(#|$)' "$map")
  jobname="JVC-GDBI-w${w}"
  logfile="$LOG_DIR/jaguar-w${w}-%j.out"

  wrap="env GENDBI_READER_THREADS=${WAVE_RTHREADS[$w]} GENDBI_BATCH_SIZE=${WAVE_BATCH[$w]} \
bash '$S04' '$map' '$OUTPUT_DIR' '$CHROM' '$action'"

  set -- --parsable --job-name="$jobname" \
         --nodes=1 --ntasks=1 --cpus-per-task="${WAVE_CPUS[$w]}" \
         --mem="${WAVE_MEM[$w]}" --time="${WAVE_TIME[$w]}" \
         --output="$logfile"
  [[ -n "$dep" ]] && set -- "$@" --dependency="afterok:$dep"

  echo "[*]  wave $w: $action   ($nsamp samples)   ${WAVE_CPUS[$w]} CPU / ${WAVE_MEM[$w]} / ${WAVE_TIME[$w]}   rthreads=${WAVE_RTHREADS[$w]} batch=${WAVE_BATCH[$w]}${dep:+   [afterok:$dep]}"

  if $DRY; then
    echo "       sbatch $* --wrap \"$wrap\""
    jid="DRYRUN-w${w}"
  else
    jid=$(sbatch "$@" --wrap "$wrap")
    echo "       submitted: job $jid   log: ${logfile/\%j/$jid}"
  fi

  echo -e "wave${w}\t$action\t$jid\t${WAVE_CPUS[$w]}\t${WAVE_MEM[$w]}\t${WAVE_TIME[$w]}\t${WAVE_RTHREADS[$w]}\t${WAVE_BATCH[$w]}\t$map\t${logfile/\%j/$jid}" >> "$MANIFEST"
  SUBMITTED+=("$jid")
  dep="$jid"
done

echo
if $DRY; then
  echo "[i]  dry run — nothing submitted. Manifest: $MANIFEST"
else
  echo "[i]  Monitor:   squeue -u \$USER --name=JVC-GDBI-w1,JVC-GDBI-w2,JVC-GDBI-w3"
  echo "[i]           squeue -j $(IFS=,; echo "${SUBMITTED[*]}")"
  echo "[i]  Logs:      $LOG_DIR/jaguar-w*-<jobid>.out"
  echo "[i]  Timings:   bash $(basename "$0") report $MANIFEST      # after all 3 finish"
  echo "[i]  If a wave FAILS, later waves stay PENDING (DependencyNeverSatisfied) —"
  echo "[i]  fix, then resubmit from that wave with action=update."
fi
