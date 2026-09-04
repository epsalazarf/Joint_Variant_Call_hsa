#!/usr/bin/env bash
# =============================================================================
# Title : JAGUAR — sequential N-wave launcher for Step 04 (GenomicsDBImport)
# About : Submits N chained SLURM jobs for one chromosome:
#           wave 1        create  (fresh GenomicsDB)   + consolidate (trivial)
#           waves 2..N    update  (afterok previous)    append a fragment only
#         Consolidation on `update` merges fragments and its cost scales with
#         the whole DB (~13.5 h for 93 samples on chr1 — worse than a full
#         rebuild), so updates never consolidate. GenotypeGVCFs reads
#         multi-fragment DBs fine. All logs -> LOG_DIR + a manifest.
#
# Usage :
#   bash run_jaguar_waves.sh [options] [cohort_root] [output_dir]      # submit
#   bash run_jaguar_waves.sh --dry-run [options] [cohort_root] [output_dir]
#   bash run_jaguar_waves.sh report [manifest_file]                    # post-run timings
#
#   -c, --chrom CHR    chromosome (default: chr22)
#   -w, --waves N      number of waves (default: 3)
#   -t, --tag   NAME   cohort tag — namespaces maps/jobs/logs/manifest so two
#                      cohorts tested in this dir don't collide (default: jaguar)
#       --cpus  N      --cpus-per-task per wave (default: 2 — cores don't help)
#       --mem   SIZE   --mem per wave    (default: chromosome-scaled, see below)
#       --hours H      --time in hours   (default: chromosome-scaled, see below)
#       --batch N      GenomicsDBImport --batch-size (default 50; drop to 25 to
#                      roughly halve import memory on the big chromosomes)
#
#   Memory scales with CHROMOSOME SIZE, not just sample count (the per-sample
#   TileDB fragment data is ~5x bigger on chr1 than chr22). test02: chr1 create
#   of 47 samples used 27 GB; the update+consolidate hit 32 GB (100%). So the
#   --mem / --hours defaults below are keyed to the contig:
#     chr1, chr2                 -> 64G / 20h
#     chr3-8, chrX               -> 32G / 12h
#     chr9-22, chrY, chrM        -> 16G /  6h
#   Override with --mem / --hours for anything unusual.
#
#   cohort_root  dir containing bam/<SAMPLE>/chrom_gvcf/...
#                (default: /mnt/data/amedina/mramirezc/JAGUAR_JVC)
#   output_dir   where genomicsdb/<CHR> is built. MUST be group storage, NOT
#                $HOME (FENIX home has a hard ~5 GB quota; the copy-back would
#                fail after the import). Default: the JVCdev test area.
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S04="$REPO/bin/04_gatk_GenomicsDB_import.sh"
MAPS_DIR="$HERE"
LOG_DIR="$HERE/logs"
DEFAULT_COHORT="/mnt/data/amedina/mramirezc/JAGUAR_JVC"
DEFAULT_OUTPUT_DIR="/mnt/data/amedina/${USER:-esalazarf}/JVCdev/test_gendbi_jag22"

# ---------- report mode (before option parsing) -----------------------------
if [[ "${1:-}" == "report" ]]; then
  manifest="${2:-$(ls -t "$LOG_DIR"/*_run_*.manifest 2>/dev/null | head -1 || true)}"
  [[ -s "${manifest:-}" ]] || { echo "[X]  no manifest found (looked in $LOG_DIR)"; exit 1; }
  echo "[i]  manifest: $manifest"; echo
  jids=$(awk -F'\t' '/^wave/{print $3}' "$manifest" | paste -sd, -)
  [[ -n "$jids" ]] || { echo "[X]  no job ids in manifest"; exit 1; }
  if command -v sacct >/dev/null; then
    sacct -j "$jids" --units=G \
      --format=JobID%14,JobName%18,State%12,Elapsed%12,TotalCPU%12,MaxRSS%10,AllocCPUS%9,ReqMem%8,Start%20,End%20
  else
    echo "[!]  sacct not available; raw manifest:"; cat "$manifest"
  fi
  exit 0
fi

# ---------- options --------------------------------------------------------
# --mem / --hours default to a size class keyed on the contig (hg38 lengths):
#   chr1 chr2  -> the two biggest;  chr3-8 + chrX  -> mid;  rest -> small.
chrom_class() {
  case "$1" in
    chr1|chr2)                                      echo big ;;
    chr3|chr4|chr5|chr6|chr7|chr8|chrX)             echo mid ;;
    *)                                              echo small ;;
  esac
}
default_mem()   { case "$(chrom_class "$1")" in big) echo 64G ;; mid) echo 32G ;; *) echo 16G ;; esac; }
default_hours() { case "$(chrom_class "$1")" in big) echo 20  ;; mid) echo 12  ;; *) echo 6  ;; esac; }

CHROM="chr22" ; N_WAVES=3 ; TAG="jaguar" ; CPUS=2 ; MEM="" ; HOURS="" ; BATCH=50 ; DRY=false
pos=()
while (( $# )); do
  case "$1" in
    --dry-run)     DRY=true; shift ;;
    -c|--chrom)    CHROM="${2:?}"; shift 2 ;;
    -w|--waves)    N_WAVES="${2:?}"; shift 2 ;;
    -t|--tag)      TAG="${2:?}"; shift 2 ;;
    --cpus)        CPUS="${2:?}"; shift 2 ;;
    --mem)         MEM="${2:?}"; shift 2 ;;
    --hours)       HOURS="${2:?}"; shift 2 ;;
    --batch)       BATCH="${2:?}"; shift 2 ;;
    -h|--help)     sed -n '3,42p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
    -*)            echo "[X]  unknown option: $1  (--help)"; exit 2 ;;
    *)             pos+=("$1"); shift ;;
  esac
done
COHORT_ROOT="${pos[0]:-$DEFAULT_COHORT}"
OUTPUT_DIR="${pos[1]:-$DEFAULT_OUTPUT_DIR}"
(( N_WAVES >= 1 )) || { echo "[X]  --waves must be >= 1"; exit 2; }
MEM="${MEM:-$(default_mem "$CHROM")}"
HOURS="${HOURS:-$(default_hours "$CHROM")}"

fmt_time() { printf '%d:00:00' "$1"; }   # H -> H:00:00

# ---------- checks --------------------------------------------------------
[[ -f "$S04" ]] || { echo "[X]  Step 04 script not found: $S04"; exit 1; }
[[ -d "$COHORT_ROOT" ]] || { echo "[X]  cohort_root not a directory: $COHORT_ROOT"; exit 1; }
COHORT_ROOT="$(readlink -f "$COHORT_ROOT")"
$DRY || command -v sbatch >/dev/null || { echo "[X]  sbatch not found"; exit 1; }
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
OUTPUT_DIR="$(readlink -f "$OUTPUT_DIR")"
case "$OUTPUT_DIR/" in
  "${HOME}"/*)
    echo "[X]  output_dir is under \$HOME ($HOME) — FENIX home has a hard quota."
    echo "[i]  Pass a group-storage path, or edit DEFAULT_OUTPUT_DIR."
    exit 1 ;;
esac

# ---------- wave plan ----------------------------------------------------
#   action:      wave 1 = create, rest = update
#   consolidate: wave 1 ONLY. On `create` it is a ~4 s no-op (single fragment).
#                On `update` it merges fragments and the cost scales with the
#                whole DB — measured at ~13.5 h for 93 samples on chr1 (test02),
#                which is worse than a full `create` rebuild. So updates just
#                append a fragment; GenotypeGVCFs reads multi-fragment DBs fine.
#                Force it on a wave with GENDBI_CONSOLIDATE=true if you must.
wave_action()  { (( $1 == 1 )) && echo create || echo update; }
wave_consol()  { (( $1 == 1 )) && echo true || echo false; }
wave_time()    { fmt_time "$HOURS"; }

# ---------- maps: build if absent --------------------------------------------
need_maps=false
for (( w=1; w<=N_WAVES; w++ )); do
  [[ -s "$MAPS_DIR/${TAG}_${CHROM}_wave${w}.sample_map.tsv" ]] || need_maps=true
done
if $need_maps; then
  echo "[i]  wave maps missing for ${TAG}/${CHROM} — generating with make_wave_maps.sh"
  bash "$HERE/make_wave_maps.sh" --chrom "$CHROM" --waves "$N_WAVES" --tag "$TAG" --cohort "$COHORT_ROOT" --out "$MAPS_DIR"
  echo
fi

# ---------- submit -----------------------------------------------------------
TS="$(date +%Y%m%d-%H%M%S)"
MANIFEST="$LOG_DIR/${TAG}_run_${CHROM}_${TS}.manifest"
{
  echo -e "# ${TAG} ${CHROM} Step-04 wave run\t$TS"
  echo -e "# cohort_root\t$COHORT_ROOT"
  echo -e "# output_dir\t$OUTPUT_DIR"
  echo -e "# s04\t$S04"
  echo -e "# columns: wave\taction\tjobid\tcpus\tmem\ttime\tconsolidate\tmap\tlogfile"
} > "$MANIFEST"

echo "[&]  ${TAG} ${CHROM} — Step 04 wave launcher   ($TS)"
echo "[i]  tag         : $TAG"
echo "[i]  chromosome  : $CHROM     waves: $N_WAVES"
echo "[i]  cohort_root : $COHORT_ROOT"
echo "[i]  output_dir  : $OUTPUT_DIR"
echo "[i]  per wave    : ${CPUS} CPU / ${MEM} / ${HOURS}h  batch=${BATCH}  ($(chrom_class "$CHROM") chromosome)"
echo "[i]  logs        : $LOG_DIR"
echo "[i]  manifest    : $MANIFEST"
echo "[i]  dry-run     : $DRY"
echo

dep="" ; SUBMITTED=()
for (( w=1; w<=N_WAVES; w++ )); do
  action="$(wave_action "$w")" ; consol="$(wave_consol "$w")" ; wtime="$(wave_time "$w")"
  map="$MAPS_DIR/${TAG}_${CHROM}_wave${w}.sample_map.tsv"
  [[ -s "$map" ]] || { echo "[X]  missing map: $map"; exit 1; }
  nsamp=$(grep -cvE '^[[:space:]]*(#|$)' "$map")
  jobname="JVC-GDBI-${TAG}-${CHROM}-w${w}"
  logfile="$LOG_DIR/${TAG}-${CHROM}-w${w}-%j.out"

  wrap="env GENDBI_READER_THREADS=${CPUS} GENDBI_BATCH_SIZE=${BATCH} GENDBI_CONSOLIDATE=${consol} \
bash '$S04' '$map' '$OUTPUT_DIR' '$CHROM' '$action'"

  set -- --parsable --job-name="$jobname" \
         --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --mem="$MEM" --time="$wtime" \
         --output="$logfile"
  [[ -n "$dep" ]] && set -- "$@" --dependency="afterok:$dep"

  echo "[*]  wave $w: $action  ($nsamp samples)  ${CPUS} CPU / ${MEM} / ${wtime}  consolidate=${consol}${dep:+   [afterok:$dep]}"

  if $DRY; then
    echo "       sbatch $* --wrap \"$wrap\""
    jid="DRYRUN-w${w}"
  else
    jid=$(sbatch "$@" --wrap "$wrap")
    echo "       submitted: job $jid   log: ${logfile/\%j/$jid}"
  fi

  echo -e "wave${w}\t$action\t$jid\t${CPUS}\t${MEM}\t${wtime}\t${consol}\t$map\t${logfile/\%j/$jid}" >> "$MANIFEST"
  SUBMITTED+=("$jid")
  dep="$jid"
done

echo
if $DRY; then
  echo "[i]  dry run — nothing submitted. Manifest: $MANIFEST"
else
  echo "[i]  Monitor:  squeue -u \$USER | grep JVC-GDBI-${TAG}-${CHROM}"
  echo "[i]            squeue -j $(IFS=,; echo "${SUBMITTED[*]}")"
  echo "[i]  Logs:     $LOG_DIR/${TAG}-${CHROM}-w*-<jobid>.out"
  echo "[i]  Timings:  bash $(basename "$0") report $MANIFEST      # after all $N_WAVES finish"
  echo "[i]  A failed wave leaves the rest PENDING (DependencyNeverSatisfied);"
  echo "[i]  fix it, then resubmit from that wave with action=update."
fi
