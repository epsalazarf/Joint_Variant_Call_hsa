#!/usr/bin/env bash
# =============================================================================
# Title : Step 05 (GenotypeGVCFs) per-chromosome SLURM launcher
# About : Submits one sbatch job per requested chromosome, each running
#         05_gatk_GenotypeGVCFs.sh with GENO_PHASE=scatter (genotype only, no
#         gather — see that script's header for why running the combined
#         GENO_PHASE=all per chromosome from parallel jobs would corrupt the
#         cohort VCF). Once every scatter job succeeds, one final job runs
#         GENO_PHASE=gather (afterok all scatter jobs) to concatenate the
#         per-chromosome VCFs into the cohort-level callset. Chromosomes run
#         independently — there is no wave/incremental concept here, unlike
#         Step 04's JAGUAR launcher (test/jaguar/run_jaguar_waves.sh), which
#         this script's option/manifest/report conventions otherwise mirror.
#
# Usage :
#   bash GENOTYPE.seq_batch-slurmer.sh [options] <genomicsdb_path> <output_path> <cohort_name>
#   bash GENOTYPE.seq_batch-slurmer.sh --dry-run [options] <genomicsdb_path> <output_path> <cohort_name>
#   bash GENOTYPE.seq_batch-slurmer.sh report <manifest_file>
#
#   -c, --chroms LIST   chr1..chr22 | chrX | chrY | chrM | autosomes | all | a
#                       comma-separated list of the above (default: all)
#       --cpus  N       --cpus-per-task per scatter job (default: 2 — GenotypeGVCFs
#                       is effectively single-threaded per contig)
#       --mem   SIZE    --mem per scatter job (default: chromosome-scaled, see below)
#       --hours H       --time in hours per scatter job (default: chromosome-scaled)
#       --gather-mem SIZE    --mem for the final gather job (default: 8G)
#       --gather-hours H     --time in hours for the final gather job (default: 2)
#
#   Memory/time class boundaries REUSE Step 04's measured GenomicsDBImport
#   contig-size reasoning (chr1/chr2 biggest, chr3-8+chrX mid, rest small) as a
#   starting point — GenotypeGVCFs' own footprint is UNMEASURED (see Step 05's
#   STATUS block). Treat these as placeholders and override with --mem/--hours
#   once a real run gives actual numbers to work from:
#     chr1, chr2                 -> 64G / 20h
#     chr3-8, chrX               -> 32G / 12h
#     chr9-22, chrY, chrM        -> 16G /  6h
#   The gather job (bcftools concat + a sample-set consistency check) is
#   expected to be far lighter than any scatter job; its 8G/2h default is also
#   unmeasured.
#
#   genomicsdb_path   the Step 04 output_path (workspaces under .../genomicsdb/<CHR>)
#   output_path       where Step 05 writes chrom_vcf/, logs/, and the gathered VCF
#   cohort_name       label used in output filenames
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S05="$REPO/bin/05_gatk_GenotypeGVCFs.sh"

# ---------- report mode (before option parsing) ------------------------------
if [[ "${1:-}" == "report" ]]; then
  manifest="${2:?Usage: $(basename "$0") report <manifest_file>}"
  [[ -s "$manifest" ]] || { echo "[X]  manifest not found: $manifest"; exit 1; }
  echo "[i]  manifest: $manifest"; echo
  jids=$(awk -F'\t' '/^(scatter|gather)/{print $3}' "$manifest" | paste -sd, -)
  [[ -n "$jids" ]] || { echo "[X]  no job ids in manifest"; exit 1; }
  if command -v sacct >/dev/null; then
    sacct -j "$jids" --units=G \
      --format=JobID%14,JobName%24,State%12,Elapsed%12,TotalCPU%12,MaxRSS%10,AllocCPUS%9,ReqMem%8,Start%20,End%20
  else
    echo "[!]  sacct not available; raw manifest:"; cat "$manifest"
  fi
  exit 0
fi

# ---------- chromosome-class sizing (borrowed from Step 04's JAGUAR launcher; UNMEASURED for this step) ----
chrom_class() {
  case "$1" in
    chr1|chr2)                          echo big ;;
    chr3|chr4|chr5|chr6|chr7|chr8|chrX) echo mid ;;
    *)                                  echo small ;;
  esac
}
default_mem()   { case "$(chrom_class "$1")" in big) echo 64G ;; mid) echo 32G ;; *) echo 16G ;; esac; }
default_hours() { case "$(chrom_class "$1")" in big) echo 20  ;; mid) echo 12  ;; *) echo 6  ;; esac; }
fmt_time()      { printf '%d:00:00' "$1"; }   # H -> H:00:00

# ---------- options ------------------------------------------------------------
CHROM_SEL="all" ; CPUS=2 ; MEM="" ; HOURS="" ; GATHER_MEM="8G" ; GATHER_HOURS=2 ; DRY=false
pos=()
while (( $# )); do
  case "$1" in
    --dry-run)       DRY=true; shift ;;
    -c|--chroms)     CHROM_SEL="${2:?}"; shift 2 ;;
    --cpus)          CPUS="${2:?}"; shift 2 ;;
    --mem)           MEM="${2:?}"; shift 2 ;;
    --hours)         HOURS="${2:?}"; shift 2 ;;
    --gather-mem)    GATHER_MEM="${2:?}"; shift 2 ;;
    --gather-hours)  GATHER_HOURS="${2:?}"; shift 2 ;;
    -h|--help)       sed -n '3,43p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
    -*)              echo "[X]  unknown option: $1  (--help)"; exit 2 ;;
    *)               pos+=("$1"); shift ;;
  esac
done
GDB_PATH="${pos[0]:?Usage: $(basename "$0") [options] <genomicsdb_path> <output_path> <cohort_name>}"
OUTPUT_PATH="${pos[1]:?Usage: $(basename "$0") [options] <genomicsdb_path> <output_path> <cohort_name>}"
COHORT="${pos[2]:?Usage: $(basename "$0") [options] <genomicsdb_path> <output_path> <cohort_name>}"

# ---------- checks --------------------------------------------------------
[[ -f "$S05" ]] || { echo "[X]  Step 05 script not found: $S05"; exit 1; }
[[ -d "$GDB_PATH/genomicsdb" ]] || { echo "[X]  No genomicsdb/ under: $GDB_PATH (run Step 04 first)"; exit 1; }
GDB_PATH="$(readlink -f "$GDB_PATH")"
$DRY || command -v sbatch >/dev/null || { echo "[X]  sbatch not found"; exit 1; }
mkdir -p "$OUTPUT_PATH"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"
LOG_DIR="${OUTPUT_PATH}/logs"
mkdir -p "$LOG_DIR"

# ---------- resolve the chromosome list (bash-side mirror of Step 05's own selector) ----
AUTOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
ALLCHROMS=("${AUTOSOMES[@]}" chrX chrY chrM)
case "$CHROM_SEL" in
  autosomes) CHROMS=("${AUTOSOMES[@]}") ;;
  all)       CHROMS=("${ALLCHROMS[@]}") ;;
  *,*)       IFS=',' read -r -a CHROMS <<< "$CHROM_SEL" ;;
  *)         CHROMS=("$CHROM_SEL") ;;
esac
for c in "${CHROMS[@]}"; do
  case "$c" in
    chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY|chrM|chrMT) ;;
    *) echo "[X]  Invalid chromosome: '$c'"; exit 1 ;;
  esac
done
(( ${#CHROMS[@]} >= 1 )) || { echo "[X]  Empty chromosome selection: '$CHROM_SEL'"; exit 1; }
CHROM_CSV=$(IFS=','; echo "${CHROMS[*]}")   # exact set handed to the gather-phase invocation

# ---------- submit: one scatter job per chromosome ----------------------------
TS="$(date +%Y%m%d-%H%M%S)"
MANIFEST="$LOG_DIR/genotype_${COHORT}_run_${TS}.manifest"
{
  echo -e "# Step 05 (GenotypeGVCFs) scatter+gather run\t$TS"
  echo -e "# genomicsdb_path\t$GDB_PATH"
  echo -e "# output_path\t$OUTPUT_PATH"
  echo -e "# cohort\t$COHORT"
  echo -e "# chroms\t$CHROM_CSV"
  echo -e "# s05\t$S05"
  echo -e "# columns: phase\tchrom\tjobid\tcpus\tmem\ttime\tlogfile"
} > "$MANIFEST"

echo "[&]  ${COHORT} — Step 05 scatter+gather launcher   ($TS)"
echo "[i]  genomicsdb_path : $GDB_PATH"
echo "[i]  output_path     : $OUTPUT_PATH"
echo "[i]  cohort          : $COHORT"
echo "[i]  chromosomes     : ${#CHROMS[@]}  ($CHROM_CSV)"
echo "[i]  per-scatter job : ${CPUS} CPU; mem/time chromosome-scaled unless --mem/--hours given"
echo "[i]  gather job      : ${GATHER_MEM} / ${GATHER_HOURS}h"
echo "[i]  logs            : $LOG_DIR"
echo "[i]  manifest        : $MANIFEST"
echo "[i]  dry-run         : $DRY"
echo

SCATTER_JIDS=()
for chrom in "${CHROMS[@]}"; do
  mem="${MEM:-$(default_mem "$chrom")}"
  hours="${HOURS:-$(default_hours "$chrom")}"
  wtime="$(fmt_time "$hours")"
  jobname="JVC-GENO-${COHORT}-${chrom}"
  logfile="$LOG_DIR/genotype-${chrom}-%j.out"

  wrap="env GENO_PHASE=scatter bash '$S05' '$GDB_PATH' '$OUTPUT_PATH' '$chrom' '$COHORT'"

  set -- --parsable --job-name="$jobname" \
         --nodes=1 --ntasks=1 --cpus-per-task="$CPUS" --mem="$mem" --time="$wtime" \
         --output="$logfile"

  echo "[*]  scatter $chrom  (${CPUS} CPU / ${mem} / ${wtime})"

  if $DRY; then
    echo "       sbatch $* --wrap \"$wrap\""
    jid="DRYRUN-${chrom}"
  else
    jid=$(sbatch "$@" --wrap "$wrap")
    echo "       submitted: job $jid   log: ${logfile/\%j/$jid}"
  fi

  echo -e "scatter\t$chrom\t$jid\t${CPUS}\t${mem}\t${wtime}\t${logfile/\%j/$jid}" >> "$MANIFEST"
  SCATTER_JIDS+=("$jid")
done

# ---------- submit: one gather job, afterok every scatter job -----------------
dep=$(IFS=':'; echo "${SCATTER_JIDS[*]}")
jobname="JVC-GENO-${COHORT}-gather"
logfile="$LOG_DIR/genotype-gather-%j.out"
wrap="env GENO_PHASE=gather bash '$S05' '$GDB_PATH' '$OUTPUT_PATH' '$CHROM_CSV' '$COHORT'"

set -- --parsable --job-name="$jobname" \
       --nodes=1 --ntasks=1 --cpus-per-task=2 --mem="$GATHER_MEM" --time="$(fmt_time "$GATHER_HOURS")" \
       --output="$logfile" --dependency="afterok:$dep"

echo
echo "[*]  gather  (afterok: ${SCATTER_JIDS[*]})"

if $DRY; then
  echo "       sbatch $* --wrap \"$wrap\""
  gather_jid="DRYRUN-gather"
else
  gather_jid=$(sbatch "$@" --wrap "$wrap")
  echo "       submitted: job $gather_jid   log: ${logfile/\%j/$gather_jid}"
fi
echo -e "gather\t-\t$gather_jid\t2\t${GATHER_MEM}\t$(fmt_time "$GATHER_HOURS")\t${logfile/\%j/$gather_jid}" >> "$MANIFEST"

echo
echo "[i]  manifest : $MANIFEST"
echo "[i]  report   : bash '$0' report '$MANIFEST'"

#EOF
