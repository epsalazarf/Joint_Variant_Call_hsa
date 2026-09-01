#!/usr/bin/env bash

# =============================================================================
# Title       : GATK GenomicsDBImport [FENIX]
# Description : Consolidates per-sample, per-chromosome GVCFs (Step 03a output)
#               into GenomicsDB workspaces — one workspace per chromosome — as
#               the input for joint genotyping (Step 05). Adapted for LAVIS-FENIX.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-08-31
# Version     : 1.0
# Usage       : 04_gatk_GenomicsDB_import.sh <sample_map> <output_path> <chrom> [action]
#             : sample_map  — TSV, 2 columns: <sample_label> <TAB> <chrom_gvcf_dir>
#             : output_path — workspaces are created at <output_path>/genomicsdb/<chrom>
#             :               MUST be group storage (/mnt/data/...), never $HOME
#             :               (FENIX home quota ~5 GB); the script refuses a $HOME path
#             : chrom       — chr1..chr22 | chrX | chrY | chrM | autosomes | all
#             : action      — create (default) | update
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
# =============================================================================
#
# ---------------------------------------------------------------------------
#  SAMPLE MAP — what it is and how to generate it
# ---------------------------------------------------------------------------
#
# A plain text file, one line per sample, TWO columns separated by a TAB:
#
#     LUPUS001<TAB>/mnt/data/amedina/esalazarf/lupus/LUPUS001/chrom_gvcf
#     LUPUS002<TAB>/mnt/data/amedina/esalazarf/lupus/LUPUS002/chrom_gvcf
#     ...
#
#   Column 1  A readable label for the sample. Used for cohort membership and
#             error messages. The REAL sample name stored in the database is
#             read from each GVCF header (`bcftools query -l`); if the label and
#             the header disagree, the script warns and keeps the header name.
#   Column 2  The per-sample `chrom_gvcf/` directory produced by Step 03a
#             (03_gatk_haplotype_caller.sh). For each chromosome the script
#             picks the file matching  *.<chrom>.g.vcf.gz  inside that folder.
#
# Generate it from a cohort directory that holds one sub-directory per sample,
# each containing a `chrom_gvcf/` folder:
#
#     cd /path/to/cohort
#     for d in */chrom_gvcf; do
#       s=$(basename "$(dirname "$d")")
#       printf '%s\t%s\n' "$s" "$(readlink -f "$d")"
#     done > cohort.sample_map.tsv
#
# Lines that are blank or start with '#' are ignored. Always eyeball the file
# before launching.
#
# ---------------------------------------------------------------------------
#  INCREMENTAL IMPORT — samples arriving in waves (e.g. lupus cohort)
# ---------------------------------------------------------------------------
#
#   First wave :  ... <map_wave1.tsv>  <out>  chr22  create
#   Later waves:  ... <map_wave2.tsv>  <out>  chr22  update
#
#   `update` calls GATK --genomicsdb-update-workspace-path and adds ONLY the
#   samples listed in the new map to the existing per-chromosome workspaces,
#   in place. Do NOT re-list samples that are already imported. The set of
#   intervals (chromosomes) cannot change on update — a workspace that does
#   not exist yet must be built with `create` first.
#
# ---------------------------------------------------------------------------
#  FENIX resourcing (per-chromosome job)
# ---------------------------------------------------------------------------
#
#   MEASURED (seff 276894: 31 samples, chr22, 8 cores): CPU efficiency 9.53%,
#   6.68 GB RAM, 72 min wall. For one interval GenomicsDBImport is effectively
#   serial (the tail `Consolidating` step) + I/O-bound — extra cores/RAM do not
#   help. Request small; scale --time with contig size:
#
#     --cpus-per-task=2   --mem=16G   (-Xmx6g; ~1.5 h for chr22, the smallest)
#
#   The TileDB workspace does many small random reads/writes — bad on NFS — so
#   this script builds it on /scratch when available and copies the finished
#   workspace back to <output_path> (see Sept-2026 scratch notes in
#   INSTRUCTIONS.md: unlike Step 02, this step genuinely benefits from scratch).
# =============================================================================

set -euo pipefail

# <ARGUMENTS> ---------------------------------------------------------------

SAMPLE_MAP="${1:?Usage: $(basename "$0") <sample_map> <output_path> <chrom> [create|update]}"
OUTPUT_PATH="${2:?Usage: $(basename "$0") <sample_map> <output_path> <chrom> [create|update]}"
CHROM_ARG="${3:?Usage: $(basename "$0") <sample_map> <output_path> <chrom> [create|update]}"
ACTION="${4:-create}"

case "$ACTION" in
  create|update) ;;
  *) echo "[X]  Invalid action: '$ACTION' (expected: create | update)"; exit 1 ;;
esac

# <\ARGS> -----------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GATK GenomicsDBImport [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking inputs..."

[ -s "$SAMPLE_MAP" ] || { echo "[X]  CANCELLED. Sample map not found or empty: $SAMPLE_MAP"; exit 1; }
SAMPLE_MAP="$(readlink -f "$SAMPLE_MAP")"
mkdir -p "$OUTPUT_PATH"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"

n_samples=$(grep -cvE '^[[:space:]]*(#|$)' "$SAMPLE_MAP" || true)
[ "${n_samples:-0}" -ge 1 ] || { echo "[X]  CANCELLED. No usable sample lines in: $SAMPLE_MAP"; exit 1; }

echo "[<]  Sample map : $SAMPLE_MAP  ($n_samples samples)"
echo "[i]  Output     : ${OUTPUT_PATH}/genomicsdb/<chrom>"
echo "[i]  Action     : ${ACTION}"

# Guard: the finished workspace is copied back to OUTPUT_PATH. A home directory
# on FENIX has a hard ~5 GB quota — a per-chromosome workspace for a real cohort
# will blow through it and the copy-back fails AFTER the (long) import. Refuse a
# home-dir output, and do a real write test on the target filesystem.
case "$OUTPUT_PATH/" in
  "${HOME}"/*)
    echo "[X]  CANCELLED. Output is under \$HOME (${HOME}) — home has a hard quota on FENIX."
    echo "[i]  Point it at group storage, e.g. /mnt/data/amedina/${USER}/JVCdev/<name>"
    exit 1 ;;
esac
if ! { mkdir -p "$OUTPUT_PATH/genomicsdb" && _wt="$OUTPUT_PATH/genomicsdb/.writetest.$$" \
       && dd if=/dev/zero of="$_wt" bs=1M count=64 status=none 2>/dev/null; }; then
  rm -f "${_wt:-}" 2>/dev/null || true
  echo "[X]  CANCELLED. Cannot write ~64 MB to ${OUTPUT_PATH}/genomicsdb (quota or permissions)."
  df -h "$OUTPUT_PATH" 2>/dev/null | sed 's/^/[i]    /'
  exit 1
fi
rm -f "$_wt"
avail_kb=$(df -Pk "$OUTPUT_PATH" 2>/dev/null | awk 'NR==2{print $4}')
[ -n "${avail_kb:-}" ] && echo "[i]  Output free : $(( avail_kb / 1024 )) MB on $(df -Ph "$OUTPUT_PATH" 2>/dev/null | awk 'NR==2{print $6}')"

# Options (env-overridable so a launcher can vary resources per run)
#   GENDBI_READER_THREADS  — GenomicsDBImport --reader-threads
#                            (default: $SLURM_CPUS_PER_TASK, else 2)
#   GENDBI_BATCH_SIZE      — GenomicsDBImport --batch-size            (default 50)
#   GENDBI_JAVA_MEM        — Java -Xms/-Xmx                (default 6G remote / 4G local)
#   GENDBI_STRICT_GVCF     — true: abort on a GVCF that looks truncated;
#                            false (default): warn and keep it in the cohort
njobs="${GENDBI_READER_THREADS:-${SLURM_CPUS_PER_TASK:-2}}"
BATCH_SIZE="${GENDBI_BATCH_SIZE:-50}"
STRICT_GVCF="${GENDBI_STRICT_GVCF:-false}"
TRUNCATED_SEEN=""   # labels of samples whose GVCF looked truncated (warn mode)
SMALL_SEEN=""       # labels of samples whose GVCF is a valid but tiny outlier

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
  MEM="${GENDBI_JAVA_MEM:-6G}"   # seff 276894: whole job used 6.68 GB RSS; heap need is small
else
  env_type="local"
  MEM="${GENDBI_JAVA_MEM:-4G}"
fi
echo "[i]  Environment: $env_type"
echo "[i]  Knobs      : reader-threads=${njobs}  batch-size=${BATCH_SIZE}  java-mem=${MEM}  strict-gvcf=${STRICT_GVCF}"

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

# Load modules on remote (work-around due to faulty parser [ARC02])
if [[ "$env_type" == "remote" ]]; then
  echo "[i]  Loading modules..."
  module unload oracle-java 2>/dev/null || true
  module load oracle-java/25.0.2
  module load gatk
  module load bcftools
fi

# Guard: reference genome (GenomicsDBImport uses it to validate contigs)
if [ -z "${ref_gnm:-}" ]; then
  echo "[X]  Missing ref_gnm in config: $CONFIG_FILE"
  exit 1
fi
[ -f "${ref_gnm}" ] || { echo "[X]  Reference genome not found: ${ref_gnm}"; exit 1; }
echo "[i]  Reference  : ${ref_gnm}"

# Scratch: this step DOES benefit from scratch (random-access TileDB I/O).
# Build the workspace on scratch, copy the finished product back to OUTPUT_PATH.
if [[ -n "${scratch_root:-}" && "$env_type" == "remote" ]]; then
  USE_SCRATCH=true
  SCRATCH_BASE="${scratch_root}/${USER}/gendbi_${SLURM_JOB_ID:-$$}_${EPOCHSECONDS}"
else
  USE_SCRATCH=false
  SCRATCH_BASE="${OUTPUT_PATH}/.gendbi_work_${SLURM_JOB_ID:-$$}_${EPOCHSECONDS}"
fi
mkdir -p "$SCRATCH_BASE/tmp" "$SCRATCH_BASE/gdb" "$SCRATCH_BASE/maps"
export TMPDIR="$SCRATCH_BASE/tmp"
echo "[i]  Scratch    : ${USE_SCRATCH}  (${SCRATCH_BASE})"

# Cleanup: wipe scratch on success. On failure, KEEP it only if GenomicsDBImport
# actually produced a workspace there (so a copy-back failure — e.g. quota — does
# not throw away the expensive compute; the finisher prints how to recover).
GDB_BUILT=false
cleanup() {
  local rc=$?
  if (( rc != 0 )) && [[ "$GDB_BUILT" == true && -d "$SCRATCH_BASE/gdb" ]] \
       && find "$SCRATCH_BASE/gdb" -name callset.json -print -quit 2>/dev/null | grep -q .; then
    echo
    echo "[!]  Exited $rc AFTER GenomicsDBImport succeeded — scratch workspace PRESERVED:"
    echo "[!]    $SCRATCH_BASE/gdb"
    echo "[!]  Recover WITHOUT recomputing, once the destination has room:"
    echo "[!]    mkdir -p '${OUTPUT_PATH}/genomicsdb'"
    echo "[!]    cp -r '$SCRATCH_BASE/gdb/'* '${OUTPUT_PATH}/genomicsdb/'"
    echo "[!]    bash '$0' '$SAMPLE_MAP' '$OUTPUT_PATH' '$CHROM_ARG' create   # verifies + finishes"
    echo "[!]  (scratch auto-purges after a few days — copy it out soon)"
  else
    rm -rf "$SCRATCH_BASE"
  fi
}
trap cleanup EXIT

# <\ENV> --------------------------------------------------------------------

# <FUNCTIONS> ---------------------------------------------------------------

AUTOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
           chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)
ALLCHROMS=("${AUTOSOMES[@]}" chrX chrY chrM)

## Resolve the requested selector into a list of chromosomes.
resolve_chrom_list() {
  case "$CHROM_ARG" in
    autosomes) CHROMS=("${AUTOSOMES[@]}") ;;
    all)       CHROMS=("${ALLCHROMS[@]}") ;;
    chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY|chrM|chrMT)
               CHROMS=("$CHROM_ARG") ;;
    *) echo "[X]  Invalid chrom selector: '$CHROM_ARG'"
       echo "[i]  Expected: chr1..chr22 | chrX | chrY | chrM | autosomes | all"
       exit 1 ;;
  esac
}

## Pick the mitochondrial contig notation actually present in a chrom_gvcf dir.
## hg38 (UCSC) uses 'chrM'; some call sets use 'chrMT'. Returns the effective
## contig name for the workspace + the -L interval.
effective_chrom() {
  local want="$1" dir="$2"
  if [[ "$want" != "chrM" && "$want" != "chrMT" ]]; then
    echo "$want"; return 0
  fi
  shopt -s nullglob
  local m=( "$dir"/*.chrM.g.vcf.gz )
  local mt=( "$dir"/*.chrMT.g.vcf.gz )
  shopt -u nullglob
  if (( ${#m[@]} )); then
    echo "chrM"
  elif (( ${#mt[@]} )); then
    echo "[!]  WARNING: data uses 'chrMT' notation for the mitochondrion" >&2
    echo "[!]           workspace + interval will use 'chrMT' to match" >&2
    echo "chrMT"
  else
    echo "$want"   # let the downstream guard produce the missing-file error
  fi
}

## Cheap integrity check: a well-formed BGZF (bgzip) file ends with a fixed
## 28-byte empty-block EOF marker. Its absence means the file was truncated
## (interrupted copy/transfer, killed writer, full disk). Same test htslib uses.
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
gvcf_looks_truncated() {
  local f="$1" tail_hex
  tail_hex="$(tail -c 28 "$f" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')"
  [[ "$tail_hex" != "$BGZF_EOF" ]]
}

## Build the per-chromosome sample-name-map GATK consumes.
##   $1 effective chrom   $2 destination map file
build_chrom_sample_map() {
  local chr="$1" dest="$2"
  : > "$dest"
  local seen=" "        # space-delimited list of sample names already added
  local count=0

  while IFS=$'\t' read -r label dir _rest || [[ -n "${label:-}" ]]; do
    [[ -z "${label// }" || "$label" == \#* ]] && continue
    [[ -z "${dir:-}" ]] && { echo "[X]  Malformed sample-map line (need 2 tab-separated columns): '$label'"; exit 1; }

    dir="${dir%/}"
    [[ -d "$dir" ]] || { echo "[X]  chrom_gvcf dir not found for '$label': $dir"; exit 1; }

    shopt -s nullglob
    local hits=( "$dir"/*."$chr".g.vcf.gz )
    shopt -u nullglob

    (( ${#hits[@]} == 0 )) && { echo "[X]  No '*.${chr}.g.vcf.gz' for sample '$label' in: $dir"; exit 1; }
    (( ${#hits[@]} > 1 ))  && { echo "[X]  Ambiguous — ${#hits[@]} files match '*.${chr}.g.vcf.gz' for '$label' in: $dir"; exit 1; }

    local gvcf="${hits[0]}"
    [[ -s "${gvcf}.tbi" || -s "${gvcf}.csi" ]] || { echo "[X]  Missing index (.tbi/.csi) for: $gvcf"; exit 1; }

    # Truncation guard
    if gvcf_looks_truncated "$gvcf"; then
      if [[ "$STRICT_GVCF" == "true" ]]; then
        echo "[X]  Truncated GVCF (no BGZF EOF marker) for sample '$label': $gvcf"
        echo "[i]  Re-run Step 03a for this sample, or drop it from the map. (GENDBI_STRICT_GVCF=true)"
        exit 1
      fi
      echo "[!]  WARNING: '$label' GVCF has no BGZF EOF marker — likely TRUNCATED: $gvcf"
      echo "[!]           $(ls -l "$gvcf" | awk '{print $5" bytes"}'); kept in cohort (GENDBI_STRICT_GVCF=false)."
      echo "[!]           GenomicsDBImport will most likely fail on this file."
      TRUNCATED_SEEN="${TRUNCATED_SEEN}${label} "
    fi

    local sm
    sm="$(bcftools query -l "$gvcf" 2>/dev/null)"
    [[ "$(wc -l <<< "$sm")" -eq 1 && -n "$sm" ]] || { echo "[X]  Expected exactly one sample in GVCF header: $gvcf"; exit 1; }

    [[ "$sm" != "$label" ]] && echo "[!]  WARNING: map label '$label' != GVCF sample '$sm' — using '$sm'"
    case "$seen" in *" ${sm} "*) echo "[X]  Duplicate sample '$sm' in cohort (last seen as label '$label')"; exit 1 ;; esac
    seen="${seen}${sm} "

    printf '%s\t%s\n' "$sm" "$gvcf" >> "$dest"
    (( ++count ))
  done < "$SAMPLE_MAP"

  (( count >= 1 )) || { echo "[X]  Empty sample map generated for $chr"; exit 1; }

  # Size-outlier check: a valid-but-tiny GVCF (passes the BGZF EOF test) still
  # points at a near-empty sample — usually a bad BAM upstream. Warn, don't block.
  local sizes=() gv gs med cutoff
  while IFS=$'\t' read -r _ gv; do
    gs=$(stat -c%s "$gv" 2>/dev/null || stat -f%z "$gv" 2>/dev/null || echo 0)
    sizes+=("$gs")
  done < "$dest"
  med=$(printf '%s\n' "${sizes[@]}" | sort -n | awk '{v[NR]=$1} END{ if(NR) print (NR%2)? v[(NR+1)/2] : int((v[NR/2]+v[NR/2+1])/2) }')
  if [[ -n "${med:-}" ]] && (( med > 0 )); then
    cutoff=$(( med / 2 ))
    while IFS=$'\t' read -r sm gv; do
      gs=$(stat -c%s "$gv" 2>/dev/null || stat -f%z "$gv" 2>/dev/null || echo 0)
      if (( gs < cutoff )); then
        echo "[!]  WARNING: '$sm' GVCF unusually small ($gs b vs ${med} b median) — near-empty sample? check its ${chr} coverage."
        SMALL_SEEN="${SMALL_SEEN}${sm} "
      fi
    done < "$dest"
  fi

  echo "[i]    $chr: $count samples -> $dest"
}

## Import one chromosome into its own GenomicsDB workspace.
import_one_chrom() {
  local req="$1"
  local step_timestamp=$EPOCHSECONDS

  # Effective contig notation is decided from the first sample's dir.
  local first_dir
  first_dir="$(grep -vE '^[[:space:]]*(#|$)' "$SAMPLE_MAP" | head -1 | cut -f2 || true)"
  first_dir="${first_dir%/}"
  local chr
  chr="$(effective_chrom "$req" "$first_dir")"

  local ws_final="${OUTPUT_PATH}/genomicsdb/${chr}"
  local ws_work="${SCRATCH_BASE}/gdb/${chr}"
  local map="${SCRATCH_BASE}/maps/${chr}.sample_map.tsv"

  echo
  echo "[*]  Chromosome ${req}$( [[ "$chr" != "$req" ]] && echo "  (effective: ${chr})" )"
  echo "[&]  $(date +%Y%m%d-%H%M)"

  # --- skip / precondition logic ---
  if [[ "$ACTION" == "create" ]]; then
    if [[ -s "${ws_final}/callset.json" && -s "${ws_final}/vidmap.json" ]]; then
      echo "[i]  Already completed (workspace exists: ${ws_final})"
      return 0
    fi
  else # update
    if [[ ! -s "${ws_final}/callset.json" ]]; then
      echo "[X]  action=update but no workspace to update: ${ws_final}"
      echo "[i]  Build it first with:  ... <map> ${OUTPUT_PATH} ${req} create"
      exit 1
    fi
  fi

  # --- build the per-chromosome sample map ---
  build_chrom_sample_map "$chr" "$map"

  # --- stage workspace (update: copy existing to scratch) ---
  rm -rf "$ws_work"
  local mode_opts
  if [[ "$ACTION" == "create" ]]; then
    mode_opts=(--genomicsdb-workspace-path "$ws_work" -L "$chr")
  else
    echo "[i]  Staging existing workspace to scratch for update..."
    cp -r "$ws_final" "$ws_work"
    mode_opts=(--genomicsdb-update-workspace-path "$ws_work")
  fi

  # --- run GenomicsDBImport ---
  # (JDK-25 "restricted method" / sun.misc.Unsafe warnings from GATK 4.6's native
  #  libs are harmless — GATK 4.6 targets JDK 17. Not worth extra --java-options.)
  set -o xtrace
  gatk --java-options "-Xms${MEM} -Xmx${MEM}" GenomicsDBImport \
    --sample-name-map "$map" \
    --reference "$ref_gnm" \
    --tmp-dir "$TMPDIR" \
    --batch-size "$BATCH_SIZE" \
    --reader-threads "$njobs" \
    --genomicsdb-shared-posixfs-optimizations true \
    --consolidate \
    --verbosity ERROR \
    "${mode_opts[@]}"
  set +o xtrace

  [[ -s "${ws_work}/callset.json" && -s "${ws_work}/vidmap.json" ]] \
    || { echo "[X]  CANCELLED: GenomicsDBImport failed for ${chr}, workspace incomplete: ${ws_work}"; exit 1; }
  GDB_BUILT=true   # from here on, a failure must not throw away the scratch workspace

  # --- copy finished workspace back (stage via .new, then swap) ---
  local ws_bytes ws_avail_kb
  ws_bytes=$(du -sk "$ws_work" 2>/dev/null | awk '{print $1}')
  mkdir -p "${OUTPUT_PATH}/genomicsdb"
  ws_avail_kb=$(df -Pk "${OUTPUT_PATH}/genomicsdb" 2>/dev/null | awk 'NR==2{print $4}')
  echo "[i]  Workspace ${chr}: ~$(( ${ws_bytes:-0} / 1024 )) MB to copy back; ~$(( ${ws_avail_kb:-0} / 1024 )) MB free at destination"
  if [[ -n "${ws_bytes:-}" && -n "${ws_avail_kb:-}" ]] && (( ws_bytes + ws_bytes/10 > ws_avail_kb )); then
    echo "[X]  CANCELLED: not enough room at ${OUTPUT_PATH}/genomicsdb for the ${chr} workspace."
    echo "[!]  GenomicsDBImport already succeeded — recover from scratch (see below), do NOT recompute."
    exit 1
  fi

  local ws_new="${OUTPUT_PATH}/genomicsdb/.${chr}.new.$$"
  rm -rf "$ws_new"
  if ! cp -r "$ws_work" "$ws_new"; then
    rm -rf "$ws_new"
    echo "[X]  CANCELLED: copy-back failed for ${chr} (quota/space/permissions)."
    exit 1
  fi
  rm -rf "$ws_final"
  mv "$ws_new" "$ws_final"
  rm -rf "$ws_work"

  [[ -s "${ws_final}/callset.json" ]] || { echo "[X]  CANCELLED: copy-back verification failed for ${chr}: ${ws_final}"; exit 1; }
  local n_this; n_this="$(grep -cvE '^[[:space:]]*$' "$map" || true)"
  echo "[>]  ${ws_final}"
  echo "[!]  Chromosome ${chr} — ${ACTION} done (${n_this} samples $( [[ "$ACTION" == update ]] && echo "added" || echo "imported" ))"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Finisher: verify every requested workspace and report.
finisher() {
  echo
  local missing=0
  for req in "${CHROMS[@]}"; do
    local d="${OUTPUT_PATH}/genomicsdb/${req}"
    # MT notation fallback — only for the mitochondrion
    if [[ ! -d "$d" && ( "$req" == chrM || "$req" == chrMT ) ]]; then
      [[ -d "${OUTPUT_PATH}/genomicsdb/chrMT" ]] && d="${OUTPUT_PATH}/genomicsdb/chrMT"
      [[ -d "${OUTPUT_PATH}/genomicsdb/chrM"  ]] && d="${OUTPUT_PATH}/genomicsdb/chrM"
    fi
    if [[ -s "${d}/callset.json" ]]; then
      echo "[>]  ${d}"
    else
      echo "[X]  MISSING workspace for ${req}"
      missing=1
    fi
  done

  if [[ -n "$TRUNCATED_SEEN" ]]; then
    echo
    echo "[!]  Truncated-looking GVCF(s) were kept in the cohort: ${TRUNCATED_SEEN}"
    echo "[!]  Re-run Step 03a for those samples and re-import with action=update,"
    echo "[!]  or rebuild the workspace without them."
  fi
  if [[ -n "$SMALL_SEEN" ]]; then
    echo
    echo "[!]  Unusually small (but valid) GVCF(s) imported: ${SMALL_SEEN}"
    echo "[!]  Likely near-empty samples — verify their coverage before Step 05."
  fi

  echo
  if (( missing == 0 )); then
    echo "[$] GATK GenomicsDBImport [FENIX] completed successfully!"
    echo "[i]  Next: Step 05 (GenotypeGVCFs) reads  gendb://<workspace>  per chromosome."
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 0
  else
    echo "[X]  GATK GenomicsDBImport [FENIX] INCOMPLETE — see MISSING entries above."
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 1
  fi
}

# <\FUNCTIONS> --------------------------------------------------------------

# <MAIN> ------------------------------------------------------------------------

main() {
  resolve_chrom_list

  echo
  echo "[i]  Chromosomes to process (${#CHROMS[@]}): ${CHROMS[*]}"
  [[ ${#CHROMS[@]} -gt 1 ]] && echo "[i]  Note: chromosomes run serially in this process — use the per-chromosome" \
                            && echo "[i]        SLURM launcher to scatter them across jobs for large cohorts."

  for req in "${CHROMS[@]}"; do
    import_one_chrom "$req"
  done

  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
