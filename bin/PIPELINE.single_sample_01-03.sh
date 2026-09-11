#!/bin/bash
# =============================================================================
# Title       : Single-Sample Sequential Pipeline Launcher [FENIX]
# Description : Submits Steps 01→02→03 as chained SLURM jobs for one sample.
#               Each step is held until the previous one succeeds (afterok).
#               When USE_SCRATCH=true, inputs are read directly from NFS while
#               TMPDIR + intermediate/output files use /scratch, then finals are
#               copied back (verified) to SAMPLE_DIR. This keeps bulky transient
#               intermediates off the group's persistent quota.
#
#               Scratch root is derived at RUN TIME from the running user's
#               PRIMARY group:  <scratch_base>/<primary-group>/<user>/job_<id>
#               so it works for amedina members and guests alike, with no
#               dependency on per-node supplementary-group propagation
#               (2026-09 incident: /scratch/groups/amedina was untraversable
#               from node09 for a guest → 36 jobs died instantly).
#
#               NOTE: benchmarking showed scratch gives no wall-time speedup for
#               these sequential GATK steps (CPU-bound) — its value is storage
#               hygiene and cluster citizenship, not speed. See INSTRUCTIONS.md.
# Usage       : PIPELINE.single_sample_01-03.sh [sample_dir]
#               sample_dir — directory named as the sample ID, containing FASTQ
#                            files. Defaults to $PWD.
# =============================================================================

set -euo pipefail

# <TOGGLES> -------------------------------------------------------------------
# Set USE_SCRATCH=true on FENIX (recommended). Set false only for testing or
# when the scratch filesystem is unavailable.
USE_SCRATCH=true

# SCATTER_S03=true : run Step 03 as a 25-way SLURM array (one HaplotypeCaller
#   task per canonical chromosome) writing straight into chrom_gvcf/, then a
#   gather job builds the whole-genome canon_chr GVCF. Per-sample wall time
#   ~2 h (chr1) instead of ~15 h, and a failure costs one chromosome, not the
#   sample. Set false for the legacy single whole-genome HaplotypeCaller job.
SCATTER_S03=true
ARRAY_CONC=6          # max chromosome tasks running at once per sample
# <\TOGGLES> ------------------------------------------------------------------

# Canonical chromosomes, in reference (.dict) order — array index N = Nth entry.
CANON_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
CANON_N=25

# <ARGUMENTS> -----------------------------------------------------------------

SAMPLE_DIR="$(realpath "${1:-$PWD}")"
SAMPLE_ID="$(basename "$SAMPLE_DIR")"

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
S01="${SCRIPT_DIR}/01_bwa_map_fastq_reads.sh"
S02="${SCRIPT_DIR}/02_gatk_bam_qc_workflow.sh"
S03="${SCRIPT_DIR}/03_gatk_haplotype_caller.sh"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

CONFIG_FILE="${SCRIPT_DIR}/../config/config.yaml"

if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
else
  env_type="local"
fi

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

: "${scratch_base:=}"
: "${sbatch_exclude:=}"
: "${min_free_gb:=0}"

# --exclude flag for every sbatch (broken-scratch nodes etc.); empty = no flag
EXCLUDE_ARG=()
[[ -n "$sbatch_exclude" ]] && EXCLUDE_ARG=(--exclude="$sbatch_exclude")

# <\ENVIRONMENT> --------------------------------------------------------------

# Output-presence helpers (needed by the CHECKS gate below, and again for the
# per-step resume decisions further down).
have_output() { find "$SAMPLE_DIR" -maxdepth 1 -name "$1" | grep -q .; }

# S01 output check: accepts both *.sort.bam (current convention) and *.sorted.bam
# (legacy naming seen on pre-existing SLE cohort dirs — same pair the S02 stage
# already treats as equivalent when it discovers inputs; keep this consistent
# with that, or a legacy dir with only .sorted.bam looks "unmapped" and gets
# S01 wrongly resubmitted even though S02's output already exists).
have_sort_bam() {
  find "$SAMPLE_DIR" -maxdepth 1 \( -name '*.sort.bam' -o -name '*.sorted.bam' \) | grep -q .
}

# <CHECKS> --------------------------------------------------------------------

for f in "$S01" "$S02" "$S03"; do
  [[ -f "$f" ]] || { echo "<ERROR> Script not found: $f"; exit 1; }
done

[[ -d "$SAMPLE_DIR" ]] || { echo "<ERROR> Sample directory not found: $SAMPLE_DIR"; exit 1; }

# FASTQs are only required when Step 01 would actually run. A resume-only call
# on a dir that already has a sorted BAM (or further downstream output) — e.g.
# archived/processed cohorts where the raw FASTQs live elsewhere entirely —
# must not be blocked by their absence.
if ! have_sort_bam && ! have_output "*.rmdup.mqfilt.bqsr.bam"; then
  if ! find "$SAMPLE_DIR" -maxdepth 1 -name "*.f*q.gz" | grep -qE '_R?1\.f[^.]*q\.gz$'; then
    echo "<ERROR> No FASTQ R1 files found in: $SAMPLE_DIR (and no existing S01/S02 output to resume from)"
    exit 1
  fi
fi

command -v sbatch &>/dev/null || { echo "<ERROR> sbatch not found — SLURM environment required."; exit 1; }

if [[ "$USE_SCRATCH" == true && -z "$scratch_base" ]]; then
  echo "<ERROR> USE_SCRATCH=true but scratch_base is not set in config (${env_type} section)."
  exit 1
fi

# Output-filesystem free-space guard (GNU df). Abort before submitting if the
# persistent store is low — a full FS silently truncates copy-backs and makes
# SLURM unable to even create job logs.
if [[ "${min_free_gb:-0}" -gt 0 ]]; then
  free_gb=$(df -BG --output=avail "$SAMPLE_DIR" 2>/dev/null | tail -1 | tr -dc '0-9')
  if [[ -n "$free_gb" && "$free_gb" -lt "$min_free_gb" ]]; then
    echo "<ERROR> Only ${free_gb}G free on the output filesystem (need ${min_free_gb}G). Aborting — free space or lower min_free_gb."
    exit 1
  fi
  echo "[i]  Output FS free : ${free_gb:-?}G  (min ${min_free_gb}G)"
fi

# <\CHECKS> -------------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

LOG_DIR="${SAMPLE_DIR}/log"
mkdir -p "$LOG_DIR"

echo
echo "[&] Single-Sample Pipeline Launcher [FENIX]"
echo "[i]  Sample     : $SAMPLE_ID"
echo "[i]  Dir        : $SAMPLE_DIR"
echo "[i]  Logs       : $LOG_DIR"
echo "[i]  Scratch    : ${USE_SCRATCH} $([ "$USE_SCRATCH" = true ] && echo "(${scratch_base}/<primary-group>/\$USER/job_\$SLURM_JOB_ID)" || true)"
[[ ${#EXCLUDE_ARG[@]} -gt 0 ]] && echo "[i]  Excluding   : $sbatch_exclude"
echo

# Helper: build a scratch-enabled wrap for a pipeline step.
#
# Strategy (see Sept-2026 benchmark in INSTRUCTIONS.md): inputs are read
# DIRECTLY from NFS (SAMPLE_DIR) — large, read once, sequential — so staging
# them to scratch is pure overhead. Scratch holds TMPDIR (tool spill — the real
# IOPS traffic) and intermediate + output files; only verified finals are
# copied back. Kill-safe: an EXIT trap wipes scratch even if the job is
# SIGTERM'd, so nothing leaks; the one exception is a FAILED copy-back
# verification, which preserves scratch so the data can be recovered.
#
#   $1  — pipeline script path (expands at submission time)
#   $2  — input mode: "fastq" | "sort_bams" | "bqsr_bam"
#   $3  — space-separated globs of output files to copy back (under scratch root)
#
# Submission-time expansion (from this script): $1, SAMPLE_DIR, scratch_base
# Job-runtime expansion (escaped \$): SLURM_JOB_ID, USER, SCRATCH_JOB, _grp ...
#
build_scratch_wrap() {
  local script="$1"
  local input_mode="$2"
  local copy_back_globs="$3"

  # -- resolve inputs (read from NFS) and point output at scratch --
  local stage_inputs=""
  case "$input_mode" in
    fastq)
      stage_inputs="cd '${SAMPLE_DIR}'
SCRIPT_ARGS=('${SAMPLE_DIR}' \"\$SCRATCH_JOB\")"
      ;;
    sort_bams)
      stage_inputs="
_nfs_list=\$(find '${SAMPLE_DIR}' -maxdepth 1 \( -name '*.sort.bam' -o -name '*.sorted.bam' \) | sort | paste -sd,)
[[ -z \"\$_nfs_list\" ]] && { echo '<ERROR> No *.sort.bam files in ${SAMPLE_DIR}'; exit 1; }
SCRIPT_ARGS=(\"\$_nfs_list\" \"\$SCRATCH_JOB\")"
      ;;
    bqsr_bam)
      stage_inputs="
_bqsr=\$(find '${SAMPLE_DIR}' -maxdepth 1 -name '*.rmdup.mqfilt.bqsr.bam' | head -1)
[[ -z \"\$_bqsr\" ]] && { echo '<ERROR> No *.rmdup.mqfilt.bqsr.bam in ${SAMPLE_DIR}'; exit 1; }
SCRIPT_ARGS=(\"\$_bqsr\" \"\$SCRATCH_JOB\")"
      ;;
    bqsr_bam_scatter)
      # one array task = one canonical chromosome (index = SLURM_ARRAY_TASK_ID)
      stage_inputs="
_bqsr=\$(find '${SAMPLE_DIR}' -maxdepth 1 -name '*.rmdup.mqfilt.bqsr.bam' | head -1)
[[ -z \"\$_bqsr\" ]] && { echo '<ERROR> No *.rmdup.mqfilt.bqsr.bam in ${SAMPLE_DIR}'; exit 1; }
_canon=(${CANON_CHROMS})
_chrom=\"\${_canon[\$((SLURM_ARRAY_TASK_ID - 1))]:-}\"
[[ -z \"\$_chrom\" ]] && { echo \"<ERROR> no chromosome for array index \$SLURM_ARRAY_TASK_ID\" >&2; exit 1; }
echo \"[i] array task \$SLURM_ARRAY_TASK_ID -> \$_chrom\"
SCRIPT_ARGS=(\"\$_bqsr\" \"\$SCRATCH_JOB\" \"\$_chrom\")"
      ;;
  esac

  # -- verified copy-back: one _cb call per matched file --
  local copy_back=""
  for glob in $copy_back_globs; do
    copy_back+="
for _f in \"\$SCRATCH_JOB\"/${glob}; do [ -e \"\$_f\" ] && _cb \"\$_f\" '${SAMPLE_DIR}'; done"
  done
  if [[ "$input_mode" == "fastq" ]]; then
    copy_back+="
for _d in \"\$SCRATCH_JOB\"/*-fastqc; do [ -d \"\$_d\" ] && { cp -r \"\$_d\" '${SAMPLE_DIR}/' || { echo '<ERROR> cp fastqc failed' >&2; _cb_fail=1; }; }; done"
  fi
  if [[ "$input_mode" == "bqsr_bam" || "$input_mode" == "bqsr_bam_scatter" ]]; then
    copy_back+="
if [ -d \"\$SCRATCH_JOB/chrom_gvcf\" ]; then
  mkdir -p '${SAMPLE_DIR}/chrom_gvcf'
  for _f in \"\$SCRATCH_JOB\"/chrom_gvcf/*; do [ -e \"\$_f\" ] && _cb \"\$_f\" '${SAMPLE_DIR}/chrom_gvcf'; done
fi"
  fi

  printf '%s' "
set -euo pipefail

# --- scratch root: primary group of the RUNNING user (node-independent) ---
_grp=\"\$(id -gn)\"
SCRATCH_ROOT=\"${scratch_base}/\${_grp}/\${USER}\"
SCRATCH_JOB=\"\${SCRATCH_ROOT}/job_\${SLURM_JOB_ID}\"
if ! ( mkdir -p \"\$SCRATCH_JOB/tmp\" && : > \"\$SCRATCH_JOB/tmp/.w\" && rm -f \"\$SCRATCH_JOB/tmp/.w\" ) 2>/dev/null; then
  echo \"<ERROR> scratch unavailable on \$(hostname -s): \$SCRATCH_ROOT (primary group \$_grp)\" >&2
  echo \"<ERROR> add \$(hostname -s) to sbatch_exclude in config/config.yaml, then resubmit.\" >&2
  exit 1
fi
trap 'rm -rf \"\$SCRATCH_JOB\"' EXIT HUP TERM
export TMPDIR=\"\$SCRATCH_JOB/tmp\"
# Keep crash artefacts on scratch (auto-wiped), not in the repo / CWD: a JVM
# SIGSEGV on a bad node otherwise drops core.<pid> + hs_err_pid<pid>.log next to
# the launcher (2026-09: node15 crashed 7 L066 shards, littered the repo root).
ulimit -c 0 2>/dev/null || true
cd \"\$SCRATCH_JOB\"
module load samtools >/dev/null 2>&1 || true
echo \"[i] Scratch: \$SCRATCH_JOB  (group \$_grp, node \$(hostname -s))\"
${stage_inputs}

bash '${script}' \"\${SCRIPT_ARGS[@]}\"

# --- verified copy-back --------------------------------------------------------
_cb_fail=0
_cb() {   # \$1 = source file (under scratch)   \$2 = destination directory
  local s=\"\$1\" d=\"\$2/\$(basename \"\$1\")\"
  cp -f \"\$s\" \"\$d\" || { echo \"<ERROR> cp failed: \$s\" >&2; _cb_fail=1; return; }
  [ \"\$(stat -c%s \"\$s\" 2>/dev/null)\" = \"\$(stat -c%s \"\$d\" 2>/dev/null)\" ] || {
    echo \"<ERROR> size mismatch after copy: \$d\" >&2; _cb_fail=1; return; }
  case \"\$s\" in
    *.bam)                       command -v samtools >/dev/null && { samtools quickcheck \"\$d\" 2>/dev/null || { echo \"<ERROR> quickcheck failed: \$d\" >&2; _cb_fail=1; }; } ;;
    *.vcf.gz|*.g.vcf.gz|*.tbi)   gzip -t \"\$d\" 2>/dev/null || { echo \"<ERROR> gzip -t failed: \$d\" >&2; _cb_fail=1; } ;;
  esac
}
${copy_back}

if [ \"\$_cb_fail\" -ne 0 ]; then
  echo \"<ERROR> copy-back verification FAILED — scratch preserved for recovery: \$SCRATCH_JOB\" >&2
  trap - EXIT
  exit 1
fi
echo \"[i] Copy-back verified; scratch will be removed.\"
"
}

# --- Build wraps -------------------------------------------------------------

if [[ "$USE_SCRATCH" == true ]]; then

  WRAP_S01="$(build_scratch_wrap "$S01" "fastq" \
    "*.sort.bam *.sort.bam.bai *.sort.stats.txt")"

  WRAP_S02="$(build_scratch_wrap "$S02" "sort_bams" \
    "*.rmdup.mqfilt.bqsr.bam *.rmdup.mqfilt.bqsr.bam.bai \
     *.bqsr_table.txt *-dups.txt *.mosdepth.* *.metrics.txt *.pdf")"

  # S03 whole-genome: keep only the canon_chr GVCF + per-chrom split
  # (chrom_gvcf/ is the Step 04 input). all-contigs raw_variants.g.vcf.gz is NOT
  # copied back. S03 scatter: each array task writes one chrom_gvcf/ entry.
  WRAP_S03="$(build_scratch_wrap "$S03" "bqsr_bam" \
    "*.raw_variants.canon_chr.g.vcf.gz *.raw_variants.canon_chr.g.vcf.gz.tbi")"
  WRAP_S03_SCATTER="$(build_scratch_wrap "$S03" "bqsr_bam_scatter" "")"

else

  WRAP_S01="set -euo pipefail; cd '${SAMPLE_DIR}'; bash '${S01}' '${SAMPLE_DIR}' '${SAMPLE_DIR}'"

  WRAP_S02="
    set -euo pipefail
    bam_list=\$(find '${SAMPLE_DIR}' -maxdepth 1 \
                 \\( -name '*.sort.bam' -o -name '*.sorted.bam' \\) | sort | paste -sd,)
    [[ -z \"\$bam_list\" ]] && { echo '<ERROR> No *.sort.bam in ${SAMPLE_DIR}'; exit 1; }
    bash '${S02}' \"\$bam_list\" '${SAMPLE_DIR}'"

  WRAP_S03="
    set -euo pipefail
    bqsr_bam=\$(find '${SAMPLE_DIR}' -maxdepth 1 -name '*.rmdup.mqfilt.bqsr.bam' | head -1)
    [[ -z \"\$bqsr_bam\" ]] && { echo '<ERROR> No *.rmdup.mqfilt.bqsr.bam in ${SAMPLE_DIR}'; exit 1; }
    bash '${S03}' \"\$bqsr_bam\" '${SAMPLE_DIR}'"

  WRAP_S03_SCATTER="
    set -euo pipefail
    bqsr_bam=\$(find '${SAMPLE_DIR}' -maxdepth 1 -name '*.rmdup.mqfilt.bqsr.bam' | head -1)
    [[ -z \"\$bqsr_bam\" ]] && { echo '<ERROR> No *.rmdup.mqfilt.bqsr.bam in ${SAMPLE_DIR}'; exit 1; }
    _canon=(${CANON_CHROMS})
    _chrom=\"\${_canon[\$((SLURM_ARRAY_TASK_ID - 1))]:-}\"
    [[ -z \"\$_chrom\" ]] && { echo \"<ERROR> no chromosome for array index \$SLURM_ARRAY_TASK_ID\" >&2; exit 1; }
    bash '${S03}' \"\$bqsr_bam\" '${SAMPLE_DIR}' \"\$_chrom\""

fi

# S03 gather (scatter mode only): concat the 25 per-chrom GVCFs → canon_chr GVCF.
# The pipeline does NOT consume this file — Step 04 (GenomicsDBImport) reads
# chrom_gvcf/ directly. It is produced purely for ARCHIVAL: one whole-genome
# GVCF per sample on the backup volume is far easier to handle than 25 shards.
# ~13 min with a full bcftools concat (recompress) — deliberately not --naive;
# the cost is trivial next to the ~3 h scatter and a clean single-stream file is
# worth more for long-term storage. Writes straight to SAMPLE_DIR via .tmp rename.
WRAP_S03_GATHER="
set -euo pipefail
command -v bcftools >/dev/null || module load bcftools >/dev/null 2>&1 || true
cd '${SAMPLE_DIR}'
_canon=(${CANON_CHROMS})
_files=()
for _c in \"\${_canon[@]}\"; do
  _f=\$(find chrom_gvcf -maxdepth 1 -name \"*.raw_vars.\${_c}.g.vcf.gz\" 2>/dev/null | head -1)
  [[ -z \"\$_f\" ]] && { echo \"<ERROR> gather: missing chrom_gvcf for \$_c\" >&2; exit 1; }
  _files+=(\"\$_f\")
done
_pref=\$(basename \"\${_files[0]}\"); _pref=\${_pref%%.raw_vars.*}
_out=\"\${_pref}.raw_variants.canon_chr.g.vcf.gz\"
echo \"[i] gather \${#_files[@]} chromosomes -> \$_out\"
bcftools concat --output-type z --output \"\${_out}.tmp\" \"\${_files[@]}\"
bcftools index --tbi --force \"\${_out}.tmp\"
mv \"\${_out}.tmp\" \"\$_out\"
mv \"\${_out}.tmp.tbi\" \"\${_out}.tbi\"
echo \"[i] gather done: \$_out\"
"

# Resumability: under scratch each step starts in a fresh empty dir, so the
# steps' own skip logic never sees prior outputs in SAMPLE_DIR. We decide HERE
# whether to submit each step, from the presence of its final output:
#
#   S01 → *.sort.bam
#   S02 → *.rmdup.mqfilt.bqsr.bam
#   S03 → *.raw_variants.canon_chr.g.vcf.gz
#
# Skipped steps drop out of the afterok chain; the next submitted step depends
# on the most recent job actually submitted (if any).

# One chrom_gvcf/ entry present and non-empty (with its .tbi)?
have_chrom_gvcf() {
  local f
  f=$(find "$SAMPLE_DIR/chrom_gvcf" -maxdepth 1 -name "*.raw_vars.$1.g.vcf.gz" 2>/dev/null | head -1)
  [ -n "$f" ] && [ -s "$f" ] && [ -s "${f}.tbi" ]
}
# Comma list of 1-based array indices whose chromosome GVCF is missing ("" = none)
chrom_missing_indices() {
  local i=0 c out=""
  for c in $CANON_CHROMS; do
    i=$((i + 1))
    have_chrom_gvcf "$c" || out="${out:+$out,}$i"
  done
  printf '%s' "$out"
}

DEP=""          # afterok dependency for the next job to submit
SUBMITTED=()    # job ids actually submitted (for the monitor hint)

# --- Step 01: BWA Alignment -------------------------------------------------
# 8 CPUs / 20G / 48h — bwa mem -t8; repair.sh pinned -Xmx8g; sort TMPDIR→scratch.
# Walltimes are 2x the observed worst case, rounded to a 24h grid: high cluster
# load stretches these jobs and defq has no time cap, so headroom beats reruns.

if have_sort_bam; then
  echo "[SKIP] Step 01 — sorted BAM(s) already present in ${SAMPLE_DIR}"
else
  JOB01=$(sbatch \
    --job-name="${SAMPLE_ID}-S01-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=8 \
    --mem=20G --time=48:00:00 \
    "${EXCLUDE_ARG[@]}" \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S01" \
    | awk '{print $4}')
  echo "[>] Step 01 submitted  — Job ${JOB01}  (8 CPUs / 20G / 48h)"
  DEP="afterok:${JOB01}"
  SUBMITTED+=("$JOB01")
fi

# --- Step 02: BAM QC + BQSR -----------------------------------------------
# 4 CPUs / 28G / 48h — MarkDuplicates -Xmx16g, BQSR -Xmx8g; BaseRecalibrator is
# the single-threaded tail (~7h at 8GB, ~12h at 9.5GB observed — 48h wall).

if have_output "*.rmdup.mqfilt.bqsr.bam"; then
  echo "[SKIP] Step 02 — analysis-ready BAM already present in ${SAMPLE_DIR}"
  # Don't let a lingering S01 dependency (e.g. S01 was redundantly resubmitted
  # for a legacy dir, or is simply unrelated) leak forward onto S03 — S02's
  # output already exists independent of whatever S01 is/was doing.
  DEP=""
else
  JOB02=$(sbatch \
    --job-name="${SAMPLE_ID}-S02-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=4 \
    --mem=28G --time=48:00:00 \
    "${EXCLUDE_ARG[@]}" \
    ${DEP:+--dependency="$DEP"} \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S02" \
    | awk '{print $4}')
  echo "[>] Step 02 submitted  — Job ${JOB02}  (4 CPUs / 28G / 48h)${DEP:+  [${DEP}]}"
  DEP="afterok:${JOB02}"
  SUBMITTED+=("$JOB02")
fi

# --- Step 03: HaplotypeCaller --------------------------------------------
# SCATTER_S03=true  : 25-way array (one chrom each), 2 CPU / 13G / 48h per task,
#                     %${ARRAY_CONC} concurrent; then a gather job for canon_chr.
# SCATTER_S03=false : one whole-genome job, 4 CPU / 32G / 48h.
# Either way "done" = all 25 chrom_gvcf/ entries present; canon_chr GVCF is the
# archival roll-up.

canon_present=false
have_output "*.raw_variants.canon_chr.g.vcf.gz" && canon_present=true

if [[ "$SCATTER_S03" == true ]]; then

  miss="$(chrom_missing_indices)"

  if [[ -z "$miss" && "$canon_present" == true ]]; then
    echo "[SKIP] Step 03 — 25/25 chrom GVCFs + canon roll-up already present"
  else
    ARR_DEP=""
    if [[ -n "$miss" ]]; then
      nmiss=$(awk -F, '{print NF}' <<< "$miss")
      [[ "$nmiss" -eq "$CANON_N" ]] && miss="1-${CANON_N}"   # tidy the all-missing case
      JOB03=$(sbatch \
        --job-name="${SAMPLE_ID}-S03-${EPOCHSECONDS}" \
        --array="${miss}%${ARRAY_CONC}" \
        --nodes=1 --ntasks=1 --cpus-per-task=2 \
        --mem=13G --time=48:00:00 \
        "${EXCLUDE_ARG[@]}" \
        ${DEP:+--dependency="$DEP"} \
        --output="${LOG_DIR}/%x.%A_%a.log" \
        --wrap "$WRAP_S03_SCATTER" \
        | awk '{print $4}')
      echo "[>] Step 03 array submitted — Job ${JOB03}  (${nmiss} chrom(s): ${miss}; 2 CPU / 13G / 48h; %${ARRAY_CONC})${DEP:+  [${DEP}]}"
      ARR_DEP="afterok:${JOB03}"
      SUBMITTED+=("$JOB03")
    fi
    gdep="${ARR_DEP:-$DEP}"
    JOB03G=$(sbatch \
      --job-name="${SAMPLE_ID}-S03g-${EPOCHSECONDS}" \
      --nodes=1 --ntasks=1 --cpus-per-task=2 \
      --mem=8G --time=24:00:00 \
      "${EXCLUDE_ARG[@]}" \
      ${gdep:+--dependency="$gdep"} \
      --output="${LOG_DIR}/%x.%j.log" \
      --wrap "$WRAP_S03_GATHER" \
      | awk '{print $4}')
    echo "[>] Step 03 gather submitted — Job ${JOB03G}${gdep:+  [${gdep}]}"
    SUBMITTED+=("$JOB03G")
  fi

else

  if [[ -z "$(chrom_missing_indices)" && "$canon_present" == true ]]; then
    echo "[SKIP] Step 03 — final GVCF + chrom_gvcf/ already present in ${SAMPLE_DIR}"
  else
    JOB03=$(sbatch \
      --job-name="${SAMPLE_ID}-S03-${EPOCHSECONDS}" \
      --nodes=1 --ntasks=1 --cpus-per-task=4 \
      --mem=32G --time=48:00:00 \
      "${EXCLUDE_ARG[@]}" \
      ${DEP:+--dependency="$DEP"} \
      --output="${LOG_DIR}/%x.%j.log" \
      --wrap "$WRAP_S03" \
      | awk '{print $4}')
    echo "[>] Step 03 submitted  — Job ${JOB03}  (4 CPUs / 32G / 48h)${DEP:+  [${DEP}]}"
    SUBMITTED+=("$JOB03")
  fi

fi

# <\MAIN> ---------------------------------------------------------------------

echo
if [[ ${#SUBMITTED[@]} -eq 0 ]]; then
  echo "[i] Nothing to submit — all step outputs already present in ${SAMPLE_DIR}."
else
  echo "[i] Monitor jobs:"
  echo "    squeue -u \$USER | grep ${SAMPLE_ID}"
  echo "    squeue -j $(IFS=,; echo "${SUBMITTED[*]}")"
  echo "[i] Logs: ${LOG_DIR}/"
fi
echo
