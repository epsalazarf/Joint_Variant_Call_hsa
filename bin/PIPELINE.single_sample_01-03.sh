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
# <\TOGGLES> ------------------------------------------------------------------

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

# <CHECKS> --------------------------------------------------------------------

for f in "$S01" "$S02" "$S03"; do
  [[ -f "$f" ]] || { echo "<ERROR> Script not found: $f"; exit 1; }
done

[[ -d "$SAMPLE_DIR" ]] || { echo "<ERROR> Sample directory not found: $SAMPLE_DIR"; exit 1; }

if ! find "$SAMPLE_DIR" -maxdepth 1 -name "*.f*q.gz" | grep -qE '_R?1\.f[^.]*q\.gz$'; then
  echo "<ERROR> No FASTQ R1 files found in: $SAMPLE_DIR"
  exit 1
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
  if [[ "$input_mode" == "bqsr_bam" ]]; then
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

  # S03: keep only the canonical-chromosome GVCF + the per-chrom split
  # (chrom_gvcf/ is the Step 04 input). The all-contigs raw_variants.g.vcf.gz is
  # redundant and NOT copied back.
  WRAP_S03="$(build_scratch_wrap "$S03" "bqsr_bam" \
    "*.raw_variants.canon_chr.g.vcf.gz *.raw_variants.canon_chr.g.vcf.gz.tbi")"

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

fi

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

have_output() { find "$SAMPLE_DIR" -maxdepth 1 -name "$1" | grep -q .; }

DEP=""          # afterok dependency for the next job to submit
SUBMITTED=()    # job ids actually submitted (for the monitor hint)

# --- Step 01: BWA Alignment -------------------------------------------------
# 8 CPUs / 20G / 14h — bwa mem -t8; repair.sh pinned -Xmx8g; sort TMPDIR→scratch

if have_output "*.sort.bam"; then
  echo "[SKIP] Step 01 — sorted BAM(s) already present in ${SAMPLE_DIR}"
else
  JOB01=$(sbatch \
    --job-name="${SAMPLE_ID}-S01-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=8 \
    --mem=20G --time=14:00:00 \
    "${EXCLUDE_ARG[@]}" \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S01" \
    | awk '{print $4}')
  echo "[>] Step 01 submitted  — Job ${JOB01}  (8 CPUs / 20G / 14h)"
  DEP="afterok:${JOB01}"
  SUBMITTED+=("$JOB01")
fi

# --- Step 02: BAM QC + BQSR -----------------------------------------------
# 4 CPUs / 28G / 20h — MarkDuplicates -Xmx16g, BQSR -Xmx8g; BaseRecalibrator is
# the single-threaded tail (~7h at 8GB, ~10h at 9.5GB — hence 20h wall).

if have_output "*.rmdup.mqfilt.bqsr.bam"; then
  echo "[SKIP] Step 02 — analysis-ready BAM already present in ${SAMPLE_DIR}"
else
  JOB02=$(sbatch \
    --job-name="${SAMPLE_ID}-S02-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=4 \
    --mem=28G --time=20:00:00 \
    "${EXCLUDE_ARG[@]}" \
    ${DEP:+--dependency="$DEP"} \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S02" \
    | awk '{print $4}')
  echo "[>] Step 02 submitted  — Job ${JOB02}  (4 CPUs / 28G / 20h)${DEP:+  [${DEP}]}"
  DEP="afterok:${JOB02}"
  SUBMITTED+=("$JOB02")
fi

# --- Step 03: HaplotypeCaller --------------------------------------------
# 4 CPUs / 32G / 24h — GVCF-mode HaplotypeCaller is ~10h at 4.8GB and 14-18h at
# 8-11GB (single sample, whole genome). A per-chromosome scatter is the real fix
# for throughput — see docs; until then the wall has to cover the largest sample.

if have_output "*.raw_variants.canon_chr.g.vcf.gz" \
   && [ -d "${SAMPLE_DIR}/chrom_gvcf" ] \
   && [ -n "$(ls -A "${SAMPLE_DIR}/chrom_gvcf" 2>/dev/null)" ]; then
  echo "[SKIP] Step 03 — final GVCF + chrom_gvcf/ already present in ${SAMPLE_DIR}"
else
  JOB03=$(sbatch \
    --job-name="${SAMPLE_ID}-S03-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=4 \
    --mem=32G --time=24:00:00 \
    "${EXCLUDE_ARG[@]}" \
    ${DEP:+--dependency="$DEP"} \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S03" \
    | awk '{print $4}')
  echo "[>] Step 03 submitted  — Job ${JOB03}  (4 CPUs / 32G / 24h)${DEP:+  [${DEP}]}"
  SUBMITTED+=("$JOB03")
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
