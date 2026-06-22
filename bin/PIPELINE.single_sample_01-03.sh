#!/bin/bash
# =============================================================================
# Title       : Single-Sample Sequential Pipeline Launcher [FENIX]
# Description : Submits Steps 01→02→03 as chained SLURM jobs for one sample.
#               Each step is held until the previous one succeeds (afterok).
#               When USE_SCRATCH=true, inputs are read directly from NFS while
#               TMPDIR + intermediate/output files use /scratch, then finals are
#               copied back. This keeps bulky transient intermediates off the
#               group's persistent quota. NOTE: benchmarking showed scratch
#               gives no wall-time speedup for these sequential GATK steps (they
#               are CPU-bound) — its value here is storage/quota and cluster
#               citizenship, not speed. See INSTRUCTIONS.md.
# Usage       : PIPELINE.single_sample.sh [sample_dir]
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

if [[ "$USE_SCRATCH" == true ]]; then
  if [[ -z "${scratch_root:-}" ]]; then
    echo "<ERROR> USE_SCRATCH=true but scratch_root is not set in config (${env_type} section)."
    exit 1
  fi
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
echo "[i]  Scratch    : ${USE_SCRATCH} $([ "$USE_SCRATCH" = true ] && echo "(${scratch_root}/\$USER/job_\$SLURM_JOB_ID)" || true)"
echo

# Helper: build a scratch-enabled wrap for a pipeline step.
#
# Strategy (see Sept-2026 benchmark in INSTRUCTIONS.md): inputs are read
# DIRECTLY from NFS (SAMPLE_DIR) — they are large and read once sequentially,
# so staging them to scratch is pure overhead with no speedup. Scratch is used
# only for TMPDIR (tool temp/spill — the genuinely IOPS-heavy traffic) and for
# intermediate + output files, which are then copied back to SAMPLE_DIR. This
# keeps bulky transient intermediates off the group's persistent quota without
# paying a useless copy-in tax.
#
#   $1  — pipeline script path (expands at submission time)
#   $2  — "inputs" descriptor, one of:
#           "fastq"      : run S01 from SAMPLE_DIR (FASTQs read in place)
#           "sort_bams"  : read *.sort.bam from SAMPLE_DIR (NFS), output to scratch
#           "bqsr_bam"   : read *.rmdup.mqfilt.bqsr.bam from SAMPLE_DIR (NFS)
#   $3  — glob(s) of output files/dirs to copy back (space-separated, relative to scratch)
#
# Variables that expand at submission time (from this script): $1, SAMPLE_DIR, scratch_root
# Variables that expand at job runtime (escaped with \$):       SLURM_JOB_ID, USER, SCRATCH_JOB
#
build_scratch_wrap() {
  local script="$1"
  local input_mode="$2"
  local copy_back_globs="$3"

  # -- resolve inputs (read from NFS) and point output at scratch --
  local stage_inputs=""
  case "$input_mode" in
    fastq)
      # Script 01 uses `find .` for FASTQ discovery, so it must run from SAMPLE_DIR.
      # FASTQs read in place; output (BAMs) goes to scratch, copied back below.
      stage_inputs="cd '${SAMPLE_DIR}'
SCRIPT_ARGS=('${SAMPLE_DIR}' \"\$SCRATCH_JOB\")"
      ;;
    sort_bams)
      stage_inputs="
# Read *.sort.bam directly from NFS; output to scratch
_nfs_list=\$(find '${SAMPLE_DIR}' -maxdepth 1 \( -name '*.sort.bam' -o -name '*.sorted.bam' \) | sort | paste -sd,)
[[ -z \"\$_nfs_list\" ]] && { echo '<ERROR> No *.sort.bam files in ${SAMPLE_DIR}'; exit 1; }
SCRIPT_ARGS=(\"\$_nfs_list\" \"\$SCRATCH_JOB\")"
      ;;
    bqsr_bam)
      stage_inputs="
# Read *.bqsr.bam directly from NFS; output to scratch
_bqsr=\$(find '${SAMPLE_DIR}' -maxdepth 1 -name '*.rmdup.mqfilt.bqsr.bam' | head -1)
[[ -z \"\$_bqsr\" ]] && { echo '<ERROR> No *.rmdup.mqfilt.bqsr.bam in ${SAMPLE_DIR}'; exit 1; }
SCRIPT_ARGS=(\"\$_bqsr\" \"\$SCRATCH_JOB\")"
      ;;
  esac

  # -- copy-back: each glob is applied under SCRATCH_JOB --
  local copy_back=""
  for glob in $copy_back_globs; do
    copy_back+="
find \"\$SCRATCH_JOB\" -maxdepth 1 -name '${glob}' -exec cp -r {} '${SAMPLE_DIR}/' \\; 2>/dev/null || true"
  done
  # Always copy the chrom_gvcf directory if present (step 03)
  if [[ "$input_mode" == "bqsr_bam" ]]; then
    copy_back+="
[[ -d \"\$SCRATCH_JOB/chrom_gvcf\" ]] && cp -r \"\$SCRATCH_JOB/chrom_gvcf\" '${SAMPLE_DIR}/' 2>/dev/null || true"
  fi

  printf '%s' "
set -euo pipefail
SCRATCH_JOB=\"${scratch_root}/\${USER}/job_\${SLURM_JOB_ID}\"
mkdir -p \"\$SCRATCH_JOB/tmp\"
export TMPDIR=\"\$SCRATCH_JOB/tmp\"
echo \"[i] Scratch: \$SCRATCH_JOB\"
${stage_inputs}
bash '${script}' \"\${SCRIPT_ARGS[@]}\"
${copy_back}
echo \"[i] Copying back done. Cleaning scratch...\"
rm -rf \"\$SCRATCH_JOB\"
echo \"[i] Scratch removed.\"
"
}

# --- Build wraps -------------------------------------------------------------

if [[ "$USE_SCRATCH" == true ]]; then

  WRAP_S01="$(build_scratch_wrap "$S01" "fastq" \
    "*.sort.bam *.sort.bam.bai *.sort.stats.txt")"

  WRAP_S02="$(build_scratch_wrap "$S02" "sort_bams" \
    "*.rmdup.mqfilt.bqsr.bam *.rmdup.mqfilt.bqsr.bam.bai \
     *.bqsr_table.txt *-dups.txt *.mosdepth.* *.metrics.txt *.pdf")"

  WRAP_S03="$(build_scratch_wrap "$S03" "bqsr_bam" \
    "*.g.vcf.gz *.g.vcf.gz.tbi")"

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

# Resumability: each step's internal skip logic checks its OUTPUT_PATH, but with
# scratch that is a fresh empty dir every run — so prior outputs in SAMPLE_DIR are
# never seen and steps redo all work. We therefore decide HERE whether to submit
# each step, based on whether its final output already exists in SAMPLE_DIR.
#
#   S01 → *.sort.bam                        (per-lane sorted BAMs)
#   S02 → *.rmdup.mqfilt.bqsr.bam           (analysis-ready BAM)
#   S03 → *.raw_variants.canon_chr.g.vcf.gz (final GVCF)
#
# Skipped steps drop out of the afterok chain; the next submitted step depends on
# the most recent job actually submitted (if any).

have_output() { find "$SAMPLE_DIR" -maxdepth 1 -name "$1" | grep -q .; }

DEP=""          # afterok dependency for the next job to submit
SUBMITTED=()    # job ids actually submitted (for the monitor hint)

# --- Step 01: BWA Alignment --------------------------------------------------
# 8 CPUs / 16G — matches njobs=8 FENIX default; samtools sort TMPDIR → scratch

if have_output "*.sort.bam"; then
  echo "[SKIP] Step 01 — sorted BAM(s) already present in ${SAMPLE_DIR}"
else
  JOB01=$(sbatch \
    --job-name="${SAMPLE_ID}-S01-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=8 \
    --mem=16G \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S01" \
    | awk '{print $4}')
  echo "[>] Step 01 submitted  — Job ${JOB01}  (8 CPUs / 16G)"
  DEP="afterok:${JOB01}"
  SUBMITTED+=("$JOB01")
fi

# --- Step 02: BAM QC + BQSR --------------------------------------------------
# 8 CPUs / 32G — MarkDuplicatesSpark -Xmx24G + GC threads; biggest scratch win

if have_output "*.rmdup.mqfilt.bqsr.bam"; then
  echo "[SKIP] Step 02 — analysis-ready BAM already present in ${SAMPLE_DIR}"
else
  JOB02=$(sbatch \
    --job-name="${SAMPLE_ID}-S02-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=8 \
    --mem=32G \
    ${DEP:+--dependency="$DEP"} \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S02" \
    | awk '{print $4}')
  echo "[>] Step 02 submitted  — Job ${JOB02}  (8 CPUs / 32G)${DEP:+  [${DEP}]}"
  DEP="afterok:${JOB02}"
  SUBMITTED+=("$JOB02")
fi

# --- Step 03: HaplotypeCaller ------------------------------------------------
# 4 CPUs / 32G — --native-pair-hmm-threads 4, -Xmx20G on FENIX

if have_output "*.raw_variants.canon_chr.g.vcf.gz"; then
  echo "[SKIP] Step 03 — final GVCF already present in ${SAMPLE_DIR}"
else
  JOB03=$(sbatch \
    --job-name="${SAMPLE_ID}-S03-${EPOCHSECONDS}" \
    --nodes=1 --ntasks=1 --cpus-per-task=4 \
    --mem=32G \
    ${DEP:+--dependency="$DEP"} \
    --output="${LOG_DIR}/%x.%j.log" \
    --wrap "$WRAP_S03" \
    | awk '{print $4}')
  echo "[>] Step 03 submitted  — Job ${JOB03}  (4 CPUs / 32G)${DEP:+  [${DEP}]}"
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
