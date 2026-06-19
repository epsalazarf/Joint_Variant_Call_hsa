#!/bin/bash
# =============================================================================
# Title       : S02 Scratch vs NFS Benchmark [FENIX]
# Description : Submits N parallel SLURM jobs of Step 02 (BAM QC) for the same
#               input, half using /scratch and half running directly on NFS,
#               to measure the real-world benefit of scratch staging.
#               The scratch arm times copy-in / run / copy-out separately so
#               staging overhead is counted against it (fair comparison).
# Usage       : S02_scratch_benchmark.sh [sample_dir] [n_iters]
#                 sample_dir — dir with *.sort.bam inputs (default: gamma/L66)
#                 n_iters    — reps per arm (default: 3)
# Notes       : Both arms run the identical 02_gatk_bam_qc_workflow.sh, so any
#               script-internal toggles (RUN_METRICS, HOUSEKEEP, ...) cancel out.
#               NFS arm leaves TMPDIR at the node default (/tmp), matching the
#               original "normal" runs; scratch arm points TMPDIR at scratch,
#               matching the pipeline launcher.
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

SAMPLE_DIR="$(realpath "${1:-/mnt/data/amedina/esalazarf/JVCdev/test_gamma/L66}")"
N_ITERS="${2:-3}"

# After timing is recorded, delete the large BAM outputs to save NFS space
# (set to true to keep full outputs for inspection).
KEEP_OUTPUTS=false

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
S02="${SCRIPT_DIR}/../02_gatk_bam_qc_workflow.sh"
CONFIG_FILE="${SCRIPT_DIR}/../../config/config.yaml"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

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

[[ -f "$S02" ]] || { echo "<ERROR> S02 script not found: $S02"; exit 1; }
[[ -d "$SAMPLE_DIR" ]] || { echo "<ERROR> Sample dir not found: $SAMPLE_DIR"; exit 1; }
command -v sbatch &>/dev/null || { echo "<ERROR> sbatch not found — SLURM required."; exit 1; }
[[ -n "${scratch_root:-}" ]] || { echo "<ERROR> scratch_root not set in config (${env_type})."; exit 1; }

# Discover input sort BAMs (absolute paths, comma-joined)
mapfile -t BAMS < <(find "$SAMPLE_DIR" -maxdepth 1 -name '*.sort.bam' | sort)
[[ ${#BAMS[@]} -gt 0 ]] || { echo "<ERROR> No *.sort.bam in $SAMPLE_DIR"; exit 1; }
BAM_CSV_NFS=$(IFS=','; echo "${BAMS[*]}")

# <\CHECKS> -------------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

SAMPLE_ID="$(basename "$SAMPLE_DIR")"
STAMP="$EPOCHSECONDS"
BENCH_ROOT="$(dirname "$SAMPLE_DIR")/S02_bench_${SAMPLE_ID}_${STAMP}"
RESULTS="${BENCH_ROOT}/results"
LOGS="${BENCH_ROOT}/log"
mkdir -p "$RESULTS" "$LOGS"

echo
echo "[&] S02 Scratch vs NFS Benchmark [FENIX]"
echo "[i]  Sample      : $SAMPLE_ID  (${#BAMS[@]} BAM input(s))"
echo "[i]  Iterations  : $N_ITERS per arm  ($((N_ITERS*2)) jobs total)"
echo "[i]  Bench root  : $BENCH_ROOT"
echo "[i]  Keep outputs: $KEEP_OUTPUTS"
echo

# Build the line that trims large outputs once timing is captured
if [[ "$KEEP_OUTPUTS" == true ]]; then
  TRIM=":"   # no-op
else
  TRIM="rm -f \"\$OUT\"/*.bam \"\$OUT\"/*.bai 2>/dev/null || true"
fi

SUBMITTED=()

for ((i=1; i<=N_ITERS; i++)); do

  # -------- NFS arm (no scratch) --------
  RES_NFS="${RESULTS}/noscratch_iter${i}.\${SLURM_JOB_ID}.txt"
  WRAP_NFS="
    set -euo pipefail
    OUT='${BENCH_ROOT}/noscratch_iter${i}'
    mkdir -p \"\$OUT\"
    echo \"[bench] NFS arm iter ${i} — reads+writes on NFS, TMPDIR=node default\"
    t0=\$(date +%s)
    bash '${S02}' '${BAM_CSV_NFS}' \"\$OUT\"
    t1=\$(date +%s)
    run_s=\$((t1-t0))
    printf 'arm=noscratch iter=%d jobid=%s copy_in_s=0 run_s=%d copy_out_s=0 total_s=%d\n' \
      ${i} \"\$SLURM_JOB_ID\" \"\$run_s\" \"\$run_s\" > '${RES_NFS}'
    ${TRIM}
    echo \"[bench] NFS arm iter ${i} done: run=\${run_s}s\"
  "

  JID=$(sbatch --parsable \
    --job-name="BENCH-S02nfs-${SAMPLE_ID}-i${i}" \
    --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G \
    --output="${LOGS}/noscratch_iter${i}.%j.log" \
    --wrap "$WRAP_NFS")
  echo "[>] NFS  arm iter ${i} submitted — Job ${JID}"
  SUBMITTED+=("$JID")

  # -------- Scratch arm --------
  RES_SCR="${RESULTS}/scratch_iter${i}.\${SLURM_JOB_ID}.txt"
  WRAP_SCR="
    set -euo pipefail
    OUT='${BENCH_ROOT}/scratch_iter${i}'
    mkdir -p \"\$OUT\"
    SJ=\"${scratch_root}/\${USER}/bench_\${SLURM_JOB_ID}\"
    mkdir -p \"\$SJ/tmp\"
    export TMPDIR=\"\$SJ/tmp\"
    echo \"[bench] Scratch arm iter ${i} — SJ=\$SJ\"

    # copy-in (timed)
    t0=\$(date +%s)
    for b in ${BAMS[*]}; do
      cp \"\$b\" \"\$SJ/\"
      [[ -f \"\${b}.bai\" ]] && cp \"\${b}.bai\" \"\$SJ/\"
    done
    t1=\$(date +%s); copy_in_s=\$((t1-t0))

    # build scratch input list
    scratch_csv=\$(find \"\$SJ\" -maxdepth 1 -name '*.sort.bam' | sort | paste -sd,)

    # run (timed)
    bash '${S02}' \"\$scratch_csv\" \"\$SJ\"
    t2=\$(date +%s); run_s=\$((t2-t1))

    # copy-out (timed): final BAM + index + small metrics/tables
    cp \"\$SJ\"/*.rmdup.mqfilt.bqsr.bam  \"\$OUT/\" 2>/dev/null || true
    cp \"\$SJ\"/*.rmdup.mqfilt.bqsr.bam.bai \"\$OUT/\" 2>/dev/null || true
    cp \"\$SJ\"/*.txt \"\$SJ\"/*.pdf \"\$OUT/\" 2>/dev/null || true
    t3=\$(date +%s); copy_out_s=\$((t3-t2))

    total_s=\$((t3-t0))
    printf 'arm=scratch iter=%d jobid=%s copy_in_s=%d run_s=%d copy_out_s=%d total_s=%d\n' \
      ${i} \"\$SLURM_JOB_ID\" \"\$copy_in_s\" \"\$run_s\" \"\$copy_out_s\" \"\$total_s\" > '${RES_SCR}'

    rm -rf \"\$SJ\"
    ${TRIM}
    echo \"[bench] Scratch arm iter ${i} done: copy_in=\${copy_in_s}s run=\${run_s}s copy_out=\${copy_out_s}s total=\${total_s}s\"
  "

  JID=$(sbatch --parsable \
    --job-name="BENCH-S02scr-${SAMPLE_ID}-i${i}" \
    --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G \
    --output="${LOGS}/scratch_iter${i}.%j.log" \
    --wrap "$WRAP_SCR")
  echo "[>] Scr  arm iter ${i} submitted — Job ${JID}"
  SUBMITTED+=("$JID")

done

# <\MAIN> ---------------------------------------------------------------------

# Write a self-contained summarizer next to the results
cat > "${BENCH_ROOT}/summarize.sh" <<'SUMM'
#!/bin/bash
# Tabulate benchmark results. Usage: ./summarize.sh   (run from bench root)
set -euo pipefail
cd "$(dirname "$(realpath "$0")")/results"
echo "Per-job results:"
cat *.txt 2>/dev/null | sort || { echo "(no results yet)"; exit 0; }
echo
echo "Averages by arm (seconds; H:MM:SS):"
awk '
  { for (j=1;j<=NF;j++){ split($j,a,"="); v[a[1]]=a[2] }
    arm=v["arm"]; n[arm]++;
    ci[arm]+=v["copy_in_s"]; r[arm]+=v["run_s"]; co[arm]+=v["copy_out_s"]; t[arm]+=v["total_s"] }
  function hms(s,  h,m){ h=int(s/3600); m=int((s%3600)/60); return sprintf("%d:%02d:%02d",h,m,s%60) }
  END{
    printf "%-10s %4s  %12s %12s %12s %12s\n","arm","n","copy_in","run","copy_out","total"
    for (k in n) printf "%-10s %4d  %12s %12s %12s %12s\n", k, n[k], \
       hms(ci[k]/n[k]), hms(r[k]/n[k]), hms(co[k]/n[k]), hms(t[k]/n[k])
  }' *.txt
SUMM
chmod +x "${BENCH_ROOT}/summarize.sh"

echo
echo "[i] Submitted $((N_ITERS*2)) jobs."
echo "[i] Monitor : squeue -u \$USER | grep BENCH-S02"
echo "[i] Results : ${RESULTS}/"
echo "[i] When all jobs finish, summarize with:"
echo "      bash ${BENCH_ROOT}/summarize.sh"
echo
