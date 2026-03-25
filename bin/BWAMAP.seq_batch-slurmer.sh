#!/bin/bash
# Title: BWAMAP Batch Slurmer (sequential)
# Usage: ./BWAMAP.seq_batch-slurmer.sh <01_bwa_map_fastq_reads_WIP.sh> [parent_dir]
# Auto-discovers sample subdirectories containing FASTQ files and submits one SLURM job per sample.
# Jobs are chained sequentially (--dependency=afterany) to avoid overloading the cluster.

BATCH_SCRIPT="$(realpath "$1")"
DIR_PATH="${2:-$PWD}"

if [[ -z "$BATCH_SCRIPT" ]]; then
  echo "Usage: $0 <01_bwa_map_fastq_reads_WIP.sh> [parent_dir]"
  exit 1
fi

if [ ! -f "$BATCH_SCRIPT" ]; then
  echo "<ERROR> Script not found: $BATCH_SCRIPT"
  exit 1
fi

if [ ! -d "$DIR_PATH" ]; then
  echo "<ERROR> Directory not found: $DIR_PATH"
  exit 1
fi

cd "$DIR_PATH" || exit 1
echo "[INFO] Working directory: $DIR_PATH"

# Auto-discover sample directories: subdirs containing at least one FASTQ R1 file
# Each subdirectory is treated as one sample (directory name = sample name)
sample_dirs=$(find . -mindepth 2 -maxdepth 2 -name "*.f*q.gz" \
  | grep -E '_R?1\.f[^.]*q\.gz$' \
  | sed 's|/[^/]*$||' | sort -u)

if [ -z "$sample_dirs" ]; then
  echo "<ERROR> No FASTQ R1 files found in subdirectories of: $DIR_PATH"
  exit 1
fi

echo "[INFO] Samples found: $(echo "$sample_dirs" | wc -l | tr -d ' ')"

# Function to submit job with optional dependency
submit_job() {
  local dependency="$1"
  local script="$2"
  local abs_sample_dir="$3"
  local jobname="BWAMAP-${sample//[^A-Za-z0-9_]/_}-$EPOCHSECONDS"

  if [[ -n "$dependency" ]]; then
    jobid=$(sbatch --dependency=afterany:"$dependency" \
                   --job-name="${jobname}" \
                   --nodes=1 --ntasks=1 --cpus-per-task=4 \
                   --mem=16G \
                   --output="${abs_sample_dir}/%x.%j.log" \
                   --wrap "cd '${abs_sample_dir}' && bash '${script}' '${abs_sample_dir}'" \
            | awk '{print $4}')
  else
    jobid=$(sbatch --job-name="${jobname}" \
                   --nodes=1 --ntasks=1 --cpus-per-task=4 \
                   --mem=16G \
                   --output="${abs_sample_dir}/%x.%j.log" \
                   --wrap "cd '${abs_sample_dir}' && bash '${script}' '${abs_sample_dir}'" \
            | awk '{print $4}')
  fi

  echo "$jobid"
}

prev_jobid=""

while IFS= read -r sample_dir; do
  [[ -z "$sample_dir" ]] && continue

  sample=$(basename "$sample_dir")
  abs_sample_dir="$(realpath "$sample_dir")"

  echo
  echo "Submitting job for: $sample"
  echo "  Directory: $abs_sample_dir"

  jobid=$(submit_job "$prev_jobid" "$BATCH_SCRIPT" "$abs_sample_dir")
  echo "  Job ID: $jobid"

  prev_jobid="$jobid"
  sleep 10

done <<< "$sample_dirs"

echo
echo "All jobs submitted. Last job in chain: $prev_jobid"
