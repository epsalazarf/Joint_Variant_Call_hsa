#!/bin/bash
# Title: BWAMAP Batch Slurmer
# Usage: ./BWAMAP.seq_batch-slurmer.sh <01_bwa_map_fastq_reads_WIP.sh> [parent_dir]
# Auto-discovers sample subdirectories containing FASTQ files and submits one SLURM job per sample.
#
# NOTE: Each sample is submitted independently; samples run in parallel.
# Intra-pipeline sequencing (01->02->03) is enforced by the operator, not SLURM.
# Run BAMQC slurmer only after all jobs submitted here have completed.

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

# Function to submit job independently (no inter-sample dependency)
submit_job() {
  local script="$1"
  local abs_sample_dir="$2"
  local jobname="BWAMAP-${sample//[^A-Za-z0-9_]/_}-$EPOCHSECONDS"

  jobid=$(sbatch --job-name="${jobname}" \
                 --nodes=1 --ntasks=1 --cpus-per-task=4 \
                 --mem=16G \
                 --output="${abs_sample_dir}/%x.%j.log" \
                 --wrap "cd '${abs_sample_dir}' && bash '${script}' '${abs_sample_dir}'" \
          | awk '{print $4}')

  echo "$jobid"
}

while IFS= read -r sample_dir; do
  [[ -z "$sample_dir" ]] && continue

  sample=$(basename "$sample_dir")
  abs_sample_dir="$(realpath "$sample_dir")"

  echo
  echo "Submitting job for: $sample"
  echo "  Directory: $abs_sample_dir"

  jobid=$(submit_job "$BATCH_SCRIPT" "$abs_sample_dir")
  echo "  Job ID: $jobid"

  sleep 1

done <<< "$sample_dirs"

echo
echo "All jobs submitted."
