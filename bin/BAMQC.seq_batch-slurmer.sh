#!/bin/bash
# Title: BAMQC Batch Slurmer (sequential, tolerant to failed jobs)
# Usage: ./BAMQC.batch-slurmer batch-list.txt 02_gatk_bam_qc_workflow.sh files-directory
# Each line in batch list: sample_name    sample_file    readgroup_string (with literal \t)

BATCH_LIST="$1"
BATCH_SCRIPT="$2"
DIR_PATH="$3"



if [[ -z "$BATCH_LIST" || -z "$BATCH_SCRIPT" || -z "$TARGET_DIR" ]]; then
  echo "Usage: $0 <batch-list.txt> <02_gatk_bam_qc_workflow.sh> <files directory>"
  exit 1
fi

if [[ "${DIR_PATH: -1}" != "/" ]]; then
  DIR_PATH="${DIR_PATH}/"
fi

# Function to submit job with optional dependency
submit_job() {
  local dependency="$1"
  shift
  local script="$1"
  shift
  local args=("$@")

  local jobname="BAMQC-${sample}-$EPOCHSECONDS"

  if [[ -n "$dependency" ]]; then
    jobid=$(sbatch --dependency=afterany:"$dependency" \
                   --job-name="${jobname}" \
                   --nodes=1 --ntasks=1 --cpus-per-task=4 \
                   --output=%x.%j.log \
                   --wrap "bash '$script' '${args[0]}' '${args[1]}' '${args[2]}'" \
            | awk '{print $4}')
  else
    jobid=$(sbatch --job-name="${jobname}" \
                   --nodes=1 --ntasks=1 --cpus-per-task=4 \
                   --output=%x.%j.log \
                   --wrap "bash '$script' '${args[0]}' '${args[1]}' '${args[2]}'" \
            | awk '{print $4}')
  fi

  echo "$jobid"
}

prev_jobid=""

while IFS=$'\t' read -r sample bamfile rgstring; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue

  echo
  echo "Preparing: $sample ($bamfile)"
  mkdir -p "$sample"
  cd "$sample" || exit
  ln -sf "${DIR_PATH}${bamfile}" .

  echo "Submitting job for: $sample"
  echo "  RG string: $rgstring"

  # Submit job and chain after previous (continue even on fail)
  jobid=$(submit_job "$prev_jobid" "$BATCH_SCRIPT" "$bamfile" "." "$rgstring")
  echo "  Job ID: $jobid"

  prev_jobid="$jobid"
  cd ..
  sleep 10
done < "$BATCH_LIST"

echo
echo "All jobs submitted. Last job in chain: $prev_jobid"
