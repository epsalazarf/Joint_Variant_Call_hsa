#!/bin/bash
# Title: BAMQC Batch Slurmer
# Usage: ./BAMQC.batch-slurmer batch-list.txt 02_gatk_bam_qc_workflow.sh files-directory
# Each line in batch list: sample_name    sample_file    readgroup_string (with literal \t)
#
# NOTE: Run this slurmer only after all Script 01 (BWAMAP) jobs have completed.
# Intra-pipeline sequencing (01->02->03) is enforced by the operator, not SLURM.
# Each sample is submitted independently; samples run in parallel.

BATCH_LIST="$1"
BATCH_SCRIPT="$2"
DIR_PATH="$3"



if [[ -z "$BATCH_LIST" || -z "$BATCH_SCRIPT" || -z "$DIR_PATH" ]]; then
  echo "Usage: $0 <batch-list.txt> <02_gatk_bam_qc_workflow.sh> <files directory>"
  exit 1
fi

if [[ "${DIR_PATH: -1}" != "/" ]]; then
  DIR_PATH="${DIR_PATH}/"
fi

# Function to submit job independently (no inter-sample dependency)
submit_job() {
  local script="$1"
  shift
  local args=("$@")

  local jobname="BAMQC-${sample}-$EPOCHSECONDS"

  jobid=$(sbatch --job-name="${jobname}" \
                 --nodes=1 --ntasks=1 --cpus-per-task=8 \
                 --mem=32G \
                 --output=%x.%j.log \
                 --wrap "bash '$script' '${args[0]}' '${args[1]}' '${args[2]}'" \
          | awk '{print $4}')

  echo "$jobid"
}

while IFS=$'\t' read -r sample bamfile rgstring; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue

  echo
  echo "Preparing: $sample ($bamfile)"
  mkdir -p "$sample"
  cd "$sample" || exit
  ln -sf "${DIR_PATH}${bamfile}" .

  echo "Submitting job for: $sample"
  echo "  RG string: $rgstring"

  jobid=$(submit_job "$BATCH_SCRIPT" "$bamfile" "." "$rgstring")
  echo "  Job ID: $jobid"

  cd ..
  sleep 1
done < "$BATCH_LIST"

echo
echo "All jobs submitted."
