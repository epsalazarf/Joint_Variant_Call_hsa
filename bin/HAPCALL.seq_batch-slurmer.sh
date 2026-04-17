#!/bin/bash
# Title: HAPCALL Batch Slurmer
# Usage: ./HAPCALL.batch-slurmer 03_gatk_haplotype_caller.sh [input directory]
# Each line in batch list: sample_name    sample_file
#
# NOTE: Run this slurmer only after all Script 02 (BAMQC) jobs have completed.
# Intra-pipeline sequencing (01->02->03) is enforced by the operator, not SLURM.
# Each sample is submitted independently; samples run in parallel.

BATCH_SCRIPT="$(realpath $1)"
DIR_PATH="${2:-$PWD}"
BAM_SUFFIX="*.rmdup.mqfilt.bqsr.bam"

# CHECK: Script parameters
if [[ -z "$BATCH_SCRIPT" || -z "$DIR_PATH" ]]; then
  echo "Usage: $0 <03_gatk_haplotype_caller.sh> <dirpath>"
  exit 1
fi

#CHECK: Directory exists
DIR_PATH=$(realpath $DIR_PATH)
if [ -d "$DIR_PATH" ]; then
  echo "[INFO] Moving to: $DIR_PATH"
  cd "$DIR_PATH"
else
  echo "<ERROR> Directory does not exist. Exiting..."
  exit 1
fi

#CHECK: Directory name format
if [[ "${DIR_PATH: -1}" != "/" ]]; then
  DIR_PATH="${DIR_PATH}/"
fi


#PREPARATIONS: Find & Create BAM Files Table
BATCH_LIST="$PWD/$(basename $PWD).hapcall_input_files.txt"
if [ -s "$BATCH_LIST" ]; then
  echo "[INFO] Reading files from pre-existing batch list: $BATCH_LIST"
else
  find  . -name ${BAM_SUFFIX} | awk 'BEGIN{FS="/";OFS="\t"}{print $(NF-1),$0}' > "$BATCH_LIST"
  [[ -s "$BATCH_LIST" ]] && echo "[INFO] Creating batch list from found files: $BATCH_LIST"
fi

# Function to submit job independently (no inter-sample dependency)
submit_job() {
  local script="$1"
  shift
  local args=("$@")
  local jobname="HAPCALL-${sample//[^A-Za-z0-9_]/_}-$EPOCHSECONDS"

  jobid=$(sbatch --job-name="${jobname}" \
                 --nodes=1 --ntasks=1 --cpus-per-task=4 \
                 --mem=32G \
                 --output=%x.%j.log \
                 --wrap "bash '$script' '${args[0]}' '$(dirname "${args[0]}")'" \
          | awk '{print $4}')

  echo "$jobid"
}

while IFS=$'\t' read -r sample bamfile ; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue

  cd "$(dirname "$bamfile")" || exit 1
  echo
  echo "Submitting job for: $sample ($bamfile)"

  jobid=$(submit_job "$BATCH_SCRIPT" "$(basename $bamfile)")
  echo "  Job ID: $jobid"

  cd "$DIR_PATH"
  sleep 1

done < "$BATCH_LIST"

echo
echo "All jobs submitted."
