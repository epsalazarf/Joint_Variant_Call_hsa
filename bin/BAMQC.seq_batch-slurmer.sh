#!/bin/bash
# Title: BAMQC Batch Slurmer
# Usage: ./BAMQC.batch-slurmer 02_gatk_bam_qc_workflow.sh [input directory]
# Each line in batch list: sample_name    sample_file    readgroup_string (with literal \t)
#
# NOTE: Run this slurmer only after all Script 01 (BWAMAP) jobs have completed.
# Intra-pipeline sequencing (01->02->03) is enforced by the operator, not SLURM.
# Each sample is submitted independently; samples run in parallel.

BATCH_SCRIPT="$(realpath $1)"
DIR_PATH="${2:-$PWD}"
BAM_PATTERN_A="*.sort.bam"
BAM_PATTERN_B="*.sorted.bam"

# CHECK: Script parameters
if [[ -z "$BATCH_SCRIPT" || -z "$DIR_PATH" ]]; then
  echo "Usage: $0 <02_gatk_bam_qc_workflow.sh> [dirpath]"
  exit 1
fi

# CHECK: Directory exists
DIR_PATH=$(realpath "$DIR_PATH")
if [ -d "$DIR_PATH" ]; then
  echo "[INFO] Moving to: $DIR_PATH"
  cd "$DIR_PATH"
else
  echo "<ERROR> Directory does not exist. Exiting..."
  exit 1
fi

# CHECK: Directory name format
if [[ "${DIR_PATH: -1}" != "/" ]]; then
  DIR_PATH="${DIR_PATH}/"
fi

# PREPARATIONS: Find & Create BAM Files Table
BATCH_LIST="$PWD/$(basename $PWD).bamqc_input_files.txt"
if [ -s "$BATCH_LIST" ]; then
  echo "[INFO] Reading files from pre-existing batch list: $BATCH_LIST"
else
  find . \( -name "$BAM_PATTERN_A" -o -name "$BAM_PATTERN_B" \) | \
    awk 'BEGIN{FS="/";OFS="\t"}{
      f=$(NF);
      gsub(/\.sort\.bam$|\.sorted\.bam$/, "", f);
      print f, $0, ""
    }' > "$BATCH_LIST"
  if [[ -s "$BATCH_LIST" ]]; then
    echo "[INFO] Creating batch list from found files: $BATCH_LIST"
  else
    echo "<ERROR> No BAM files found and no batch list provided. Exiting..."
    exit 1
  fi
fi

# Function to submit job independently (no inter-sample dependency)
submit_job() {
  local script="$1"
  shift
  local args=("$@")

  local jobname="BAMQC-${sample//[^A-Za-z0-9_]/_}-$EPOCHSECONDS"

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

  abs_bamfile="$(realpath "$bamfile")"
  bam_basename="$(basename "$bamfile")"

  echo
  echo "Preparing: $sample ($abs_bamfile)"
  mkdir -p "$sample"
  cd "$sample" || exit
  ln -sf "$abs_bamfile" .

  echo "Submitting job for: $sample"
  [[ -n "$rgstring" ]] && echo "  RG string: $rgstring"

  jobid=$(submit_job "$BATCH_SCRIPT" "$bam_basename" "." "$rgstring")
  echo "  Job ID: $jobid"

  cd "$DIR_PATH"
  sleep 1
done < "$BATCH_LIST"

echo
echo "All jobs submitted."
