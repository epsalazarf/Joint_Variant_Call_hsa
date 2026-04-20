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
  find . \( -name "$BAM_PATTERN_A" -o -name "$BAM_PATTERN_B" \) | sort | \
    awk 'BEGIN{FS="/";OFS="\t"}{
      f=$(NF);
      gsub(/\.sort\.bam$|\.sorted\.bam$/, "", f);
      split(f, parts, "_"); sid = parts[1];
      if (sid in seen) { files[sid] = files[sid] "," $0 }
      else             { seen[sid]=1; files[sid]=$0; order[++n]=sid }
    }
    END{ for(i=1;i<=n;i++){ sid=order[i]; print sid, files[sid], "" } }' > "$BATCH_LIST"
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

  # Resolve each BAM (supports comma-separated split-lane files)
  IFS=',' read -r -a bam_array <<< "$bamfile"
  abs_bams=(); basenames=()
  for bam in "${bam_array[@]}"; do
    abs_bams+=("$(realpath "$bam")")
    basenames+=("$(basename "$bam")")
  done
  bam_list=$(IFS=','; echo "${basenames[*]}")

  echo
  echo "Preparing: $sample (${#abs_bams[@]} file(s))"
  mkdir -p "$sample"
  cd "$sample" || exit
  for abs_bam in "${abs_bams[@]}"; do ln -sf "$abs_bam" .; done

  echo "Submitting job for: $sample"
  [[ -n "$rgstring" ]] && echo "  RG string: $rgstring"

  jobid=$(submit_job "$BATCH_SCRIPT" "$bam_list" "." "$rgstring")
  echo "  Job ID: $jobid"

  cd "$DIR_PATH"
  sleep 1
done < "$BATCH_LIST"

echo
echo "All jobs submitted."
