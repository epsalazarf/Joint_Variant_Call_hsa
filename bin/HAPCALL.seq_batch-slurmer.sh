#!/bin/bash
# Title: HAPCALL Batch Slurmer (sequential)
# Usage: ./BAMQC.batch-slurmer 03_gatk_haplotype_caller.sh
# Each line in batch list: sample_name    sample_file

BATCH_SCRIPT="$(realpath $1)"
DIR_PATH="$2"
BAM_SUFFIX="*.rmdup.mqfilt.bqsr.bam"

# CHECK: Script parameters
if [[ -z "$BATCH_SCRIPT" || -z "$DIR_PATH" ]]; then
  echo "Usage: $0 <03_gatk_haplotype_caller.sh> <dirpath>"
  exit 1
fi

#CHECK: Directory exists
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
find  . -name "*.rmdup.mqfilt.bqsr.bam" | awk 'BEGIN{FS="/";OFS="\t"}{print $(NF-1),$0}' > "$BATCH_LIST"

[[ -s "$BATCH_LIST" ]] && echo "[INFO] Batch list created: $BATCH_LIST"

# Function to submit job with optional dependency
submit_job() {
  local dependency="$1"
  shift
  local script="$1"
  shift
  local args=("$@")
  local jobname="HAPCALL-${sample//[^A-Za-z0-9_]/_}-$EPOCHSECONDS"


  if [[ -n "$dependency" ]]; then
    jobid=$(sbatch --dependency=afterany:"$dependency" \
                   --job-name="${jobname}" \
                   --nodes=1 --ntasks=1 --cpus-per-task=2 \
                   --mem=20G \
                   --output=%x.%j.log \
                   --wrap "bash '$script' '${args[0]}'" \
            | awk '{print $4}')
  else
    jobid=$(sbatch --job-name="${jobname}" \
                   --nodes=1 --ntasks=1 --cpus-per-task=2 \
                   --mem=20G \
                   --output=%x.%j.log \
                   --wrap "bash '$script' '${args[0]}' '$(dirname "${args[0]}")'" \
            | awk '{print $4}')
  fi

  echo "$jobid"
}

prev_jobid=""

while IFS=$'\t' read -r sample bamfile ; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue

  cd "$(dirname "$bamfile")" || exit 1
  echo
  echo "Submitting job for: $sample ($bamfile)"

  # Submit job and chain after previous (continue even on fail)
  jobid=$(submit_job "$prev_jobid" "$BATCH_SCRIPT" "$(basename $bamfile)")
  echo "  Job ID: $jobid"

  prev_jobid="$jobid"
  cd "$DIR_PATH"
  sleep 10

done < "$BATCH_LIST"

echo
echo "All jobs submitted. Last job in chain: $prev_jobid"
