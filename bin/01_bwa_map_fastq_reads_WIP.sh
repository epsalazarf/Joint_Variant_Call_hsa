#!/bin/bash

# Title: BWA FASTQ READS MAPPER [FENIX]
# About: Quality control for mapped BAM files to ready. Adapted for LAVIS-FENIX.
# Usage: 01_bwa_map_fastq_reads.sh [INPUT BAM] [OUTPUT PATH]
# Authors: Pavel Salazar-Fernandez (this version), AH, EA, FASQ, MCAA
# Source: GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Version notes: [+ALT] Reworked the steps in encapsulated functions for improved workflow management.

#<DEBUG>
set -e

#<START>
echo;echo -e "<START> GATK4 BAM QC +ALT [FENIX]"
echo "Started: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "> SETUP: Check Input Files >>"
SAMPLE_NAME=$1
OUTPUT_PATH=$2

#TODO
if ls ./${SAMPLE_NAME}*fq.gz 1> /dev/null 2>&1; then
    echo "> Sample name: ${SAMPLE_NAME}"
    BAM_base=${BAM_FILE%.*am}
    BAM_prefix=${BAM_FILE%%.*am}
    OUTPUT_PATH="${OUTPUT_PATH:-$(pwd)}"
    echo "> Output folder: ${OUTPUT_PATH}"
else
    echo -e "<ERROR> BWA MX-MAPPER CANCELLED. FASTA Files starting with ${SAMPLE_NAME}* not found."
    exit 1
fi
##</INPUT>

##<ENVIROMENT>
# Threads (no significant benefit when >4)
njobs=4

# Options
BQSR_COV=true  # Run Step 3d (auto-off if remote)
RUN_METRICS=false   # Run Step 5

# Path to config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [ -n "$SSH_CLIENT" ] || [ -n "$SSH_TTY" ] || [ -n "$SSH_CONNECTION" ]; then
    env_type="remote"
else
    env_type="local"
fi

echo "> Running on $env_type environment (SSH session detected)."

# Parse YAML-ish config into Bash variables
# This strips indentation, removes quotes, and exports as key=value
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

# Load modules if run on remote server (work-around due to faulty  parser [ARC02])
if [ $env_type == "remote" ] ; then
    echo "- Loading required modules "
    module load fastqc
    module load bwa
fi 

#GUARD: Reference files paths
if [ -z "$ref_gnm" ] || [ -z "$ref_vars" ]; then
    echo "<ERROR> Missing required reference paths. Please check your config file: $CONFIG_FILE"
    exit 1
fi

#GUARD: Reference files availability
echo "> References used: "
[ -f "${ref_gnm}" ] || ( echo " ERROR: Reference Genome not found (${ref_gnm})." ; exit 1 )
echo "  - Genome: ${ref_gnm}"
[ -f "${ref_vars}" ] || ( echo " ERROR: Reference Variants not found (${ref_vars})." ; exit 1 )
echo "  - Variants: ${ref_vars}" 

#</ENVIROMENT>

#<FUNCTIONS>

# QC STEP 0: Add ReadGroup
step1_reads_quality() {
    local step_name="Step 0: ReadGroup Assignation"
    local infile="$BAM_FILE"
    local outfile="${OUTPUT_PATH}/${BAM_base}.rg.bam"

    echo;echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || ( echo "<ERROR> Missing input: $infile" ; exit 1 )

    #SKIP
    [ -s "$outfile" ] && ( echo " [!] $step_name already completed ($outfile exists)" ; return )
    
    #COMMAND
    gatk AddOrReplaceReadGroups \
        --INPUT "$infile" \
        --RGID "$BAM_prefix" \
        --RGSM "$BAM_prefix" \
        --RGPL "ILLUMINA" \
        --RGPU "$BAM_prefix" \
        --RGLB "$BAM_prefix" \
        --VERBOSITY ERROR \
        --OUTPUT "${outfile}.tmp"

        #DETEMP
    mv "${outfile}.tmp" "$outfile"

    #CHECK
    [ -s "$outfile" ] || ( echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile" ; exit 1 )
    echo " [DONE] $step_name"
}

# -------------------
# STEP 1: QC - Run fastqc
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming required, quality looks okay.



# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference
bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

# QC STEP 2: Add ReadGroup
step2_add_readgroup() {
    local step_name="Step 0: ReadGroup Assignation"
    local infile="$BAM_FILE"
    local outfile="${OUTPUT_PATH}/${BAM_base}.rg.bam"

    echo;echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || ( echo "<ERROR> Missing input: $infile" ; exit 1 )

    #SKIP
    [ -s "$outfile" ] && ( echo " [!] $step_name already completed ($outfile exists)" ; return )
    
    #COMMAND
    gatk AddOrReplaceReadGroups \
        --INPUT "$infile" \
        --RGID "$BAM_prefix" \
        --RGSM "$BAM_prefix" \
        --RGPL "ILLUMINA" \
        --RGPU "$BAM_prefix" \
        --RGLB "$BAM_prefix" \
        --VERBOSITY ERROR \
        --OUTPUT "${outfile}.tmp"

        #DETEMP
    mv "${outfile}.tmp" "$outfile"

    #CHECK
    [ -s "$outfile" ] || ( echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile" ; exit 1 )
    echo

#<SANDBOX>

#!/bin/bash

map_reads_per_sample() {
    SAMPLE_NAME="${1:-$(pwd)}"
    REF="ref.fa"

    # find all unique lane IDs for this sample prefix
    # expects files like SAMPLE_L2_1.fq.gz / SAMPLE_L2_2.fq.gz
    lanes=$(ls ${SAMPLE_NAME}*_1.fq.gz | sed -E 's/.*(L[0-9]+)_1.fq.gz/\1/' | sort -u)

    for lane in $lanes; do
        R1="${SAMPLE_NAME}_${lane}_1.fq.gz"
        R2="${SAMPLE_NAME}_${lane}_2.fq.gz"

        # read group ID = sample + lane
        RGID="$(basename ${SAMPLE_NAME})_${lane}"
        RGLB="$(basename ${SAMPLE_NAME})"
        RGSM="$(basename ${SAMPLE_NAME})"
        RGPL="ILLUMINA"

        if [[ -f "$R2" ]]; then
            # paired-end
            bwa mem -t 8 -R "@RG\tID:${RGID}\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}" \
                "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "${RGID}.bam"
        else
            # single-end
            bwa mem -t 8 -R "@RG\tID:${RGID}\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}" \
                "$REF" "$R1" | samtools sort -@ 4 -o "${RGID}.bam"
        fi
    done

    echo "Generated per-lane BAMs for ${SAMPLE_NAME}. Next: merge & mark duplicates per sample."
}

