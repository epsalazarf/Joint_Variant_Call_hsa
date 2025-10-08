#!/bin/bash

# Title: BWA FASTQ READS MAPPER [FENIX]
# About: Maps FASTQ read files to Reference Genomes to produce mapped BAMs. Adapted for running on LAVIS-FENIX.
# Usage: 01_bwa_map_fastq_reads.sh [SAMPLE PREFIX] [OUTPUT PATH]
# Authors: Pavel Salazar-Fernandez (this version), AH, EA, FASQ, MCAA
# Source 1: GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Source 2: https://gatk.broadinstitute.org/hc/en-us/articles/360035889471-How-should-I-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs


#<DEBUG>
set -e

#<START>
echo;echo -e "<START> BWA FASTQ READS MAPPER [FENIX]"
echo "Started: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "> SETUP: Check Input Files >>"
SAMPLE_NAME=${1:-$(basename $PWD)}
OUTPUT_PATH=${2:-$PWD}

#GUARD
if ls ./${SAMPLE_NAME}*.fq.gz > /dev/null 2>&1; then
    echo "> Sample name: ${SAMPLE_NAME}"
    echo "> Output folder: ${OUTPUT_PATH}"
else
    echo -e "<ERROR> BWA FASTQ READS MAPPER cancelled. FASTA Files starting with ${SAMPLE_NAME}* not found."
    exit 1
fi
##</INPUT>

##<ENVIROMENT>
# Threads (no significant benefit when >4)
njobs=4

# Options
RUN_FASTQC=false

# Path to config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [ -n "$SSH_CLIENT" ] || [ -n "$SSH_TTY" ] || [ -n "$SSH_CONNECTION" ]; then
    env_type="remote"
    echo "> Running on $env_type environment (SSH session detected)."
else
    env_type="local"
    njobs=8
    echo "> Running on $env_type environment (no SSH session detected)."
fi

# Parse config into Bash variables
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
    module load samtools
fi 

#GUARD: Reference files paths
if [ -z "$ref_gnm" ] ; then
    echo "<ERROR> Missing required reference paths. Please check your config file: $CONFIG_FILE"
    exit 1
fi

#GUARD: Reference files availability
echo "> References used: "
[ -f "${ref_gnm}" ] || ( echo " ERROR: Reference Genome not found (${ref_gnm})." ; exit 1 )
echo "  - Genome: ${ref_gnm}"

#</ENVIROMENT>

#<FUNCTIONS>

# MAP STEP 0: Find input files
step0_map_reads_per_sample() {
    local step_name="Step 0: Map input files"
    local infile="${OUTPUT_PATH}/${BAM_base}.rg.bam"
    local outfile="${OUTPUT_PATH}/${BAM_base}.rmdup.bam"

    echo;echo -e ">> [BWAMAP] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"
    
    file_list=$(find . -name "${SAMPLE_NAME}*.f*q.gz")
    lanes=$(find . -name "${SAMPLE_NAME}*_1.f*q.gz" | sed -E 's/.*(L[0-9]+)_1.f.*q.gz/\1/' | sort -u)

        #GUARD
    [ -z "$lanes" ] && { echo "<ERROR> No files found with pattern: ${SAMPLE_NAME}*_1.f.*q.gz"; exit 1; }

    echo "  Lanes found:" 

    for lane in $lanes; do
        echo "  > ${lane}:"
        R1=$(find . -type f -name "${SAMPLE_NAME}*_${lane}_1.f*q.gz")
        echo "  $R1"
        R2=$(find . -type f -name "${SAMPLE_NAME}*_${lane}_2.f*q.gz")
        echo "  $R2"

        # read group ID = sample + lane
        RGID="$(basename ${SAMPLE_NAME})_${lane}"
        RGLB="$(basename ${SAMPLE_NAME})"
        RGSM="$(basename ${SAMPLE_NAME})"
        RGPL="ILLUMINA"

        echo -n "@RG\tID:${RGID}\tSM:${RGSM}\tLB:${RGLB}\tPL:${RGPL}" > "${OUTPUT_PATH}/${SAMPLE_NAME}_${lane}_bwa_inputs.txt"
        echo -e "\n$R1\n$R2" >> "${OUTPUT_PATH}/${SAMPLE_NAME}_${lane}_bwa_inputs.txt"

    done
    echo;echo " [DONE] $step_name"
}

# MAP STEP 1: Check Quality
step1_reads_quality() {
    local step_name="Step 1: Check Reads Quality"
    local outdir="${OUTPUT_PATH}/${SAMPLE_NAME}-fastqc/"

    echo;echo -e ">> [BWAMAP] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #TOGGLE
    [ "$RUN_FASTQC" =  true ] || { echo " [SKIP] $step_name disabled by user toggle"; return 0; }
    
    mkdir -p "$outdir"
    echo "> Folder created: ${outdir}"

    #COMMAND
    echo " [!] Running: fastqc for all $SAMPLE_NAME reads."
    fastqc \
        --outdir "$outdir" \
        --threads  "$njobs" \
        --noextract \
        --quiet  \
        $(find . -name "${SAMPLE_NAME}*.f*q.gz")

    echo " [DONE] $step_name"
}

# QC STEP 2: Map Reads to Reference
step2_bwa_mapping_per_lane() {
    local step_name="Step 2: BWA Mapping to Reference"

    echo;echo -e "> [BWAMAP] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    for lane in $lanes; do
        echo;echo "> Processing: $SAMPLE_NAME [$lane]..."
        local infile="${OUTPUT_PATH}/${SAMPLE_NAME}_${lane}_bwa_inputs.txt"
        local outfile="${OUTPUT_PATH}/${SAMPLE_NAME}.${lane}.paired.bam"

        #GUARD
        [ -f "$infile" ] || ( echo "<ERROR> Missing input: $infile" ; exit 1 )

        # Inputs array
        declare -a inputs_array
        while IFS= read -r line; do
            inputs_array+=( "$line" )
        done < "${SAMPLE_NAME}_${lane}_bwa_inputs.txt"

        #SKIP
        [ -f "$outfile" ] && ( echo " [!] $step_name already completed ($outfile exists)" ; continue )

        echo "> Mapping reads to create: ${outfile}"
        
        set -o xtrace
        #COMMAND
        bwa mem \
            -t "$njobs" \
            -R "${inputs_array[0]}" \
            -v 1 \
            "$ref_gnm" \
            "${inputs_array[1]}" \
            "${inputs_array[2]}" |
        samtools view \
            --threads "$njobs" \
            --fai-reference "$ref_gnm".fai \
            --bam \
            --output "${outfile}.tmp" -
        set +o xtrace
        
        #DETEMP
        mv "${outfile}.tmp" "$outfile"
        
        #CHECK
        [ -s "$outfile" ] || ( echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile" ; exit 1 )

        echo "> Job finished: ${outfile}"
    
    done


    echo " [DONE] $step_name"
}

# Finisher: checks for final output and reports success and timing
finisher() {
    if [ -f "${OUTPUT_PATH}/${SAMPLE_NAME}.${lane}.paired.bam" ]; then
        echo; echo -e "<SUCCESS> GATK4 BAM QC [FENIX] COMPLETED"
        echo -e "> Last output(s): ${OUTPUT_PATH}/${SAMPLE_NAME}.${lane}.paired.bam"
        echo "> Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
        echo "> Exiting..."
        exit 0
    else 
        echo; echo -e "<ERROR> GATK4 BAM QC [FENIX] UNCOMPLETED. File not found: ${SAMPLE_NAME}.${lane}.paired.bam"
        echo "> Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
        echo "> Exiting..."
        exit 1
    fi
}
#</FUNCTIONS>


#<MAIN>
main() {
    step0_map_reads_per_sample
    step1_reads_quality
    step2_bwa_mapping_per_lane
    #step2a_bwa_mapping_per_sample
    #step3_merge_sample_bams
    finisher
}

main "$@"
#</MAIN>

#<MISC>
#samtools sort Q064_D.L2.paired.bam -O bam -o Q064_D_L2.paired-sort.bam -@ 8
#samtools index Q064_D_L2.paired-sort.bam -@ 8

#<END>
