#!/bin/bash

# Title: GATK4 BAM QC [FENIX]
# About: Quality control for mapped BAM files to ready. Adapted for LAVIS-FENIX.
# Usage: 02_gatk_bam_qc_workflow.sh [INPUT BAM] [OUTPUT PATH] [RG STRING]
# Authors: Pavel Salazar-Fernandez (this version), AH, EA, FASQ, MCAA
# Source: GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Version notes: [+Workflow] Reworked the steps in encapsulated functions for improved workflow management.

#<DEBUG>
set -e

#<START>
echo;echo -e "<START> GATK4 BAM QC +ALT [FENIX]"
echo "Started: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "> SETUP: Check Input Files >>"
BAM_FILE=$1
OUTPUT_PATH=${2:-$PWD}
rg_string=${3:-""}

if [ -f "$BAM_FILE" ]; then
    echo "> Input file: $BAM_FILE"
    echo "> Output folder: ${OUTPUT_PATH}"
    BAM_name=$(basename $1)
    BAM_base=${BAM_FILE%.*am}
    BAM_prefix=${BAM_FILE%%.*am}
    FINAL_FILE="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam"
else
    echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${BAM_FILE}"
    exit 1
fi
##</INPUT>

##<ENVIROMENT>
# Threads (no significant benefit when >4)
njobs=4

# Options
BQSR_COV=false  # Run Step 3d (auto-off if remote)
RUN_METRICS=false   # Run Step 5
HOUSEKEEP=false

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
    module load gatk
    module load samtools
    module load mosdepth
    module load r
    BQSR_COV=false
fi 

#GUARD: Reference files paths
if [ -z "$ref_gnm" ] || [ -z "$ref_vars" ]; then
    echo "<ERROR> Missing required reference paths. Please check your config file: $CONFIG_FILE"
    exit 1
fi

#GUARD: Reference files availability
echo "> References used: "
[ -f "${ref_gnm}" ] || { echo " ERROR: Reference Genome not found (${ref_gnm})."; exit 1; }
echo "  - Genome: ${ref_gnm}"
[ -f "${ref_vars}" ] || { echo " ERROR: Reference Variants not found (${ref_vars})."; exit 1; }
echo "  - Variants: ${ref_vars}" 

#</ENVIROMENT>

#<FUNCTIONS>

# QC STEP 0: Add ReadGroup
step0_add_readgroup() {
    local step_name="Step 0: ReadGroup Assignation"
    local infile="$BAM_FILE"
    local outfile="${OUTPUT_PATH}/${BAM_base}.rg.bam"
    local step_timestamp=$(date +%s)

    echo;echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }
    
    #COMMAND
    if [[ "$rg_string" == @RG* ]]
    then
        local rgID="$(sed -E 's/.*ID:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
        local rgPL="$(sed -E 's/.*PL:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
        local rgPU="$(sed -E 's/.*PU:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
        local rgLB="$(sed -E 's/.*LB:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
        local rgSM="$(sed -E 's/.*SM:([a-zA-Z0-9_.]+).*/\1/' <<< "$rg_string")"
    else
        local rgID="$BAM_prefix"
        local rgPL="ILLUMINA"
        local rgPU="$BAM_prefix"
        local rgLB="$BAM_prefix"
        local rgSM="$BAM_prefix"
    fi

    gatk AddOrReplaceReadGroups \
        --INPUT "$infile" \
        --RGID "$rgID" \
        --RGPL "$rgPL" \
        --RGPU "$rgPU" \
        --RGLB "$rgLB" \
        --RGSM "$rgSM" \
        --VERBOSITY ERROR \
        --OUTPUT "${outfile}.tmp"

    #DETEMP
    mv "${outfile}.tmp" "$outfile"

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"
    echo "  &> Step Time: $( echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

# QC STEP 1: GATK Remove duplicates
step1_mark_duplicates() {
    local step_name="Step 1: Remove Duplicates"
    local step_timestamp=$(date +%s)
    local infile="${OUTPUT_PATH}/${BAM_base}.rg.bam"
    local outfile="${OUTPUT_PATH}/${BAM_base}.rmdup.bam"
    local metrics="${OUTPUT_PATH}/${BAM_base}-dups.txt"

    echo;echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    #COMMAND
    gatk MarkDuplicatesSpark \
        --input "$infile" \
        --conf 'spark.executor.cores=8' \
        --remove-all-duplicates \
        --verbosity ERROR \
        --metrics-file "$metrics" \
        --output "${outfile}.tmp"
    
    #DETEMP
    rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"
    echo "> Step Time: $( echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

# QC STEP 2: Filter by mapping quality
step2_mapping_quality_filter() {
    local step_name="Step 2: Mapping Quality Filter"
    local step_timestamp=$(date +%s)
    local infile="${OUTPUT_PATH}/${BAM_base}.rmdup.bam"
    local outfile="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam"
    local counts="${OUTPUT_PATH}/${BAM_base}.mqfilt-counts.txt"

    echo;echo -e "> [BAMQC] $step_name >>"
    #echo -e "  > Command: samtools view -F 4 -q 30 $infile"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    set -o xtrace
    #COMMAND
    samtools view "$infile" \
        -F 4 -q 30 \
        --threads "$njobs" \
        --write-index \
        --save-counts "$counts" \
        --bam \
        --output "${outfile}.tmp"
    
    set +o xtrace

    #DETEMP
    rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"
    echo "  &> Step Time: $( echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

# QC STEP 3: Base Quality Recalibration
step3_bqsr() {
    local step_name="Step 3: Base Quality Recalibration"
    local step_timestamp=$(date +%s)
    local infile="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam"
    local table="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table.txt"
    local bam_out="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam"
    local table_recal="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table_recal.txt"
    local plot="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.AnalyzeCovariates.pdf"

    echo;echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$bam_out" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    # 3.1) Build the model for recalibration
    echo;echo "> Step 3a: Model Building  >>"
    if [ -s "$table" ]; then
        echo " [SKIP] Model table already exists: $table"
    else
        #COMMAND
        gatk BaseRecalibrator \
            --input "$infile" \
            --reference "$ref_gnm" \
            --known-sites "$ref_vars" \
            --verbosity ERROR \
            --output "${table}.tmp"
        mv "${table}.tmp" "$table"
        echo " [DONE] Recalibration model built"
    fi

    # 3.2) Apply the model to adjust the base quality scores
    echo;echo "> Step 3b: Base Recalibration  >>"
    if [ -s "$bam_out" ]; then
        echo " [SKIP] Recalibrated BAM already exists: $bam_out"
    else
        #COMMAND
        gatk ApplyBQSR \
            --input "$infile" \
            --reference "$ref_gnm" \
            --bqsr-recal-file "$table" \
            --create-output-bam-index \
            --verbosity ERROR \
            --output "${bam_out}.tmp"
        #DETEMP
        rename -v "${bam_out}.tmp" "${bam_out}" "${bam_out}.tmp"*
        echo " [DONE] Recalibrated BAM created"
    fi

    echo "  &> Step Time (3a+3b): $( echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"

    # 3.3) Evaluate and compare base quality score recalibration tables
    #COMMAND
    echo;echo "> Step 3c: Evaluate Recalibration  >>"
    if [ -s "$table_recal" ]; then
        echo " [SKIP] Post-recalibration table exists: $table_recal"
    else
        gatk BaseRecalibrator \
            --input "$bam_out" \
            --reference "$ref_gnm" \
            --known-sites "$ref_vars" \
            --verbosity ERROR \
            --output "${table_recal}.tmp"
        mv "${table_recal}.tmp" "$table_recal"
        echo " [DONE] Post-recalibration table generated"
    fi

    echo " [!] Skipping AnalyzeCovariates (FENIX)"

    # 3.4) Apply the model to adjust the base quality scores
    ## ERROR: Required R LIBRARY NOT AVAILABLE ON FENIX
    echo;echo "> Step 3d: Analyze Covariates  >>"
    
    #SKIP
    [ "$BQSR_COV" = true ] || { echo " [!] Skipping AnalyzeCovariates: disabled by user toggle"; return 0; }

    #SKIP
    if [ -s "$plot" ]; then
        echo " [SKIP] Post-recalibration plots exist: $plot"
    else
        #COMMAND
        gatk AnalyzeCovariates \
            -before "$table" \
            -after "$table_recal" \
            --verbosity ERROR \
            -plots "$plot"
    fi
}

step4_mosdepth() {
# QC STEP 4: Calculation of coverage with mosdepth from samtools
    local step_name="Step 4: Calculate coverage (mosdepth)"
    local step_timestamp=$(date +%s)
    local infile="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam"
    local prefix="${OUTPUT_PATH}/${BAM_base}.rmb"
    local outfile="${prefix}.mosdepth.summary.txt"

    echo;echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"
    echo -e "  > Command: mosdepth --fast-mode --no-per-base $prefix $infile"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }
    
    #COMMAND
    set -o xtrace
    mosdepth \
        --threads "$njobs" \
        --fast-mode \
        --no-per-base \
        --fasta "$ref_gnm" \
        "$prefix" \
        "$infile"
    set +o xtrace

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"
    echo "  &> Step Time: $( echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

# QC STEP 5: Alignment & Insert Size Metrics (optional)
step5_metrics()
{
    local step_name="Step 5: Collect Alignment & Insert Size Metrics"
    local step_timestamp=$(date +%s)
    local infile="${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam"
    local align_metrics="${OUTPUT_PATH}/${BAM_base}.alignment_metrics.txt"
    local insert_metrics="${OUTPUT_PATH}/${BAM_base}.insert_size_metrics.txt"
    local insert_hist="${OUTPUT_PATH}/${BAM_base}.insert_size_histogram.pdf"

    echo; echo -e "> [BAMQC] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"
   
    #TOGGLE
    [ "$RUN_METRICS" = true ] || { echo " [SKIP] $step_name disabled by user toggle"; return 0; }

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    # Sub-step 5a: Alignment metrics
    echo;echo "  > QC Step 5a: Alignment summary metrics >>"
    if [ -s "$align_metrics" ]; then
        echo " [SKIP] Alignment metrics already exist: $align_metrics"
    else
        #COMMAND
        gatk CollectAlignmentSummaryMetrics \
            --INPUT "$infile" \
            --REFERENCE_SEQUENCE "$ref_gnm" \
            --VERBOSITY ERROR \
            --OUTPUT "${align_metrics}.tmp"
        mv "${align_metrics}.tmp" "$align_metrics"
        echo " [DONE] Alignment metrics generated"
    fi

    # Sub-step 5b: Insert size metrics
    echo;echo "  > QC Step 5b: Insert size metrics >>"
    if [ -s "$insert_metrics" ] && [ -s "$insert_hist" ]; then
        echo " [SKIP] Insert size metrics already exist"
    else
        gatk CollectInsertSizeMetrics \
            --INPUT "$infile" \
            --VERBOSITY ERROR \
            --OUTPUT "${insert_metrics}.tmp" \
            --Histogram_FILE "${insert_hist}.tmp"
        mv "${insert_metrics}.tmp" "$insert_metrics"
        mv "${insert_hist}.tmp" "$insert_hist"
        echo " [DONE] Insert size metrics generated"
    fi

    echo "  &> Step Time: $( echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

# HK Step: Remove intermediate files (optional)
housekeeping() {
    echo ; echo -e "> [HK] Housekeeping: Remove intermediate files"

    #TOGGLE
    [ "$HOUSEKEEP" = true ] || { echo " [SKIP] Housekeeping disabled by user toggle"; return 0; }

    if [ -f "$FINAL_FILE" ]; then
        rm -v \
            "${OUTPUT_PATH}"/"${BAM_base}".rg.bam \
            "${OUTPUT_PATH}"/"${BAM_base}".rmdup.ba* \
            "${OUTPUT_PATH}"/"${BAM_base}".rmdup.mqfilt.ba*
        echo "[!] WARNING: Intermediate files removed."
    else
        echo "<ERROR> Final file not found. Intermediate files not removed."
    fi;
}

# Finisher: checks for final output and reports success and timing
finisher(){
    if [ -f "$FINAL_FILE" ]; then
        echo; echo -e "<SUCCESS> GATK4 BAM QC [FENIX] COMPLETED"
        echo -e "> Final output file: ${FINAL_FILE}"
        echo "> Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
        echo "> Exiting..."
        exit 0
    else 
        echo; echo -e "<ERROR> GATK4 BAM QC [FENIX] UNCOMPLETED. Final output file not found: ${FINAL_FILE}"
        echo "> Total Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
        echo "> Exiting..."
        exit 1
    fi
}
#</FUNCTIONS>

#<MAIN>

main() {
    step0_add_readgroup
    step1_mark_duplicates
    step2_mapping_quality_filter
    step3_bqsr
    step4_mosdepth
    step5_metrics
    housekeeping
    finisher
}

main "$@"

#</MAIN>

#<END>
