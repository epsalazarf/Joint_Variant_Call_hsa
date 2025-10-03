#!/bin/bash

# Title: GATK4 BAM QC [FENIX]
# About: Quality control for mapped BAM files to ready. Adapted for LAVIS-FENIX.
# Usage: 02_gatk_bam_qc_pipeline.sh [INPUT BAM] [OUTPUT PATH]
# Authors: Pavel Salazar-Fernandez (this version), AH, EA, FASQ, MCAA
# Source: GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

#<DEBUG>
set -e

#<START>
echo;echo -e "<START> GATK4 BAM QC [FENIX]"
echo "Started: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "> SETUP: Check Input Files >>"
BAM=$1
OUTPUT_PATH=$2

if [ -f ${BAM} ]; then
    echo "> Input file: ${BAM}"
    BAM_base=${BAM%.*am}
    BAM_prefix=${BAM%%.*am}
    OUTPUT_PATH="${OUTPUT_PATH:-$(pwd)}"
    echo "> Output folder: ${OUTPUT_PATH}"
else
    echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${BAM}"
    exit 1
fi
##</INPUT>

##<ENVIROMENT>
# Threads (no significant benefit when >4)
njobs=4

# Path to config file (relative to repo root)
CONFIG_FILE="$(dirname $(readlink -f $0))/../config/config.yaml"

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

# {ARC02}

# Load modules if run on remote server (work-around due to faulty module parser [ARC02])
if [ $env_type == "remote" ] ; then
    echo "- Loading required modules "
    module load gatk
    module load samtools
    module load mosdepth
    module load r
fi 

# Sanity checks for required references
if [ -z "$ref_gnm" ] || [ -z "$ref_vars" ]; then
    echo "[!] WARNING: Missing required reference paths. Please check your config file: $CONFIG_FILE"
    exit 1
fi

# {ARC01}


echo "> References: "
[ -f "${ref_gnm}" ] && echo "  - Genome: ${ref_gnm}" || ( echo " ERROR: Reference Genome not found (${ref_gnm})." ; exit 1 )
[ -f "${ref_vars}" ] && echo "  - Variants: ${ref_vars}" ||( echo " ERROR: Reference Variants not found (${ref_vars})." ; exit 1 )

#<MAIN>
# QC STEP 0: Add ReadGroup
echo;echo -e "> QC Step 0: ReadGroup Assignation  >>"
gatk AddOrReplaceReadGroups \
    --INPUT ${BAM} \
    --RGID ${BAM_prefix} \
    --RGSM ${BAM_prefix} \
    --RGPL "ILLUMINA" \
    --RGPU ${BAM_prefix} \
    --RGLB ${BAM_prefix} \
    --VERBOSITY ERROR \
    --OUTPUT ${OUTPUT_PATH}/${BAM_base}.rg.bam

[ ! -f ${OUTPUT_PATH}/${BAM_base}.rg.bam ] && ( echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${OUTPUT_PATH}/${BAM_base}.rg.bam"; exit 1 )

# QC STEP 1: GATK Remove duplicates
echo;echo -e "> QC Step 1: Remove Duplicates  >>"
gatk MarkDuplicatesSpark \
    --input ${BAM_base}.rg.bam \
    --conf 'spark.executor.cores=8' \
    --remove-all-duplicates \
    --verbosity ERROR \
    --metrics-file ${OUTPUT_PATH}/${BAM_base}-dups.txt \
    --output ${OUTPUT_PATH}/${BAM_base}.rmdup.bam

# Picard-based (non-recommended, requires an additional sort step, not included)
# gatk MarkDuplicates\
#     --INPUT ${BAM} \
#     --REMOVE_DUPLICATES \
#     --VERBOSITY ERROR \
#     --METRICS_FILE ${OUTPUT_PATH}/${BAM_base}-dups.txt \
#     --OUTPUT ${OUTPUT_PATH}/${BAM_base}.rmdup.bam

[ ! -f ${OUTPUT_PATH}/${BAM_base}.rmdup.bam ] && ( echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${OUTPUT_PATH}/${BAM_base}.rmdup.bam"; exit 1 )

# QC STEP 2: Filter by mapping quality
echo;echo -e "> QC Step 2: Mapping Quality Filter [samtools] >>"
echo -e "  > Command: samtools view -F 4 -q 30 ${OUTPUT_PATH}/${BAM_base}.rmdup.bam"
samtools view ${OUTPUT_PATH}/${BAM_base}.rmdup.bam \
    -F 4 -q 30 \
    -@ ${njobs} \
    --write-index \
    --save-counts ${OUTPUT_PATH}/${BAM_base}.mqfilt-counts.txt \
    -b --output ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam

[ -f ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam ] \
    && ( echo "  > Step 2 finished."; cat ${OUTPUT_PATH}/${BAM_base}.mqfilt-counts.txt ) \
    || ( echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam" ; exit 1 )

# QC STEP 3: Base Quality Recalibration
echo;echo -e "> QC Step 3: Base Quality Recalibration  >>"

# 3.1) Build the model for recalibration
echo;echo "> QC Step 3a: Model Building  >>"
gatk BaseRecalibrator \
    --input ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam \
    --reference ${ref_gnm} --known-sites ${ref_vars} \
    --verbosity ERROR \
    --output ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table.txt

# 3.2) Apply the model to adjust the base quality scores
echo;echo "> QC Step 3b: Base recalibration  >>"
gatk ApplyBQSR \
    --input ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam \
    --reference ${ref_gnm} \
    --bqsr-recal-file ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table.txt \
    --create-output-bam-index \
    --verbosity ERROR \
    --output ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam

[ ! -f ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bam ] && ( echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam" ; exit 1 )

# 3.3) Evaluate and compare base quality score recalibration tables
echo;echo "> QC Step 3c: Recalculate base quality score  >>"
gatk BaseRecalibrator \
    --input ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam \
    --reference ${ref_gnm} --known-sites ${ref_vars} \
    --verbosity ERROR \
    --output ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table_recal.txt




# 3.4) Apply the model to adjust the base quality scores
## ERROR: R LIBRARY NOT AVAILABLE ON FENIX
echo;echo "> QC Step 3d: Analyze Covariates  >>"
#<NEGATE>
if false; then
gatk AnalyzeCovariates \
    -before ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table.txt \
    -after ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr_table_recal.txt \
    --verbosity ERROR \
    -plots ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.AnalyzeCovariates.pdf \

#</NEGATE>
fi; echo; echo " [!] STEP OMITTED"

# QC STEP 4: Calculation of coverage with mosdepht from samtools
echo;echo "> QC Step 4: Calculate coverage [mosdepth] >>"
mosdepth \
    --threads 4 \
    --fast-mode \
    --no-per-base \
    --fasta ${ref_gnm} \
    ${OUTPUT_PATH}/${BAM_base}.rmb \
    ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam \

[ -f ${OUTPUT_PATH}/${BAM_base}.rmb.mosdepth.summary.txt ] \
    && ( echo "  > Step 4: 'mosdepth' finished.") \
    || ( echo -e "<ERROR> GATK4 BAM QC CANCELLED. File not found: ${OUTPUT_PATH}/${BAM_base}.rmb.mosdepth.summary.txt" ; exit 1 )

# QC STEP 5: Alignment & Insert Size Metrics (optional)
echo;echo -e "> QC Step 5: Collect Alignment & Insert Size Metrics  >>"

#<NEGATE>
if false; then

echo;echo "  > QC Step 5a: Alignment summary metrics >>"
gatk CollectAlignmentSummaryMetrics \
    --INPUT ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam \
    --REFERENCE_SEQUENCE ${ref_gnm} \
    --VERBOSITY ERROR \
    --OUTPUT ${OUTPUT_PATH}/${BAM_base}.alignment_metrics.txt

echo;echo "  > QC Step 5b: Insert size metrics >>"
gatk CollectInsertSizeMetrics \
    --INPUT ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam \
    --VERBOSITY ERROR \
    --OUTPUT ${OUTPUT_PATH}/${BAM_base}.insert_size_metrics.txt \
    --Histogram_FILE ${OUTPUT_PATH}/${BAM_base}.insert_size_histogram.pdf

#</NEGATE>
fi; echo; echo " [!] STEP OMITTED"

# EX STEP: Cleaning up files
echo ; echo -e "> EX01 Housekeeping: Remove intermediate files"
if [ -f ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam ]; then
    rm -v \
        ${OUTPUT_PATH}/${BAM_base}.rg.bam \
        ${OUTPUT_PATH}/${BAM_base}.rmdup.ba* \
        ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.ba*
    echo "[!] WARNING: Intermediate files removed."
else
    echo "[!] WARNING: Intermediate files not removed."
fi;

[ -f ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam ] \
    && ( echo; echo -e "<COMPLETED> GATK4 BAM QC [FENIX]"; 
         echo "> Processing Time: $( echo $(( $EPOCHSECONDS - $script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"; 
         echo "> Exiting...";  
         exit 0) \
    || ( echo; echo -e "<ERROR> GATK4 BAM QC [FENIX] UNCOMPLETED. File not found: ${OUTPUT_PATH}/${BAM_base}.rmdup.mqfilt.bqsr.bam"; echo "> Exiting..."; exit 1)

#<SANDBOX>
#{ARC01}
#     ##<REFERENCES>
#     # Human Genome FASTA references (with index)
#     #ref_gnm="/mnt/data/fsanchezq/esalazarf/References/hg38/ref_genome/GRCh38.p14.fa"
#     ref_gnm="/mnt/data/amedina/lupus/WGS/bam_files/hg38.fa"
#     # Human variant database (usually dbSNP)
#     ref_vars="/mnt/data/fsanchezq/esalazarf/References/hg38/dbsnp157/dbSNP157.canon_chr.vcf.gz"
# else
#     echo "> Running on a local machine (no SSH session detected)."
#     ##<REFERENCES>
#     # Human Genome FASTA references (with index)
#     ref_gnm="$HOME/Data/REFx/hg38.fa"
#     # Human variant database (usually dbSNP)
#     ref_vars="$HOME/Data/REFx/dbSNP157.canon_chr.vcf.gz"
# fi

#{ARC02}

# # Load modules if specified
# if [ -n "$modules" ]; then
#     for m in $modules; do
#         echo "- Loading module: ${m}"
#         module load "$m"
#     done
# fi

#<END>
