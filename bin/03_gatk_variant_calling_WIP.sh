#!/bin/bash

# Title: GATK HAPLOTYPE CALLER [FENIX]
# About: Maps FASTQ read files to Reference Genomes to produce mapped BAMs. Adapted for running on LAVIS-FENIX.
# Usage: 003_gatk_variant_calling.sh [Mapped BAM] [OUTPUT PATH]
# Authors: Pavel Salazar-Fernandez (this version), AH, EA, FASQ, MCAA
# Source: GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

#<DEBUG>
set -e

#<START>
echo;echo -e "<START> GATK HAPLOTYPE CALLER [FENIX]"
echo "Started: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "> SETUP: Check Input Files >>"
BAM_FILE=$(readlink -f $1)
OUTPUT_PATH=${2:-$PWD}

#GUARD
if [ -f "$BAM_FILE" ]; then
    echo "> Input file: $BAM_FILE"
    echo "> Output folder: ${OUTPUT_PATH}"
    BAM_name=$(basename $1)
    #BAM_base=${BAM_FILE%.*bam}
    BAM_prefix=${BAM_FILE%%.*bam}
    FINAL_FILE="${OUTPUT_PATH}/${BAM_prefix}.snp.canon.g.vcf.gz"
else
    echo -e "<ERROR> GATK HAPLOTYPE CALLER cancelled. File not found: ${BAM_FILE}"
    exit 1
##</INPUT>

##<ENVIROMENT>
# Threads (no significant benefit when >4)
njobs=4

# Options
SPLIT_VAR_TYPE=true
HOUSEKEEP=false

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
    module load gatk
    module load bcftools
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
[ -f "${ref_vars}" ] || { echo " ERROR: Reference Variants not found (${ref_vars})."; exit 1; }
echo "  - Variants: ${ref_vars}" 

#</ENVIROMENT>

#<FUNCTIONS>

# VARCALL STEP 1: Haplotype Caller
step1_run_haplotype_caller() {
    local step_name="Step 1: Haplotype Caller"
    local infile="$BAM_FILE"
    local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"

    echo;echo -e ">> [VARCALL] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    #COMMAND
    gatk HaplotypeCaller \
            --input "$infile" \
            --reference "$ref_gnm" \
            --dbsnp "$ref_vars" \
            --emit-ref-confidence GVCF \
            --verbosity ERROR \
            --create-output-variant-index \
            --output "${outfile}.tmp"
    #DETEMP
    mv "${outfile}.tmp" "$outfile"

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"
}

# TODO: finish step 2 & 3
# QC STEP 2: Map Reads to Reference
step2_filter_only_snps() {
    local step_name="Step 2: GVCF SNP Filter"
    local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"
    local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_snps.g.vcf.gz"

    echo;echo -e ">> [VARCALL] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    set -o xtrace
    #COMMAND
    bcftools view \
        --exclude-types indels,mnps,bnd,other \
        --threads "$njobs" \
        --write-index \
        --output-type b \
        --output "$outfile"
    
    set +o xtrace

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"

    echo " [DONE] $step_name"
}

step3_extract_canon_chroms(){
    
}

# HK Step: Remove intermediate files (optional)
housekeeping() {
    echo ; echo -e "> [HK] Housekeeping: Remove intermediate files"

    #TOGGLE
    [ "$HOUSEKEEP" ] || { echo " [SKIP] Housekeeping disabled by user toggle"; return 0; }

    if [ -f "$FINAL_FILE}" ]; then
        rm -v \
            asdfg.txt
        echo "[!] WARNING: Intermediate files removed."
    else
        echo "<ERROR> Final file not found. Intermediate files not removed."
    fi;
}

# Finisher: checks for final output and reports success and timing
finisher(){
    if [ -f "$FINAL_FILE}" ]; then
        echo; echo -e "<SUCCESS> GATK HAPLOTYPE CALLER [FENIX] COMPLETED"
        echo -e "> Final output: ${FINAL_FILE}"
        echo "> Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
        echo "> Exiting..."
        exit 0
    else 
        echo; echo -e "<ERROR> GATK HAPLOTYPE CALLER [FENIX] UNCOMPLETED. File not found: ${FINAL_FILE}"
        echo "> Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
        echo "> Exiting..."
        exit 1
    fi
}

#</FUNCTIONS>

#<MAIN>
main() {
    #step0_map_reads_per_sample
    step1_run_haplotype_caller
    #step2_filter_only_snps
    #step3_extract_canon_chroms
    #housekeeping
    finisher
}

main "$@"
#</MAIN>


#<END>
