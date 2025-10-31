#!/bin/bash

# Title: GVCF SINGLE SAMPLE FILTERING [FENIX]
# About: Performs subsampling of a single sample GVCF.
# Usage: 04alt_single_sample_filtering.sh [Raw GVCF] [OUTPUT PATH]
# Authors: Pavel Salazar-Fernandez (this version)
# Source: GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

#<DEBUG>
set -e

#<START>
echo;echo -e "<START> GVCF SINGLE SAMPLE FILTERING [FENIX]"
echo "Started: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "> SETUP: Check Input Files >>"
GVCF_FILE=$1
OUTPUT_PATH=${2:-$PWD}

#GUARD
if [ -f "$GVCF_FILE" ]; then
    echo "> Input file: $GVCF_FILE"
    echo "> Output folder: ${OUTPUT_PATH}"
    GVCF_name=$(basename $1)
    #BAM_base=${GVCF_FILE%.*bam}
    GVCF_prefix=${GVCF_name%%.*.g.vcf.gz}
    FINAL_FILE="${OUTPUT_PATH}/${GVCF_prefix}.filt_snps.canon.g.vcf.gz"
else
    echo -e "<ERROR> GVCF SINGLE SAMPLE FILTERING cancelled. File not found: ${GVCF_FILE}"
    exit 1
fi
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

# Load modules if run on remote server (work-around due to faulty  parser [ARC02])
if [ $env_type == "remote" ] ; then
    echo "- Loading required modules "
    #module load gatk
    module load bcftools
fi 

#</ENVIROMENT>

#<FUNCTIONS>

step0_index_gvcf() {
    local step_name="Extra Step 1: Index GVCF + Counts"
    local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"
    local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz.tbi"

    echo;echo -e ">> [VARCALL] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }
    
    #COMMAND 1: Indexing
    bcftools index "$infile" \
        --threads "$njobs" \
        --force
    
    #COMMAND 2: Print record counts
    bcftools index "$infile" \
        --threads "$njobs" \
        --stats > ${outfile%.vcf.gz.csi}.counts.txt

    #CHECK
    [ -s "$outfile" ] && echo " [DONE] $step_name"
}

# QC STEP 2: SNP FILTER
step1_filter_only_snps() {
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
    bcftools view "$infile" \
        --exclude-types indels,mnps,bnd,other \
        --threads "$njobs" \
        --write-index \
        --output-type b \
        --output "${outfile}.tmp"
    
    set +o xtrace

    #DETEMP
    rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"
}

step1x2_extract_canon_chroms_snps(){
    local step_name="Step 2+3: GVCF Canon Chromosomes SNPs"
    local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_variants.g.vcf.gz"
    local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_snps.canon.g.vcf.gz"

    echo;echo -e ">> [VARCALL] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    set -o xtrace
    #COMMAND
    bcftools view "$infile" \
        --regions "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM" \
        --exclude-types indels,mnps,bnd,other \
        --threads "$njobs" \
        --write-index \
        --output-type b \
        --output "${outfile}.tmp"
    set +o xtrace

    #DETEMP
    rename -v "${outfile}.tmp" "${outfile}" "${outfile}.tmp"*

    #CHECK
    [ -s "$outfile" ] || { echo "<ERROR> CANCELLED: $step_name failed, output missing: $outfile"; exit 1; }
    echo " [DONE] $step_name"

}

step3_sort_chroms() {
    local step_name="Extra Step 2: Sort Chromosomes GVCF"
    local infile="${OUTPUT_PATH}/${BAM_prefix}.raw_snps.canon.g.vcf.gz"
    local outfile="${OUTPUT_PATH}/${BAM_prefix}.raw_snps.canon-sort.g.vcf.gz"

    echo;echo -e ">> [VARCALL] $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    #GUARD
    [ -f "$infile" ] || { echo "<ERROR> Missing input: $infile"; exit 1; }

    #SKIP
    [ -s "$outfile" ] && { echo " [!] $step_name already completed ($outfile exists)"; return 0; }

    set -o xtrace
    
    #COMMAND 1: Indexing
    bcftools view "$infile" \
        --regions "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM" \
        --threads "$njobs" \
        --write-index \
        --output-type b \
        --output "${outfile}.tmp"
    
    set +o xtrace

    #CHECK
    [ -s "$outfile" ] && echo " [DONE] $step_name"
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
finisher() {
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
    step0_index_gvcf
    #step1_filter_only_snps
    #step2_extract_canon_chroms
    step1x2_extract_canon_chroms_snps
    step3_sort_chroms
    #housekeeping
    finisher
}

main "$@"
#</MAIN>

#<END>
