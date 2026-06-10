#!/bin/bash

# Title: GVCF CHROM SPLIT (Reprise) [FENIX]
# About: Reads chromosome list from (G)VCF file and splits it by chromosome. [Reprise]
# Usage: gvcf_chrom_split.sh [*.raw_variants.canon_chr.g.vcf.gz] [opt: output folder location]

#<START>
echo;echo -e "<START> GVCF CHROMOSOME SPLIT [FENIX]"
echo "Time: $(date)"
script_timestamp=$(date +%s)

##<INPUT>
echo;echo -e "<SETUP> Check Input Files >>"
GVCF_FILE=$1
GVCF_PATH=$(readlink -f $GVCF_FILE)
OUTPUT_PATH=${2:-$PWD}

#GUARD
if [ -f "$GVCF_FILE" ]; then
    echo "> Input file: $GVCF_FILE"
    echo "> Output folder: ${OUTPUT_PATH}"
    GVCF_name=$(basename $1)
    GVCF_prefix=${GVCF_name%%.*.g.vcf.gz}
    GVCF_prefix=${GVCF_prefix%%.vcf*}
else
    echo -e "<ERROR> GVCF CHROM SPLIT (Reprise). File not found: ${GVCF_FILE}"
    exit 1
fi
##</INPUT>

# Threads (no significant benefit when >2)
njobs=2

# Load modules if run on remote server (work-around due to faulty  parser [ARC02])
echo "- Loading required module: bcftools "
module load bcftools

split_chroms_gvcf() {
    local step_name="Split Chromosomes (G)VCFs"
    local infile="${GVCF_PATH}"
    local outdir="${OUTPUT_PATH}/chrom_gvcf"
    local outfile="${outdir}/${GVCF_prefix}.raw_vars"

    echo
    echo -e "<EXEC> $step_name >>"
    echo "  &> $(date +%Y%m%d-%H%M)"

    # GUARD
    if [[ ! -f "$infile" ]]; then
        echo "<ERROR> Missing input: $infile" >&2
        return 1
    fi

    echo "> Creating chrosomosome directory..."
    mkdir -pv $outdir || {
        echo "<ERROR> Cannot create output directory: $outdir" >&2
        return 1
    }

    echo "> Splitting ${infile}..."

    local chroms
    if ! chroms=$(bcftools index -s "$infile" | cut -f1); then 
        echo "<ERROR> Failed to extract chromosome list from index" >&2
        return 1
    fi

    if [[ -z "$chroms" ]]; then
        echo "<ERROR> No chromosomes found in $infile" >&2
        return 1
    fi

    local count=0

    echo "> Creating chrosomosome directory..."
    mkdir -pv "${OUTPUT_PATH}/chrom_gvcf"
    
    
    #COMMAND
    while read -r C; do

        [[ -z "$C" ]] && continue
        
        echo
        echo "> Extracting chromosome $C..."
        
        if ! bcftools view "$infile" \
            --regions "$C" \
            --threads "$njobs" \
            --write-index=tbi \
            --output-type z \
            --output "${outfile}.${C}.g.vcf.gz"
        then
            echo "<ERROR> bcftools view failed for chromosome: $C" >&2
            return 1
        fi

        ((count++))

    done <<< "$chroms"

    if (( count == 0 )); then
        echo "<ERROR> No chromosome files were produced" >&2
        return 1
    fi


    #CHECK
    echo; echo "[DONE] $step_name ($count chromosomes)"
}

#<MAIN>
split_chroms_gvcf

#FINISH
echo "> Processing Time: $( echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
echo; echo "<EXIT>"
exit 0

#<END>
