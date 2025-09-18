#!/bin/bash

# NOTES
# This script adapts to a project where the BAM files are already worked the duplicates
# No forget to do this step before run this script
# The MarkDuplicates step will be commented, change the output-input flow if you want to incorporate it

# USAGE
# sh BP_GATK_BSQR.sh file.bam 6 /path/output true ".hs37d5.fa.dedup" 

BAM=$1;
CHROM=$2;
OUTPUT_PATH=$3;
AnalyzeCovariates=$4; #bool
ALT_NAME=$5; # It can be nothing

sample_name=${BAM##*/};
sample_name=${sample_name%%.*am}; # Remove all after the first dot
#sample_name=${sample_name%%.[sbc]am}; # Remove just the extension

# Paths
PATH_GATK=/software/gatk;
PATH_ASSEMBLY="/data/reference_genomes/hs37d5.fa";
PATH_DBSNP="/data/Project_Eduardo/references/All.vcf";

# Data pre-processing
## Mark Duplicates
# [ ! -d ${OUTPUT_PATH}/marked_duplicates ] && mkdir ${OUTPUT_PATH}/marked_duplicates
# ${PATH_GATK} MarkDuplicatesSpark -I ${BAM} -L ${CHROM} -R ${PATH_ASSEMBLY} -O ${OUTPUT_PATH}/marked_duplicates/${sample_name}${ALT_NAME}_marked_duplicates.bam

# Base (Quality Score) Recalibration
## BaseRecalibrator
[ ! -d ${OUTPUT_PATH}/tables ] && mkdir ${OUTPUT_PATH}/tables
${PATH_GATK} BaseRecalibrator -I ${BAM} -R ${PATH_ASSEMBLY} -L ${CHROM} --known-sites ${PATH_DBSNP} \
															-O ${OUTPUT_PATH}/tables/${sample_name}${ALT_NAME}_chr${CHROM}_recalibration.table

## Apply base quality score recalibration
[ ! -d ${OUTPUT_PATH}/BAMs ] && mkdir ${OUTPUT_PATH}/BAMs
${PATH_GATK} ApplyBQSR -I ${BAM} -R ${PATH_ASSEMBLY} -L ${CHROM} --bqsr-recal-file ${OUTPUT_PATH}/tables/${sample_name}${ALT_NAME}_chr${CHROM}_recalibration.table \
											-O ${OUTPUT_PATH}/BAMs/${sample_name}${ALT_NAME}_chr${CHROM}.bsqr.bam

## Evaluate and compare base quality score recalibration tables
if [ ${AnalyzeCovariates} = true ]
then
	${PATH_GATK} BaseRecalibrator -I ${OUTPUT_PATH}/BAMs/${sample_name}${ALT_NAME}_chr${CHROM}.bam -R ${PATH_ASSEMBLY} -L ${CHROM} \
															--known-sites ${PATH_DBSNP} -O ${OUTPUT_PATH}/tables/${sample_name}${ALT_NAME}_chr${CHROM}_recalibration_after.table

	[ ! -d ${OUTPUT_PATH}/plots ] && mkdir ${OUTPUT_PATH}/plots
	${PATH_GATK} AnalyzeCovariates -before ${OUTPUT_PATH}/tables/${sample_name}${ALT_NAME}_chr${CHROM}_recalibration.table \
																	-after ${OUTPUT_PATH}/tables/${sample_name}${ALT_NAME}_chr${CHROM}_recalibration_after.table \
																	-plots ${OUTPUT_PATH}/plots/${sample_name}${ALT_NAME}_chr${CHROM}_AnalyzeCovariates.pdf
fi

#rm -r ${OUTPUT_PATH}/marked_duplicates
#rm -r ${OUTPUT_PATH}/tables
