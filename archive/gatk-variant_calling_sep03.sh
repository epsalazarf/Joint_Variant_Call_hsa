#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# This script is for demonstration purposes only

echo "START> GATK4 Variant Calling"

#<SAMPLE>
# download data
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz
#</SAMPLE>

#echo "Run Prep files..."

################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################
echo "STEP 0: Prepare References"
#<REFERENCES>

# download reference files
#wget -P ~/Desktop/demo/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#gunzip ~/Desktop/demo/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
#samtools faidx ~/Delta/Rho/hg38.fa


# ref dict - .dict file before running haplotype caller
#gatk CreateSequenceDictionary -R ~/Delta/Rho/hg38/hg38.fa -O ~/Delta/Rho/hg38/hg38.dict


# download known sites files for BQSR from GATK resource bundle
#wget -P ~/Delta/Rho/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
#wget -P ~/Delta/Rho/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

#</REFERENCES>

###################################################### VARIANT CALLING STEPS ####################################################################

# Directories
ref="/Users/epsal/Data/Rho/hg38.fa"
known_sites="/Users/epsal/Data/Rho/dbSNP157.vcf.gz"
aligned_reads="/Users/epsal/Data/GATK_varcall_test/aligned_reads"
reads="/Users/epsal/Data/GATK_varcall_test/reads"
results="/Users/epsal/Data/GATK_varcall_test/results"
data="/Users/epsal/Data/GATK_varcall_test/data"

#<NEGATE>
if false; then


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


# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

# Directories (in-Docker)
# ref="/gatk/Delta/Rho/hg38.fa"
# known_sites="/gatk/Delta/Rho/dbSNP157.vcf.gz"
# aligned_reads="/gatk/Delta/GATK_varcall_test/aligned_reads"
# reads="/gatk/Delta/GATK_varcall_test/reads"
# results="/gatk/Delta/GATK_varcall_test/results"
# data="/gatk/Delta/GATK_varcall_test/data"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
#gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data_table.txt


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data_table.txt -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam




# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf




# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf.gz

# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf.gz --select-type SNP -O ${results}/raw_snps.vcf.gz
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf.gz --select-type INDEL -O ${results}/raw_indels.vcf.gz

#</NEGATE>
fi; echo "^ CHAIN NEGATED"

# ----------------------------------------------
# STEP 7: Filter and rename canonical chroms
# ----------------------------------------------

echo "STEP 7: Filter and rename canonical chroms"

canon_chr="/Users/epsal/Data/Rho/canon_hgvs_contigs.hg38.txt"
rename_chr="/Users/epsal/Data/Rho/canon_hgvs_to_chroms.hg38.txt"


#bcftools filter ${results}/raw_snps.vcf.gz -i "CHROM=@${canon_chr}" -Oz -o ${results}/raw_snps.canon_hgvs.vcf.gz  --threads 8

bcftools annotate ${results}/raw_snps.canon_hgvs.vcf.gz --rename-chrs ${rename_chr} --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o ${results}/raw_snps.canon_chr.vcf.gz -W --threads 8


#bcftools annotate $1 --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o ${1%.vcf*}.cpra.vcf.gz --threads "$mana" -W
#<BREAK>
exit 0
