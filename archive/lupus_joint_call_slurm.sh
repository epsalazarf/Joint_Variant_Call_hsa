#!/usr/bin/env bash

# ============================================================================
# Lupus Joint Variant Calling Pipeline (for SLURM)
# ----------------------------------------------------------------------------
# Author: Adapted from colleague's notes, reformatted by ChatGPT
# Purpose: gVCF calling, joint genotyping, VQSR, LiftOver, and concat for 
#          control and lupus human WGS cohorts.
# Requirements: GATK, bcftools, qsub (or adapted to sbatch), hg38/hg19 references
# Usage: Run individual steps as modules (e.g. bash lupus_joint_call_pipeline.sh call_controls)
# ============================================================================

set -euo pipefail
IFS=$'\n\t'

# ----------------------------------------------------------------------------
# CONFIGURATION
# ----------------------------------------------------------------------------

PROJECT_DIR="/mnt/atgc-d3/fsanchezq/earrieta/lupus"
CONTROL_LIST_DIR="/data/Project_Eduardo/lupus"
CONTROL_SCRIPT="/data/Project_Eduardo/scripts/BP_GATK_HaploCall_gVCFs_updatedSites.sh"
LUPUS_BAM_DIR="/mnt/Citosina/amedina/lupus/WGS/bam_files"
LUPUS_RG_DIR="$LUPUS_BAM_DIR/RG_bamfiles"
SCRIPTS_DIR="/home/earrieta/earrieta/scripts"
ASSEMBLY_REF="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38.fa"
CHAIN_FILE="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38ToHg19.over.chain"
HG19_REF="/mnt/atgc-d3/fsanchezq/earrieta/references/hg19.fa"

CHROMS=$(seq 1 22)
THREADS=4

log() { echo "[$(date +'%F %T')] $*"; }

# ----------------------------------------------------------------------------
# STEP 1 - gVCF Calling for Control Samples (MXL)
# ----------------------------------------------------------------------------

call_controls() {
  log "Calling gVCFs for control (MXL) samples"

  for CHROM in $CHROMS; do
    for LIST in MXL_path_bams__00 MXL_path_bams__01 MXL_path_bams__02 MXL_path_bams__03; do
      while read -r BAM; do
        sample=$(basename "${BAM}" .bam)
        sample=${sample%%.*}

        log "Submitting control gVCF call for ${sample}, chr${CHROM}"
        sh "${CONTROL_SCRIPT}" "$BAM" "$CHROM" "$CONTROL_LIST_DIR" &
      done < "${CONTROL_LIST_DIR}/${LIST}"
    done
    wait

    log "Archiving gVCFs for chr${CHROM}"
    cd "${CONTROL_LIST_DIR}"
    tar -czf "gVCFs_chr${CHROM}.tgz" gVCFs/*_chr${CHROM}.g.vcf.gz
    rm -f gVCFs/*_chr${CHROM}.g.vcf.gz*
  done
}

# ----------------------------------------------------------------------------
# STEP 2 - Add Read Groups to Lupus Patient BAMs
# ----------------------------------------------------------------------------

prepare_lupus_bams() {
  log "Preparing lupus patient BAMs with read groups"

  list="${PROJECT_DIR}/lupus_patients_bams.list"
  find "${LUPUS_BAM_DIR}" -name "*.sort.rmdup.qc.bam" > "$list"

  while read -r BAM; do
    sample=$(basename "${BAM}" .sort.rmdup.qc.bam)
    qsub -cwd -N "${sample}_RG" -o "${sample}_RG.output" -e "${sample}_RG.error" \
         "${PROJECT_DIR}/../scripts/add_ReadGroups.sh" \
         "$BAM" "$sample" "$LUPUS_RG_DIR"
  done < "$list"
}

# ----------------------------------------------------------------------------
# STEP 3 - gVCF Calling for Lupus Samples
# ----------------------------------------------------------------------------

call_lupus_gvcfs() {
  log "Calling gVCFs for lupus patient samples"

  list="${PROJECT_DIR}/lupus_patients_RG_bams.list"
  find "${LUPUS_RG_DIR}" -name "*.bam" > "$list"

  for CHROM in $CHROMS; do
    mkdir -p "${PROJECT_DIR}/here_chr${CHROM}"

    while read -r BAM; do
      sample=$(basename "${BAM}" .sort.rmdup.qc.RG.bam)
      qsub -cwd -N "${sample}_chr${CHROM}_gVCF" \
           -o "${sample}_chr${CHROM}_gVCF.output" \
           -e "${sample}_chr${CHROM}_gVCF.error" \
           -l virtual_free=30G,h_vmem=30G \
           "${SCRIPTS_DIR}/call_gVCFs.sh" \
           "$BAM" "$ASSEMBLY_REF" "chr${CHROM}" "${PROJECT_DIR}/here_chr${CHROM}"
    done < "$list"
  done
}

# ----------------------------------------------------------------------------
# STEP 4 - Generate Sample Map Files for Joint Calling
# ----------------------------------------------------------------------------

generate_sample_maps() {
  log "Generating sample maps for each chromosome"

  for CHROM in $CHROMS; do
    outmap="${PROJECT_DIR}/here_chr${CHROM}/cohort_lupus_patients_chr${CHROM}.sample_map"
    > "$outmap"

    for VCF in "${PROJECT_DIR}/here_chr${CHROM}/gVCFs/"*.vcf.gz; do
      sample=$(basename "$VCF" .sort.rmdup.qc.RG.g.vcf.gz)
      echo -e "${sample}\t${VCF}" >> "$outmap"
    done
  done
}

# ----------------------------------------------------------------------------
# STEP 5 - Joint Calling (GATK GenomicsDB + GenotypeGVCFs)
# ----------------------------------------------------------------------------

joint_calling() {
  log "Performing joint calling for all chromosomes"

  for CHROM in $CHROMS; do
    qsub -cwd -N "jointCall_chr${CHROM}" \
         -o "jointCall_chr${CHROM}.output" \
         -e "jointCall_chr${CHROM}.error" \
         -l virtual_free=10G \
         "${SCRIPTS_DIR}/jointCall.sh" \
         "$ASSEMBLY_REF" "chr${CHROM}" \
         "${PROJECT_DIR}/here_chr${CHROM}" \
         "${PROJECT_DIR}/here_chr${CHROM}/cohort_lupus_patients_chr${CHROM}.sample_map" \
         "lupus_patients_jointCall_chr${CHROM}"
  done
}

# ----------------------------------------------------------------------------
# STEP 6 - Variant Quality Score Recalibration (VQSR)
# ----------------------------------------------------------------------------

run_vqsr() {
  log "Running VQSR for joint-called VCFs"

  for CHROM in $CHROMS; do
    qsub -cwd -N "VQSR_chr${CHROM}" \
         -o "VQSR_chr${CHROM}.output" \
         -e "VQSR_chr${CHROM}.error" \
         -l virtual_free=6G \
         "${SCRIPTS_DIR}/VQSR.sh" \
         "${PROJECT_DIR}/here_chr${CHROM}/lupus_patients_jointCall_chr${CHROM}.vcf.gz" \
         "${PROJECT_DIR}/here_chr${CHROM}" \
         "chr${CHROM}"
  done
}

# ----------------------------------------------------------------------------
# STEP 7 - LiftOver to hg19 coordinates
# ----------------------------------------------------------------------------

run_liftover() {
  log "Running LiftOver to hg19"

  for CHROM in $CHROMS; do
    qsub -cwd -N "LO_chr${CHROM}" \
         -o "LO_chr${CHROM}.output" \
         -e "LO_chr${CHROM}.error" \
         -l virtual_free=15G \
         "${SCRIPTS_DIR}/LiftOver.sh" \
         "${PROJECT_DIR}/here_chr${CHROM}/VQSR/lupus_patients_jointCall_chr${CHROM}.vqsr.vcf.gz" \
         "${PROJECT_DIR}/here_chr${CHROM}" \
         "$CHAIN_FILE" "$HG19_REF"
  done
}

# ----------------------------------------------------------------------------
# STEP 8 - Concatenate All VCFs Across Chromosomes
# ----------------------------------------------------------------------------

concat_vcfs() {
  log "Concatenating all VQSR VCFs into one multisample VCF"

  listfile="${PROJECT_DIR}/VCFs_order_hg38.list"
  > "$listfile"

  for CHROM in $CHROMS; do
    echo "${PROJECT_DIR}/here_chr${CHROM}/VQSR/lupus_patients_jointCall_chr${CHROM}.vqsr.vcf.gz" >> "$listfile"
  done

  bcftools concat -a -f "$listfile" -Oz -o "${PROJECT_DIR}/lupus_patients_all-chrs.vcf.gz"
  log "Final concatenated VCF written to lupus_patients_all-chrs.vcf.gz"
}

# ----------------------------------------------------------------------------
# MAIN - Dispatch Function Calls
# ----------------------------------------------------------------------------

main() {
  case "${1:-}" in
    call_controls) call_controls ;;
    prepare_lupus_bams) prepare_lupus_bams ;;
    call_lupus_gvcfs) call_lupus_gvcfs ;;
    generate_sample_maps) generate_sample_maps ;;
    joint_calling) joint_calling ;;
    run_vqsr) run_vqsr ;;
    run_liftover) run_liftover ;;
    concat_vcfs) concat_vcfs ;;
    *) 
      echo "Usage: $0 {call_controls|prepare_lupus_bams|call_lupus_gvcfs|generate_sample_maps|joint_calling|run_vqsr|run_liftover|concat_vcfs}"
      exit 1
      ;;
  esac
}

main "$@"
