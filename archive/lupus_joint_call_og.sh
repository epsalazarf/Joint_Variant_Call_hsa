#!/bin/bash

# ============================================================================
# Script: lupus_joint_call_pipeline.sh
# Purpose: Joint variant calling pipeline for lupus and control cohorts
# Author: Adapted from raw commands by colleague
# ============================================================================

# --- CONFIGURATION ---

# Remote servers
UNAM_SSH_PORT=37533
USER=earrieta
GODZILLA=godzilla.inmegen.gob.mx
DNA_UNAM=dna.lavis.unam.mx

# Local directories
PROJECT_DIR="/mnt/atgc-d3/fsanchezq/earrieta/lupus"
CONTROL_LIST_DIR="/data/Project_Eduardo/lupus"
CONTROL_SCRIPT="/data/Project_Eduardo/scripts/BP_GATK_HaploCall_gVCFs_updatedSites.sh"
LUPUS_BAM_DIR="/mnt/Citosina/amedina/lupus/WGS/bam_files"
LUPUS_RG_DIR="$LUPUS_BAM_DIR/RG_bamfiles"
ASSEMBLY_REF="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38.fa"
CHAIN_FILE="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38ToHg19.over.chain"
HG19_REF="/mnt/atgc-d3/fsanchezq/earrieta/references/hg19.fa"

# --- CONNECTION STEPS (Manual) ---
# Step 1: Connect to INMEGEN
# ssh -p $UNAM_SSH_PORT $USER@$GODZILLA

# Step 2: VPN to UNAM
# ssh -Y $USER@$DNA_UNAM

# --- STEP 1: Call gVCFs for MXL control cohort ---

call_controls() {
  for CHROM in {1..22}; do
    for LIST in MXL_path_bams__00 MXL_path_bams__01 MXL_path_bams__02 MXL_path_bams__03; do
      while read -r BAM; do
        sample=$(basename "$BAM" .bam)
        sample=${sample%%.*}

        echo "Calling gVCF for $sample on chr$CHROM"
        sh $CONTROL_SCRIPT "$BAM" "$CHROM" "$CONTROL_LIST_DIR" &

      done < "$CONTROL_LIST_DIR/$LIST"
    done
    wait

    echo "Compressing and cleaning gVCFs for chr$CHROM"
    cd "$CONTROL_LIST_DIR"
    tar -cvzf "gVCFs_chr${CHROM}.tgz" gVCFs/*_chr${CHROM}.g.vcf.gz
    rm -f gVCFs/*_chr${CHROM}.g.vcf.gz*
  done
}

# --- STEP 2: Add read groups to Lupus BAM files ---

prepare_lupus_bams() {
  list_path="$PROJECT_DIR/lupus_patients_bams.list"
  find "$LUPUS_BAM_DIR" -name "*.sort.rmdup.qc.bam" > "$list_path"

  while read -r BAM; do
    sample=$(basename "$BAM" .sort.rmdup.qc.bam)
    qsub -cwd -N "${sample}_RG" \
         -o "${sample}_RG.output" \
         -e "${sample}_RG.error" \
         "$PROJECT_DIR/../scripts/add_ReadGroups.sh" \
         "$BAM" "$sample" "$LUPUS_RG_DIR"
  done < "$list_path"

  # Additional missing BAMs (manual override)
  for BAM in "$LUPUS_BAM_DIR"/Q096*.bam "$LUPUS_BAM_DIR"/Q047*.bam "$LUPUS_BAM_DIR"/Q048*.bam; do
    sample=$(basename "$BAM" .sort.rmdup.qc.bam)
    qsub -cwd -N "${sample}_RG" \
         -o "${sample}_RG.output" \
         -e "${sample}_RG.error" \
         "$PROJECT_DIR/../scripts/add_ReadGroups.sh" \
         "$BAM" "$sample" "$LUPUS_RG_DIR"
  done
}

# --- STEP 3: Call gVCFs for Lupus cohort ---

call_lupus_gvcfs() {
  list_path="$PROJECT_DIR/lupus_patients_RG_bams.list"
  find "$LUPUS_RG_DIR" -name "*.bam" > "$list_path"

  for CHROM in 15; do
    mkdir -p "$PROJECT_DIR/here_chr${CHROM}"

    while read -r BAM; do
      sample=$(basename "$BAM" .sort.rmdup.qc.RG.bam)
      qsub -cwd -N "${sample}_chr${CHROM}_gVCF" \
           -o "${sample}_chr${CHROM}_gVCF.output" \
           -e "${sample}_chr${CHROM}_gVCF.error" \
           -l virtual_free=30G,h_vmem=30G \
           "$HOME/earrieta/scripts/call_gVCFs.sh" \
           "$BAM" "$ASSEMBLY_REF" "chr${CHROM}" \
           "$PROJECT_DIR/here_chr${CHROM}"
    done < "$list_path"
  done
}

# --- STEP 4: Create sample maps for Joint Call ---

generate_sample_maps() {
  for CHROM in {1..22}; do
    sample_map="$PROJECT_DIR/here_chr${CHROM}/cohort_lupus_patients_chr${CHROM}.sample_map"
    mkdir -p "$(dirname "$sample_map")"
    > "$sample_map"
    
    for VCF in "$PROJECT_DIR/here_chr${CHROM}/gVCFs/"*.vcf.gz; do
      sample=$(basename "$VCF" .sort.rmdup.qc.RG.g.vcf.gz)
      echo -e "${sample}\t${VCF}" >> "$sample_map"
    done
  done
}

# --- STEP 5: Perform Joint Calling ---

joint_calling() {
  for CHROM in 15 9 8; do
    qsub -cwd -N "jointCall_chr${CHROM}" \
         -o "jointCall_chr${CHROM}.output" \
         -e "jointCall_chr${CHROM}.error" \
         -l virtual_free=10G \
         "$HOME/earrieta/scripts/jointCall.sh" \
         "$ASSEMBLY_REF" \
         "chr${CHROM}" \
         "$PROJECT_DIR/here_chr${CHROM}" \
         "$PROJECT_DIR/here_chr${CHROM}/cohort_lupus_patients_chr${CHROM}.sample_map" \
         "lupus_patients_jointCall_chr${CHROM}"
  done
}

# --- STEP 6: Variant Quality Score Recalibration (VQSR) ---

run_vqsr() {
  for CHROM in 14 8; do
    qsub -cwd -N "VQSR_chr${CHROM}" \
         -o "VQSR_chr${CHROM}.output" \
         -e "VQSR_chr${CHROM}.error" \
         -l virtual_free=6G \
         "$HOME/earrieta/scripts/VQSR.sh" \
         "$PROJECT_DIR/here_chr${CHROM}/lupus_patients_jointCall_chr${CHROM}.vcf.gz" \
         "$PROJECT_DIR/here_chr${CHROM}" \
         "chr${CHROM}"
  done
}

# --- STEP 7: LiftOver to hg19 ---

run_liftover() {
  for CHROM in 8 14; do
    qsub -cwd -N "LO_chr${CHROM}" \
         -o "LO_chr${CHROM}.output" \
         -e "LO_chr${CHROM}.error" \
         -l virtual_free=15G \
         "$HOME/earrieta/scripts/LiftOver.sh" \
         "$PROJECT_DIR/here_chr${CHROM}/VQSR/lupus_patients_jointCall_chr${CHROM}.vqsr.vcf.gz" \
         "$PROJECT_DIR/here_chr${CHROM}" \
         "$CHAIN_FILE" "$HG19_REF"
  done
}

# --- STEP 8: Concatenate Final VCFs ---

concat_vcfs() {
  vcf_list="$PROJECT_DIR/VCFs_order_hg38.list"
  > "$vcf_list"
  for CHROM in {1..22}; do
    echo "$PROJECT_DIR/here_chr${CHROM}/VQSR/lupus_patients_jointCall_chr${CHROM}.vqsr.vcf.gz" >> "$vcf_list"
  done

  bcftools concat -a -f "$vcf_list" -Oz -o "$PROJECT_DIR/lupus_patientes_all-chrs.vcf.gz"
}


# --- MAIN EXECUTION ---

# Uncomment and run the steps you need:
# call_controls
# prepare_lupus_bams
# call_lupus_gvcfs
# generate_sample_maps
# joint_calling
# run_vqsr
# run_liftover
# concat_vcfs

