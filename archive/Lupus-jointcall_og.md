# LUPUS Joint Call

```bash
ssh -p 37533 earrieta@godzilla.inmegen.gob.mx
```

```bash
# VPN a UNAM
ssh -Y earrieta@dna.lavis.unam.mx
```

Path to work

```bash
cd /mnt/atgc-d3/fsanchezq/earrieta/lupus
```

### gVCFs

### MXL

```bash
for CHROM in
do
  for list in MXL_path_bams__02 MXL_path_bams__03 MXL_path_bams__00 MXL_path_bams__01
  do
    for BAM in `cat /data/Project_Eduardo/lupus/${list}`
    do
      sample=`basename ${BAM} .bam`
      sample=${sample%%.*}

      sh /data/Project_Eduardo/scripts/BP_GATK_HaploCall_gVCFs_updatedSites.sh ${BAM} \
                                                                               ${CHROM} \
                                                                               /data/Project_Eduardo/lupus \
                                                                               &
    done
  done
done
```

```bash

cd /data/Project_Eduardo/lupus
tar -cvzf gVCFs_chr${CHROM}.tgz gVCFs/*_chr${CHROM}.g.vcf.gz
rm /data/Project_Eduardo/lupus/gVCFs/*_chr${CHROM}.g.vcf.gz*
```

### Lupus patients

```bash
ls /mnt/Citosina/amedina/lupus/WGS/bam_files/*.sort.rmdup.qc.bam > /mnt/atgc-d3/fsanchezq/earrieta/lupus/lupus_patients_bams.list
```

```bash
for BAM in `cat /mnt/atgc-d3/fsanchezq/earrieta/lupus/lupus_patients_bams.list`
do
  sample=`basename ${BAM} .sort.rmdup.qc.bam`
  
  qsub -cwd -N ${sample}_RG -o ${sample}_RG.output -e ${sample}_RG.error \
         /mnt/atgc-d3/fsanchezq/earrieta/scripts/add_ReadGroups.sh \
         ${BAM} \
         ${sample} \
         /mnt/Citosina/amedina/lupus/WGS/bam_files/RG_bamfiles
done
```

- Faltantes
    
    ```bash
    for BAM in `ls /mnt/Citosina/amedina/lupus/WGS/bam_files/Q096.sort.rmdup.qc.bam /mnt/Citosina/amedina/lupus/WGS/bam_files/Q047.sort.rmdup.qc.bam /mnt/Citosina/amedina/lupus/WGS/bam_files/Q048.sort.rmdup.qc.bam`
    do
      sample=`basename ${BAM} .sort.rmdup.qc.bam`
      
      qsub -cwd -N ${sample}_RG -o ${sample}_RG.output -e ${sample}_RG.error \
             /mnt/atgc-d3/fsanchezq/earrieta/scripts/add_ReadGroups.sh \
             ${BAM} \
             ${sample} \
             /mnt/Citosina/amedina/lupus/WGS/bam_files/RG_bamfiles
    done
    ```
    

```bash
ls /mnt/Citosina/amedina/lupus/WGS/bam_files/RG_bamfiles/*.bam > /mnt/atgc-d3/fsanchezq/earrieta/lupus/lupus_patients_RG_bams.list
```

```bash
cd /mnt/atgc-d3/fsanchezq/earrieta/lupus

ASSEMBLY="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38.fa"

for CHROM in 15
do
  mkdir /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}
  for BAM in `cat /mnt/atgc-d3/fsanchezq/earrieta/lupus/lupus_patients_RG_bams.list`
  do
    sample=`basename ${BAM} .sort.rmdup.qc.RG.bam`

    qsub -cwd -N ${sample}_chr${CHROM}_gVCF -o ${sample}_chr${CHROM}_gVCF.output -e ${sample}_chr${CHROM}_gVCF.error -l virtual_free=30G,h_vmem=30g \
         /home/earrieta/earrieta/scripts/call_gVCFs.sh \
         ${BAM} \
         ${ASSEMBLY} \
         chr${CHROM} \
         /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}

  done
done
```

- Faltantes
    
    ```bash
    cd /mnt/atgc-d3/fsanchezq/earrieta/lupus
    path=/mnt/Citosina/amedina/lupus/WGS/bam_files/RG_bamfiles
    
    ASSEMBLY="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38.fa"
    
    for CHROM in 15
    do
      for BAM in `ls ${path}/Q096.sort.rmdup.qc.RG.bam ${path}/Q047.sort.rmdup.qc.RG.bam ${path}/Q048.sort.rmdup.qc.RG.bam`
      do
        sample=`basename ${BAM} .sort.rmdup.qc.RG.bam`
    
        qsub -cwd -N ${sample}_chr${CHROM}_gVCF -o ${sample}_chr${CHROM}_gVCF.output -e ${sample}_chr${CHROM}_gVCF.error -l virtual_free=40G,h_vmem=40g \
             /home/earrieta/earrieta/scripts/call_gVCFs.sh \
             ${BAM} \
             ${ASSEMBLY} \
             chr${CHROM} \
             /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}
    
      done
    done
    ```
    

### JointCall

- Cohort
    
    ```bash
    for CHROM in 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
    do
      for VCF in `ls /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}/gVCFs/*.vcf.gz`
      do
        sample=`basename ${VCF} .sort.rmdup.qc.RG.g.vcf.gz`
        echo -e "${sample}\t${VCF}" >> /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}/cohort_lupus_patients_chr${CHROM}.sample_map
      done
    done
    ```
    

```bash
ASSEMBLY="/mnt/atgc-d3/fsanchezq/earrieta/references/hg38.fa"

for CHROM in 15 9 8
do
  qsub -cwd -N jointCall_chr${CHROM} -o jointCall_chr${CHROM}.output -e jointCall_chr${CHROM}.error -l virtual_free=10G \
       /home/earrieta/earrieta/scripts/jointCall.sh \
       ${ASSEMBLY} \
       chr${CHROM} \
       /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM} \
       /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}/cohort_lupus_patients_chr${CHROM}.sample_map \
       lupus_patients_jointCall_chr${CHROM}

done
```

### VQSR

```bash
for CHROM in 14 8
do
  qsub -cwd -N VQSR_chr${CHROM} -o VQSR_chr${CHROM}.output -e VQSR_chr${CHROM}.error -l virtual_free=6G \
      /home/earrieta/earrieta/scripts/VQSR.sh \
      /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}/lupus_patients_jointCall_chr${CHROM}.vcf.gz \
      /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM} \
      chr${CHROM}
done
```

### Lift Over

```bash
for CHROM in 8 14
do
  qsub -cwd -N LO_chr${CHROM} -o LO_chr${CHROM}.output -e LO_chr${CHROM}.error -l virtual_free=15G \
       /home/earrieta/earrieta/scripts/LiftOver.sh \
       /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}/VQSR/lupus_patients_jointCall_chr${CHROM}.vqsr.vcf.gz \
       /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM} \
       /mnt/atgc-d3/fsanchezq/earrieta/references/hg38ToHg19.over.chain \
       /mnt/atgc-d3/fsanchezq/earrieta/references/hg19.fa
done
```

### Concat

```bash
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
  ls /mnt/atgc-d3/fsanchezq/earrieta/lupus/here_chr${CHROM}/VQSR/lupus_patients_jointCall_chr${CHROM}.vqsr.vcf.gz >> \
     /mnt/atgc-d3/fsanchezq/earrieta/lupus/VCFs_order_hg38.list
done
```

```bash
bcftools concat -a -f /mnt/atgc-d3/fsanchezq/earrieta/lupus/VCFs_order.list -Oz -o /mnt/atgc-d3/fsanchezq/earrieta/lupus/lupus_patientes_all-chrs.vcf.gz
```