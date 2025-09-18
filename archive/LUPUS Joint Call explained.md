# LUPUS Joint Call explained

---

## Context

This script is part of a genomic pipeline focused on **joint variant calling** using GATK's HaplotypeCaller + GenotypeGVCFs approach — a standard for large cohort variant analysis, here applied to a **Lupus study** with multiple samples from patients and controls.

---

## Step-by-Step Breakdown

### ### Step 1: Generate gVCFs for control samples

```bash
bash /data/Project_Eduardo/scripts/BP_GATK_HaploCall_gVCFs_updatedSites.sh
```

**What it does:**
Runs a pre-existing script to generate per-sample **gVCF files** for the *control cohort* using **GATK HaplotypeCaller**.

**Why:**
The GATK best-practice pipeline separates variant calling into two stages:

* First, it emits per-sample *genomic VCFs (gVCFs)* that contain variant and non-variant information, preserving uncertainty.
* Later, these are jointly genotyped across samples.

**Expected output:**
A folder or set of files containing `.g.vcf.gz` files — one per control sample — typically indexed with `.g.vcf.gz.tbi`.

> *Note:* This step is modularized into another script — a good practice if you're running control and case samples independently.

---

### ### Step 2: Add read groups to case (Lupus) BAMs

```bash
for i in *.bam
do
  sample=$(basename "$i" .bam)
  gatk AddOrReplaceReadGroups \
    -I "$i" \
    -O "$sample"_RG.bam \
    -RGID "$sample" \
    -RGLB lib1 \
    -RGPL ILLUMINA \
    -RGPU unit1 \
    -RGSM "$sample"
done
```

**What it does:**
Iterates over all `.bam` files (Lupus samples), adding or standardizing **Read Group (RG)** information using GATK’s `AddOrReplaceReadGroups`.

**Why:**
Read Groups are **mandatory for GATK** tools. They label each BAM with metadata (e.g. sample name, platform, library) needed to distinguish technical vs. biological replicates, especially during joint calling.

**Expected output:**
A new BAM file for each sample with `_RG.bam` suffix, containing updated metadata.

> **Improvement:** You could also index these files afterward using `samtools index`.

---

### ### Step 3: Call gVCFs on Lupus (case) samples

```bash
for i in *_RG.bam
do
  sample=$(basename "$i" _RG.bam)
  gatk HaplotypeCaller \
    -R "$ref" \
    -I "$i" \
    -ERC GVCF \
    -O "$sample".g.vcf.gz
done
```

**What it does:**
Runs **HaplotypeCaller** in GVCF mode (`-ERC GVCF`) for each Lupus BAM file.

**Why:**
This produces one `.g.vcf.gz` per sample — storing not only variant sites but also confident non-variant sites, enabling **joint genotyping** in the next step.

**Expected output:**
Compressed `.g.vcf.gz` files for all Lupus samples.

> *Note:* It’s common to run HaplotypeCaller in parallel across samples via Slurm job arrays for efficiency.

---

### ### Step 4: Prepare the sample map

```bash
for i in *.g.vcf.gz
do
  echo -e "$(basename "$i" .g.vcf.gz)\t$PWD/$i" >> cohort.sample_map
done
```

**What it does:**
Generates a tab-delimited file required by GATK’s `GenomicsDBImport`. Each line maps a sample ID to the full path of its gVCF.

**Why:**
GATK tools like `GenomicsDBImport` need to know:

* The sample name (must match `-RGSM` from earlier)
* The path to the gVCF file

**Expected output:**
A `cohort.sample_map` file like:

```
sample1    /path/to/sample1.g.vcf.gz
sample2    /path/to/sample2.g.vcf.gz
...
```

> 🔎 *Potential improvement:* Sort the file if sample ordering is critical for downstream consistency or comparison.

---

### ### Step 5: Create chromosome list

```bash
for chr in {1..22} X Y
do
  echo -e "chr$chr" >> chr.list
done
```

**What it does:**
Generates a plain text file listing chromosomes `chr1` to `chr22`, then `chrX` and `chrY`.

**Why:**
This file is used to split work per chromosome — a common optimization for joint calling since chromosomes are independent for most variant analyses.

**Expected output:**
A `chr.list` text file with:

```
chr1
chr2
...
chr22
chrX
chrY
```

> *Gotcha:* This assumes your reference FASTA uses `chr` prefixes. If your reference is just `1, 2, ..., X`, this will cause mismatches.

---

### ### Step 6: Run GenomicsDBImport and GenotypeGVCFs per chromosome

```bash
for chr in $(cat chr.list)
do
  gatk GenomicsDBImport \
    --genomicsdb-workspace-path "$chr"_workspace \
    --intervals "$chr" \
    --sample-name-map cohort.sample_map \
    --reader-threads 4 \
    --batch-size 50

  gatk GenotypeGVCFs \
    -R "$ref" \
    -V gendb://"$chr"_workspace \
    -O "$chr"_cohort.vcf.gz
done
```

**What it does:**

1. **`GenomicsDBImport`**:
   
   * Merges gVCFs into a GenomicsDB workspace per chromosome.
   * Batching and multithreading are enabled to reduce memory/time usage.

2. **`GenotypeGVCFs`**:
   
   * Converts the gVCF data into **jointly genotyped** VCFs — one per chromosome.

**Why:**
Joint calling ensures that **all samples are genotyped together**, improving variant quality and enabling cohort-wide filtering.

**Expected output:**

* `chrN_workspace/` directories for each chromosome.
* Joint-called, compressed VCFs: `chr1_cohort.vcf.gz`, etc.

> *Suggestion:* Consider adding `--tmp-dir` and `--java-options` if running on shared systems to manage disk and memory.

---

## Summary Table

| Step | Action                    | Tool                                  | Purpose                           | Output               |
| ---- | ------------------------- | ------------------------------------- | --------------------------------- | -------------------- |
| 1    | Run control gVCF pipeline | GATK HaplotypeCaller                  | Prepare control gVCFs             | `.g.vcf.gz` files    |
| 2    | Add read groups           | GATK AddOrReplaceReadGroups           | Metadata for joint calling        | `_RG.bam`            |
| 3    | Call gVCFs on Lupus cases | GATK HaplotypeCaller                  | Case gVCFs                        | `.g.vcf.gz`          |
| 4    | Create sample map         | bash                                  | Needed for GenomicsDBImport       | `cohort.sample_map`  |
| 5    | Chromosome list           | bash                                  | Enables per-chromosome processing | `chr.list`           |
| 6    | Joint calling             | GATK GenomicsDBImport + GenotypeGVCFs | Final VCFs per chromosome         | `chrX_cohort.vcf.gz` |

---

## Final Thoughts and Potential Enhancements

* **Parallelization**: Split Step 6 across Slurm jobs per chromosome using job arrays.
* **Compression and indexing**: Make sure all VCFs are bgzipped and tabix-indexed before downstream use.
* **QC metrics**: Consider adding a step for `VariantEval`, `VQSR`, or hard filtering.
* **Logging**: Implement `tee` or `>> log.txt` on key commands to track progress.

---