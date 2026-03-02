# Joint Variant Calling Pipeline [FENIX]

This repository contains a modular workflow for performing **joint variant discovery** following [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).  
The pipeline is organized into **three sequential steps**, each encapsulated in its own script, plus a **wrapper** to run the full process end-to-end.

> **Note:** The scripts are under active development. Logic, style, and usage will converge as the modules mature.

---

## Pipeline Overview

1. **Step 01 — FASTQ Alignment** *(in development)*
  - Input: raw FASTQ reads
  - Process: align reads to a reference genome (using `bwa mem`), produce raw BAM/CRAM files, generate alignment statistics.
  - Output: coordinate-sorted, indexed BAM files.

2. **Step 02 — Mapped BAM QC (this script)**  
  - Input: mapped BAM from Step 01
  - Process: perform read group assignment, duplicate removal, mapping quality filtering, base quality recalibration, and coverage estimation.
  - Output: quality-controlled, analysis-ready BAM files.

3. **Step 03 — GVCF generation** *(in testing)*  
  - Input: QC’d BAMs from Step 02
  - Process: call SNPs/indels, generate GVCFs, splits by chromosome.
  - Output: single-sample VCF ready for merge with other samples.

4. **Step 04 — Joint Variant Calling** *(in development)*
  - Process: perform joint genotyping, and filter variants.

---

## Step 02: GATK4 BAM QC [FENIX]

This module prepares mapped BAMs for downstream variant calling, ensuring high-quality alignments and base quality scores.

### Usage

```bash
02_gatk_bam_qc_workflow.sh [INPUT BAM] [OUTPUT PATH] [RG STRING]
```

### Input

- `INPUT_BAM` : BAM file aligned to reference genome (output of Step 01).
- `OUTPUT_PATH` : destination folder for QC outputs (default: current working directory).

### Output

- `${BAM}.rg.bam` — BAM with read groups assigned
- `${BAM}.rmdup.bam` — duplicate-removed BAM
- `${BAM}.rmdup.mqfilt.bam` — mapping-quality filtered BAM
- `${BAM}.rmdup.mqfilt.bqsr.bam` — recalibrated BAM ready for variant calling
- Coverage and quality metrics (`.txt`, `.pdf`, `.summary.txt`)

### Quality Control Steps

1. **Assign Read Groups**  
  GATK `AddOrReplaceReadGroups` for downstream compatibility.
  
2. **Remove Duplicates**  
  GATK `MarkDuplicatesSpark` to discard PCR duplicates.
  
3. **Mapping Quality Filtering**  
  `samtools view` with mapping quality threshold (`-q 30`).
  
4. **Base Quality Recalibration (BQSR)**  
  GATK `BaseRecalibrator` and `ApplyBQSR` using known variant sites.
  
5. **Coverage Calculation**  
  `mosdepth` summarizing per-sample coverage.
  
6. **Optional Metrics**  
  GATK `CollectAlignmentSummaryMetrics` and `CollectInsertSizeMetrics`.

---

## Step 03: GATK4 HaplotypeCaller — Per-Sample GVCF Generation [FENIX]

This module performs per-sample variant calling using GATK4 HaplotypeCaller in GVCF mode and prepares chromosome-level subsets for downstream joint genotyping.

### Usage

```bash
03_gatk_haplotype_caller.sh [MAPPED_BAM] [OUTPUT_PATH]
```

### Input

* `MAPPED_BAM` : Analysis-ready BAM file (BQSR-processed; output of Step 02).
* `OUTPUT_PATH` : Destination directory for GVCF outputs (default: current working directory).
* `config/config.yaml` : Environment-specific configuration file containing:
  * `ref_gnm` — reference genome FASTA
  * `ref_vars` — known variant sites (e.g., dbSNP)

### Output

* `${BAM}.raw_variants.g.vcf.gz` — raw per-sample GVCF
* `${BAM}.raw_variants.canon_chr.g.vcf.gz` — canonical-chromosome GVCF (chr1–22, X, Y, M)
* `chrom_gvcf/${BAM}.raw_vars.<CHR>.g.vcf.gz` — per-chromosome GVCFs

### Processing Steps

1. **HaplotypeCaller (GVCF mode)**
   Produces a single-sample GVCF suitable for joint genotyping.

2. **Extract Canonical Chromosomes**
   Subsets to: chr1–chr22, chrX, chrY, chrM using `bcftools view`, producing a cleaned GVCF focused on primary chromosomes.

3. **Split by Chromosome**
   Automatically detects chromosomes from the GVCF index and generates one GVCF per chromosome inside a subfolder.

---

## Dependencies

The pipeline expects the following tools available in `$PATH` (modules will be loaded automatically if run on the FENIX remote HPC system):

- `gatk` (v4.x)
- `samtools` (>=1.15)
- `mosdepth` (>=0.3.x)
- `fastqc` (>=0.12.1)
- `R` (for optional plots)

Reference data required (configured via `config/config.yaml`):

- Reference genome FASTA (`ref_gnm`) with index
- Known variant sites VCF (`ref_vars`)

---

## Development Notes

- **Consistency:** Each script follows the same input-output convention and sanity checks.
- **Extensibility:** Steps are modular; users can run only the portions relevant to their analysis.
- **Portability:** Scripts detect whether they are running on a local machine or a remote HPC environment and adapt accordingly.

---

## Authors

- Pavel Salazar-Fernandez (current maintainer): epsalazarf@gmail.com
- Dr. Federico Sanchez-Quinto (project leader)
