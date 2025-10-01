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
3. **Step 03 — Joint Variant Calling** *(in development)*
  
  - Input: QC’d BAMs from Step 02
  - Process: call SNPs/indels, generate GVCFs, perform joint genotyping, and filter variants.
  - Output: multi-sample VCF/BCF ready for downstream analysis.
4. **Wrapper Script** *(in development)*
  
  - Runs steps 01 → 03 in tandem, managing intermediate outputs and sanity checks.

---

## Step 02: GATK4 BAM QC [FENIX]

This module prepares mapped BAMs for downstream variant calling, ensuring high-quality alignments and base quality scores.

### Usage

```bash
gatk_mapped_bam_qc.sh [INPUT_BAM] [OUTPUT_PATH]
```

### Input

- `INPUT_BAM` : BAM file aligned to reference genome (output of Step 01).
- `OUTPUT_PATH` : destination folder for QC outputs (default: current working directory).

### Output

The script produces:

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

## Example Workflow

```bash
# Step 01 - Align FASTQs (placeholder)
bash align_fastq_reads.sh sample_R1.fastq.gz sample_R2.fastq.gz ./results/align/

# Step 02 - QC mapped BAM
bash gatk_mapped_bam_qc.sh ./results/align/sample.bam ./results/qc/

# Step 03 - Joint variant calling (placeholder)
bash joint_variant_calling.sh ./results/qc/ ./results/variants/

# Wrapper - Full pipeline run
bash run_joint_calling_pipeline.sh sample_config.yaml
```

---

## Development Notes

- **Consistency:** Each script follows the same input-output convention and sanity checks.
- **Extensibility:** Steps are modular; users can run only the portions relevant to their analysis.
- **Portability:** Scripts detect whether they are running on a local machine or a remote HPC environment and adapt accordingly.

---

## Authors

- Pavel Salazar-Fernandez (current maintainer): epsalazarf@gmial.com
- Dr. Federico Sanchez-Quinto (project leader)
