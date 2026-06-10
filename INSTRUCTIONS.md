# Pipeline Instructions

Step-by-step usage guide for the Joint Variant Calling Pipeline.

> For a high-level overview, see [README.md](README.md).

---

## Prerequisites

1. **Configure paths** in `config/config.yaml` — set `ref_gnm`, `ref_vars`, and any other keys for your environment (`remote:` for FENIX, `local:` for workstation).
2. **Load modules** (FENIX only): scripts call `module load` automatically based on the `modules` key in the config.
3. All scripts must be run from the repo root or with explicit paths.

---

## Step 00 — Scan FASTQ Pairs (optional)

Inspect a sample directory before mapping to verify FASTQ pairs are correctly detected.

```bash
bash bin/00_scan_fastq_pairs.sh [input_dir]
```

Defaults to `$PWD`. Prints detected pairs without running anything.

---

## Step 01 — FASTQ Alignment

Maps raw reads to hg38, embeds read groups, and produces coordinate-sorted indexed BAM files.

### Single sample

```bash
# Run from the directory containing the sample's FASTQ files
bash bin/01_bwa_map_fastq_reads.sh [input_dir] [output_path]
```

Both arguments default to `$PWD`.

### Batch (SLURM)

```bash
bash bin/BWAMAP.seq_batch-slurmer.sh bin/01_bwa_map_fastq_reads.sh /path/to/samples/
```

Auto-discovers sample subdirectories and submits one job per sample (parallel). Wait for all jobs to finish before running Step 02.

```bash
squeue -u $USER | grep BWAMAP
```

### Supported FASTQ naming conventions

| Pattern | Example |
|---------|---------|
| `{SAMPLE}_{BARCODE}-1A_{PLATE}_{LANE}_{R}.fq.gz` | `L23_CKDN...-1A_227CC2LT4_L5_1.fq.gz` |
| `{SAMPLE}_{R}.fq.gz` | `HL078_1.fq.gz` |
| `{SAMPLE}_R{R}.fastq.gz` | `EGAN00004552350_R1.fastq.gz` |

R0 (index read) files are ignored automatically.

### Output

| File | Description |
|------|-------------|
| `{SAMPLE}.sort.bam` + `.bai` | Sorted, indexed BAM (single-pair samples) |
| `{SAMPLE}_{PLATE_LANE}.sort.bam` + `.bai` | BAM per read group (multi-lane samples) |

---

## Step 02 — BAM QC and Preprocessing

Deduplicates, optionally filters by mapping quality, and applies BQSR to produce analysis-ready BAMs.

### Single sample

```bash
# Standard — single BAM from Step 01
bash bin/02_gatk_bam_qc_workflow.sh SAMPLE.sort.bam [output_path]

# Multiplexed — multiple BAMs from the same sample (merged during deduplication)
bash bin/02_gatk_bam_qc_workflow.sh "SAMPLE_plateA.sort.bam,SAMPLE_plateB.sort.bam" [output_path]

# Legacy — BAM without read groups (RG string triggers step 0)
bash bin/02_gatk_bam_qc_workflow.sh SAMPLE.bam [output_path] "@RG\tID:...\tSM:...\tPL:ILLUMINA\t..."
```

The fourth argument (`add_rg`) controls read group assignment: `auto` (default) | `true` | `false`.

### Batch (SLURM)

```bash
bash bin/BAMQC.seq_batch-slurmer.sh batch-list.txt bin/02_gatk_bam_qc_workflow.sh /path/to/bams/
```

`batch-list.txt` is tab-separated: `sample_name  sample_file  readgroup_string`.  
Each sample is submitted as an independent job (8 CPUs / 32 GB).

```bash
squeue -u $USER | grep BAMQC
```

### Output

| File | Description |
|------|-------------|
| `SAMPLE.rmdup.bam` | Duplicates marked (intermediate) |
| `SAMPLE-dups.txt` | Duplicate metrics |
| `SAMPLE.rmdup.mqfilt.bam` | MQ ≥ 30 filtered (intermediate, if `MQ_FILTER=true`) |
| `SAMPLE.rmdup.mqfilt.bqsr.bam` | **Final analysis-ready BAM** |
| `SAMPLE.rmb.mosdepth.*` | Coverage summary |

### Toggles (edit inside the script)

| Variable | Default | Effect |
|----------|---------|--------|
| `HOUSEKEEP` | `true` | Remove intermediate BAMs on success |
| `MQ_FILTER` | `false` | Apply MQ ≥ 30 filter (originally for aDNA) |
| `RUN_METRICS` | `false` | Collect alignment + insert-size metrics |
| `BQSR_EVAL` | `false` | Run post-BQSR evaluation inline (use 02a instead) |
| `REMOVE_DUPS` | `false` | Remove duplicates instead of marking |

---

## Step 02a — Retroactive BQSR Evaluation (optional)

Generates before/after BQSR covariate plots without re-running the full Step 02. Use when Step 02 was run with `BQSR_EVAL=false`.

**Prerequisite:** `SAMPLE.bqsr_table.txt` must exist in the output directory (produced by Step 02).

```bash
bash bin/02a_bqsr_evaluate.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
```

### Output

| File | Description |
|------|-------------|
| `SAMPLE.bqsr_table_recal.txt` | Post-BQSR recalibration table |
| `SAMPLE.bqsr_covariates.pdf` | Before/after covariate plots |
| `SAMPLE.bqsr_covariates.csv` | Intermediate covariate data |

---

## Step 03a — HaplotypeCaller (standard coverage)

Runs GATK HaplotypeCaller in GVCF mode. Use for samples with normal WGS depth (≥ 15×).

### Single sample

```bash
bash bin/03_gatk_haplotype_caller.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
```

### Batch (SLURM)

```bash
bash bin/HAPCALL.seq_batch-slurmer.sh bin/03_gatk_haplotype_caller.sh /path/to/bqsr/bams/
```

Auto-discovers `*.rmdup.mqfilt.bqsr.bam` files and submits one job per sample (4 CPUs / 32 GB).

```bash
squeue -u $USER | grep HAPCALL
```

### Output

| File | Description |
|------|-------------|
| `*.raw_variants.g.vcf.gz` | Full raw GVCF |
| `*.raw_variants.canon_chr.g.vcf.gz` | Canonical chromosomes only |
| `chrom_gvcf/*.raw_vars.{CHR}.g.vcf.gz` | Per-chromosome GVCFs (input for Step 04) |

---

## Step 03b — GLIMPSE2 Imputation (low-coverage) — PAUSED

> **Status: Paused.** Reference panel chunks have not yet been generated. This step cannot be run until `GLIMPSE2_chunk` and `GLIMPSE2_split_reference` outputs are available and the `ref_panel` / `ref_gmap` config keys are populated.

For low-pass WGS (~3.5–8×). Uses GLIMPSE2 for reference-panel-based imputation and phasing. Produces the same per-chromosome VCF layout as Step 03a and feeds into Step 04.

```bash
# Do NOT run — blocked on reference chunk generation
bash bin/03_glimpse2_imputation.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
```

Config keys required (`config/config.yaml`):

| Key | Description |
|-----|-------------|
| `ref_panel` | Directory of per-chromosome reference BCF files |
| `ref_gmap` | Directory of per-chromosome genetic map files |

---

## Steps 04–06 — Joint Genotyping and Filtering (stubs)

These scripts are scaffolded but not yet implemented.

| Step | Script | Purpose |
|------|--------|---------|
| 04 | `04_gatk_GenomicsDB_import.sh` | Import per-sample GVCFs into GenomicsDB |
| 05 | `05_gatk_GenotypeGVCFs.sh` | Joint genotyping across all samples |
| 06 | `06_gatk_vqsr.sh` | VQSR filtering of joint VCF |

Do not use these for production runs. Watch for status updates in the pipeline overview table in [README.md](README.md).

---

## Monitoring Jobs

```bash
squeue -u $USER              # all your jobs
squeue -u $USER | grep BAMQC
squeue -u $USER | grep HAPCALL
```

## Logs

Each script writes a timestamped log to `log/`. Check there first if a job fails.

---

## Typical End-to-End Run (single sample, FENIX)

```bash
# 1. Align
bash bin/01_bwa_map_fastq_reads.sh /data/sample01/ /output/bams/

# 2. QC + BQSR
bash bin/02_gatk_bam_qc_workflow.sh /output/bams/sample01.sort.bam /output/bqsr/

# 3. Variant calling
bash bin/03_gatk_haplotype_caller.sh /output/bqsr/sample01.rmdup.mqfilt.bqsr.bam /output/gvcf/
```

For batches: use the `*_batch-slurmer.sh` wrappers in sequence, waiting for each stage to complete before submitting the next.
