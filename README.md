# Joint Variant Calling Pipeline [FENIX]

This repository contains a modular workflow for performing **joint germline variant discovery** following [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).
The pipeline is designed for human WGS data (hg38) and adapted for the **LAVIS-FENIX HPC** environment.

> **Note:** Scripts are under active development. Steps 01–03 are functional; Steps 04–06 are stubs pending implementation. Steps 01 and 02 must be used together for new samples — script 01 embeds proper read groups so script 02's RG correction step is skipped automatically.
>
> For low-coverage WGS (~3.5–8X), use `03_glimpse2_imputation.sh` instead of `03_gatk_haplotype_caller.sh`. Both scripts accept the same input (`.rmdup.mqfilt.bqsr.bam`) and produce per-chromosome VCFs under `chrom_gvcf/`.

---

## Pipeline Overview

| Step | Script | Status | Input → Output |
|------|--------|--------|----------------|
| 01 | `01_bwa_map_fastq_reads.sh` | Functional | FASTQ → sorted BAM per read group |
| 02 | `02_gatk_bam_qc_workflow.sh` | Functional | raw BAM → analysis-ready BAM |
| 03a | `03_gatk_haplotype_caller.sh` | Functional | BAM → per-chromosome GVCFs *(standard coverage)* |
| 03b | `03_glimpse2_imputation.sh` | Functional | BAM → per-chromosome imputed VCFs *(low-coverage)* |
| 04 | `04_gatk_GenomicsDB_import.sh` | Stub | GVCFs → GenomicsDB |
| 05 | `05_gatk_GenotypeGVCFs.sh` | Stub | GenomicsDB → joint VCF |
| 06 | `06_gatk_vqsr.sh` | Stub | joint VCF → filtered VCF |

Each step is run sequentially; the output of one step is the input for the next. Steps 03a and 03b are alternatives — choose based on coverage depth.

---

## Step 01 — FASTQ Alignment

Maps raw FASTQ reads to the reference genome per read group and produces coordinate-sorted, indexed BAM files.

### Usage

```bash
bash 01_bwa_map_fastq_reads.sh [input_dir] [output_path]
```

Both arguments default to the current directory. The script must be run from the directory containing the FASTQ files; `input_dir` is used only to derive the sample name.

### Cluster submission

```bash
bash BWAMAP.seq_batch-slurmer.sh bin/01_bwa_map_fastq_reads.sh /path/to/samples/
```

Auto-discovers sample subdirectories containing FASTQ files and chains jobs sequentially on SLURM.

### Supported FASTQ naming conventions

| Provider | Pattern | Example |
|----------|---------|---------|
| Multi-lane (old) | `{SAMPLE}_{BARCODE}-1A_{PLATE}_{LANE}_{R}.fq.gz` | `L23_CKDN...-1A_227CC2LT4_L5_1.fq.gz` |
| Single-pair | `{SAMPLE}_{R}.fq.gz` | `HL078_1.fq.gz` |
| EGAN-style | `{SAMPLE}_R{R}.fastq.gz` | `EGAN00004552350_R1.fastq.gz` |

R0 (index read) files are automatically ignored.

### Output

- `{SAMPLE}.bam` — mapped BAM (single pair)
- `{SAMPLE}_{PLATE_LANE}.bam` — mapped BAM per read group (multi-lane)
- `{SAMPLE}[-{RG}].sort.bam` + `.bai` — coordinate-sorted, indexed BAM (final output)

### Processing steps

1. **Map Input Files** — detect FASTQ pairs, write inputs manifest per read group
2. **Annotate Read Groups** *(optional, `BUILD_RG=false`)* — prepend `@RG` tag to manifest; only meaningful for multi-lane samples
3. **Read Quality Check** *(optional, `RUN_FASTQC=false`)* — FastQC report for all sample reads
4. **BWA Mapping** — `bwa mem` per read group, piped to `samtools view` (unsorted BAM)
5. **Sort and Index** — `samtools sort` + `samtools index` per BAM

---

## Step 02 — BAM QC and Preprocessing

Prepares mapped BAMs for variant calling: optional read group correction (legacy only), deduplication with merge support for multiplexed samples, optional MQ filtering, base quality recalibration, and coverage estimation.

### Usage

```bash
# New workflow — single sorted BAM from script 01 (RGs already embedded)
bash 02_gatk_bam_qc_workflow.sh SAMPLE.sort.bam [output_path]

# New workflow — multiplexed sample (multiple BAMs merged during deduplication)
bash 02_gatk_bam_qc_workflow.sh "SAMPLE_plateA_lane1.sort.bam,SAMPLE_plateB_lane3.sort.bam" [output_path]

# Legacy workflow — BAM with missing RGs; explicit RG string triggers step 0
bash 02_gatk_bam_qc_workflow.sh SAMPLE.bam [output_path] "@RG\tID:...\tSM:...\tPL:ILLUMINA\t..."
```

Argument 4 (`add_rg`): `auto` (default — step 0 runs only if RG string is provided) | `true` | `false`

### Cluster submission

```bash
bash BAMQC.seq_batch-slurmer.sh batch-list.txt bin/02_gatk_bam_qc_workflow.sh /path/to/bams/
```

Batch list (legacy): three tab-separated columns — `sample_name  sample_file  readgroup_string`

### Output

All output files are named using the sample name as base (e.g. `SAMPLE01`).

| File | Description |
|------|-------------|
| `SAMPLE.rg.bam` | Read groups corrected *(step 0, legacy workflow only)* |
| `SAMPLE.rmdup.bam` | Duplicates marked *(intermediate)* |
| `SAMPLE-dups.txt` | Duplicate metrics |
| `SAMPLE.rmdup.mqfilt.bam` | MQ ≥ 30 filtered *(intermediate, if `MQ_FILTER=true`)* |
| `SAMPLE.rmdup.mqfilt.bqsr.bam` | BQSR-recalibrated — **final analysis-ready BAM** |
| `SAMPLE.rmb.mosdepth.*` | Coverage summary |
| `SAMPLE.*.metrics.txt`, `*.pdf` | Alignment and insert-size metrics *(if `RUN_METRICS=true`)* |

### Processing steps

0. **Assign Read Groups** *(optional — legacy only; skipped automatically for BAMs from script 01)* — GATK `AddOrReplaceReadGroups`; controlled by `add_rg` argument
1. **Mark Duplicates** — GATK `MarkDuplicatesSpark`; accepts one or more input BAMs (merged on the fly); set `REMOVE_DUPS=true` to remove instead of mark
2. **MQ Filter** *(optional, `MQ_FILTER=false`)* — `samtools view -q 30`; disabled by default; originally for aDNA
3. **BQSR** — GATK `BaseRecalibrator` + `ApplyBQSR` using known variant sites
4. **Coverage** — `mosdepth` fast-mode genome-wide summary
5. **Alignment Metrics** *(optional, `RUN_METRICS=true`)* — GATK `CollectAlignmentSummaryMetrics` + `CollectInsertSizeMetrics`

---

## Step 03 — Per-Sample GVCF Generation

Runs GATK HaplotypeCaller in GVCF mode and prepares chromosome-level subsets for joint genotyping.

### Usage

```bash
bash 03_gatk_haplotype_caller.sh <bqsr.bam> [output_path]
```

### Cluster submission

```bash
bash HAPCALL.seq_batch-slurmer.sh bin/03_gatk_haplotype_caller.sh /path/to/bqsr/bams/
```

Auto-discovers `*.rmdup.mqfilt.bqsr.bam` files and chains jobs sequentially.

### Output

| File | Description |
|------|-------------|
| `*.raw_variants.g.vcf.gz` | Full raw GVCF |
| `*.raw_variants.canon_chr.g.vcf.gz` | Canonical chromosomes (chr1–22, X, Y, M) |
| `chrom_gvcf/*.raw_vars.{CHR}.g.vcf.gz` | Per-chromosome GVCFs |

### Processing steps

1. **HaplotypeCaller** — GVCF mode with dbSNP annotation
2. **Index GVCF** — `bcftools index --tbi`
3. **Extract Canonical Chromosomes** — `bcftools view` filtering to primary chromosomes
4. **Split by Chromosome** — one GVCF per chromosome for parallel joint genotyping

---

## Step 03b — Per-Sample Imputation and Phasing (Low-Coverage WGS)

An alternative to Step 03a for low-pass WGS data (~3.5–8X coverage). Uses [GLIMPSE2](https://odelaneau.github.io/GLIMPSE/) for reference-panel-based imputation and phasing, producing phased genotype calls from shallow sequencing without a separate genotype likelihood step.

Reference panel preparation (Steps 1–2) is cached in a shared directory and reused across samples, so the first sample incurs the setup cost and subsequent samples skip directly to imputation.

### Usage

```bash
bash 03_glimpse2_imputation.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
```

### Output

| File | Description |
|------|-------------|
| `chrom_gvcf/SAMPLE.imputed.chrN.vcf.gz` + `.tbi` | Per-chromosome phased, imputed VCF |
| `SAMPLE.imputed.canon_chr.vcf.gz` + `.tbi` | Whole-genome merged VCF (chr1–22, X, Y, M) |
| `imputed_chunks/` | Intermediate chunk BCFs and ligate lists *(removed if `HOUSEKEEP=true`)* |

### Processing steps

1. **Chunk Chromosomes** *(cached)* — `GLIMPSE2_chunk` defines imputation windows per chromosome using the reference panel and genetic map
2. **Split Reference** *(cached)* — `GLIMPSE2_split_reference` builds binary reference panel files per chunk
3. **Phase and Impute** — `GLIMPSE2_phase` computes genotype likelihoods from the BAM and imputes genotypes per chunk
4. **Ligate** — `GLIMPSE2_ligate` merges chunks into a single chromosome-wide phased BCF, converted to VCF.gz
5. **Collect** — `bcftools concat` merges all per-chromosome VCFs into the final whole-genome output

### Configuration keys (in `config/config.yaml`)

| Key | Description |
|-----|-------------|
| `ref_panel` | Directory of per-chromosome reference panel BCF files (`reference_panel.chrN.bcf`) |
| `ref_gmap` | Directory of per-chromosome genetic map files (`chrN.b38.gmap.gz`); chrY/chrM skipped if absent |

Chromosomes for which the reference panel BCF or genetic map is missing are skipped gracefully with an `[i]` message.

---

## Configuration

All environment-specific paths are set in `config/config.yaml` under `remote:` and `local:` keys. Scripts auto-detect the environment via SSH session variables and parse the config with an embedded `awk` snippet (no external YAML parser required).

| Key | Description |
|-----|-------------|
| `ref_gnm` | hg38 reference FASTA (must have `.fai` and BWA index) |
| `ref_vars` | dbSNP VCF (canonical chromosomes, bgzipped + tabixed) |
| `ref_panel` | Directory of GLIMPSE2 reference panel BCF files (Step 03b only) |
| `ref_gmap` | Directory of GLIMPSE2 genetic map files, one per chromosome (Step 03b only) |

---

## Dependencies

Tools are loaded automatically via `module load` on FENIX; must be in `$PATH` for local use.

| Tool | Version | Used in |
|------|---------|---------|
| `bwa` | any | Step 01 |
| `samtools` | ≥ 1.15 | Steps 01, 02 |
| `fastqc` | ≥ 0.12.1 | Step 01 (optional) |
| `bbtools` (`repair.sh`) | ≥ 38.00 | Step 01 (single-pair samples) |
| `gatk` | ≥ 4.x | Steps 02, 03a |
| `bcftools` | any | Steps 02, 03a, 03b |
| `mosdepth` | ≥ 0.3.x | Step 02 (optional) |
| `R` | ≥ 4.4.1 | Step 02 (optional plots) |
| `GLIMPSE2_chunk` / `GLIMPSE2_split_reference` / `GLIMPSE2_phase` / `GLIMPSE2_ligate` | ≥ 2.0 | Step 03b |

---

## Authors

- Pavel Salazar-Fernandez (current maintainer): epsalazarf@gmail.com
- Dr. Federico Sanchez-Quinto (project leader)
