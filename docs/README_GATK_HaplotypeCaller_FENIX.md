# GATK HaplotypeCaller [FENIX] — Step 03a

**Script:** `bin/03_gatk_haplotype_caller.sh` (v1.2)
**Author:** Pavel Salazar-Fernandez et al.
**Source:** [GATK4 Best Practices – Germline Short Variant Discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932)

---

## Overview

Per-sample germline calling (SNPs + indels) with **GATK4 HaplotypeCaller in GVCF mode** (`-ERC GVCF`, dbSNP annotation) from the analysis-ready BAM produced by Step 02. Its main product is `chrom_gvcf/`: one GVCF per canonical chromosome, which Step 04 (GenomicsDBImport) reads.

The script has two modes:

| Mode | Call | Use |
|------|------|-----|
| **scatter** (one chromosome) | `03_gatk_haplotype_caller.sh <bam> <out> chrN` | Default in the launcher: 25 of these run as a SLURM array |
| **whole-genome** | `03_gatk_haplotype_caller.sh <bam> <out>` | One long job (~15 h), then split by chromosome |

`chrN` is one of `chr1`..`chr22`, `chrX`, `chrY`, `chrM`. `chrM` is called with `--sample-ploidy 1`.

---

## Usage

### Recommended — via the launcher

```bash
bash bin/PIPELINE.single_sample_01-03.sh /path/to/SAMPLE_ID/
```

With `SCATTER_S03=true` (default), Step 03 is submitted as a **25-task array** (2 CPU / 13 GB / 48 h each, `ARRAY_CONC=6` at a time). Once every chromosome succeeds, a **gather** job concatenates them into the whole-genome `canon_chr` GVCF, which is kept for archiving only. Per-sample wall time is about 2 h (the chr1 task). A failed task only costs that chromosome, and relaunching resubmits only the chromosomes whose GVCF + `.tbi` is missing.

### Manual

```bash
# one chromosome → <out>/chrom_gvcf/SAMPLE.raw_vars.chr7.g.vcf.gz
bash bin/03_gatk_haplotype_caller.sh SAMPLE.rmdup.mqfilt.bqsr.bam <out> chr7

# whole genome
bash bin/03_gatk_haplotype_caller.sh SAMPLE.rmdup.mqfilt.bqsr.bam <out>
```

### Legacy batch

```bash
bash bin/supp/HAPCALL.seq_batch-slurmer.sh bin/03_gatk_haplotype_caller.sh /path/to/bqsr/bams/
```

One whole-genome job per sample (4 CPU / 32 GB). No scatter, no scratch, no verified copy-back.

---

## Workflow

**Scatter mode**
1. HaplotypeCaller `--intervals chrN` → hidden `.tmp` GVCF + index in `chrom_gvcf/`
2. Rename to the final name only on success (skip if the final GVCF + `.tbi` already exist)

**Whole-genome mode**
1. HaplotypeCaller → `SAMPLE.raw_variants.g.vcf.gz`
2. Index (`bcftools index --tbi`)
3. Keep canonical chromosomes → `SAMPLE.raw_variants.canon_chr.g.vcf.gz`
4. Split by chromosome → `chrom_gvcf/SAMPLE.raw_vars.<CHR>.g.vcf.gz`

The sample prefix is the BAM filename up to the first `.` (so `SAMPLE.rmdup.mqfilt.bqsr.bam` → `SAMPLE`).

---

## Resources

| Mode | HC heap | pair-HMM threads |
|------|---------|------------------|
| scatter (FENIX) | `-Xms2G -Xmx10G` | 2 |
| whole-genome (FENIX) | `-Xms20G -Xmx20G` | 4 |
| local | 8G / 12G max | 2 / 4 |

The small initial heap in scatter mode stops short tasks (chrY, chrM) from reserving memory they never use, which makes them easier for SLURM to place.

---

## Output

| File | Description |
|------|-------------|
| `chrom_gvcf/SAMPLE.raw_vars.<CHR>.g.vcf.gz` + `.tbi` | Per-chromosome GVCFs — **Step 04 input** |
| `SAMPLE.raw_variants.canon_chr.g.vcf.gz` + `.tbi` | Whole-genome canonical GVCF (gather / whole-genome mode) — archival |
| `SAMPLE.raw_variants.g.vcf.gz` | All-contigs GVCF — whole-genome mode only; not copied back from scratch by the launcher |

To check the depth of finished GVCFs, run `bin/supp/SUMMARY.sample_stats_S01-S03.sh <batch> --gvcf-depth`.

---

## Dependencies

- `gatk` ≥ 4.x (on FENIX: `oracle-java/25.0.2` + `gatk` modules, loaded by the script)
- `bcftools`, `samtools`
- `ref_gnm` and `ref_vars` in `config/config.yaml`
- bash ≥ 4 (`EPOCHSECONDS`)

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| `CANCELLED. File not found` | Wrong BAM path | Check the path / symlink |
| `Not a canonical chromosome` | Bad 3rd argument | Use `chr1`..`chr22`, `chrX`, `chrY`, `chrM` |
| `Missing reference paths` | Config keys empty | Set `ref_gnm` / `ref_vars` for the environment |
| Array task dies instantly with `scratch unavailable on <node>` | Broken scratch mount on that node | Add the node to `sbatch_exclude` in the config, then relaunch |
| JVM `SIGSEGV` on one node | Bad hardware | Exclude the node (`sbatch_exclude`) and relaunch; only failed chromosomes rerun |
| `chrom_gvcfs` < 25/25 in the summary table | Some array tasks failed | Check `log/SAMPLE-S03-*_<task>.log`, then relaunch the sample |
