# Joint Variant Calling Pipeline [FENIX]

Modular workflow for **joint germline variant discovery** following [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). Designed for human WGS data (hg38) on the **LAVIS-FENIX HPC**.

> See [INSTRUCTIONS.md](INSTRUCTIONS.md) for step-by-step usage and [docs/PIPELINE_STATUS.md](docs/PIPELINE_STATUS.md) for the current status report.

---

## Quick Run (FENIX)

### One sample

Submit Steps 01→02→03 as chained SLURM jobs with one command. The sample directory must be named after the sample ID and contain its FASTQ pair(s) — or, when resuming, the output of an earlier step.

```bash
# From inside the sample directory
bash bin/PIPELINE.single_sample_01-03.sh

# Or pointing to the sample directory explicitly
bash bin/PIPELINE.single_sample_01-03.sh /path/to/SAMPLE_ID/
```

### A whole batch

```bash
bash bin/supp/BATCH.submit_S01-S03.sh /path/to/BatchNN/
```

Runs the single-sample launcher on every sample subdirectory of the batch. One sample failing to submit doesn't stop the others.

### What the launcher does

- **Chains the steps** with `--dependency=afterok`, so each one waits for the previous to succeed.
- **Resumes**: a step whose final output already exists is skipped (S01 → `*.sort.bam` / `*.sorted.bam`, S02 → `*.rmdup.mqfilt.bqsr.bam`, S03 → all 25 `chrom_gvcf/` files + the `canon_chr` GVCF). Re-running the launcher on a partly-done or finished sample is always safe.
- **Scatters Step 03** (`SCATTER_S03=true`, default): a 25-task SLURM array, one HaplotypeCaller job per canonical chromosome, then a small gather job. That's ~2 h per sample instead of ~15 h, and a failure costs one chromosome. A resubmit only re-runs the missing chromosomes.
- **Uses scratch** (`USE_SCRATCH=true`, default): inputs are read directly from NFS, while `TMPDIR`, intermediates and outputs go to `<scratch_base>/<your primary group>/$USER/job_<jobid>`. Only the final outputs are copied back, and each copy is checked (size, `samtools quickcheck` / `gzip -t`). If a check fails, the scratch dir is kept so the data can be recovered. Scratch doesn't make these steps faster (see [INSTRUCTIONS.md](INSTRUCTIONS.md#scratch-storage-scratch--what-it-is-and-isnt-for)); it keeps bulky intermediates off the group quota.
- **Checks before submitting**: it refuses to start if the output filesystem has less than `min_free_gb` free, and keeps jobs off the nodes listed in `sbatch_exclude` (both in `config/config.yaml`).

Logs go to `SAMPLE_ID/log/`:

| Log | Step |
|-----|------|
| `SAMPLE-S01-<epoch>.<jobid>.log` | 01 alignment |
| `SAMPLE-S02-<epoch>.<jobid>.log` | 02 BAM QC + BQSR |
| `SAMPLE-S03-<epoch>.<arrayjob>_<task>.log` | 03 HaplotypeCaller, one per chromosome (task N = Nth of chr1..chr22, chrX, chrY, chrM) |
| `SAMPLE-S03g-<epoch>.<jobid>.log` | 03 gather (whole-genome `canon_chr` GVCF) |

```bash
squeue -u $USER | grep SAMPLE_ID
```

### Checking results

```bash
bash bin/supp/SUMMARY.sample_stats_S01-S03.sh /path/to/BatchNN/ [--gvcf-depth[=fast|full]]
```

Writes one TSV row per sample: stage reached, file sizes, time per step, reads / % mapped, % duplication, insert size, mean depth (mosdepth), and optionally depth re-derived from the GVCFs. It only reads files. See [INSTRUCTIONS.md](INSTRUCTIONS.md#per-sample-summary-table).

For more control over single steps across many samples, the per-step `*.seq_batch-slurmer.sh` wrappers are still available — see [INSTRUCTIONS.md](INSTRUCTIONS.md).

---

## Pipeline Overview

```
FASTQ reads
    │
    ▼
[01] bwa mem alignment          → sorted BAM per read group
    │
    ▼
[02] BAM QC & preprocessing     → analysis-ready BAM (.rmdup.mqfilt.bqsr.bam)
    │
    ├──► [02a] Retroactive BQSR evaluation (optional)
    │
    ├──► [03a] GATK HaplotypeCaller   → per-chromosome GVCFs  [standard coverage]
    │
    └──► [03b] GLIMPSE2 imputation    → per-chromosome VCFs   [low-coverage, PAUSED]
              
[03a] per-sample GVCFs (whole cohort)
              │
              ▼
         [04] GenomicsDB import       → one workspace per chromosome
              │
              ▼
         [05] GenotypeGVCFs           → joint-genotyped VCF
              │
              ▼
         [06] VQSR filtering          → filtered VCF
```

| Step | Script | Status |
|------|--------|--------|
| 01 | `01_bwa_map_fastq_reads.sh` | Functional |
| 02 | `02_gatk_bam_qc_workflow.sh` | Functional |
| 02a | `supp/02a_bqsr_evaluate.sh` | Functional |
| 03a | `03_gatk_haplotype_caller.sh` | Functional — per-chromosome scatter validated in production |
| 03b | `03_glimpse2_imputation.sh` | **Paused** — reference chunks not yet generated |
| 04 | `04_gatk_GenomicsDB_import.sh` | **Functional** — validated on FENIX (JAGUAR chr22, 93 samples, 3-wave incremental); cohort-specific all-chromosome launcher in `supp/JAGUAR.S04_run_all_chroms.sh` |
| 05 | `05_gatk_GenotypeGVCFs.sh` | Stub (on hold) |
| 06 | `06_gatk_vqsr.sh` | Stub (on hold) |

See [docs/S04_GenomicsDBImport_design.md](docs/S04_GenomicsDBImport_design.md) for Step 04 design notes.

Steps 03a and 03b are alternatives — choose based on coverage depth. Steps 04–06 are **cohort-level** (run once over all samples) and feed from 03a output.

---

## Configuration

Edit `config/config.yaml`, which has environment-specific values under `remote:` (FENIX) and `local:` keys. Scripts auto-detect the environment (SSH session → `remote`) and parse the config without an external YAML parser.

| Key | Used in | Description |
|-----|---------|-------------|
| `ref_gnm` | 01, 02, 03a, 04 | hg38 reference FASTA (indexed; BWA index for Step 01) |
| `ref_vars` | 02, 03a | dbSNP VCF (canonical chromosomes, bgzipped + tabixed) |
| `ref_panel` | 03b | GLIMPSE2 reference panel directory |
| `ref_gmap` | 03b | Genetic maps directory |
| `scratch_base` | S01–03 launcher | Job scratch = `<scratch_base>/<primary group>/<user>/job_<id>` |
| `scratch_root` | 04, supp scripts | Fixed scratch tree (group-pinned) |
| `sbatch_exclude` | S01–03 launcher | Comma-separated nodes to keep jobs off (`""` = none) |
| `min_free_gb` | S01–03 launcher | Abort submission if the output filesystem has less free space (`0` = no check) |
| `modules` | — | Module list (informational; scripts load their own modules) |

---

## Dependencies

Loaded with `module load` on FENIX by each script.

| Tool | Steps |
|------|-------|
| `bwa`, `samtools` ≥ 1.15 | 01, 02 |
| `fastqc`, `bbtools` | 01 (QC, pair repair) |
| `gatk` ≥ 4.x (with `oracle-java/25.0.2` on FENIX) | 02, 03a, 04 |
| `bcftools` | 03a, 03b, summary (`--gvcf-depth`) |
| `mosdepth` ≥ 0.3 | 02 |
| `R` ≥ 4.4.1 | 02 / 02a (optional BQSR plots) |
| `GLIMPSE2_*` | 03b |

---

## Authors

- Pavel Salazar-Fernandez (maintainer): epsalazarf@gmail.com
- Dr. Federico Sanchez-Quinto (project leader)
