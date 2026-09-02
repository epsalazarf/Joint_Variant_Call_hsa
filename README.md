# Joint Variant Calling Pipeline [FENIX]

Modular workflow for **joint germline variant discovery** following [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-). Designed for human WGS data (hg38) on the **LAVIS-FENIX HPC**.

> See [INSTRUCTIONS.md](INSTRUCTIONS.md) for step-by-step usage.

---

## Quick Run (single sample, FENIX)

Submit Steps 01→02→03 as chained SLURM jobs with one command. The sample directory must be named after the sample ID and contain at least one FASTQ pair.

```bash
# From inside the sample directory
bash bin/PIPELINE.single_sample.sh

# Or pointing to the sample directory explicitly
bash bin/PIPELINE.single_sample.sh /path/to/SAMPLE_ID/
```

Each step is held until the previous one succeeds (`--dependency=afterok`). Logs are written to `SAMPLE_ID/log/S0{1,2,3}.<jobid>.log`.

```bash
# Monitor all three jobs
squeue -u $USER | grep PIPELINE
```

**Scratch I/O (`USE_SCRATCH=true`, default):** each job stages its working files through `/scratch` instead of writing directly to NFS. Input BAMs are copied to a per-job scratch directory (`/scratch/groups/amedina/$USER/job_$SLURM_JOB_ID`), all intermediate and output files are written there, and only the final outputs are copied back to the sample directory on completion. This is especially beneficial for Step 02 (`MarkDuplicatesSpark`), which generates 2–3× BAM-size temp files during Spark local-mode shuffle. Set `USE_SCRATCH=false` at the top of `PIPELINE.single_sample.sh` to bypass scratch (e.g. for local testing).

For batches or more control over each step, use the individual `*_batch-slurmer.sh` wrappers described in [INSTRUCTIONS.md](INSTRUCTIONS.md).

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
              │
              ▼
         [04] GenomicsDB import       → GenomicsDB
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
| 02a | `02a_bqsr_evaluate.sh` | Functional |
| 03a | `03_gatk_haplotype_caller.sh` | Functional |
| 03b | `03_glimpse2_imputation.sh` | **Paused** — reference chunks not yet generated |
| 04 | `04_gatk_GenomicsDB_import.sh` | **Functional** — validated on FENIX (JAGUAR chr22, 93 samples, 3-wave incremental); whole-genome scatter pending |
| 05 | `05_gatk_GenotypeGVCFs.sh` | Stub |
| 06 | `06_gatk_vqsr.sh` | Stub |

See [docs/PIPELINE_STATUS.md](docs/PIPELINE_STATUS.md) for the full status report and [docs/S04_GenomicsDBImport_design.md](docs/S04_GenomicsDBImport_design.md) for Step 04 design notes.

Steps 03a and 03b are alternatives — choose based on coverage depth. Steps 04–06 feed from 03a output.

---

## Configuration

Edit `config/config.yaml` with environment-specific paths under `remote:` and `local:` keys. Scripts auto-detect the environment and parse the config without an external YAML parser.

| Key | Used in |
|-----|---------|
| `ref_gnm` | Steps 01, 02, 03a |
| `ref_vars` | Step 02 (BQSR) |
| `ref_panel` | Step 03b |
| `ref_gmap` | Step 03b |

---

## Dependencies

| Tool | Steps |
|------|-------|
| `bwa`, `samtools` | 01, 02 |
| `fastqc`, `bbtools` | 01 (optional) |
| `gatk` ≥ 4.x | 02, 03a |
| `bcftools` | 02, 03a, 03b |
| `mosdepth` | 02 (optional) |
| `R` ≥ 4.4.1 | 02 (optional plots) |
| `GLIMPSE2_*` | 03b |

---

## Authors

- Pavel Salazar-Fernandez (maintainer): epsalazarf@gmail.com
- Dr. Federico Sanchez-Quinto (project leader)
