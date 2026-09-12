# Pipeline Instructions

Step-by-step usage guide for the Joint Variant Calling Pipeline.

> For a high-level overview, see [README.md](README.md).

---

## Script Layout

```
bin/
├── 01_bwa_map_fastq_reads.sh       # Step 01 — FASTQ alignment
├── 02_gatk_bam_qc_workflow.sh      # Step 02 — BAM QC + BQSR
├── 03_gatk_haplotype_caller.sh     # Step 03a — HaplotypeCaller (standard coverage)
├── 03_glimpse2_imputation.sh       # Step 03b — GLIMPSE2 imputation (low-coverage)
├── 04_gatk_GenomicsDB_import.sh    # Step 04 — GenomicsDB import (beta)
├── 05_gatk_GenotypeGVCFs.sh        # Step 05 — Joint genotyping (draft)
├── 06_gatk_vqsr.sh                 # Step 06 — VQSR / hard-filtering (draft)
└── supp/
    ├── 00_scan_fastq_pairs.sh              # Inspect FASTQ pairs before mapping
    ├── 00_glimpse2_ref_panel_prep.sh       # Prepare GLIMPSE2 binary reference panel
    ├── 02a_bqsr_evaluate.sh                # Retroactive BQSR covariate plots
    ├── 03-s4_gvcf_chrom_split.sh           # Split GVCFs by chromosome
    ├── BWAMAP.seq_batch-slurmer.sh         # Batch SLURM launcher for Step 01
    ├── BAMQC.seq_batch-slurmer.sh          # Batch SLURM launcher for Step 02
    ├── HAPCALL.seq_batch-slurmer.sh        # Batch SLURM launcher for Step 03a
    ├── GENOTYPE.seq_batch-slurmer.sh       # Per-chromosome scatter+gather launcher for Step 05
    ├── VQSR.seq_batch-slurmer.sh           # SLURM launcher for Step 06
    ├── PIPELINE.single_sample.sh           # End-to-end single-sample launcher
    └── run_pipeline.sh                     # Full pipeline wrapper (stub)
```

---

## Prerequisites

1. **Configure paths** in `config/config.yaml` — set `ref_gnm`, `ref_vars`, and any other keys for your environment (`remote:` for FENIX, `local:` for workstation).
2. **Load modules** (FENIX only): scripts call `module load` automatically based on the `modules` key in the config.
3. All scripts must be run from the repo root or with explicit paths.

---

## Step 00 — Scan FASTQ Pairs (optional)

Inspect a sample directory before mapping to verify FASTQ pairs are correctly detected.

```bash
bash bin/supp/00_scan_fastq_pairs.sh [input_dir]
```

Defaults to `$PWD`. Prints detected pairs without running anything.

---

## Step 00 — GLIMPSE2 Reference Panel Prep (one-time setup)

Converts a phased reference panel VCF into the binary format required by GLIMPSE2.
Run once per chromosome before using Step 03b. Jobs can be submitted in parallel.

```bash
# With QC (normalize + filter to biallelic SNPs)
sbatch bin/supp/00_glimpse2_ref_panel_prep.sh <panel.chr1.vcf.gz>

# Without QC (panel already normalized — e.g. H1K2 ARPnorm)
sbatch bin/supp/00_glimpse2_ref_panel_prep.sh <panel.chr1.vcf.gz> "" false

# All chromosomes in parallel
for vcf in /path/to/panel/panel.chr*.vcf.gz; do
  sbatch bin/supp/00_glimpse2_ref_panel_prep.sh "$vcf" "" false
done
```

Arguments: `<panel.vcf.gz>  [output_panel_dir]  [run_qc: true|false]`  
`output_panel_dir` defaults to the `ref_panel` path in `config/config.yaml`.

### Output

| Path | Description |
|------|-------------|
| `<ref_panel>/reference_panel.<CHR>.bcf` | QC'd BCF (indexed) |
| `<ref_panel>/../glimpse2_cache/chunks.<CHR>.txt` | Imputation window definitions |
| `<ref_panel>/../glimpse2_cache/ref_panel.<CHR>.chunk*.bin` | Binary panels for `GLIMPSE2_phase` |

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
bash bin/supp/BWAMAP.seq_batch-slurmer.sh bin/01_bwa_map_fastq_reads.sh /path/to/samples/
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
bash bin/supp/BAMQC.seq_batch-slurmer.sh batch-list.txt bin/02_gatk_bam_qc_workflow.sh /path/to/bams/
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
bash bin/supp/02a_bqsr_evaluate.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
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
# whole genome (one long call, then split by chromosome)
bash bin/03_gatk_haplotype_caller.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]

# one chromosome only — writes straight into chrom_gvcf/ (used by the scatter)
bash bin/03_gatk_haplotype_caller.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path] chr7
```

### Recommended: per-chromosome scatter via the launcher

`bin/PIPELINE.single_sample_01-03.sh` (`SCATTER_S03=true`, default) submits Step 03
as a **25-way SLURM array** — one `HaplotypeCaller -L <chrom>` task per canonical
chromosome (2 CPU / 13 GB each, `%ARRAY_CONC` concurrent), followed by a small
gather job that concatenates them into the whole-genome `canon_chr` GVCF.

Per-sample wall time drops from ~15 h to ~2 h (the chr1 task), and a failure
costs one chromosome, not the sample. `chrom_gvcf/` — the Step 04 input — is
identical either way. Resubmitting re-runs only the chromosomes whose GVCF is
missing. Set `SCATTER_S03=false` for the legacy single whole-genome job.

### Batch (SLURM, legacy)

```bash
bash bin/supp/HAPCALL.seq_batch-slurmer.sh bin/03_gatk_haplotype_caller.sh /path/to/bqsr/bams/
```

Auto-discovers `*.rmdup.mqfilt.bqsr.bam` files and submits one whole-genome job
per sample (4 CPUs / 32 GB). No scatter, no verified copy-back.

```bash
squeue -u $USER | grep -E 'S03|HAPCALL'
```

### Output

| File | Description |
|------|-------------|
| `chrom_gvcf/*.raw_vars.{CHR}.g.vcf.gz` | Per-chromosome GVCFs (**input for Step 04**) |
| `*.raw_variants.canon_chr.g.vcf.gz` | Canonical-chromosome roll-up (archival / QC) |
| `*.raw_variants.g.vcf.gz` | Full all-contigs GVCF — whole-genome mode only |

---

## Step 03b — GLIMPSE2 Imputation (low-coverage)

For low-pass WGS (~3.5–8×). Uses GLIMPSE2 for reference-panel-based imputation and phasing.
Produces the same per-chromosome VCF layout as Step 03a and feeds into Step 04.

> **Prerequisite:** run `bin/supp/00_glimpse2_ref_panel_prep.sh` for all chromosomes first
> and confirm `ref_panel` / `ref_gmap` keys are populated in `config/config.yaml`.

```bash
bash bin/03_glimpse2_imputation.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
```

Config keys required (`config/config.yaml`):

| Key | Description |
|-----|-------------|
| `ref_panel` | Directory of per-chromosome reference BCF files |
| `ref_gmap` | Directory of per-chromosome genetic map files |

---

## Step 04 — GenomicsDB Import (cohort-level)

Consolidates the per-chromosome GVCFs from Step 03a into GenomicsDB workspaces,
**one per chromosome**, for joint genotyping. This is a cohort step: run it once
for the whole set of samples you want genotyped together.

> **Beta** — local tests pass; not yet validated on FENIX with real data.
> See [docs/S04_GenomicsDBImport_design.md](docs/S04_GenomicsDBImport_design.md).

### 1. Build the sample map

A 2-column, TAB-separated file — one line per sample:

```
LUPUS001	/mnt/data/amedina/esalazarf/lupus/LUPUS001/chrom_gvcf
LUPUS002	/mnt/data/amedina/esalazarf/lupus/LUPUS002/chrom_gvcf
```

Column 2 is the `chrom_gvcf/` directory produced by Step 03a. Generate it from a
cohort directory that holds one sub-directory per sample:

```bash
cd /path/to/cohort
for d in */chrom_gvcf; do
  s=$(basename "$(dirname "$d")")
  printf '%s\t%s\n' "$s" "$(readlink -f "$d")"
done > cohort.sample_map.tsv
```

Always eyeball the file before launching. Lines that are blank or start with `#`
are ignored.

### 2. Run (single chromosome — start here for testing)

```bash
bash bin/04_gatk_GenomicsDB_import.sh cohort.sample_map.tsv <output_path> chr22 create
```

| Argument | Meaning |
|----------|---------|
| `sample_map` | the TSV from step 1 |
| `output_path` | workspaces created at `<output_path>/genomicsdb/<chrom>` — **use group storage (`/mnt/data/...`), never `$HOME`** (FENIX home has a ~5 GB quota; the script refuses a `$HOME` path) |
| `chrom` | `chr1`..`chr22` \| `chrX` \| `chrY` \| `chrM` \| `autosomes` \| `all` |
| `action` | `create` (default) \| `update` |

`autosomes` / `all` process every chromosome **serially in one process** — fine
for tests and small cohorts. A per-chromosome SLURM launcher will come later.

### 3. Adding samples later (waves)

```bash
# first wave
bash bin/04_gatk_GenomicsDB_import.sh wave1.sample_map.tsv <out> chr22 create
# later wave — only the NEW samples in this map
bash bin/04_gatk_GenomicsDB_import.sh wave2.sample_map.tsv <out> chr22 update
```

### Output

| Path | Description |
|------|-------------|
| `<output_path>/genomicsdb/<CHR>/` | GenomicsDB workspace (input for Step 05 as `gendb://...`) |

The workspace is built on `/scratch` (when available) and copied back on
success.

---

## Step 05 — Joint Genotyping (cohort-level)

Runs GenotypeGVCFs against each per-chromosome GenomicsDB workspace from
Step 04, then gathers the results into one cohort-level VCF.

> **Draft** — written against Step 04's design, not yet run on FENIX or
> against a real workspace. See the STATUS block in the script header before
> using it for anything but testing.

### Run (single chromosome — start here for testing)

```bash
bash bin/05_gatk_GenotypeGVCFs.sh <genomicsdb_path> <output_path> chr22 my_cohort
```

| Argument | Meaning |
|----------|---------|
| `genomicsdb_path` | the `output_path` given to Step 04 — workspaces are read from `<genomicsdb_path>/genomicsdb/<chrom>` |
| `output_path` | where this step writes its own outputs (see below) |
| `chrom` | `chr1`..`chr22` \| `chrX` \| `chrY` \| `chrM` \| `autosomes` \| `all` \| a comma-separated list of the above |
| `cohort_name` | optional label used in output filenames (default: `cohort`) |

`autosomes` / `all` genotype every chromosome **serially in one process** —
fine for tests and small cohorts.

### Scatter across SLURM jobs (real cohorts)

```bash
bash bin/supp/GENOTYPE.seq_batch-slurmer.sh [options] <genomicsdb_path> <output_path> <cohort_name>
```

Submits one job per chromosome (`GENO_PHASE=scatter`), then a final job
(`GENO_PHASE=gather`, `--dependency=afterok` on every scatter job) that
concatenates them into the cohort VCF. See `--help` for `--chroms` /
`--cpus` / `--mem` / `--hours` / `--gather-mem` / `--gather-hours`, and
`report <manifest_file>` for post-run `sacct` timings. Memory/time defaults
are chromosome-scaled the way Step 04's are, borrowed as a starting point —
see the STATUS block in `bin/05_gatk_GenotypeGVCFs.sh` for why that's a
placeholder, not a measurement, for this step.

Do **not** call `05_gatk_GenotypeGVCFs.sh` directly with `GENO_PHASE=all` for
a chromosome set that a scatter run already split across jobs — each
`GENO_PHASE=all` (or `scatter`) job only knows its own chromosome(s), so
running the gather step per-job would each independently (and wrongly)
overwrite the cohort VCF. Follow scatter with `GENO_PHASE=gather` (same
chrom selector), which the launcher already does for you.

### Output

| Path | Description |
|------|-------------|
| `<output_path>/chrom_vcf/<cohort>.joint.<CHR>.vcf.gz` | per-chromosome joint-genotyped VCF |
| `<output_path>/<cohort>.joint.vcf.gz` | gathered cohort-level callset (input for Step 06) |

Before gathering, the script checks that every per-chromosome VCF carries the
same sample set — an uneven Step 04 wave import (a sample added to some
chromosomes' workspaces but not others) is caught here rather than silently
producing a malformed callset.

---

## Step 06 — VQSR / Hard-Filtering (cohort-level)

Filters the gathered joint VCF from Step 05. Two paths, chosen automatically
from cohort size or forced explicitly:

- `vqsr` — VariantRecalibrator + ApplyVQSR (SNP and INDEL modes, chained).
  Needs the GATK hg38 resource bundle (`ref_hapmap`, `ref_omni`,
  `ref_1kg_snp`, `ref_mills` in `config/config.yaml`).
- `hard-filter` — GATK's standard hard-filter expressions (SNP and INDEL,
  merged). No extra resource files; the right choice while a cohort is too
  small for VQSR's Gaussian mixture model to converge.

> **Draft** — written against GATK4's published Best Practices recipes, not
> yet run on FENIX or against a real joint VCF, and the VQSR resource bundle
> is not staged yet (`ref_hapmap`/`ref_omni`/`ref_1kg_snp`/`ref_mills` are
> `EDIT_THIS` placeholders in `config/config.yaml`). See the STATUS block in
> the script header before using it for anything but testing.

### Run

```bash
bash bin/06_gatk_vqsr.sh <joint_vcf> <output_path> [vqsr|hard-filter|auto] my_cohort
```

| Argument | Meaning |
|----------|---------|
| `joint_vcf` | the gathered `<cohort>.joint.vcf.gz` from Step 05 |
| `output_path` | where the filtered VCF and intermediates are written |
| `mode` | `vqsr` \| `hard-filter` \| `auto` (default: `auto`) |
| `cohort_name` | optional label for output filenames (default: derived from `joint_vcf`'s own name) |

`auto` counts samples in `joint_vcf` and picks `vqsr` at `VQSR_MIN_SAMPLES`
(default 30 — a commonly cited GATK rule of thumb, not a project-specific
measurement) or more, else `hard-filter`. Pass the mode explicitly once the
cohort-size decision in [docs/PIPELINE_STATUS.md](docs/PIPELINE_STATUS.md) is
settled, rather than relying on the heuristic.

### Submit as a SLURM job

```bash
bash bin/supp/VQSR.seq_batch-slurmer.sh [options] <joint_vcf> <output_path> [mode] [cohort_name]
```

Unlike Steps 04/05, Step 06 is per-cohort rather than per-chromosome, so this
launcher submits a single job (no scatter/gather). See `--help` for
`--cpus`/`--mem`/`--hours`, and `report <manifest_file>` for post-run `sacct`
timings. Its `--mem`/`--hours` defaults (16G/24h) are unmeasured placeholders,
same caveat as the standalone script's STATUS block.

### Output

| Path | Description |
|------|-------------|
| `<output_path>/<cohort>.filtered.vcf.gz` | final callset — `FILTER` is annotated (`PASS` or a reason), records are not removed |
| `<output_path>/vqsr_work/` or `<output_path>/hardfilter_work/` | per-mode intermediates (recal/tranches files, or per-type selected/filtered VCFs) |

Extract `PASS`-only variants downstream with `bcftools view -f PASS`.

---

## Scratch Storage (`/scratch`) — what it is and isn't for

The single-sample launcher (`bin/supp/PIPELINE.single_sample.sh`) has a `USE_SCRATCH` toggle (default `true` on FENIX). When enabled it:

- reads **inputs directly from NFS** (`/mnt/data`),
- points **`TMPDIR` and all intermediate/output files at `/scratch`**,
- copies only the **final outputs back** to the sample directory, then wipes scratch.

### Benchmark finding (Sept 2026)

A controlled benchmark of Step 02 (the most I/O-heavy step — MarkDuplicates + BQSR), 3 iterations per arm, scratch vs. direct-NFS, run in parallel on the same day:

| Arm | copy-in | run (compute) | copy-out | **total** |
|-----|---------|---------------|----------|-----------|
| NFS | — | 15:41:58 | — | **15:41:58** |
| scratch | 0:02:29 | 15:36:02 | 0:02:48 | **15:41:19** |

**Scratch gave essentially no wall-time benefit (~39 s out of ~15h42m, within noise).** An earlier single-run comparison that *appeared* to show a 5.5h gain was node/scheduling variance, not scratch.

**Why:** these GATK steps are **CPU-bound with sequential, streaming I/O**, not bandwidth- or IOPS-bound. They demand only ~1–3 MB/s from storage, which any filesystem supplies trivially; sequential reads are also NFS's best case (readahead + page cache). The bottleneck is the single-threaded CPU work, which is identical on either filesystem. Staging a large read-once BAM to scratch is therefore pure overhead.

### So why still use scratch?

Its value here is **storage hygiene and cluster citizenship, not speed**:

- Bulky transient intermediates (e.g. `*.rmdup.bam`) never touch the group's **persistent quota** on `/mnt/data` — they live and die on the auto-wiped scratch volume.
- Routing **`TMPDIR`** to scratch offloads tool temp/spill traffic (the genuinely IOPS-heavy part) from shared NFS.

### Guidance for future cluster scripts

1. **Don't stage large, read-once inputs to scratch** — read them directly from NFS. Copy-in only pays off for random-access, multi-pass, or many-small-file workloads.
2. **Always point `TMPDIR` (and tool tmp flags: `samtools sort -T`, `gatk --tmp-dir`, Spark tmp) at scratch** — cheap, and it captures exactly the IOPS/spill traffic scratch is good for.
3. **Write intermediates + outputs to scratch, copy only finals back** — the quota win, speed-neutral.
4. **Don't expect wall-time speedups for sequential GATK steps** — judge scratch on quota and shared-resource health, not job runtime.

To reproduce or extend the benchmark: `bin/supp/S02_scratch_benchmark.sh [sample_dir] [n_iters]`.

---

## Monitoring Jobs

```bash
squeue -u $USER              # all your jobs
squeue -u $USER | grep BWAMAP
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

For batches: use the `*_batch-slurmer.sh` wrappers in `bin/supp/` in sequence, waiting for each stage to complete before submitting the next.
