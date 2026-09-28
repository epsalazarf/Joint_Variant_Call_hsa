# Pipeline Instructions

Step-by-step usage guide for the Joint Variant Calling Pipeline.

> For a high-level overview, see [README.md](README.md).

---

## Script Layout

```
bin/
├── PIPELINE.single_sample_01-03.sh # Recommended: S01→S02→S03 chained SLURM launcher (one sample)
├── 01_bwa_map_fastq_reads.sh       # Step 01 — FASTQ alignment
├── 02_gatk_bam_qc_workflow.sh      # Step 02 — BAM QC + BQSR
├── 03_gatk_haplotype_caller.sh     # Step 03a — HaplotypeCaller (whole-genome or one chromosome)
├── 03_glimpse2_imputation.sh       # Step 03b — GLIMPSE2 imputation (low-coverage, paused)
├── 04_gatk_GenomicsDB_import.sh    # Step 04 — GenomicsDB import (cohort-level)
├── 05_gatk_GenotypeGVCFs.sh        # Step 05 — Joint genotyping (stub)
├── 06_gatk_vqsr.sh                 # Step 06 — VQSR filtering (stub)
└── supp/
    ├── BATCH.submit_S01-S03.sh             # Run the launcher on every sample of a batch dir
    ├── SUMMARY.sample_stats_S01-S03.sh     # Per-sample TSV: stage, sizes, times, QC, depth
    ├── 00_scan_fastq_pairs.sh              # Inspect FASTQ pairs before mapping
    ├── 00_glimpse2_ref_panel_prep.sh       # Prepare GLIMPSE2 binary reference panel
    ├── 02a_bqsr_evaluate.sh                # Retroactive BQSR covariate plots
    ├── 03-s4_gvcf_chrom_split.sh           # Split GVCFs by chromosome
    ├── BWAMAP.seq_batch-slurmer.sh         # Per-step batch launcher for Step 01
    ├── BAMQC.seq_batch-slurmer.sh          # Per-step batch launcher for Step 02
    ├── HAPCALL.seq_batch-slurmer.sh        # Per-step batch launcher for Step 03a (whole-genome, legacy)
    ├── S02_scratch_benchmark.sh            # Scratch-vs-NFS benchmark for Step 02
    ├── JAGUAR.S04_run_all_chroms.sh        # One-shot: JAGUAR full-cohort Step 04, all chromosomes
    ├── LUPUS25.setup_and_submit_S01-S03.sh # One-shot: stage + submit pending Lupus2025 samples
    ├── SLEmx_legacy.submit_S03.sh          # One-shot: S03 for the legacy SLE cohort (mbravog)
    └── run_pipeline.sh                     # Full pipeline wrapper (stub)
```

The `JAGUAR.*`, `LUPUS25.*` and `SLEmx_legacy.*` scripts are cohort-specific one-offs with hard-coded paths. Keep them for reference, and copy one as a template rather than reusing it on another cohort.

---

## Prerequisites

1. **Configure paths** in `config/config.yaml` — set `ref_gnm`, `ref_vars`, and any other keys for your environment (`remote:` for FENIX, `local:` for workstation).
2. **Load modules** (FENIX only): scripts call `module load` automatically based on the `modules` key in the config.
3. All scripts must be run from the repo root or with explicit paths.
4. On FENIX, keep outputs on group storage (`/mnt/data/...`), never under `$HOME` (hard 5 GB quota).

---

## Recommended: One-Command Launchers (Steps 01→03)

For new samples, use the launchers instead of running the steps by hand. They run everything in the sample directory, so no hand-off between steps is needed.

```bash
# One sample (dir named after the sample ID, holding its FASTQs)
bash bin/PIPELINE.single_sample_01-03.sh /path/to/SAMPLE_ID/

# Every sample subdirectory of a batch
bash bin/supp/BATCH.submit_S01-S03.sh /path/to/BatchNN/
```

| Behaviour | Details |
|-----------|---------|
| Chaining | S01 → S02 → S03 via `--dependency=afterok` |
| Resume | A step is skipped when its final output exists (S01 `*.sort.bam`/`*.sorted.bam`, S02 `*.rmdup.mqfilt.bqsr.bam`, S03 all 25 `chrom_gvcf/` + `canon_chr` GVCF). FASTQs are only required if S01 still has to run. |
| S03 scatter | `SCATTER_S03=true` (default): 25-task array (one chromosome each, `%ARRAY_CONC` concurrent, default 6), then a gather job that builds the `canon_chr` GVCF. Resubmitting re-runs only the missing chromosomes. |
| Scratch | `USE_SCRATCH=true` (default) — see [Scratch Storage](#scratch-storage-scratch--what-it-is-and-isnt-for) |
| Guards | Aborts if the output FS has < `min_free_gb` free; `sbatch_exclude` nodes are avoided (both in `config/config.yaml`) |
| Logs | `SAMPLE_ID/log/SAMPLE-S0N-<epoch>.<jobid>.log` (S03 array: `SAMPLE-S03-<epoch>.<arrayjob>_<task>.log`, gather: `SAMPLE-S03g-...`) |

Resources requested per job:

| Job | CPUs | Memory | Walltime |
|-----|------|--------|----------|
| S01 | 8 | 20 GB | 48 h |
| S02 | 4 | 28 GB | 48 h |
| S03 per chromosome (scatter) | 2 | 13 GB | 48 h |
| S03 gather | 2 | 8 GB | 24 h |
| S03 whole-genome (`SCATTER_S03=false`) | 4 | 32 GB | 48 h |

Heaps are pinned inside the scripts (e.g. MarkDuplicates `-Xmx16g`, BQSR `-Xmx8g`, repair.sh `-Xmx8g`). Tools otherwise size their heap to the node's RAM rather than the SLURM limit.

**Re-running a step on purpose:** "resume" never redoes finished work. To force a step, remove its final output from the sample dir first, then relaunch.

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
bash bin/supp/BWAMAP.seq_batch-slurmer.sh bin/01_bwa_map_fastq_reads.sh /path/to/fastq_samples/ /path/to/bams/
```

Auto-discovers sample subdirectories and submits one job per sample (parallel). Wait for all jobs to finish before running Step 02.

**Output location (3rd argument, optional):**
- Given → BAMs are written to `/path/to/bams/<SAMPLE>/` (created if missing). Pass the same `/path/to/bams/` to the Step 02 batch slurmer.
- Omitted → BAMs are written **inside each FASTQ sample directory** (legacy behaviour). In that case point the Step 02 slurmer at the FASTQ parent directory, or symlink the `*.sort.bam` files into your BAM directory first — the Step 02 slurmer only searches the directory it is given.

> For new samples, the one-command launchers (`PIPELINE.single_sample_01-03.sh` / `BATCH.submit_S01-S03.sh`, see top of this file / README) are simpler: they chain 01→02→03 in the sample directory and avoid this hand-off entirely.

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
bash bin/supp/BAMQC.seq_batch-slurmer.sh bin/02_gatk_bam_qc_workflow.sh /path/to/bams/
```

Recursively finds every `*.sort.bam` / `*.sorted.bam` under `/path/to/bams/` (use the Step 01 output directory), groups them by sample ID (text before the first `_`, so per-lane BAMs are merged), symlinks them into `/path/to/bams/<SAMPLE>/`, and submits one independent job per sample (8 CPUs / 32 GB, parallel). Step 02 outputs land in `/path/to/bams/<SAMPLE>/`.

The discovered list is saved as `<dir>.bamqc_input_files.txt` (tab-separated: `sample_name  sample_file(s)  readgroup_string`). If that file already exists it is reused instead of re-scanning — edit it to add an RG string for legacy BAMs, or delete it to force a fresh scan.

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
| `SAMPLE.rmdup.mqfilt.bqsr_table.txt` | BQSR recalibration table |
| `SAMPLE.rmb.mosdepth.*` | Coverage summary (`--fast-mode`, no per-base) |
| `SAMPLE.alignment_metrics.txt` | Picard alignment summary (if `RUN_METRICS=true`) |
| `SAMPLE.insert_size_metrics.txt` + `SAMPLE.insert_size_histogram.pdf` | Insert-size metrics (if `RUN_METRICS=true`) |

Files are named after the sample (`SAMPLE.rmdup.*`), not the input BAM (`SAMPLE.sort.rmdup.*` before v1.2). Output from before v1.2 isn't recognised by the skip logic; rename those files or rerun.

### Toggles (edit inside the script)

| Variable | Default | Effect |
|----------|---------|--------|
| `HOUSEKEEP` | `true` | Remove intermediate BAMs on success |
| `MQ_FILTER` | `false` | Apply MQ ≥ 30 filter (originally for aDNA) |
| `RUN_METRICS` | `true` | Collect alignment + insert-size metrics |
| `BQSR_EVAL` | `false` | Run post-BQSR evaluation inline (use 02a instead) |
| `BQSR_COV` | `true` | AnalyzeCovariates plots (forced off on FENIX; needs `BQSR_EVAL=true`) |
| `REMOVE_DUPS` | `false` | Remove duplicates instead of marking |

---

## Step 02a — Retroactive BQSR Evaluation (optional)

Generates before/after BQSR covariate plots without re-running the full Step 02. Use when Step 02 was run with `BQSR_EVAL=false`.

**Prerequisite:** the pre-BQSR table `SAMPLE.rmdup.mqfilt.bqsr_table.txt` (written by Step 02) in the output directory. The table is looked up next to the BAM as `<BAM name minus .bam>_table.txt`, so legacy runs (`L46_1.sorted.rmdup.mqfilt.bqsr_table.txt`) and older `SAMPLE.bqsr_table.txt` outputs are picked up automatically.

```bash
bash bin/supp/02a_bqsr_evaluate.sh <sample.rmdup.mqfilt.bqsr.bam> [output_path]
```

### Output

| File | Description |
|------|-------------|
| `SAMPLE.rmdup.mqfilt.bqsr_table_recal.txt` | Post-BQSR recalibration table (reused if Step 02 already wrote it) |
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

> **Validated on FENIX** — JAGUAR chr22, 93 samples, 3-wave incremental import.
> Design notes, resourcing and known limitations:
> [docs/S04_GenomicsDBImport_design.md](docs/S04_GenomicsDBImport_design.md).
> Test harness (any chromosome / cohort): [test/jaguar/README.md](test/jaguar/README.md).
> Coordinate with the maintainer before running it; it is a cohort-level step, run once after **all** samples finish Step 03.

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

`autosomes` / `all` process every chromosome **serially in one process**. That's fine
for tests and small cohorts, but a full cohort takes days that way. To run
all chromosomes in parallel, submit one SLURM job per chromosome.
`bin/supp/JAGUAR.S04_run_all_chroms.sh` does this for JAGUAR (arrays by
chromosome-size class, `--dry-run` / `--status` modes). Copy it for a new
cohort; a generic launcher is not built yet.

### 3. Adding samples later (waves)

```bash
# first wave
bash bin/04_gatk_GenomicsDB_import.sh wave1.sample_map.tsv <out> chr22 create
# later wave — only the NEW samples in this map
bash bin/04_gatk_GenomicsDB_import.sh wave2.sample_map.tsv <out> chr22 update
```

Never consolidate on `update`: it took 13.5 h on chr1 in testing. `--consolidate` is opt-in and off by default.

### Output

| Path | Description |
|------|-------------|
| `<output_path>/genomicsdb/<CHR>/` | GenomicsDB workspace (input for Step 05 as `gendb://...`) |

The workspace is built on `/scratch` (when available) and copied back on
success.

---

## Steps 05–06 — Joint Genotyping and Filtering (stubs)

Not yet implemented; on hold for now (Step 04 is validated, so this is the next step to build).

| Step | Script | Purpose |
|------|--------|---------|
| 05 | `bin/05_gatk_GenotypeGVCFs.sh` | Joint genotyping across all samples (`gendb://` per chromosome) |
| 06 | `bin/06_gatk_vqsr.sh` | VQSR (or hard-filter) of the joint VCF — **needs new config keys for VQSR resources** |

Do not use these for production runs. Watch the pipeline overview table in
[README.md](README.md) and [docs/PIPELINE_STATUS.md](docs/PIPELINE_STATUS.md).

---

## Scratch Storage (`/scratch`) — what it is and isn't for

The single-sample launcher (`bin/PIPELINE.single_sample_01-03.sh`) has a `USE_SCRATCH` toggle (default `true` on FENIX). When enabled it:

- reads **inputs directly from NFS** (`/mnt/data`),
- points **`TMPDIR` and all intermediate/output files at `/scratch`**, under `<scratch_base>/<primary group of the running user>/$USER/job_<jobid>`. That works for group members and guests alike, whatever node the job lands on.
- copies only the **final outputs back** to the sample directory, **verifying each one** (size match, then `samtools quickcheck` for BAMs or `gzip -t` for GVCFs), then wipes scratch,
- **preserves** scratch if any check fails, so the data can be recovered; the log prints the path. The scratch dir is also removed if the job is killed (SIGTERM/HUP).

If a node's scratch mount is broken, jobs on it fail immediately with `scratch unavailable on <node>`. Add that node to `sbatch_exclude` in `config/config.yaml` and resubmit.

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
squeue -u $USER                  # all your jobs
squeue -u $USER | grep SAMPLE_ID # one sample's chain (launcher job names: SAMPLE-S01-..., SAMPLE-S03g-...)
squeue -u $USER | grep -E 'BWAMAP|BAMQC|HAPCALL'   # per-step slurmers
```

## Logs

- Launcher jobs: `SAMPLE_ID/log/` (see [the launcher table](#recommended-one-command-launchers-steps-0103)).
- Per-step slurmers: `%x.%j.log` in the sample/output dir (BWAMAP, BAMQC) or the submission dir (HAPCALL).

Each step script prints `[*]` step headers, `[&]  Step time:` lines, a final `[&]  Total time:`, and `completed successfully!` on success. Search a log for `[X]` or `<ERROR>` first.

---

## Per-sample Summary Table

Scans the outputs and logs of Steps 01–03 and writes one TSV row per sample. It only reads files, so it's safe to run while jobs are still going. It follows symlinked files and folders (legacy BAMs, staged FASTQs, read-only batches from other users) and understands both current and legacy log wording.

```bash
bash bin/supp/SUMMARY.sample_stats_S01-S03.sh <batch_or_sample_dir> [more_dirs...] [-o out.tsv]
    [--gvcf-depth[=fast|full]] [-j N]
```

Default output: `./<first_dir>.S01-S03_summary.<date>.tsv`.

| Columns | Source |
|---------|--------|
| `batch` | parent folder — tells apart the same sample ID in several batches |
| `stage` | furthest step whose final output exists (S00–S03) |
| `fastq_gib`, `sort_bam_gib`, `bqsr_bam_gib`, `chrom_gvcf_gib`, `canon_gvcf_gib` | file sizes |
| `raw_reads`, `pct_mapped`, `pct_proper_pair` | `*.sort.stats.txt` (samtools stats, summed over read groups) |
| `pct_dup` | `*-dups.txt` (Picard, pooled across libraries) |
| `median_insert` | `*.insert_size_metrics.txt` |
| `mean_depth_auto`, `mean_depth_chrX/chrY/chrM` | `*.mosdepth.summary.txt` (autosomes = Σbases / Σlength over chr1–22) |
| `s01_time`, `s02_time` | `Total time` of the most recent **successful** log for that step |
| `s03_time_max`, `s03_time_sum` | scatter: slowest chromosome (≈ wall time) and total compute |
| `*_runs` | logs found per step (> 1 = failed attempts or resumes) |
| `chrom_gvcfs` | `n/25` per-chromosome GVCFs present with index |
| `chrom_gvcf_fmt`, `canon_gvcf_fmt` | real container: `vcf`, `bcf` or `mixed(...)`. Legacy S03 ran `bcftools -O b` with a `.g.vcf.gz` name, and GATK can't read bcftools BCF, so anything other than `vcf` must be fixed before Step 04 |
| `gvcf_depth_region`, `gvcf_depth`, `gvcf_vs_mosdepth` | only with `--gvcf-depth` (see below) |

**`--gvcf-depth`**: re-derives mean depth from the GVCFs as a length-weighted mean of `FORMAT/DP`. Each reference block counts for the bases it spans (ref-block DP is GATK's block median), and each variant site counts once. That makes it valid on GVCFs, unlike `bcftools stats`, whose DP histogram counts a 5 kb block the same as one SNP.
- `fast` (default when you pass only `--gvcf-depth`): chr20 only, about a minute per sample. Use it as a quick check on finished samples.
- `full`: chr1–22, `-j` chromosomes in parallel (default 4). Use it for the final report.
- `gvcf_vs_mosdepth` compares against mosdepth over the **same** region. HaplotypeCaller's DP only counts filtered reads (MAPQ ≥ 20, downsampled), so a ratio a little below 1 is expected. Look for samples that stand out from the rest of the batch.

---

## Typical End-to-End Run (single sample, FENIX)

Recommended: `bash bin/PIPELINE.single_sample_01-03.sh /data/sample01/` (see above). To run the steps by hand instead:

```bash
# 1. Align
bash bin/01_bwa_map_fastq_reads.sh /data/sample01/ /output/bams/

# 2. QC + BQSR
bash bin/02_gatk_bam_qc_workflow.sh /output/bams/sample01.sort.bam /output/bqsr/

# 3. Variant calling
bash bin/03_gatk_haplotype_caller.sh /output/bqsr/sample01.rmdup.mqfilt.bqsr.bam /output/gvcf/
```

For batches: `bin/supp/BATCH.submit_S01-S03.sh`. Alternatively, run the per-step `*.seq_batch-slurmer.sh` wrappers in `bin/supp/` in sequence, waiting for each stage to finish before submitting the next.
