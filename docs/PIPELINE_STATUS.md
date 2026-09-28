# Pipeline Status Report

**Project:** Joint Variant Calling Pipeline [FENIX] — GATK4 germline short-variant discovery (hg38 WGS)
**Date:** 2026-09-28 (previous: 2026-08-31)
**Reference workflow:** [GATK4 Best Practices — Germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

---

## 1. Step status at a glance

| Step | Script | Status | Runs on | Notes |
|------|--------|--------|---------|-------|
| 01 | `bin/01_bwa_map_fastq_reads.sh` | ✅ Functional | per read-group | bwa mem, RG embedded, coordinate-sorted BAM |
| 02 | `bin/02_gatk_bam_qc_workflow.sh` | ✅ Functional | per sample | MarkDuplicates (Picard) → MQ filter → BQSR → mosdepth |
| 02a | `bin/supp/02a_bqsr_evaluate.sh` | ✅ Functional | per sample | retroactive BQSR covariate plots |
| 03a | `bin/03_gatk_haplotype_caller.sh` | ✅ Functional | per sample, **per chromosome** | 25-way scatter (launcher default) → `chrom_gvcf/*.raw_vars.<CHR>.g.vcf.gz` + gathered canon-chr GVCF; validated on Batch16 |
| 03b | `bin/03_glimpse2_imputation.sh` | ⏸ Paused | per sample | low-coverage imputation; ref chunks not generated; **not** a GVCF producer |
| **04** | **`bin/04_gatk_GenomicsDB_import.sh`** | ✅ **Validated on FENIX** | **per cohort, per chromosome** | **JAGUAR chr22, 93 samples, 3-wave incremental (test01, 2026-09-01); chr1 timing (test02); SLE cohort (test03). All-chromosome parallel launcher exists for JAGUAR only (`supp/JAGUAR.S04_run_all_chroms.sh`)** |
| 05 | `bin/05_gatk_GenotypeGVCFs.sh` | ⛔ Stub (empty) | per cohort, per chromosome | on hold until S04 proven |
| 06 | `bin/06_gatk_vqsr.sh` | ⛔ Stub (empty) | per cohort | on hold; **config gap — VQSR resources missing** |
| — | `bin/PIPELINE.single_sample_01-03.sh` | ✅ Functional | per sample | chained S01→S03 SLURM launcher: resume, S03 scatter, verified scratch copy-back, free-space guard, node exclusion |
| — | `bin/supp/BATCH.submit_S01-S03.sh` | ✅ Functional | per batch dir | loops the launcher over every sample subdir |
| — | `bin/supp/SUMMARY.sample_stats_S01-S03.sh` | 🟡 New (2026-09-28) | per batch dir | per-sample TSV of stage / sizes / times / QC / depth; optional GVCF depth. Tested on synthetic data only |
| — | `bin/supp/run_pipeline.sh` | ⛔ Stub (2 lines) | — | end-to-end wrapper (superseded in practice by the launcher) |

Legend: ✅ working · 🟡 built, not yet run on real cluster data · ⏸ deliberately paused · ⛔ not implemented

---

## 2. Data flow and I/O contracts

```
FASTQ (per sample dir)
   │  01_bwa_map_fastq_reads.sh
   ▼
<SAMPLE>.sort.bam (+ .bai)                              [per read-group]
   │  02_gatk_bam_qc_workflow.sh
   ▼
<SAMPLE>.rmdup.mqfilt.bqsr.bam (+ .bai)                 [analysis-ready BAM]
   │  03_gatk_haplotype_caller.sh
   ▼
<SAMPLE>.raw_variants.canon_chr.g.vcf.gz                [per-sample GVCF]
chrom_gvcf/<SAMPLE>.raw_vars.<CHR>.g.vcf.gz (+ .tbi)    [per-chromosome GVCFs]  ◄── input to S04
   │  04_gatk_GenomicsDB_import.sh   (cohort-level; one workspace per chromosome)
   ▼
<out>/genomicsdb/<CHR>/                                 [GenomicsDB workspace per chromosome]
   │  05_gatk_GenotypeGVCFs.sh       (gendb://<CHR>, per chromosome, then gather)   ── NOT BUILT
   ▼
<cohort>.joint.vcf.gz                                   [joint-genotyped multi-sample VCF]
   │  06_gatk_vqsr.sh                (SNP + INDEL passes, or hard-filter)          ── NOT BUILT
   ▼
<cohort>.filtered.vcf.gz                                [final callset]
```

### S04 input contract (implemented)

- **Sample map** — a 2-column TSV, one line per sample:
  `<sample_label> <TAB> <path to that sample's chrom_gvcf/ directory>`
  Blank lines and `#` comments ignored. The real sample name written into the
  database is read from each GVCF header (`bcftools query -l`); a label/header
  mismatch produces a warning, not an error.
- **Output path** — workspaces are created at `<output_path>/genomicsdb/<CHR>`.
- **Chromosome selector** — `chr1`..`chr22` | `chrX` | `chrY` | `chrM` | `autosomes` | `all`.
- **Action** — `create` (default) or `update` (adds new samples to existing
  per-chromosome workspaces via `--genomicsdb-update-workspace-path`).

### S04 output contract

- `<output_path>/genomicsdb/<CHR>/` — a GenomicsDB (TileDB) workspace containing
  `callset.json`, `vidmap.json`, `vcfheader.vcf`, `__tiledb_workspace.tdb`, and
  the fragment directory `<CHR>$1$<len>`.
- Consumed by Step 05 as `gendb://<output_path>/genomicsdb/<CHR>`.

---

## 3. What blocks the next steps

| Blocker | Affects | Owner action |
|---------|---------|--------------|
| No **generic** per-chromosome SLURM launcher for S04 | S04 on new cohorts | generalise `supp/JAGUAR.S04_run_all_chroms.sh` (hard-coded paths) into `GENDBI.seq_batch-slurmer.sh` |
| VQSR resource files absent from `config/config.yaml` | S05→S06 | add `ref_hapmap` / `ref_omni` / `ref_1kg_snp` / `ref_mills`; stage files on FENIX |
| Cohort size for VQSR vs hard-filter undecided | S06 design | confirm expected N; small early cohorts need hard-filtering path |
| S03b produces phased VCF, not GVCF | pipeline diagram accuracy | when un-pausing 03b, add a `bcftools merge` path — do **not** route through S04/S05 |

---

## 3b. Operational changes since the 2026-08-31 report

- **S03 per-chromosome scatter** (`SCATTER_S03=true`): 25-task array plus a gather job, about 2 h per sample instead of ~15 h. A resubmit only re-runs missing chromosomes. Validated on Batch16.
- **Scratch hardening:** the scratch root now comes from the running user's primary group (`scratch_base` in config), so guests work on any node. Each copy-back is checked (size + `samtools quickcheck` / `gzip -t`), and scratch is kept if a check fails. The EXIT/TERM/HUP trap wipes scratch, and `ulimit -c 0` keeps JVM crash dumps out of the repo.
- **Memory pinning:** `-Xmx` is set explicitly for MarkDuplicates, BQSR, HaplotypeCaller and repair.sh. Without it these tools size their heap to the node's RAM rather than the SLURM limit.
- **Walltimes:** every step is at 48 h (gather 24 h), about 2× the worst case observed; `defq` has no time cap.
- **Launcher resume fixes:** legacy `*.sorted.bam` dirs are recognised as S01-done, and a skipped S02 no longer leaves S03 waiting on an unrelated S01 job.
- **Config guards:** `sbatch_exclude` (node09 re-enabled 2026-09-22; node15 excluded 2026-09-11 for bad RAM and re-enabled 2026-09-28 after the IT fix, now `""`) and `min_free_gb` (the launcher aborts if the output FS is low).
- **Batch tooling:** `BATCH.submit_S01-S03.sh` submits a whole batch, and `SUMMARY.sample_stats_S01-S03.sh` reports on it. The BWAMAP slurmer takes an optional `output_dir`.
- **Cohort one-offs:** `LUPUS25.setup_and_submit_S01-S03.sh` (39 pending Lupus2025 samples, B16–B20), `SLEmx_legacy.submit_S03.sh` (32 legacy SLE samples, L23/L24 lane merge), `JAGUAR.S04_run_all_chroms.sh` (full-cohort S04).
- **S04:** `--consolidate` is opt-in and never used on `update` (chr1 test02: 13.5 h). Resources were right-sized from `seff`, and a `GENDBI_PROFILE` sampler was added.

---

## 4. Best-practices conformance (steps 01–04)

| GATK BP expectation | This pipeline | Verdict |
|---------------------|---------------|---------|
| Map with bwa mem, mark duplicates, BQSR | S01 + S02 | ✅ |
| Per-sample GVCF via HaplotypeCaller `-ERC GVCF` | S03a | ✅ |
| Consolidate GVCFs with GenomicsDBImport, one workspace **per interval** | S04 — one per chromosome | ✅ |
| Use `--sample-name-map` for scalable cohorts | S04 builds it per chromosome | ✅ |
| Keep GenomicsDB off networked filesystems during import | S04 builds on `/scratch`, copies back | ✅ |
| Incremental import for growing cohorts | S04 `update` action | ✅ |
| Joint genotyping with GenotypeGVCFs (`gendb://`) | S05 — not built | ⛔ pending |
| Filter with VQSR (or hard-filter for small cohorts) | S06 — not built | ⛔ pending |

---

## 5. Local test evidence for S04 (2026-08-31)

Ran against synthetic 3-sample data (chr22 + chrM, tiny reference), GATK 4.6.2.0:

| Scenario | Result |
|----------|--------|
| `create` chr22 (2 samples) | ✅ workspace built, copied back, finisher OK |
| `update` chr22 (+1 sample) | ✅ 3 samples present in DB (verified via `SelectVariants`) |
| `create` chrM (data uses `chrM`) | ✅ |
| `create` chrM where data uses `chrMT` | ✅ warns, builds `genomicsdb/chrMT`, finisher finds it |
| rerun `create` on existing workspace | ✅ `[SKIP]` |
| `update` with no existing workspace | ✅ errors with "build it first" hint, exit 1 |
| invalid selector `chr23` | ✅ error, exit 1 |
| sample missing that chromosome's GVCF | ✅ error names the sample, exit 1 |
| duplicate sample in cohort | ✅ error, exit 1 |
| missing / empty sample map | ✅ error, exit 1 |
| truncated GVCF (no BGZF EOF marker) | ✅ `[!] WARNING` then GATK fails "Premature end of file"; prior waves' workspace intact |
| `GENDBI_STRICT_GVCF=true` on truncated GVCF | ✅ aborts before calling GATK, exit 1 |
| 3-wave create→update→update via `test/jaguar/` helpers (synthetic) | ✅ DB grows 2→4; wave 3 fails only on the bad sample |
| `run_jaguar_waves.sh --dry-run` | ✅ correct chained sbatch cmds + manifest |

Since exercised on FENIX (2026-09): real module load, `/scratch` build +
copy-back, the real 93-sample JAGUAR cohort (chr22, 3 waves) and chr1 timing.
See `test/jaguar/README.md`. Still not exercised: the `autosomes`/`all` serial loop.
