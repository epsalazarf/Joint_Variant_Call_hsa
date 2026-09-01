# Pipeline Status Report

**Project:** Joint Variant Calling Pipeline [FENIX] — GATK4 germline short-variant discovery (hg38 WGS)
**Date:** 2026-08-31
**Reference workflow:** [GATK4 Best Practices — Germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

---

## 1. Step status at a glance

| Step | Script | Status | Runs on | Notes |
|------|--------|--------|---------|-------|
| 01 | `bin/01_bwa_map_fastq_reads.sh` | ✅ Functional | per read-group | bwa mem, RG embedded, coordinate-sorted BAM |
| 02 | `bin/02_gatk_bam_qc_workflow.sh` | ✅ Functional | per sample | MarkDuplicates (Picard) → MQ filter → BQSR → mosdepth |
| 02a | `bin/supp/02a_bqsr_evaluate.sh` | ✅ Functional | per sample | retroactive BQSR covariate plots |
| 03a | `bin/03_gatk_haplotype_caller.sh` | ✅ Functional | per sample | GVCF mode → canon-chr GVCF → `chrom_gvcf/*.raw_vars.<CHR>.g.vcf.gz` |
| 03b | `bin/03_glimpse2_imputation.sh` | ⏸ Paused | per sample | low-coverage imputation; ref chunks not generated; **not** a GVCF producer |
| **04** | **`bin/04_gatk_GenomicsDB_import.sh`** | 🟡 **Implemented — needs FENIX validation** | **per cohort, per chromosome** | **new (2026-08-31); tested locally on synthetic data** |
| 05 | `bin/05_gatk_GenotypeGVCFs.sh` | ⛔ Stub (empty) | per cohort, per chromosome | on hold until S04 proven |
| 06 | `bin/06_gatk_vqsr.sh` | ⛔ Stub (empty) | per cohort | on hold; **config gap — VQSR resources missing** |
| — | `bin/supp/run_pipeline.sh` | ⛔ Stub (2 lines) | — | end-to-end wrapper |

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
| S04 not yet run on FENIX with real GVCFs | S04 → prod | **run `test/jaguar_chr22/`** — 93 samples, 3 waves (create+update+update), chr22; confirm scratch copy-back + module load + truncation flag |
| No per-chromosome SLURM launcher for S04 | S04 batch use | build `GENDBI.seq_batch-slurmer.sh` after S04 validated (S04 can loop chroms serially meanwhile) |
| VQSR resource files absent from `config/config.yaml` | S05→S06 | add `ref_hapmap` / `ref_omni` / `ref_1kg_snp` / `ref_mills`; stage files on FENIX |
| Cohort size for VQSR vs hard-filter undecided | S06 design | confirm expected N; small early cohorts need hard-filtering path |
| S03b produces phased VCF, not GVCF | pipeline diagram accuracy | when un-pausing 03b, add a `bcftools merge` path — do **not** route through S04/S05 |

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
| 3-wave create→update→update via `test/jaguar_chr22/` helpers (synthetic) | ✅ DB grows 2→4; wave 3 fails only on the bad sample |
| `run_jaguar_waves.sh --dry-run` | ✅ correct chained sbatch cmds + manifest |

Not yet exercised: real FENIX module load, `/scratch` path, the real 93-sample
cohort, `autosomes`/`all` multi-chromosome serial loop, `sacct` report mode.
