# JAGUAR — Step 04 (GenomicsDBImport) test harness

Wave-based validation of `bin/04_gatk_GenomicsDB_import.sh` on FENIX, using the
JAGUAR cohort (93 WGS samples, `/mnt/data/amedina/mramirezc/JAGUAR_JVC`). One
chromosome at a time; the cohort is split into N sequential "waves" so the
incremental (`update`) path gets exercised.

Everything is chromosome-parameterised — `--chrom chrN`. Output for every run
goes to `/mnt/data/amedina/$USER/JVCdev/test_gendbi_jag22/genomicsdb/<chrN>/`
(group storage — **never** the repo / `$HOME`; per-chromosome workspaces sit
side by side there).

## Files

| File | Purpose |
|------|---------|
| `make_wave_maps.sh` | scan the cohort for `*.<chrom>.g.vcf.gz` and write N wave sample-maps |
| `check_setup.sh` | read-only pre-flight sweep — run before launching |
| `run_jaguar_waves.sh` | submit the N chained SLURM jobs (`--dry-run` / `report` modes) |
| `JAGUAR_JVC_chr22_list.txt` | the original chr22 `ls -l` listing (reference; scan mode doesn't need it) |

Generated at run time, git-ignored: `jaguar_<chrom>_wave*.sample_map.tsv`, `logs/`.

## Run it (on FENIX)

```bash
cd /path/to/repo && git pull

CHR=chr1        # or chr22, chrX, ...
WAVES=2         # 93 samples -> 47 + 46

# 1. wave maps (scans the cohort directly)
bash test/jaguar/make_wave_maps.sh --chrom $CHR --waves $WAVES

# 2. pre-flight — resolve any FAIL (a known tiny/truncated GVCF is an OK WARN)
bash test/jaguar/check_setup.sh --chrom $CHR

# 3. see the exact sbatch commands
bash test/jaguar/run_jaguar_waves.sh --dry-run --chrom $CHR --waves $WAVES

# 4. submit  (add --cpus / --mem / --hours to push the envelope)
bash test/jaguar/run_jaguar_waves.sh --chrom $CHR --waves $WAVES --cpus 4 --mem 32G --hours 12

# monitor / timings
squeue -u $USER | grep JVC-GDBI-$CHR
bash test/jaguar/run_jaguar_waves.sh report            # after all waves finish
```

**Wave plan** (built automatically): wave 1 = `create` (+ `--consolidate`, a
~4 s no-op on a fresh import); waves 2…N = `update`, **no consolidate** — each
just appends a fragment. GenotypeGVCFs reads a multi-fragment DB fine.
Consolidation on `update` is pathologically slow (test02) — never do it.

## Runs

### test01 — chr22, 3 waves of 31  (jobs 278359-61, 2026-09-01) — PASS

93 samples imported. Findings folded into the code:

| wave | db total | wall | |
|------|----------|------|--|
| 1 create | 31 | 54 min | consolidate ~3 s |
| 2 update | 62 | 147 min | consolidate ~96 min (this run still did it) |
| 3 update | 93 | 190 min | + tiny-sample warning |

- `--consolidate` cost scales with **total** DB size (rewrites the whole array),
  not the wave → now on for `create` only (test02 showed why).
- more cores/RAM don't help (serial + I/O-bound). `seff`: CPU eff 9.53% at 8
  cores, RAM 6.68 GB. Default 2 CPU / 16 G.
- `EGAN00004696518` (24 KB, near-empty but valid) flagged by pre-build warning +
  finisher, imported anyway.
- workspace files `0700` / group `fsanchezq` — other amedina users can't read
  it; make the output dir setgid the data group (or `umask 0027`). `amgrp` in
  `~/bin` does this.
- benign `.__consolidation_lock` (0 B) left in the array dir.

### test02 — chr1, 2 waves  (jobs 279778-79, 2026-09-03) — PASS

Upper-limit run: chr1 (biggest chromosome), 47 `create` + 46 `update` → 93
samples, ~30 GB, one fragment. Both waves succeeded.

| wave | action | db total | import | consolidate | step time |
|------|--------|----------|--------|-------------|-----------|
| 1 | create | 47 | 460 min | ~4 s | 7 h 41 m |
| 2 | update | 93 | 445 min | **~13.5 h** | **21 h 01 m** |

**Findings:**

- Import scales ~linearly: ~9.7 min/sample on chr1 (vs ~1.7 on chr22).
  Predictable from `chrom_length × samples`.
- **Consolidation on `update` is the bottleneck** — 13.5 h to merge fragments
  across 93 samples of chr1, *longer than a full `create` rebuild* (~15 h). It's
  a ~4 s no-op on `create`. → **updates no longer consolidate**; `wave_consol`
  is wave 1 only. GenotypeGVCFs reads multi-fragment DBs. Rebuild with `create`
  if you ever need one fragment.
- Wave 2 used 21 h of its 24 h walltime — an update+consolidate would time out
  on a bigger cohort. Another reason to drop it.
- **Memory scales with contig size** (`seff`): chr22 6.7 GB, chr1 `create` 27 GB,
  chr1 `update`+consolidate **32 GB / 32 GB (100%)**. Native TileDB, not heap.
  → launcher now sets `--mem`/`--time` per contig class (chr1-2 64G/20h ·
  chr3-8,X 32G/12h · rest 16G/6h) + a `--batch` knob (50→25 halves import RAM).
- CPU efficiency 13–18% at 4 cores — cores useless, confirmed again.
- Storage: ~332 MB/sample on chr1 → whole-genome 93-sample DB ≈ **375–400 GB**.
- setgid on the output dir set group=`amedina` on the new files, but GenomicsDB
  still writes them `0700` — run `amgrp <workspace>` after each import so
  labmates can read it.

### First run (2026-08-31) — failed on quota, fixed

Wave 1's import **succeeded** (72 min on scratch) but the copy-back went to a
`$HOME` path (hard ~5 GB quota) → `Disk quota exceeded`, and the old trap wiped
scratch. Fixed: output defaults to group storage; `run_jaguar_waves.sh`,
`check_setup.sh`, and S04 all refuse a `$HOME` output + write-test the
destination; S04 now **preserves the scratch workspace** if the import
succeeded but copy-back fails, and prints the `cp` + re-run to recover.

The JDK-25 warnings in the logs (`System::load`, `sun.misc.Unsafe`) are harmless
— GATK 4.6 targets JDK 17.
