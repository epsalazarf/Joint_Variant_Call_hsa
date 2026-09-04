# Step 04 (GenomicsDBImport) test harness

Wave-based validation of `bin/04_gatk_GenomicsDB_import.sh` on FENIX. Started
with the JAGUAR cohort (93 WGS samples), now also validated against the SLE
cohort (much larger and structurally messier — see test03). Any cohort works:
sample discovery **scans** the given `--cohort` root recursively for
`.../chrom_gvcf/*.<chrom>.g.vcf.gz`, so it doesn't care about directory depth
or GVCF filename conventions (see "Cohort assumptions" below).

Everything is parameterised — `--chrom chrN`, `--waves N`, `--tag NAME`. `--tag`
namespaces the wave-map / job / log / manifest filenames so more than one
cohort's test artifacts can live in this directory at once without colliding
(default `jaguar`; use e.g. `--tag sle` for a different cohort). Output for
every run goes to `/mnt/data/amedina/$USER/JVCdev/test_gendbi_jag22/genomicsdb/<chrN>/`
(group storage — **never** the repo / `$HOME`; per-chromosome workspaces from
different test runs sit side by side there — mind `--output` if you don't want
that).

## Files

| File | Purpose |
|------|---------|
| `make_wave_maps.sh` | scan a cohort for `*.<chrom>.g.vcf.gz` and write N wave sample-maps |
| `check_setup.sh` | read-only pre-flight sweep — run before launching |
| `run_jaguar_waves.sh` | submit the N chained SLURM jobs (`--dry-run` / `report` modes) |
| `JAGUAR_JVC_chr22_list.txt` | the original JAGUAR chr22 `ls -l` listing (reference only; scan mode doesn't need it) |

Generated at run time, git-ignored: `<tag>_<chrom>_wave*.sample_map.tsv`, `logs/`.

## Cohort assumptions (what the scan needs, and what it doesn't)

- Needs: somewhere under `--cohort`, a directory named `chrom_gvcf/` per sample,
  containing exactly one `*.<chrom>.g.vcf.gz` (+ `.tbi`) for that sample.
- Does **not** need: a fixed depth (JAGUAR: `bam/<S>/chrom_gvcf`; SLE:
  `BAMQC/BatchNN/<S>/chrom_gvcf`), or the GVCF filename to start with the
  sample/directory name (SLE prefixes vary — `L46_1`, `Q005_D_2`,
  `Q00a_D_1`, ...). The **label** used in the sample map is always the
  `chrom_gvcf/`'s parent directory name; Step 04 re-derives the real sample
  name from the GVCF header regardless, so a label/header mismatch is only
  ever a warning, never an error.
- A sample directory with no GVCF for the requested chromosome yet (pipeline
  still running) is silently skipped, not an error — the scan only finds what
  exists.
- Two different `chrom_gvcf/` dirs with the **same** parent directory name
  (e.g. reused across two batches) is a hard error — check for a genuine ID
  collision before renaming/merging.

## Run it (on FENIX)

```bash
cd /path/to/repo && git pull

CHR=chr1        # or chr22, chrX, ...
WAVES=2
TAG=jaguar      # or e.g. sle — keeps artifacts from different cohorts apart
COHORT=/mnt/data/amedina/mramirezc/JAGUAR_JVC   # cohort root to scan

# 1. wave maps (scans the cohort directly)
bash test/jaguar/make_wave_maps.sh --chrom $CHR --waves $WAVES --tag $TAG --cohort $COHORT

# 2. pre-flight — resolve any FAIL (a known tiny/truncated GVCF is an OK WARN)
bash test/jaguar/check_setup.sh --chrom $CHR --tag $TAG

# 3. see the exact sbatch commands
bash test/jaguar/run_jaguar_waves.sh --dry-run --chrom $CHR --waves $WAVES --tag $TAG $COHORT

# 4. submit  (add --cpus / --mem / --hours to push the envelope)
bash test/jaguar/run_jaguar_waves.sh --chrom $CHR --waves $WAVES --tag $TAG $COHORT

# monitor / timings
squeue -u $USER | grep JVC-GDBI-$TAG-$CHR
bash test/jaguar/run_jaguar_waves.sh report            # latest manifest across ALL tags
bash test/jaguar/run_jaguar_waves.sh report logs/<tag>_run_<chrom>_<timestamp>.manifest  # a specific one
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
- **Read check passed** — `SelectVariants` on `chr1:1-500000` returned 169,233
  records for the 93-sample cohort. The `No valid combination operation found
  for INFO field DB / InbreedingCoeff / MLEAC / MLEAF` lines are **expected** —
  those are per-sample annotations GenomicsDB can't merge; GenotypeGVCFs
  recomputes site-level stats. The read was I/O-bound (~10× wall/CPU in the
  GenomicsDB iterator, ~8 min warm-up) → **Step 05 should stage the workspace
  to `/scratch` before reading it** (see review notes).

### test03 — SLE cohort, chr22 (heterogeneous cohort / harness-generality test)

Purpose: stress the harness itself, not just S04 — the SLE cohort is larger and
structurally messier than JAGUAR, so it's a good check that the scripts don't
silently assume JAGUAR's layout.

Cohort snapshot (`SLEmx-b38_list.txt`, 2026-09): ~127 sample directories across
13 batches, but only **52** currently have a `chr22` GVCF+index (Batches
05/06/08/09/12/15 haven't reached Step 03a yet for chr22 — expected to be
skipped by the scan, not treated as an error). Confirmed quirks the harness now
handles:

- **Different directory depth** than JAGUAR: `BAMQC/BatchNN/<sample>/chrom_gvcf/`
  instead of `bam/<sample>/chrom_gvcf/`. `make_wave_maps.sh`'s scan is now fully
  recursive (`find --cohort -path '*/chrom_gvcf/*.<chrom>.g.vcf.gz'`) instead of
  assuming a `bam/` root.
- **GVCF filename ≠ sample/directory name**: `L46/chrom_gvcf/L46_1.raw_vars...`,
  `Q005/chrom_gvcf/Q005_D_2.raw_vars...`, prefixes also seen: `_1`, `_2`, `_5`,
  `_D_1`, `_D_2`, `_D_3`. `make_wave_maps.sh` and `check_setup.sh` previously
  reconstructed the expected filename as `<label>.raw_vars.<chrom>.g.vcf.gz` —
  **wrong** for this cohort. Fixed to glob `*.<chrom>.g.vcf.gz` inside the
  directory instead (matches how Step 04 itself already resolved the file).
- Inconsistent sample-ID punctuation even within one batch: `Q005` vs `Q_034`
  vs `Q00a`. Harmless — the directory name is just used as an opaque label.
- One sample dir (`Q_034`) is empty (no GVCF yet) — scan skips it, no crash.
- No duplicate sample labels across batches, no ambiguous (>1 file) `chrom_gvcf/`
  dirs, no missing `.tbi` — checked directly against the file listing.
- Coverage is more heterogeneous than JAGUAR: chr22 GVCF sizes range
  10.1–53.3 MB (JAGUAR was a tighter 34–61 MB band) — a real test of the
  tiny-file heuristic's threshold on a noisier size distribution; no false
  positive at the current 50%-of-median cutoff.

Added `--tag` (default `jaguar`) to `make_wave_maps.sh` / `check_setup.sh` /
`run_jaguar_waves.sh` so this cohort's maps/jobs/logs/manifest
(`sle_chr22_wave*.sample_map.tsv`, `JVC-GDBI-sle-chr22-w*`, etc.) don't collide
with JAGUAR's in this same directory.

Run (needs the SLE cohort root path on FENIX — ask the maintainer if unset):

```bash
COHORT=<absolute path to the dir containing SLEmx-b38/>
bash test/jaguar/make_wave_maps.sh --chrom chr22 --waves 3 --tag sle --cohort "$COHORT/SLEmx-b38/BAMQC"
bash test/jaguar/check_setup.sh --chrom chr22 --tag sle
bash test/jaguar/run_jaguar_waves.sh --dry-run --chrom chr22 --waves 3 --tag sle "$COHORT/SLEmx-b38/BAMQC"
bash test/jaguar/run_jaguar_waves.sh --chrom chr22 --waves 3 --tag sle "$COHORT/SLEmx-b38/BAMQC"
```

52 samples / 3 waves ≈ 18 samples per wave — smaller per-wave than JAGUAR
(31) despite the larger cohort, since ~75 sample dirs aren't chr22-ready yet.
Results TBD.

### First run (2026-08-31) — failed on quota, fixed

Wave 1's import **succeeded** (72 min on scratch) but the copy-back went to a
`$HOME` path (hard ~5 GB quota) → `Disk quota exceeded`, and the old trap wiped
scratch. Fixed: output defaults to group storage; `run_jaguar_waves.sh`,
`check_setup.sh`, and S04 all refuse a `$HOME` output + write-test the
destination; S04 now **preserves the scratch workspace** if the import
succeeded but copy-back fails, and prints the `cp` + re-run to recover.

The JDK-25 warnings in the logs (`System::load`, `sun.misc.Unsafe`) are harmless
— GATK 4.6 targets JDK 17.
