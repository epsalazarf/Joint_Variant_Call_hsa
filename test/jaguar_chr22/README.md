# JAGUAR chr22 — Step 04 first test

First real-data validation of `bin/04_gatk_GenomicsDB_import.sh` on FENIX.
**93 samples**, chr22 only, imported in **three waves of 31** to exercise the
incremental (`update`) path twice.

| Wave | Action | Samples | First … last | SLURM (CPU / mem / time) | consolidate |
|------|--------|---------|--------------|--------------------------|-------------|
| 1 | `create` | 31 | EGAN00004552304 … EGAN00004552334 | 2 / 16G / 3h | yes |
| 2 | `update` | 31 | EGAN00004552335 … EGAN00004696501 | 2 / 16G / 3h | **no** |
| 3 | `update` | 31 | EGAN00004696502 … EGAN00004710984 | 2 / 16G / 6h | yes (final) |

Sizing settled by test01 (results below): GenomicsDBImport for one interval is
serial + I/O-bound — 2 CPU / 16 G is plenty. The expensive variable is
`--consolidate` (rewrites the whole array); only wave 1 and the final wave pay
it now. All 2 reader-threads / batch 50.

### Known bad sample (left in on purpose)

`EGAN00004696518` (wave 3) is **24 KB** vs a ~48 MB median. From the first
FENIX run it is a **valid bgzip that is simply near-empty** (the BGZF EOF
marker is present — `check_setup.sh` reported no truncation), not a
transfer-truncated file. Two failure modes, handled differently:

| Case | Detected by | Step 04 behaviour |
|------|-------------|-------------------|
| valid but tiny (this sample) | size vs median | `[!] WARNING: ... unusually small`; **imports anyway** — near-empty sample in the DB |
| BGZF-truncated (no EOF marker) | `gvcf_looks_truncated` | `[!] WARNING: ... TRUNCATED`; GenomicsDBImport then fails `Premature end of file`; `GENDBI_STRICT_GVCF=true` aborts first |

So wave 3 will most likely **succeed** and `EGAN00004696518` lands in the DB as
a near-empty sample. Check its chr22 coverage / call count before Step 05; if
it is genuinely bad, re-run Step 03a for it and re-import with `action=update`
(GenomicsDB keeps the newer copy). The finisher lists it under
"Unusually small (but valid) GVCF(s) imported".

## Files

| File | Purpose |
|------|---------|
| `JAGUAR_JVC_chr22_list.txt` | the `ls -l` listing of the 93 chr22 GVCFs (input) |
| `make_wave_maps.sh` | split the list into 3 wave sample-maps (absolute paths) |
| `check_setup.sh` | read-only pre-flight sweep — run before launching |
| `run_jaguar_waves.sh` | submit the 3 chained SLURM jobs (`--dry-run` supported) |

Generated at run time (git-ignored): `jaguar_chr22_wave*.sample_map.tsv`, `logs/`, `out/`.

## Run it (on FENIX)

Cohort root is `/mnt/data/amedina/mramirezc/JAGUAR_JVC` (the dir containing
`bam/`) — it is the built-in default, so no path argument is needed. Pass one
only to point at a different cohort.

```bash
cd /path/to/repo

# 1. build the three wave maps
bash test/jaguar_chr22/make_wave_maps.sh

# 2. pre-flight (resolve any FAIL before continuing; the truncated-GVCF WARN is expected)
bash test/jaguar_chr22/check_setup.sh

# 3. dry-run to see the exact sbatch commands
bash test/jaguar_chr22/run_jaguar_waves.sh --dry-run

# 4. submit the chained jobs
bash test/jaguar_chr22/run_jaguar_waves.sh

# monitor
squeue -u $USER --name=JVC-GDBI-w1,JVC-GDBI-w2,JVC-GDBI-w3

# 5. after all three finish — resource-vs-runtime comparison
bash test/jaguar_chr22/run_jaguar_waves.sh report
```

Output workspace: `/mnt/data/amedina/$USER/JVCdev/test_gendbi_jag22/genomicsdb/chr22/`
(group storage — **never** the repo/`$HOME`; consumed later by Step 05).
Wave maps + SLURM logs stay in `test/jaguar_chr22/` (tiny, git-ignored).

## Test01 result (2026-09-01, jobs 278359-61) — PASS

All three waves finished `[$] ... completed successfully!`. chr22 GenomicsDB at
`/mnt/data/amedina/esalazarf/JVCdev/test_gendbi_jag22/genomicsdb/chr22` holds
**93 samples**, ~5.3 GB.

| wave | action | db total | wall | note |
|------|--------|----------|------|------|
| 1 | create | 31 | 54 min | consolidated |
| 2 | update | 62 | 147 min | **consolidated the whole 62-sample array** |
| 3 | update | 93 | 190 min | consolidate + the tiny sample warning |

**Findings folded back into the code:**

- `--consolidate` cost scales with *total* DB size, not the wave (it rewrites the
  whole array). S04 now consolidates on `create` only; `run_jaguar_waves.sh`
  sets `GENDBI_CONSOLIDATE` true / **false** / true across the three waves so
  only wave 1 and the last wave pay it. Re-running test01 should drop wave 2 to
  ~50 min.
- More cores/RAM don't help (wave 3 had 4 CPU, still slower — bigger DB). Waves
  are now all 2 CPU / 16 G.
- `EGAN00004696518` (24 KB) was flagged by the pre-build warning *and* the
  finisher, and imported as a near-empty sample. **Check its content before
  Step 05:** `gatk SelectVariants -R <ref> -V gendb://.../genomicsdb/chr22 -sn EGAN00004696518 -L chr22 -O /tmp/x.vcf.gz` then look at the record count.
- Workspace files are `0700` / group `fsanchezq` → other amedina users can't
  read the DB. Make the output dir setgid the data group (or `umask 0027`)
  before a shared run; `check_setup.sh` now warns about this.
- `.__consolidation_lock` (0 B) left in the array dir is a benign TileDB
  leftover — does not block reads.

Verify the DB reads before Step 05 (S04's finisher prints the exact command).

### First run (2026-08-31) — failed on quota, fixed

Jobs 276894/276895/276896. Wave 1's GenomicsDBImport **succeeded** (72 min on
scratch) but the copy-back went to `<repo>/test/jaguar_chr22/out`, and the repo
is in `$HOME` on FENIX (hard ~5 GB quota) → `Disk quota exceeded`, and the old
trap wiped scratch. Fixed:

- output now defaults to `/mnt/data/amedina/$USER/JVCdev/test_gendbi_jag22`
  (group storage); `run_jaguar_waves.sh` and `check_setup.sh` both refuse a
  `$HOME` output and do a real write-test.
- S04 refuses a `$HOME` output, checks free space before the copy, and if the
  import succeeded but copy-back fails it **keeps the scratch workspace** and
  prints the `cp` + re-run commands to finish without recomputing.

Before re-running: `rm -rf ~/GitHub/lambda/Joint_Variant_Call_hsa/test/jaguar_chr22/out`
(leftover partial copy), then `git pull` and start again — wave 1 must redo the
72 min (last night's scratch was already purged).

The JDK-25 warnings in the log (`System::load` restricted method, `sun.misc.Unsafe`)
are harmless — GATK 4.6 targets JDK 17.
