# JAGUAR chr22 — Step 04 first test

First real-data validation of `bin/04_gatk_GenomicsDB_import.sh` on FENIX.
**93 samples**, chr22 only, imported in **three waves of 31** to exercise the
incremental (`update`) path twice.

| Wave | Action | Samples | First … last | SLURM (CPU / mem / time) | reader-threads / batch |
|------|--------|---------|--------------|--------------------------|------------------------|
| 1 | `create` | 31 | EGAN00004552304 … EGAN00004552334 | 2 / 16G / 6h | 2 / 50 |
| 2 | `update` | 31 | EGAN00004552335 … EGAN00004696501 | 2 / 16G / 6h | 2 / 50 |
| 3 | `update` | 31 | EGAN00004696502 … EGAN00004710984 | 4 / 16G / 8h | 4 / 100 |

**Resource experiment — mostly already answered.** `seff 276894` (the first
wave-1 attempt at 8 cores / 32 GB) reported **CPU efficiency 9.53%** and
**6.68 GB** RAM used over 72 min. GenomicsDBImport for one interval is
effectively serial (the tail `Consolidating` step) and I/O-bound, so more
cores/RAM don't shorten it. Waves are now small; wave 3 keeps 4 cores purely as
a final confirmation. `run_jaguar_waves.sh report` will show whether even that
made a difference.

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

## What success looks like

- All three waves finish `[$] ... completed successfully!`; workspace sample
  count grows 31 → 62 → 93.
- Wave 3 logs `[!] WARNING: 'EGAN00004696518' GVCF unusually small` and imports
  it anyway. Verify that sample's chr22 content afterward:
  `bcftools view -H gendb://... -s EGAN00004696518 | head` (via Step 05, or
  `gatk SelectVariants`).
- `run_jaguar_waves.sh report` shows `Elapsed` / `TotalCPU` / `MaxRSS` per wave —
  compare wave 2 vs wave 3 to judge whether the extra CPUs/RAM helped.

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
