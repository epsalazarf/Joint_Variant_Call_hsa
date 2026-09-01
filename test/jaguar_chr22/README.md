# JAGUAR chr22 — Step 04 first test

First real-data validation of `bin/04_gatk_GenomicsDB_import.sh` on FENIX.
**93 samples**, chr22 only, imported in **three waves of 31** to exercise the
incremental (`update`) path twice.

| Wave | Action | Samples | First … last | SLURM (CPU / mem / time) | reader-threads / batch |
|------|--------|---------|--------------|--------------------------|------------------------|
| 1 | `create` | 31 | EGAN00004552304 … EGAN00004552334 | 8 / 32G / 12h | 8 / 50 |
| 2 | `update` | 31 | EGAN00004552335 … EGAN00004696501 | 8 / 64G / 24h | 8 / 50 |
| 3 | `update` | 31 | EGAN00004696502 … EGAN00004710984 | 16 / 128G / 48h | 16 / 100 |

Wave 3 deliberately doubles wave 2's CPU + RAM so we can see whether more
resources actually shorten GenomicsDBImport (mostly-native, mildly parallel —
this is the experiment, not a real requirement). Big `--mem` is head-room only.

### Known bad sample (left in on purpose)

`EGAN00004696518` (wave 3) is **24 KB** vs a ~48 MB median — a truncated GVCF.
It is left in the map to confirm Step 04 flags it. Expected behaviour:

- `make_wave_maps.sh` and `check_setup.sh` print a `TINY` / `no BGZF EOF marker` warning.
- Step 04 prints `[!] WARNING: ... likely TRUNCATED` and then GenomicsDBImport
  fails with `Premature end of file`, naming the file. Wave 3 exits non-zero.
- **Waves 1–2 (62 samples) stay intact** — the failed update only touches the
  scratch copy of the workspace.
- Recovery: re-run Step 03a for that sample, then re-run wave 3 as `update`
  (or set `GENDBI_STRICT_GVCF=true` to abort wave 3 before calling GATK).

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

Output workspace: `test/jaguar_chr22/out/genomicsdb/chr22/` (consumed later by Step 05).

## What success looks like

- Waves 1 and 2 finish `[$] ... completed successfully!`, workspace sample count
  grows 31 → 62.
- Wave 3 fails on `EGAN00004696518` (expected); the other 30 wave-3 samples are
  **not** imported (GenomicsDBImport is all-or-nothing per run). After fixing the
  bad GVCF, re-running wave 3 imports all 31 and the DB reaches 93.
- `run_jaguar_waves.sh report` shows `Elapsed` / `TotalCPU` / `MaxRSS` per wave —
  compare wave 2 vs wave 3 to judge whether the extra CPUs/RAM helped.
