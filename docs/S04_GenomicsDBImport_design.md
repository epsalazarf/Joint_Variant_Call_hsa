# Step 04 — GenomicsDBImport: design notes

**Script:** `bin/04_gatk_GenomicsDB_import.sh`
**Status:** **validated on FENIX** — JAGUAR chr22, 93 samples, 3-wave incremental (test01, 2026-09-01). Harness: `test/jaguar/` (any chromosome). Whole-genome scatter + shared-perms setup still to do.
**Source:** [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) · [GenomicsDBImport tool doc](https://gatk.broadinstitute.org/hc/en-us/articles/360057439331-GenomicsDBImport)

---

## Purpose

Consolidate the per-sample, per-chromosome GVCFs from Step 03a into **GenomicsDB
workspaces — one per chromosome** — which Step 05 (`GenotypeGVCFs`) reads via
`gendb://`. This is the cohort-level fan-in point: every sample that will be
joint-genotyped together must be imported into the same set of workspaces.

---

## Decisions

### 1. Input contract: explicit sample map (not directory auto-discovery)

**Chosen:** a hand-checkable 2-column TSV — `sample_label <TAB> chrom_gvcf_dir`.

**Rejected:** pointing the script at a parent directory and letting it discover
`*/chrom_gvcf/`. Auto-discovery is one wrong path away from importing the wrong
cohort or (via a future slurmer) flooding the queue. An explicit list is the
same lesson the `*_batch-slurmer.sh` scripts still need to learn.

The script header carries a copy-paste block to generate the map:

```bash
cd /path/to/cohort
for d in */chrom_gvcf; do
  s=$(basename "$(dirname "$d")")
  printf '%s\t%s\n' "$s" "$(readlink -f "$d")"
done > cohort.sample_map.tsv
```

The real database sample name is taken from each GVCF header
(`bcftools query -l`), so a typo in column 1 is a warning, not a silent
mislabel. Duplicate sample names, missing per-chromosome GVCFs, and missing
indexes are hard errors that name the offending sample.

### 2. Incremental import for sample waves

**Chosen:** `create` by default; `update` as the 4th positional argument.

```
First wave :  04_gatk_GenomicsDB_import.sh <map_wave1> <out> chr22 create
Later waves :  04_gatk_GenomicsDB_import.sh <map_wave2> <out> chr22 update
```

`update` uses `--genomicsdb-update-workspace-path`, which appends the new
samples in place. The lupus cohort is expected to arrive in batches, so this is
a first-class feature, not an afterthought. Constraints (surfaced as errors):
intervals cannot change on update; the workspace must already exist; do not
re-list already-imported samples.

### 3. Per-chromosome workspaces + a chromosome floodgate

**Chosen:** one workspace per chromosome (`genomicsdb/<CHR>`), selector arg:

| Selector | Expands to |
|----------|------------|
| `chr1` … `chr22`, `chrX`, `chrY`, `chrM` | that one chromosome |
| `autosomes` | `chr1`…`chr22` |
| `all` | `chr1`…`chr22`, `chrX`, `chrY`, `chrM` |

Multi-chromosome selectors run **serially in one process** for now — fine for
testing and small cohorts. Real scatter (one SLURM job per chromosome) will
come from a `GENDBI.seq_batch-slurmer.sh` wrapper once S04 is proven. Testing
starts with `chr22` (smallest autosome).

### 4. Scratch: build on `/scratch`, copy the workspace back

Unlike the sequential GATK steps (02/03), GenomicsDBImport's TileDB workspace
does **small random I/O** and should not be built directly on NFS. When
`scratch_root` is set (FENIX), the script:

1. builds the workspace under `${scratch_root}/${USER}/gendbi_<jobid>_<ts>/gdb/<CHR>`,
2. points `TMPDIR` at the same scratch area,
3. copies the finished workspace back to `<out>/genomicsdb/<CHR>` via a
   `.<CHR>.new` staging dir + atomic `mv` (an interrupted copy never corrupts
   an existing workspace),
4. wipes scratch on exit (trap).

For `update`, the existing workspace is first copied **to** scratch, updated
there, then copied back — so a failed update leaves the NFS original intact.

Off FENIX (`scratch_root` empty) the same staging happens under a hidden
`.gendbi_work_*` dir inside the output path.

### 5. Mitochondrial contig: `chrM` standard, `chrMT` fallback

`chrM` (UCSC / `hg38.fa`) is the house standard. If a cohort's GVCFs use
`chrMT` instead, the script detects it from the first sample's directory,
prints a warning, and builds `genomicsdb/chrMT` with `-L chrMT` so the interval
matches the data. The finisher applies the same fallback when verifying.

### 6. Engine: GenomicsDB, not CombineGVCFs

`CombineGVCFs` is simpler for tiny cohorts but does not scale and has no
incremental-add. Cohorts here only grow, and decision 2 (waves) depends on the
GenomicsDB update path, so GenomicsDB it is — even for the first small test.

---

### 7. Truncation guard + env-overridable knobs (added for the JAGUAR test)

- Every input GVCF is checked for the 28-byte **BGZF EOF marker** (same test
  htslib uses). Missing marker ⇒ `[!] WARNING: ... likely TRUNCATED` and the
  label is listed again in the finisher. Default is *warn and keep going*
  (`GENDBI_STRICT_GVCF=false`); set it `true` to abort before calling GATK.
  Rationale: a truncated GVCF makes GenomicsDBImport fail the **whole run**
  (it is all-or-nothing per invocation), so flagging it early saves a wasted job.
- `--reader-threads` now defaults to `$SLURM_CPUS_PER_TASK` (else 4). Env
  overrides for a launcher to vary per run: `GENDBI_READER_THREADS`,
  `GENDBI_BATCH_SIZE`, `GENDBI_JAVA_MEM`. Echoed at startup as the `Knobs` line.

## FENIX resourcing

**Measured** — `seff 276894` (first wave-1 attempt: 31 samples, chr22, `create`,
8 cores / 32 GB):

| metric | value | reading |
|--------|-------|---------|
| CPU efficiency | 9.53% | ~0.76 core busy over 72 min wall — effectively serial + I/O-bound |
| Memory used | 6.68 GB / 32 GB | heap need is small; most is native TileDB |
| Wall clock | 72 min | dominated by the tail `Consolidating GenomicsDB array` step |

GenomicsDBImport for a single interval does not parallelise usefully (the
consolidate step is serial; the rest is NFS read + BeeGFS write). **Extra cores
and RAM do not shorten it.** Suggested per-chromosome job:

```
--cpus-per-task=2   --mem=16G   --time=<scale with contig size; ~1.5 h for chr22>
```

Java `-Xmx6g` (`GENDBI_JAVA_MEM`), `--reader-threads 2`, `--batch-size 50`
(imports all samples at once below 50; caps open file handles / memory above
that), `--genomicsdb-shared-posixfs-optimizations true`.

### `--consolidate` — do not run it on `update`

| test | wave | action | db total | wall | consolidate portion |
|------|------|--------|----------|------|---------------------|
| 01 chr22 | 1 | create | 31 | 54 min | ~3 s (one fragment) |
| 01 chr22 | 2 | update | 62 | 147 min | ~96 min |
| 01 chr22 | 3 | update | 93 | 190 min | ~140 min |
| 02 chr1 | 1 | create | 47 | 7 h 40 m | ~4 s (one fragment) |
| 02 chr1 | 2 | update | 93 | **21 h** | **~13.5 h** |

On `create` (a single fresh fragment) `--consolidate` is a no-op. On `update`
it merges fragments across the **whole array** — and for chr1→93 samples that
was **13.5 h**, longer than a full `create` rebuild of all 93 samples (~15 h).

**Decision:** `update` waves never consolidate (S04 default; the launcher only
sets it on wave 1). GenotypeGVCFs reads a multi-fragment DB fine. If a
single-fragment DB is genuinely wanted (e.g. read performance turns out to
matter), **rebuild the chromosome with `create` on all current samples** —
don't `GENDBI_CONSOLIDATE=true` an update.

### Storage

chr1 at 93 samples = **~30 GB** (~332 MB/sample; chr22 was ~57 MB/sample).
chr1 is ~8% of the genome, so a full 93-sample **whole-genome** GenomicsDB is
**~375–400 GB**. Size it into the plan for the real cohort.

### Shared readability

GenomicsDB writes workspace files mode `0700`, owned by your primary group
(`fsanchezq` for the maintainer). For other amedina members to run Step 05 on
the DB, the **output directory must be setgid** to the data group
(`chmod g+s`, `chgrp`) or the job must run with `umask 0027`. `check_setup.sh`
warns when the output dir is not setgid.

---

## Known limitations / follow-ups

- **No SLURM launcher yet.** `GENDBI.seq_batch-slurmer.sh` (one job per
  chromosome, with the same `autosomes` / `all` gate) is deferred until S04 is
  validated on FENIX.
- **Not in `run_pipeline.sh` / `PIPELINE.single_sample.sh`.** S04 is
  cohort-level, not single-sample; it needs its own launcher.
- **One interval == one whole chromosome.** Fine for now. Very large cohorts
  may want finer scatter (e.g. `chr1:1-125000000`) — revisit if chr1/chr2
  imports get slow.
- **Reference `-R` required.** Passed for contig validation, matching S03.
  Local `hg38.fa.gz` will not work with `-R` (GATK rejects bgzipped FASTA);
  use a plain `.fa` locally. See `docs/pipeline_review_notes.txt`.
- **bash ≥ 4** for the real timing lines (`EPOCHSECONDS`); consistent with
  S02/S03. The script otherwise avoids bash-4-only constructs.

---

## First FENIX test (pending): JAGUAR cohort, chr22

See `test/jaguar/`. 93 samples, 3 waves of 31 (create + update + update),
escalating resources per wave. Includes one deliberately-truncated GVCF
(`EGAN00004696518`) to confirm the truncation flag. Helpers:
`make_wave_maps.sh`, `check_setup.sh` (pre-flight), `run_jaguar_waves.sh`
(chained sbatch + manifest + `report` mode for `sacct` timings).

## Local validation (2026-08-31, GATK 4.6.2.0, synthetic 3-sample data)

`create` chr22 · `update` chr22 (+1, verified 3 samples in DB) · `create`
`chrM` · `chrMT` fallback · rerun skip · `update` with no workspace (error) ·
invalid selector (error) · sample missing a chromosome's GVCF (error, names
sample) · duplicate sample (error) · missing sample map (error). All pass;
all error paths exit non-zero. See `docs/PIPELINE_STATUS.md` §5.
