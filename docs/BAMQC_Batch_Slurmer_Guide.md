## BAMQC Batch Slurmer

**Script:** `bin/supp/BAMQC.seq_batch-slurmer.sh`
**Purpose:** Submit Step 02 (`02_gatk_bam_qc_workflow.sh`) for every sample in a directory, one independent SLURM job per sample (parallel).

> For new samples, prefer the one-command launchers — `bin/PIPELINE.single_sample_01-03.sh <sample_dir>` or `bin/supp/BATCH.submit_S01-S03.sh <batch_dir>` — which chain Steps 01→02→03 in the sample directory with no manual hand-off. Use the per-step slurmers below only when you need to run a single step across many samples.

### Usage

```bash
bash bin/supp/BAMQC.seq_batch-slurmer.sh bin/02_gatk_bam_qc_workflow.sh /path/to/bams/
```

`/path/to/bams/` defaults to `$PWD`.

### Getting BAMs from Step 01 into place

The slurmer only searches the directory it is given (recursively). Step 01's batch slurmer decides where BAMs land:

```bash
# Recommended: send Step 01 output to the BAM directory
bash bin/supp/BWAMAP.seq_batch-slurmer.sh bin/01_bwa_map_fastq_reads.sh /path/to/fastq_samples/ /path/to/bams/
bash bin/supp/BAMQC.seq_batch-slurmer.sh  bin/02_gatk_bam_qc_workflow.sh  /path/to/bams/
```

If Step 01 was run **without** the third argument, the `*.sort.bam` files are inside each FASTQ sample directory. Either point this slurmer at the FASTQ parent directory, or symlink the BAMs into your BAM directory first. Both work; keep the original BAMs in place if you symlink.

### What it does

1. Finds every `*.sort.bam` / `*.sorted.bam` under the directory.
2. Groups them by sample ID (the text before the first `_` in the filename), so per-lane BAMs such as `SAMPLE_L1.sort.bam,SAMPLE_L2.sort.bam` are merged in one Step 02 run.
3. Writes the list to `<dir>/<dirname>.bamqc_input_files.txt` (tab-separated: `sample_name  sample_file(s)  readgroup_string`).
   - If this file already exists it is **reused** instead of re-scanning. Edit it to add an `@RG` string for legacy BAMs without read groups, or delete it to force a fresh scan.
4. For each sample: creates `<dir>/<SAMPLE>/`, symlinks the BAM(s) there, and submits Step 02 with that folder as the output directory.

### SLURM resources

Each job requests 1 node, 1 task, 8 CPUs, 32 GB. Logs: `<dir>/<SAMPLE>/BAMQC-<SAMPLE>-<epoch>.<jobid>.log`.

```bash
squeue -u $USER | grep BAMQC
```

### Output

Final analysis-ready BAM: `<dir>/<SAMPLE>/<SAMPLE>.rmdup.mqfilt.bqsr.bam` — pass `<dir>` to the Step 03 slurmer (`HAPCALL.seq_batch-slurmer.sh`).
