## BAMQC Batch Slurmer

**Title:** BAMQC Batch Slurmer (sequential, tolerant to failed jobs)  
**Purpose:** Automatically submit a series of BAM quality control (QC) jobs to a SLURM cluster, running them sequentially (each job starts after the previous one finishes).  
Useful when processing multiple BAM files through a single QC workflow script.

### Usage

```bash
./BAMQC.batch-slurmer batch-list.txt 02_gatk_bam_qc_workflow.sh /path/to/bam/files
```

Each line in `batch-list.txt` must contain **three tab-separated fields**:

```
sample_name    sample_file    readgroup_string
```

**Example:**
```
sample1    sample1.bam    ID:sample1\tSM:sample1\tPL:ILLUMINA
sample2    sample2.bam    ID:sample2\tSM:sample2\tPL:ILLUMINA
```

### What it does

1. Reads each line of the batch list and:
   - Creates a directory named after the sample.
   - Symlinks the BAM file from the provided directory.
   - Submits a SLURM job running the specified workflow script (`02_gatk_bam_qc_workflow.sh`).
2. Chains jobs so that each one runs **after the previous finishes** (even if the previous job fails).
3. Generates a log file per job in the format:  
   ```
   BAMQC-samplename.<jobid>.log
   ```

### SLURM Resources

Each job requests:
- 1 node, 1 task  
- 4 CPUs (`--cpus-per-task=4`)  
- 32 GB memory  
- Output log written to the working directory

### Notes

- Ensure the BAM files in `/path/to/bam/files` match the names listed in `batch-list.txt`.
- Make sure `02_gatk_bam_qc_workflow.sh` is executable and tested for a single sample before running the batch.
- Jobs will be chained sequentially; check job progress with:
  ```bash
  squeue -u $USER | grep BAMQC
  ```
- The final message shows the **last job ID** in the chain.
