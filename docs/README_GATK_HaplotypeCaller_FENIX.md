# GATK HaplotypeCaller [FENIX]

**Author:** Pavel Salazar-Fernandez et al.  
**Source:** [GATK4 Best Practices – Germline Short Variant Discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932)  
**Version:** FENIX pipeline module  

---

## Overview

This script automates **GATK4 HaplotypeCaller** to perform germline variant calling (SNPs and Indels) from an analysis-ready BAM file.  
It is designed for reproducible use in both **local** and **HPC (FENIX)** environments.

---

## Usage

```bash
bash 03_gatk_variant_calling.sh [Mapped_BAM] [Output_Path]
```

**Example:**
```bash
bash 03_gatk_variant_calling.sh sample123.sorted.bam /path/to/output/
```

**Output:**  
`sample123.raw_variants.canon_sort.g.vcf.gz`

---

## Workflow Summary

1. **Setup** – Validates inputs, parses config (`config/config.yaml`), detects local or SSH environment.  
2. **Variant Calling** – Runs `gatk HaplotypeCaller` to produce a raw GVCF.  
3. **Indexing** – Indexes output using `bcftools index`.  
4. **Canonical Filtering** – Keeps only canonical chromosomes (1–22, X, Y, M).  
5. **Finisher** – Validates final output and reports runtime.

---

## Limitations

- **Single BAM per run:**  
  This script currently accepts **only one BAM file per sample**.  
  If a sample has multiple BAMs (e.g., from multiple lanes or libraries), **merge them before running** this script.

- **Future feature:**  
  Upcoming versions will support **automatic merging of multiple BAMs** per sample before variant calling.

---

## Dependencies

- `GATK >= 4.0`  
- `bcftools`  
- `awk`, `bash >= 4.0`  
- Optional: `module` environment (for remote/HPC use)

---

## Configuration

Reference files and paths are specified in:
```
config/config.yaml
```
Separate configurations for **local** and **remote** environments are supported.

---

## Output

- Final GVCF: `sample.raw_variants.canon_sort.g.vcf.gz`  
- Indexed with `.tbi`  
- Ready for joint genotyping or downstream filtering.

---

## Quick Troubleshooting

| Issue | Cause | Solution |
|-------|--------|-----------|
| `<ERROR> File not found` | Incorrect BAM or reference path | Check file paths and ensure symbolic links are valid |
| `<ERROR> Missing reference paths` | Missing or misconfigured YAML entries | Verify `ref_gnm` and `ref_vars` in `config/config.yaml` |
| GATK command not found | Module not loaded or PATH not set | Load required modules (`module load gatk`, `bcftools`) |
| Script exits early | Wrong arguments | Run as `bash 03_gatk_variant_calling.sh [BAM] [OUTPUT_PATH]` |
| Output missing or empty | Step failed silently | Re-run with verbose mode and check disk space |

---
