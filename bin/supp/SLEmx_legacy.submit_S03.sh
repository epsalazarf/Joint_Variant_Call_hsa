#!/usr/bin/env bash
# =============================================================================
# Title       : SLEmx legacy cohort — S03 submission for mbravog [FENIX]
# Description : One-shot script for pre-existing SLE samples that already have
#               a finished Step 02 output (some named with the older
#               "*_<suffix>.sorted.rmdup.mqfilt.bqsr.bam" convention from
#               before this repo's current Step 01/02). For each target
#               sample:
#                 1. best-effort fix group ownership/perms so the amedina
#                    group (mbravog + esalazarf) can read/write it
#                 2. submit Step 03 (per-chromosome scatter) via the shared
#                    launcher, which auto-skips S01/S02 since their outputs
#                    already exist
#               Two samples (Batch09/L23, Batch09/L24) have TWO independent
#               per-lane BQSR BAMs sitting in the same folder instead of one
#               merged BAM — those are archived first so the launcher's normal
#               multi-BAM discovery re-runs Step 02 as a merge, then Step 03
#               runs once on the merged result. See MERGE_SAMPLES below.
#
#               Safe to re-run: the launcher's own resume logic (per sample,
#               per chromosome) skips whatever is already done.
# Usage       : bash bin/supp/SLEmx_legacy.submit_S03.sh
#               (paths are absolute; run from anywhere after `git pull`)
# =============================================================================

set -uo pipefail   # no -e: one sample's failure must not abort the rest

R="/mnt/data/amedina/esalazarf/Joint_Variant_Call_hsa"
LAUNCHER="$R/bin/PIPELINE.single_sample_01-03.sh"

[ -f "$LAUNCHER" ] || {
  echo "<ERROR> Launcher not found: $LAUNCHER"
  echo "<ERROR> git pull the repo at $R first."
  exit 1
}

# --- Single-BQSR-BAM samples: Step 02 already complete, just need Step 03 ---
BAMS=(
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch06/Q067/Q067_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch06/Q071/Q071_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch06/Q073/Q073_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch06/Q075/Q075_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch06/Q077/Q077_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch06/Q080/Q080_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Q104/Q104_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Q105/Q105_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Q106/Q106_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Q107/Q107_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Q108/Q108_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Qro_109/Qro_109_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Qro_112/Qro_112_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch08/Qro_113/Qro_113_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L56/L56_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/Q064/Q064_D_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/Qro_128/Qro_128_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/Qro_115/Qro_115_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/Q00c/Q00c_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/Q00e/Q00e_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/Q00f/Q00f_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/L24/L24_5.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/L53/L53_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/L55/L55_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/L68/L68_6.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12/Q030/Q030_D_1.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch15/Q_001/Q_001_7.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch15/Q_003/Q_003_7.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch15/Q_007/Q_007_7.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch15/Q_008/Q_008_7.sorted.rmdup.mqfilt.bqsr.bam"
)

# --- Two-lane samples: Batch09/L23 and Batch09/L24 each have TWO independent
# per-lane BQSR BAMs instead of one merged BAM. Per maintainer decision
# (2026-09-11): archive the stale per-lane outputs, then let the launcher's
# normal multi-BAM discovery re-run Step 02 as a proper merge (both *.sorted.bam
# lane inputs -> one <sample>.rmdup.mqfilt.bqsr.bam), then Step 03 once.
# Format: "sample_dir|lane_bqsr_bam1|lane_bqsr_bam2"
MERGE_SAMPLES=(
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L23|/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L23/L23_1.sorted.rmdup.mqfilt.bqsr.bam|/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L23/L23_2.sorted.rmdup.mqfilt.bqsr.bam"
  "/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L24|/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L24/L24_1.sorted.rmdup.mqfilt.bqsr.bam|/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch09/L24/L24_2.sorted.rmdup.mqfilt.bqsr.bam"
)

echo "[i] ${#BAMS[@]} single-BAM samples + ${#MERGE_SAMPLES[@]} two-lane merge samples"
echo

# --- group-access fix: make sure amedina (mbravog + esalazarf) can read/write.
# Tries the lab's `amgrp` helper first if present; always also does the
# equivalent by hand (chgrp + g+rwX + setgid on dirs so anything written back
# — chrom_gvcf/, logs — inherits group amedina too). Failure here usually means
# the directory's owner/perms need a sysadmin fix; we report it and move on.
fix_group() {
  local d="$1" ok=true
  command -v amgrp >/dev/null 2>&1 && amgrp "$d" >/dev/null 2>&1
  chgrp -R amedina "$d" || ok=false
  chmod -R g+rwX "$d" || ok=false
  find "$d" -type d -exec chmod g+s {} + 2>/dev/null || ok=false
  $ok
}

FAILED_PERMS=()
FAILED_READ=()
SUBMITTED=()

echo "===== single-BAM samples ====="
for b in "${BAMS[@]}"; do
  d="$(dirname "$b")"
  s="$(basename "$d")"
  echo "-- $s --"

  if ! fix_group "$d"; then
    echo "   [!] group-perm fix failed on $d — needs sysadmin"
    FAILED_PERMS+=("$d")
  fi

  if [ ! -r "$b" ]; then
    echo "   [X] still can't read $b — skipping (needs sysadmin fix first)"
    FAILED_READ+=("$b")
    continue
  fi

  bash "$LAUNCHER" "$d"
  SUBMITTED+=("$s")
  echo
done

echo "===== two-lane merge samples (L23, L24) ====="
for entry in "${MERGE_SAMPLES[@]}"; do
  IFS='|' read -r d lane1 lane2 <<< "$entry"
  s="$(basename "$d")"
  echo "-- $s (merging $(basename "$lane1") + $(basename "$lane2")) --"

  if ! fix_group "$d"; then
    echo "   [!] group-perm fix failed on $d — needs sysadmin"
    FAILED_PERMS+=("$d")
  fi

  if [ ! -r "$lane1" ] || [ ! -r "$lane2" ]; then
    echo "   [X] can't read one or both lane BAMs — skipping (needs sysadmin fix first)"
    FAILED_READ+=("$lane1" "$lane2")
    continue
  fi

  # Archive (never delete) the stale per-lane BQSR outputs + sidecars so the
  # launcher's S02 skip-check doesn't see them and re-runs S02 as a merge of
  # both *.sorted.bam lane inputs instead. Explicit patterns, not a blanket
  # "${pref}.*" glob: pref (e.g. "L23_1.sorted") is itself a prefix of the
  # lane's *.sorted.bam INPUT we must keep — ".rmdup*" etc. only ever matches
  # Step 02's own outputs, never the pre-existing sorted BAM.
  archive="$d/_prelane_bqsr_archive_$(date +%Y%m%d)"
  mkdir -p "$archive"
  for lane in "$lane1" "$lane2"; do
    pref="$(basename "$lane" .rmdup.mqfilt.bqsr.bam)"
    mv -f "$d/${pref}".rmdup* "$archive/" 2>/dev/null
    mv -f "$d/${pref}"-dups.txt "$archive/" 2>/dev/null
    mv -f "$d/${pref}".rmb.mosdepth.* "$archive/" 2>/dev/null
    mv -f "$d/${pref}".alignment_metrics.txt "$archive/" 2>/dev/null
    mv -f "$d/${pref}".insert_size_* "$archive/" 2>/dev/null
    mv -f "$d/${pref}".mqfilt-counts.txt "$archive/" 2>/dev/null
    mv -f "$d/${pref}".bqsr_covariates.* "$archive/" 2>/dev/null
  done
  echo "   [i] archived stale per-lane BQSR output(s) -> $archive"

  bash "$LAUNCHER" "$d"
  SUBMITTED+=("$s")
  echo
done

echo "=============================================================="
echo "[i] Submitted S03 for ${#SUBMITTED[@]}/$((${#BAMS[@]} + ${#MERGE_SAMPLES[@]})) samples."
if [ ${#FAILED_PERMS[@]} -gt 0 ]; then
  echo "[!] Group-permission fix failed on:"
  printf '      %s\n' "${FAILED_PERMS[@]}"
fi
if [ ${#FAILED_READ[@]} -gt 0 ]; then
  echo "[X] Skipped (unreadable even after perms fix) — needs sysadmin:"
  printf '      %s\n' "${FAILED_READ[@]}"
fi
echo
echo "[i] Monitor: squeue -u \$USER"
echo "[i] A sample is done when its chrom_gvcf/ has 25 GVCFs+.tbi and"
echo "    <sample>.raw_variants.canon_chr.g.vcf.gz exists."
#EOF
