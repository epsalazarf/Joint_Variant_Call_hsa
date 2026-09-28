#!/usr/bin/env bash
# =============================================================================
# Title       : Per-sample S01→S03 Summary Table [FENIX]
# Description : Scans the outputs and logs left by Steps 01–03 in each sample
#               directory and writes one TSV row per sample: pipeline stage
#               reached, file sizes, processing time per step, and key QC stats
#               (reads, % mapped, % duplication, insert size, mean depth).
#               Optionally re-derives depth from the GVCFs (--gvcf-depth).
#               Read-only: never modifies anything in the sample dirs.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-09-28
# Version     : 1.1
# Usage       : SUMMARY.sample_stats_S01-S03.sh <dir> [dir ...] [-o out.tsv]
#                 [--gvcf-depth[=fast|full]] [-j threads]
#               dir — a sample dir (as laid out by PIPELINE.single_sample_01-03.sh)
#                     or a batch dir whose immediate subdirs are sample dirs.
#               -o  — output TSV (default: ./<first_dir>.S01-S03_summary.<date>.tsv)
#               --gvcf-depth       — same as --gvcf-depth=fast
#               --gvcf-depth=fast  — GVCF depth on chr20 only (~2% of genome, ~1 min/sample)
#               --gvcf-depth=full  — GVCF depth over chr1–22 (one pass over the GVCFs)
#               -j  — parallel chromosome jobs per sample in full mode (default 4)
#
# GVCF depth: length-weighted mean of FORMAT/DP (each record weighted by the
#   bases it spans: END-POS+1 for ref blocks, 1 for variant sites; ref-block DP
#   is GATK's block median). Read from chrom_gvcf/ (fallback: canon_chr GVCF).
#   It is HaplotypeCaller's FILTERED depth (MAPQ>=20, downsampled), so expect it
#   a bit below mosdepth; gvcf_vs_mosdepth = ratio against mosdepth's mean over
#   the SAME region (chr20 or chr1–22). Off by default: columns are NA.
#
# Sources per column (all optional — a missing source leaves "NA"):
#   S01  *.sort.bam / *.sorted.bam, *.sort.stats.txt (samtools stats SN lines)
#   S02  *.rmdup.mqfilt.bqsr.bam, *-dups.txt (Picard), *.insert_size_metrics.txt,
#        *.rmb.mosdepth.summary.txt (autosomal mean = Σbases/Σlength, chr1–22,
#        same as archive/filter_mosdepth_nuclear.sh)
#   S03  chrom_gvcf/*.raw_vars.<chrom>.g.vcf.gz, *.raw_variants.canon_chr.g.vcf.gz
#   Time "[&]  Total time:" line of the most recent SUCCESSFUL log per step.
#        Logs are searched in <sample>/log/ and <sample>/ and classified by
#        their banner, so filenames don't matter. S03 scatter: s03_time_max is
#        the slowest chromosome (≈ wall time), s03_time_sum the total compute.
#        *_runs = number of logs found for that step (>1 means resumes/reruns;
#        the time reported covers only the last successful run).
# =============================================================================

set -uo pipefail   # no -e: one unreadable sample must not abort the table

# <ARGUMENTS> -----------------------------------------------------------------

OUT_TSV=""
GVCF_DEPTH=off   # off | fast | full
NJOBS=4
INPUTS=()
while [ $# -gt 0 ]; do
  case "$1" in
    -o) OUT_TSV="${2:?-o needs a file name}"; shift 2 ;;
    -j) NJOBS="${2:?-j needs a number}"; shift 2 ;;
    --gvcf-depth|--gvcf-depth=fast) GVCF_DEPTH=fast; shift ;;
    --gvcf-depth=full)              GVCF_DEPTH=full; shift ;;
    --gvcf-depth=*) echo "<ERROR> --gvcf-depth takes fast or full"; exit 1 ;;
    -h|--help) sed -n '2,45p' "$0"; exit 0 ;;
    *)  INPUTS+=("$1"); shift ;;
  esac
done
[ ${#INPUTS[@]} -gt 0 ] || { echo "Usage: $(basename "$0") <dir> [dir ...] [-o out.tsv]"; exit 1; }

[ -n "$OUT_TSV" ] || OUT_TSV="$PWD/$(basename "$(cd "${INPUTS[0]}" && pwd)").S01-S03_summary.$(date +%F).tsv"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

CANON_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
CANON_N=25

# Banner strings printed by each step script (first lines of its log)
BANNER_S01='BWA FASTQ Reads Mapper'
BANNER_S02='GATK4 BAM QC'
BANNER_S03='GATK HaplotypeCaller'

AUTOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
FAST_CHROM="chr20"   # GATK's usual benchmark chromosome: mid-sized, typical content

# bcftools only needed for --gvcf-depth (module on FENIX)
if [ "$GVCF_DEPTH" != off ] && ! command -v bcftools >/dev/null; then
  module load bcftools >/dev/null 2>&1 || true
  command -v bcftools >/dev/null || { echo "<ERROR> --gvcf-depth needs bcftools in PATH"; exit 1; }
fi

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## Size of one file in bytes (GNU stat on FENIX, BSD stat locally)
fsize() { stat -c%s "$1" 2>/dev/null || stat -f%z "$1" 2>/dev/null || echo 0; }

## Sum sizes of the files given as args → GiB (2 dp), "NA" if none exist
sum_gib() {
  local f total=0 n=0
  for f in "$@"; do
    [ -f "$f" ] || continue
    total=$(( total + $(fsize "$f") )); n=$(( n + 1 ))
  done
  [ "$n" -gt 0 ] && awk -v b="$total" 'BEGIN{printf "%.2f", b/1024^3}' || echo NA
}

## "H:MM:SS" / "MM:SS" / "SS" → seconds
to_sec() { awk -F: '{s=0; for(i=1;i<=NF;i++) s=s*60+$i; print s}' <<< "$1"; }

## seconds → HH:MM:SS ("NA" passes through)
to_hms() {
  [ "$1" = NA ] && { echo NA; return; }
  printf '%02d:%02d:%02d' $(( $1 / 3600 )) $(( ($1 % 3600) / 60 )) $(( $1 % 60 ))
}

## Logs for one step, oldest → newest (by mtime), matched on banner text
step_logs() {   # $1 = sample dir   $2 = banner
  local d="$1"
  find "$d/log" "$d" -maxdepth 1 -type f -name '*.log' 2>/dev/null \
    | while IFS= read -r f; do
        head -30 "$f" 2>/dev/null | grep -qF "$2" && printf '%s\t%s\n' "$(fsize_mtime "$f")" "$f"
      done | sort -n | cut -f2- | awk '!seen[$0]++'
}
fsize_mtime() { stat -c%Y "$1" 2>/dev/null || stat -f%m "$1" 2>/dev/null || echo 0; }

## Total-time seconds from a log if it completed successfully, else empty
log_ok_seconds() {
  grep -q 'completed successfully' "$1" 2>/dev/null || return 0
  local t
  t=$(grep -E '^\[&\] +Total time:' "$1" | tail -1 | awk '{print $NF}')
  [ -n "$t" ] && to_sec "$t"
}

## S01/S02: "<n_runs>\t<seconds of latest successful run | NA>"
step_time() {   # $1 = sample dir   $2 = banner
  local logs n=0 sec="NA" f s
  logs=$(step_logs "$1" "$2")
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    n=$(( n + 1 ))
    s=$(log_ok_seconds "$f"); [ -n "$s" ] && sec="$s"   # newest success wins
  done <<< "$logs"
  printf '%s\t%s\n' "$n" "$sec"
}

## S03: "<n_runs>\t<max sec>\t<sum sec>" — scatter logs grouped per chromosome
## (latest success per chrom); a whole-genome log counts as a single "chrom".
s03_time() {
  local logs n f s key
  logs=$(step_logs "$1" "$BANNER_S03")
  n=$(grep -c . <<< "$logs")
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    s=$(log_ok_seconds "$f"); [ -n "$s" ] || continue
    key=$(sed -nE 's/.*Mode: scatter \((chr[^)]+)\).*/\1/p' "$f" | head -1)
    printf '%s\t%s\n' "${key:-WG}" "$s"
  done <<< "$logs" \
    | awk -F'\t' -v n="$n" '
        { last[$1] = $2 }                         # newest success per chrom
        END {
          m = 0; t = 0; k = 0
          for (c in last) { k++; t += last[c]; if (last[c] > m) m = last[c] }
          if (k == 0) print n "\tNA\tNA"; else print n "\t" m "\t" t
        }'
}

## samtools stats (summed over per-RG files) → "raw_reads\tpct_mapped\tpct_proper"
samtools_stats() {
  [ $# -gt 0 ] || { printf 'NA\tNA\tNA\n'; return; }
  awk -F'\t' '
    $1 == "SN" && $2 == "raw total sequences:"   { raw  += $3 }
    $1 == "SN" && $2 == "reads mapped:"          { map  += $3 }
    $1 == "SN" && $2 == "reads properly paired:" { prop += $3 }
    END {
      if (raw > 0) printf "%d\t%.2f\t%.2f\n", raw, 100*map/raw, 100*prop/raw
      else         print "NA\tNA\tNA"
    }' "$@"
}

## Picard duplication metrics → % duplication (pooled across libraries)
dup_pct() {
  [ -f "${1:-}" ] || { echo NA; return; }
  awk -F'\t' '
    /^## METRICS CLASS/ { hdr = 1; next }
    hdr == 1 { for (i = 1; i <= NF; i++) col[$i] = i; hdr = 2; next }
    hdr == 2 && NF < 2 { exit }
    hdr == 2 {
      ue += $col["UNPAIRED_READS_EXAMINED"];  pe += $col["READ_PAIRS_EXAMINED"]
      ud += $col["UNPAIRED_READ_DUPLICATES"]; pd += $col["READ_PAIR_DUPLICATES"]
    }
    END {
      tot = ue + 2*pe
      if (tot > 0) printf "%.2f\n", 100*(ud + 2*pd)/tot; else print "NA"
    }' "$1"
}

## Picard insert size metrics → MEDIAN_INSERT_SIZE of the first (main) orientation
median_insert() {
  [ -f "${1:-}" ] || { echo NA; return; }
  awk -F'\t' '
    /^## METRICS CLASS/ { hdr = 1; next }
    hdr == 1 { for (i = 1; i <= NF; i++) col[$i] = i; hdr = 2; next }
    hdr == 2 { print $col["MEDIAN_INSERT_SIZE"]; found = 1; exit }
    END { if (!found) print "NA" }' "$1"
}

## mosdepth summary → "auto_mean\tchrX_mean\tchrY_mean\tchrM_mean"
mosdepth_means() {
  [ -f "${1:-}" ] || { printf 'NA\tNA\tNA\tNA\n'; return; }
  awk -F'\t' '
    function v(c) { return (c in m) ? sprintf("%.2f", m[c]) : "NA" }
    $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ { len += $2; bases += $3 }
    $1 ~ /^chr[XYM]$/                 { m[$1] = $4 }
    END {
      printf "%s\t%s\t%s\t%s\n", (len > 0 ? sprintf("%.2f", bases/len) : "NA"),
             v("chrX"), v("chrY"), v("chrM")
    }' "$1"
}

## One chromosome of a GVCF → "<Σ DP×len>\t<Σ len>" (records without DP skipped)
gvcf_chrom_sums() {   # $1 = GVCF   $2 = chrom (region filter; needs .tbi)
  bcftools query -r "$2" -f '%POS\t%INFO/END\t[%DP]\n' "$1" 2>/dev/null \
    | awk -F'\t' '$3 != "." { L = ($2 == "." ? 1 : $2 - $1 + 1); b += L * $3; n += L }
                  END { printf "%.0f\t%.0f\n", b, n }'
}
export -f gvcf_chrom_sums

## GVCF depth for a sample → "<region>\t<mean depth>" (NA when no GVCF)
gvcf_depth() {   # $1 = sample dir
  local d="$1" chroms c f canon
  [ "$GVCF_DEPTH" = fast ] && chroms="$FAST_CHROM" || chroms="$AUTOSOMES"
  local region; [ "$GVCF_DEPTH" = fast ] && region="$FAST_CHROM" || region="auto"
  canon=$(first_file "$d" '*.raw_variants.canon_chr.g.vcf.gz')

  # one "<gvcf> <chrom>" job per chromosome: per-chrom shard, else canon GVCF
  for c in $chroms; do
    f=$(first_file "$d/chrom_gvcf" "*.raw_vars.${c}.g.vcf.gz")
    if   [ -n "$f" ] && [ -s "${f}.tbi" ];         then echo "$f $c"
    elif [ -n "$canon" ] && [ -s "${canon}.tbi" ]; then echo "$canon $c"
    else echo "MISSING $c"
    fi
  done \
    | xargs -P "$NJOBS" -n 2 bash -c '[ "$0" = MISSING ] && echo MISSING || gvcf_chrom_sums "$0" "$1"' \
    | awk -F'\t' -v r="$region" '
        $1 == "MISSING" || $2 == 0 { miss = 1; next }   # no shard, or chrom absent/empty
        { b += $1; n += $2 }
        END { if (miss || n == 0) print r "\tNA"; else printf "%s\t%.2f\n", r, b/n }'
}

## mosdepth mean over one named chrom, or the autosomes ("auto")
mosdepth_region_mean() {   # $1 = summary file   $2 = chrom | auto
  [ -f "${1:-}" ] || { echo NA; return; }
  awk -F'\t' -v r="$2" '
    (r == "auto" && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/) || $1 == r { len += $2; b += $3 }
    END { if (len > 0) printf "%.2f\n", b/len; else print "NA" }' "$1"
}

## First match of a glob pattern in dir (maxdepth 1), empty if none
first_file() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | sort | head -1; }

## Collect files matching one or more patterns into the global array FILES
collect() {   # $1 = dir, $2.. = patterns
  local d="$1" p f; shift
  FILES=()
  for p in "$@"; do
    while IFS= read -r f; do [ -n "$f" ] && FILES+=("$f"); done \
      < <(find "$d" -maxdepth 1 -type f -name "$p" 2>/dev/null | sort)
  done
}

## Does a dir look like a pipeline sample dir?
is_sample_dir() {
  [ -d "$1/log" ] || [ -d "$1/chrom_gvcf" ] || \
    find "$1" -maxdepth 1 -type f \( -name '*.f*q.gz' -o -name '*.bam' \) 2>/dev/null | grep -q .
}

## One TSV row for one sample dir
summarize_sample() {
  local d="$1" sample stage="S00" f
  sample=$(basename "$d")

  # --- inputs (FASTQs, excluding S01's own repaired/singleton intermediates)
  collect "$d" '*.fastq.gz' '*.fq.gz'
  local fq=()
  for f in ${FILES[@]+"${FILES[@]}"}; do
    case "$f" in *.repaired.fastq.gz|*_singletons.fastq.gz) ;; *) fq+=("$f") ;; esac
  done
  local fastq_gib; fastq_gib=$(sum_gib ${fq[@]+"${fq[@]}"})

  # --- S01
  collect "$d" '*.sort.bam' '*.sorted.bam'
  local n_sort=${#FILES[@]} sort_gib
  sort_gib=$(sum_gib ${FILES[@]+"${FILES[@]}"})
  [ "$n_sort" -gt 0 ] && stage="S01"
  collect "$d" '*.sort.stats.txt' '*.sorted.stats.txt'
  local sstats; sstats=$(samtools_stats ${FILES[@]+"${FILES[@]}"})
  local s01; s01=$(step_time "$d" "$BANNER_S01")

  # --- S02
  local bqsr; bqsr=$(first_file "$d" '*.rmdup.mqfilt.bqsr.bam')
  local bqsr_gib="NA"
  [ -n "$bqsr" ] && { bqsr_gib=$(sum_gib "$bqsr"); stage="S02"; }
  local dups;   dups=$(dup_pct "$(first_file "$d" '*-dups.txt')")
  local isize;  isize=$(median_insert "$(first_file "$d" '*.insert_size_metrics.txt')")
  local depth;  depth=$(mosdepth_means "$(first_file "$d" '*.mosdepth.summary.txt')")
  local s02; s02=$(step_time "$d" "$BANNER_S02")

  # --- S03
  local n_chr=0 c
  collect "$d/chrom_gvcf" '*.raw_vars.chr*.g.vcf.gz'
  local chr_gib; chr_gib=$(sum_gib ${FILES[@]+"${FILES[@]}"})
  for c in $CANON_CHROMS; do
    f=$(first_file "$d/chrom_gvcf" "*.raw_vars.${c}.g.vcf.gz")
    [ -n "$f" ] && [ -s "$f" ] && [ -s "${f}.tbi" ] && n_chr=$(( n_chr + 1 ))
  done
  local canon; canon=$(first_file "$d" '*.raw_variants.canon_chr.g.vcf.gz')
  local canon_gib="NA"; [ -n "$canon" ] && canon_gib=$(sum_gib "$canon")
  [ "$n_chr" -eq "$CANON_N" ] && stage="S03"
  local s03; s03=$(s03_time "$d")

  # --- optional GVCF depth (vs mosdepth over the same region)
  local g_region="NA" g_depth="NA" g_ratio="NA"
  if [ "$GVCF_DEPTH" != off ]; then
    local g; g=$(gvcf_depth "$d")   # not inline: IFS=$'\t' would leak into it
    IFS=$'\t' read -r g_region g_depth <<< "$g"
    local m; m=$(mosdepth_region_mean "$(first_file "$d" '*.mosdepth.summary.txt')" "$g_region")
    [ "$g_depth" != NA ] && [ "$m" != NA ] && \
      g_ratio=$(awk -v g="$g_depth" -v m="$m" 'BEGIN { if (m > 0) printf "%.3f", g/m; else print "NA" }')
  fi

  # --- assemble (times → HH:MM:SS)
  local s01_n s01_t s02_n s02_t s03_n s03_max s03_sum
  IFS=$'\t' read -r s01_n s01_t          <<< "$s01"
  IFS=$'\t' read -r s02_n s02_t          <<< "$s02"
  IFS=$'\t' read -r s03_n s03_max s03_sum <<< "$s03"

  printf '%s\t' \
    "$sample" "$stage" "$fastq_gib" \
    "$n_sort" "$sort_gib" "$sstats" "$(to_hms "$s01_t")" "$s01_n" \
    "$bqsr_gib" "$dups" "$isize" "$depth" "$(to_hms "$s02_t")" "$s02_n" \
    "${n_chr}/${CANON_N}" "$chr_gib" "$canon_gib" \
    "$(to_hms "$s03_max")" "$(to_hms "$s03_sum")" "$s03_n" \
    "$g_region" "$g_depth"
  printf '%s\n' "$g_ratio"
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

echo
echo "[$] Per-sample S01-S03 Summary >>"

# Resolve sample dirs: an input that is itself a sample dir is used as-is,
# otherwise its immediate subdirectories that look like sample dirs.
SAMPLE_DIRS=()
for in_dir in "${INPUTS[@]}"; do
  [ -d "$in_dir" ] || { echo "<ERROR> Not a directory: $in_dir"; exit 1; }
  in_dir="$(cd "$in_dir" && pwd)"
  if is_sample_dir "$in_dir"; then
    SAMPLE_DIRS+=("$in_dir")
  else
    for d in "$in_dir"/*/; do
      d="${d%/}"
      [ -d "$d" ] && is_sample_dir "$d" && SAMPLE_DIRS+=("$d")
    done
  fi
done
[ ${#SAMPLE_DIRS[@]} -gt 0 ] || { echo "<ERROR> No sample directories found under: ${INPUTS[*]}"; exit 1; }
echo "[i]  Samples: ${#SAMPLE_DIRS[@]}"
[ "$GVCF_DEPTH" != off ] && echo "[i]  GVCF depth: $GVCF_DEPTH ($([ "$GVCF_DEPTH" = fast ] && echo "$FAST_CHROM" || echo "chr1-22, -j $NJOBS"))"

{
  printf '%s\t' \
    sample stage fastq_gib \
    n_sort_bams sort_bam_gib raw_reads pct_mapped pct_proper_pair s01_time s01_runs \
    bqsr_bam_gib pct_dup median_insert mean_depth_auto mean_depth_chrX mean_depth_chrY mean_depth_chrM s02_time s02_runs \
    chrom_gvcfs chrom_gvcf_gib canon_gvcf_gib s03_time_max s03_time_sum s03_runs \
    gvcf_depth_region gvcf_depth
  printf '%s\n' gvcf_vs_mosdepth
  for d in "${SAMPLE_DIRS[@]}"; do
    echo "[&]  $(basename "$d")" >&2
    summarize_sample "$d"
  done
} > "${OUT_TSV}.tmp"

mv "${OUT_TSV}.tmp" "$OUT_TSV"
echo "[>]  $OUT_TSV"
echo "[i]  Stage counts: $(cut -f2 "$OUT_TSV" | tail -n +2 | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"

# <\MAIN> ---------------------------------------------------------------------

#EOF
