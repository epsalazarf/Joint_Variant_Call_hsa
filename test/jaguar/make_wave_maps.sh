#!/usr/bin/env bash
# =============================================================================
# Title : Step 04 test harness — build wave sample-maps for GenomicsDBImport
# About : Splits a cohort's per-chromosome GVCFs into N sequential "waves" and
#         writes one Step-04 sample-map (label <TAB> chrom_gvcf_dir) per wave.
#         Wave 1 is imported with `create`, waves 2..N with `update` — this
#         exercises the incremental-import path.
#
#         Sample discovery SCANS the cohort recursively for any
#         `.../chrom_gvcf/*.<chrom>.g.vcf.gz` file, so it does not assume a
#         particular directory depth (JAGUAR: bam/<S>/chrom_gvcf/..., SLE:
#         BAMQC/BatchNN/<S>/chrom_gvcf/...) or that the GVCF filename starts
#         with the sample/directory name (SLE file prefixes vary: L46_1,
#         Q005_D_2, Q00a_D_1, ...). The label is always the chrom_gvcf/'s
#         parent directory name; Step 04 itself re-derives the REAL sample name
#         from the GVCF header, so a label/header mismatch is just a warning.
#
# Usage : make_wave_maps.sh [options]
#   -c, --chrom  CHR     chromosome (default: chr22)
#   -w, --waves  N       number of waves (default: 3)
#   -t, --tag    NAME    cohort tag — namespaces the output files so different
#                        cohorts don't collide (default: jaguar)
#       --cohort DIR     dir to scan for .../chrom_gvcf/*.<chrom>.g.vcf.gz
#                        (default: /mnt/data/amedina/mramirezc/JAGUAR_JVC)
#       --out    DIR     where to write the maps (default: next to this script)
#       --list   FILE    optional `ls -l`/`find -ls` listing to take the sample
#                        order + sizes from instead of scanning live (sizes are
#                        only used for the tiny-file heuristic)
#   -h, --help
#
# Out : <out>/<tag>_<chrom>_wave<K>.sample_map.tsv     (K = 1..waves)
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

CHROM="chr22"
N_WAVES=3
TAG="jaguar"
COHORT_ROOT="/mnt/data/amedina/mramirezc/JAGUAR_JVC"
OUT_DIR="$HERE"
FILE_LIST=""
TINY_FRACTION_PCT=50   # flag any GVCF smaller than this % of the median size

die() { echo "[X]  $*" >&2; exit 1; }

while (( $# )); do
  case "$1" in
    -c|--chrom)  CHROM="${2:?}"; shift 2 ;;
    -w|--waves)  N_WAVES="${2:?}"; shift 2 ;;
    -t|--tag)    TAG="${2:?}"; shift 2 ;;
    --cohort)    COHORT_ROOT="${2:?}"; shift 2 ;;
    --out)       OUT_DIR="${2:?}"; shift 2 ;;
    --list)      FILE_LIST="${2:?}"; shift 2 ;;
    -h|--help)   sed -n '3,29p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
    *)           die "unknown argument: $1  (--help)" ;;
  esac
done

[[ -d "$COHORT_ROOT" ]] || die "cohort not a directory: $COHORT_ROOT"
COHORT_ROOT="$(readlink -f "$COHORT_ROOT")"
mkdir -p "$OUT_DIR"

# auto-pick a bundled list if none given and one exists for this tag+chromosome
[[ -z "$FILE_LIST" && -s "$HERE/${TAG}_${CHROM}_list.txt" ]] && FILE_LIST="$HERE/${TAG}_${CHROM}_list.txt"
[[ -z "$FILE_LIST" && -s "$HERE/JAGUAR_JVC_${CHROM}_list.txt" && "$TAG" == "jaguar" ]] && FILE_LIST="$HERE/JAGUAR_JVC_${CHROM}_list.txt"

# --- collect samples (ordered), their chrom_gvcf dirs, and GVCF sizes --------
SAMPLES=() ; DIRS=() ; SIZES=()
if [[ -n "$FILE_LIST" ]]; then
  [[ -s "$FILE_LIST" ]] || die "file list not found: $FILE_LIST"
  echo "[i]  source     : list  $FILE_LIST"
  while IFS=$'\t' read -r _id _dir _sz; do
    [[ -n "$_id" ]] && { SAMPLES+=("$_id"); DIRS+=("$_dir"); SIZES+=("$_sz"); }
  done < <(
    awk -v c="$CHROM" -v root="$COHORT_ROOT" '
      $0 ~ "chrom_gvcf/" && $0 ~ ("\\." c "\\.g\\.vcf\\.gz$") {
        path = $NF
        sub(/^\.\//, root "/", path)          # ./SLEmx-b38/... -> <root>/SLEmx-b38/...
        n = split(path, a, "/")
        id = a[n-1]                            # parent-of-chrom_gvcf dir name
        dir = ""
        for (k = 1; k < n; k++) dir = dir a[k] (k < n-1 ? "/" : "")
        print id "\t" dir "\t" $5
      }' "$FILE_LIST"
  )
else
  echo "[i]  source     : scan  $COHORT_ROOT/**/chrom_gvcf/*.${CHROM}.g.vcf.gz"
  while IFS= read -r g; do
    d="$(dirname "$g")"
    s="$(basename "$(dirname "$d")")"
    sz=$(stat -c%s "$g" 2>/dev/null || stat -f%z "$g" 2>/dev/null || echo 0)
    SAMPLES+=("$s"); DIRS+=("$d"); SIZES+=("$sz")
  done < <(find "$COHORT_ROOT" -type f -path '*/chrom_gvcf/*' \
               -name "*.${CHROM}.g.vcf.gz" 2>/dev/null | sort)
fi

N=${#SAMPLES[@]}
(( N > 0 )) || die "no ${CHROM} GVCFs found under $COHORT_ROOT"

# duplicate label guard (two different chrom_gvcf dirs, same parent dir name)
dupes=$(printf '%s\n' "${SAMPLES[@]}" | sort | uniq -d)
if [[ -n "$dupes" ]]; then
  die "duplicate sample label(s) found in different directories — rename or exclude:
$dupes"
fi

MEDIAN=$(printf '%s\n' "${SIZES[@]}" | sort -n | awk '{v[NR]=$1} END{print (NR%2)? v[(NR+1)/2] : int((v[NR/2]+v[NR/2+1])/2)}')
TINY_CUTOFF=$(( MEDIAN * TINY_FRACTION_PCT / 100 ))

echo "[i]  chromosome : $CHROM"
echo "[i]  cohort     : $COHORT_ROOT"
echo "[i]  tag        : $TAG"
echo "[i]  samples    : $N   waves: $N_WAVES"
echo "[i]  median GVCF : $MEDIAN bytes   tiny-flag cutoff: <$TINY_CUTOFF bytes"
echo

# --- split into waves (ceil) and write maps ---------------------------------
PER=$(( (N + N_WAVES - 1) / N_WAVES ))
missing=0 ; tiny=0

for (( w=1; w<=N_WAVES; w++ )); do
  start=$(( (w-1) * PER ))
  (( start >= N )) && { echo "[!]  wave $w empty (more waves than samples) — skipped"; continue; }
  map="$OUT_DIR/${TAG}_${CHROM}_wave${w}.sample_map.tsv"
  : > "$map"
  count=0
  for (( i=start; i<start+PER && i<N; i++ )); do
    s="${SAMPLES[i]}" ; dir="${DIRS[i]}" ; sz="${SIZES[i]:-0}"
    printf '%s\t%s\n' "$s" "$dir" >> "$map"
    (( ++count ))
    shopt -s nullglob
    hits=( "$dir"/*."$CHROM".g.vcf.gz )
    shopt -u nullglob
    if [[ ! -d "$dir" ]]; then
      echo "[!]  wave $w: MISSING dir  $dir" ; missing=1
    elif (( ${#hits[@]} == 0 )); then
      echo "[!]  wave $w: MISSING gvcf in  $dir" ; missing=1
    elif (( ${#hits[@]} > 1 )); then
      echo "[!]  wave $w: AMBIGUOUS — ${#hits[@]} matches in $dir" ; missing=1
    elif (( MEDIAN > 0 && sz > 0 && sz < TINY_CUTOFF )); then
      echo "[!]  wave $w: TINY ($sz b vs $MEDIAN b median)  $s" ; tiny=1
    fi
  done
  echo "[>]  $map   ($count samples)"
done

echo
echo "[i]  Wave maps in : $OUT_DIR"
echo "[i]  Action/wave  : wave 1 = create , waves 2+ = update"
(( missing == 0 )) || echo "[!]  Missing/ambiguous GVCFs — fix --cohort or re-run Step 03a before launching."
(( tiny == 0 ))    || echo "[!]  Tiny GVCF(s) flagged — Step 04 will warn (near-empty / truncated)."
echo "[i]  Next: bash $HERE/check_setup.sh --chrom $CHROM --tag $TAG"
