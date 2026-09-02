#!/usr/bin/env bash
# =============================================================================
# Title : JAGUAR — build wave sample-maps for Step 04 (GenomicsDBImport)
# About : Splits a cohort's per-chromosome GVCFs into N sequential "waves" and
#         writes one Step-04 sample-map (label <TAB> chrom_gvcf_dir) per wave.
#         Wave 1 is imported with `create`, waves 2..N with `update` — this
#         exercises the incremental-import path.
#
# Usage : make_wave_maps.sh [options]
#   -c, --chrom  CHR     chromosome (default: chr22)
#   -w, --waves  N       number of waves (default: 3)
#       --cohort DIR     dir containing bam/<SAMPLE>/chrom_gvcf/...
#                        (default: /mnt/data/amedina/mramirezc/JAGUAR_JVC)
#       --out    DIR     where to write the maps (default: next to this script)
#       --list   FILE    optional `ls -l` listing to take the sample order +
#                        sizes from; otherwise the cohort is scanned directly
#   -h, --help
#
# Out : <out>/jaguar_<chrom>_wave<K>.sample_map.tsv     (K = 1..waves)
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

CHROM="chr22"
N_WAVES=3
COHORT_ROOT="/mnt/data/amedina/mramirezc/JAGUAR_JVC"
OUT_DIR="$HERE"
FILE_LIST=""
TINY_FRACTION_PCT=50   # flag any GVCF smaller than this % of the median size

die() { echo "[X]  $*" >&2; exit 1; }

while (( $# )); do
  case "$1" in
    -c|--chrom)  CHROM="${2:?}"; shift 2 ;;
    -w|--waves)  N_WAVES="${2:?}"; shift 2 ;;
    --cohort)    COHORT_ROOT="${2:?}"; shift 2 ;;
    --out)       OUT_DIR="${2:?}"; shift 2 ;;
    --list)      FILE_LIST="${2:?}"; shift 2 ;;
    -h|--help)   sed -n '3,22p' "$0" | sed 's/^#\s\{0,1\}//'; exit 0 ;;
    *)           die "unknown argument: $1  (--help)" ;;
  esac
done

[[ -d "$COHORT_ROOT" ]] || die "cohort not a directory: $COHORT_ROOT"
COHORT_ROOT="$(readlink -f "$COHORT_ROOT")"
[[ -d "$COHORT_ROOT/bam" ]] || echo "[!]  WARNING: no 'bam/' under $COHORT_ROOT"
mkdir -p "$OUT_DIR"

# auto-pick a bundled list if none given and one exists for this chromosome
[[ -z "$FILE_LIST" && -s "$HERE/JAGUAR_JVC_${CHROM}_list.txt" ]] && FILE_LIST="$HERE/JAGUAR_JVC_${CHROM}_list.txt"

# --- collect samples (ordered) and their GVCF sizes --------------------------
SAMPLES=() ; SIZES=()
if [[ -n "$FILE_LIST" ]]; then
  [[ -s "$FILE_LIST" ]] || die "file list not found: $FILE_LIST"
  echo "[i]  source     : list  $FILE_LIST"
  while IFS=$'\t' read -r _id _sz; do
    [[ -n "$_id" ]] && { SAMPLES+=("$_id"); SIZES+=("$_sz"); }
  done < <(
    awk -v c="$CHROM" '
      $0 ~ "chrom_gvcf/" && $0 ~ ("\\." c "\\.g\\.vcf\\.gz$") {
        n = split($NF, a, "/"); id = a[n]; sub(/\.raw_vars\..*/, "", id)
        print id "\t" $5
      }' "$FILE_LIST"
  )
else
  echo "[i]  source     : scan  $COHORT_ROOT/bam/*/chrom_gvcf/*.${CHROM}.g.vcf.gz"
  while IFS= read -r g; do
    s="$(basename "$(dirname "$(dirname "$g")")")"
    sz=$(stat -c%s "$g" 2>/dev/null || stat -f%z "$g" 2>/dev/null || echo 0)
    SAMPLES+=("$s"); SIZES+=("$sz")
  done < <(find "$COHORT_ROOT/bam" -mindepth 3 -maxdepth 3 -type f \
               -name "*.${CHROM}.g.vcf.gz" 2>/dev/null | sort)
fi

N=${#SAMPLES[@]}
(( N > 0 )) || die "no ${CHROM} GVCFs found"

MEDIAN=$(printf '%s\n' "${SIZES[@]}" | sort -n | awk '{v[NR]=$1} END{print (NR%2)? v[(NR+1)/2] : int((v[NR/2]+v[NR/2+1])/2)}')
TINY_CUTOFF=$(( MEDIAN * TINY_FRACTION_PCT / 100 ))

echo "[i]  chromosome : $CHROM"
echo "[i]  cohort     : $COHORT_ROOT"
echo "[i]  samples    : $N   waves: $N_WAVES"
echo "[i]  median GVCF : $MEDIAN bytes   tiny-flag cutoff: <$TINY_CUTOFF bytes"
echo

# --- split into waves (ceil) and write maps ---------------------------------
PER=$(( (N + N_WAVES - 1) / N_WAVES ))
missing=0 ; tiny=0

for (( w=1; w<=N_WAVES; w++ )); do
  start=$(( (w-1) * PER ))
  (( start >= N )) && { echo "[!]  wave $w empty (more waves than samples) — skipped"; continue; }
  map="$OUT_DIR/jaguar_${CHROM}_wave${w}.sample_map.tsv"
  : > "$map"
  count=0
  for (( i=start; i<start+PER && i<N; i++ )); do
    s="${SAMPLES[i]}" ; sz="${SIZES[i]:-0}"
    dir="$COHORT_ROOT/bam/$s/chrom_gvcf"
    gvcf="$dir/$s.raw_vars.$CHROM.g.vcf.gz"
    printf '%s\t%s\n' "$s" "$dir" >> "$map"
    (( ++count ))
    if [[ ! -e "$gvcf" ]]; then
      echo "[!]  wave $w: MISSING  $gvcf" ; missing=1
    elif (( MEDIAN > 0 && sz > 0 && sz < TINY_CUTOFF )); then
      echo "[!]  wave $w: TINY ($sz b vs $MEDIAN b median)  $s" ; tiny=1
    fi
  done
  echo "[>]  $map   ($count samples)"
done

echo
echo "[i]  Wave maps in : $OUT_DIR"
echo "[i]  Action/wave  : wave 1 = create , waves 2+ = update"
(( missing == 0 )) || echo "[!]  Missing GVCFs — fix --cohort or re-run Step 03a before launching."
(( tiny == 0 ))    || echo "[!]  Tiny GVCF(s) flagged — Step 04 will warn (near-empty / truncated)."
echo "[i]  Next: bash $HERE/check_setup.sh --chrom $CHROM"
