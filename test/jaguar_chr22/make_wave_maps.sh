#!/usr/bin/env bash
# =============================================================================
# Title : JAGUAR chr22 — build wave sample-maps for Step 04 (GenomicsDBImport)
# About : Splits the JAGUAR chr22 GVCF file list into N sequential "waves" and
#         writes one Step-04 sample-map (label <TAB> chrom_gvcf_dir) per wave.
#         Wave 1 is imported with `create`, waves 2..N with `update` — this is
#         how we exercise the incremental-import path more than once.
# Usage : make_wave_maps.sh [cohort_root] [file_list] [n_waves] [out_dir]
#           cohort_root — dir that CONTAINS  bam/<SAMPLE>/chrom_gvcf/...
#                         default: /mnt/data/amedina/mramirezc/JAGUAR_JVC (FENIX)
#                         run `ls` there — you should see a `bam/` folder
#           file_list   — default: ./JAGUAR_JVC_chr22_list.txt (next to this script)
#           n_waves     — default: 3
#           out_dir     — default: next to this script
# Out   : <out_dir>/jaguar_chr22_wave<K>.sample_map.tsv   (K = 1..n_waves)
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"

COHORT_ROOT="${1:-/mnt/data/amedina/mramirezc/JAGUAR_JVC}"   # JAGUAR test default
FILE_LIST="${2:-$HERE/JAGUAR_JVC_chr22_list.txt}"
N_WAVES="${3:-3}"
OUT_DIR="${4:-$HERE}"

CHROM="chr22"
TINY_FRACTION_PCT=50   # flag any GVCF smaller than this % of the median size

# --- checks ----------------------------------------------------------------
[[ -d "$COHORT_ROOT" ]] || { echo "[X]  cohort_root not a directory: $COHORT_ROOT"; exit 1; }
COHORT_ROOT="$(readlink -f "$COHORT_ROOT")"
[[ -s "$FILE_LIST" ]]   || { echo "[X]  file list not found: $FILE_LIST"; exit 1; }
[[ -d "$COHORT_ROOT/bam" ]] || echo "[!]  WARNING: no 'bam/' under $COHORT_ROOT — check the path"
mkdir -p "$OUT_DIR"

# --- parse sample IDs from the list, preserving order --------------------------
# list lines look like:  -rw-r--r-- 1 u g SIZE DATE ./bam/<SAMPLE>/chrom_gvcf/<SAMPLE>.raw_vars.chr22.g.vcf.gz
SAMPLES=()
while IFS= read -r _id; do
  [[ -n "$_id" ]] && SAMPLES+=("$_id")
done < <(
  awk -v c="$CHROM" '
    $0 ~ "chrom_gvcf/" && $0 ~ ("\\." c "\\.g\\.vcf\\.gz") {
      n = split($NF, a, "/"); id = a[n]; sub(/\.raw_vars\..*/, "", id); print id
    }' "$FILE_LIST"
)
N=${#SAMPLES[@]}
(( N > 0 )) || { echo "[X]  no chr22 GVCF entries found in $FILE_LIST"; exit 1; }

# --- median size for the tiny-file flag -------------------------------------
MEDIAN=$(
  awk -v c="$CHROM" '$0 ~ "chrom_gvcf/" && $0 ~ ("\\." c "\\.g\\.vcf\\.gz") {print $5}' "$FILE_LIST" \
    | sort -n | awk '{v[NR]=$1} END{print (NR%2)? v[(NR+1)/2] : int((v[NR/2]+v[NR/2+1])/2)}'
)
TINY_CUTOFF=$(( MEDIAN * TINY_FRACTION_PCT / 100 ))

echo "[i]  cohort_root : $COHORT_ROOT"
echo "[i]  file list   : $FILE_LIST"
echo "[i]  samples     : $N   waves: $N_WAVES"
echo "[i]  median GVCF : $MEDIAN bytes   tiny-flag cutoff: <$TINY_CUTOFF bytes"
echo

# --- split into waves (ceil) and write maps -------------------------------------
PER=$(( (N + N_WAVES - 1) / N_WAVES ))
declare -a WAVE_IDS
missing=0 ; tiny=0

for (( w=1; w<=N_WAVES; w++ )); do
  start=$(( (w-1) * PER ))
  (( start >= N )) && { echo "[!]  wave $w is empty (fewer samples than waves) — skipping"; continue; }
  map="$OUT_DIR/jaguar_${CHROM}_wave${w}.sample_map.tsv"
  : > "$map"
  count=0
  for (( i=start; i<start+PER && i<N; i++ )); do
    s="${SAMPLES[i]}"
    dir="$COHORT_ROOT/bam/$s/chrom_gvcf"
    gvcf="$dir/$s.raw_vars.$CHROM.g.vcf.gz"
    printf '%s\t%s\n' "$s" "$dir" >> "$map"
    (( ++count ))
    if [[ ! -e "$gvcf" ]]; then
      echo "[!]  wave $w: MISSING  $gvcf"
      missing=1
    else
      sz=$(stat -c%s "$gvcf" 2>/dev/null || stat -f%z "$gvcf")
      if (( sz < TINY_CUTOFF )); then
        echo "[!]  wave $w: TINY ($sz b, likely truncated)  $s"
        tiny=1
      fi
    fi
  done
  echo "[>]  $map   ($count samples)"
  WAVE_IDS+=("$w")
done

echo
echo "[i]  Wave maps written to: $OUT_DIR"
echo "[i]  Action per wave      : wave 1 = create , waves 2+ = update"
(( missing == 0 )) || echo "[!]  Some GVCFs are missing — fix the cohort_root or re-run Step 03a before launching."
(( tiny == 0 ))    || echo "[!]  Tiny/truncated GVCF(s) detected — Step 04 will warn (or fail in strict mode). See notes."
echo "[i]  Next: bash test/jaguar_chr22/check_setup.sh $OUT_DIR"
