#!/usr/bin/env bash
# =============================================================================
# Title : JAGUAR chr22 — pre-flight check for Step 04 (GenomicsDBImport)
# About : Fast read-only sanity sweep to run on FENIX BEFORE launching the
#         wave jobs. Prints PASS / WARN / FAIL lines; exits non-zero if any FAIL.
#         Nothing here writes data or submits jobs.
# Usage : bash check_setup.sh [maps_dir] [output_dir]
#           maps_dir   — where jaguar_chr22_wave*.sample_map.tsv live
#                        (default: next to this script)
#           output_dir — where Step 04 will create genomicsdb/<CHR>
#                        (default: <repo>/test/jaguar_chr22/out)
# =============================================================================

set -uo pipefail   # NOT -e: we want to run every check and tally at the end

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
MAPS_DIR="${1:-$HERE}"
OUTPUT_DIR="${2:-/mnt/data/amedina/${USER:-esalazarf}/JVCdev/test_gendbi_jag22}"

S04="$REPO/bin/04_gatk_GenomicsDB_import.sh"
CONFIG="$REPO/config/config.yaml"
CHROM="chr22"
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"

pass=0 warn=0 fail=0
P(){ echo "  [PASS] $*"; ((pass++)); }
W(){ echo "  [WARN] $*"; ((warn++)); }
F(){ echo "  [FAIL] $*"; ((fail++)); }
hd(){ echo; echo "== $* =="; }

# --- environment ----------------------------------------------------------------
hd "Environment"
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then P "SSH session (remote/FENIX) detected"; else W "not an SSH session — Step 04 will treat this as 'local'"; fi
if command -v sbatch >/dev/null; then P "sbatch found ($(command -v sbatch))"; else F "sbatch not found — cannot submit wave jobs"; fi
if command -v module >/dev/null 2>&1 || type module >/dev/null 2>&1; then P "'module' command available"; else W "'module' not available in this shell (ok if tools already on PATH)"; fi

# --- tools --------------------------------------------------------------------
hd "Tools (loading modules like Step 04 does)"
if type module >/dev/null 2>&1; then
  module unload oracle-java 2>/dev/null || true
  module load oracle-java/25.0.2 2>/dev/null || W "could not 'module load oracle-java/25.0.2'"
  module load gatk 2>/dev/null      || W "could not 'module load gatk'"
  module load bcftools 2>/dev/null  || W "could not 'module load bcftools'"
fi
if command -v gatk >/dev/null; then
  v=$(gatk --version 2>&1 | grep -m1 -iE 'GATK|Toolkit' || echo '?')
  P "gatk on PATH — $v"
else
  F "gatk not on PATH"
fi
if command -v bcftools >/dev/null; then P "bcftools on PATH — $(bcftools --version 2>/dev/null | head -1)"; else F "bcftools not on PATH"; fi
if command -v java >/dev/null; then P "java — $(java -version 2>&1 | head -1)"; else W "java not on PATH (gatk wrapper may still find one)"; fi

# --- Step 04 script ----------------------------------------------------------
hd "Step 04 script"
if [[ -f "$S04" ]]; then
  P "found: $S04"
  [[ -x "$S04" ]] && P "executable" || W "not executable (run: chmod +x '$S04')  — 'bash $S04' still works"
  if bash -n "$S04" 2>/tmp/s04nchk.$$; then P "bash -n clean"; else F "bash -n errors: $(cat /tmp/s04nchk.$$)"; fi
  rm -f /tmp/s04nchk.$$
else
  F "not found: $S04"
fi

# --- config + reference ----------------------------------------------------------
hd "Config + reference"
if [[ -f "$CONFIG" ]]; then
  P "config: $CONFIG"
  env_type=remote; [[ -z "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]] && env_type=local
  eval "$(
    awk -v env="$env_type" '
      BEGIN{in_env=0}
      $1 ~ env":" {in_env=1;next}
      in_env && /^[^[:space:]]/ {in_env=0}
      in_env && /^[[:space:]]+[a-zA-Z0-9_]+:/ {
        gsub(":","=",$1); sub(/^[[:space:]]+/,"",$1)
        gsub(/^"/,"",$2); gsub(/"$/,"",$2); print $1 $2
      }' "$CONFIG"
  )"
  if [[ -n "${ref_gnm:-}" && -f "${ref_gnm}" ]]; then
    P "ref_gnm: $ref_gnm"
    [[ -f "${ref_gnm}.fai" ]]                 && P ".fai present"  || F "missing ${ref_gnm}.fai (samtools faidx)"
    dict="${ref_gnm%.*}.dict"
    [[ -f "$dict" || -f "${ref_gnm}.dict" ]]  && P ".dict present" || F "missing sequence dictionary (gatk CreateSequenceDictionary)"
    case "$ref_gnm" in *.gz) F "ref_gnm is bgzipped — GATK -R needs a plain .fa";; esac
  else
    F "ref_gnm missing or file not found (${ref_gnm:-unset})"
  fi
  if [[ -n "${scratch_root:-}" ]]; then
    if [[ -d "$scratch_root" ]]; then
      td="$scratch_root/${USER:-nouser}/.jvc_write_test.$$"
      if mkdir -p "$td" 2>/dev/null && : > "$td/x" 2>/dev/null; then P "scratch_root writable: $scratch_root"; rm -rf "$td"; else W "scratch_root not writable: $scratch_root"; fi
    else
      W "scratch_root set but not a directory: $scratch_root (Step 04 will fall back to output dir)"
    fi
  else
    W "scratch_root not set for '$env_type' — workspace builds under the output dir (slower TileDB I/O on NFS)"
  fi
else
  F "config not found: $CONFIG"
fi

# --- wave maps + every referenced GVCF -----------------------------------------
hd "Wave maps + GVCFs"
shopt -s nullglob
maps=( "$MAPS_DIR"/jaguar_${CHROM}_wave*.sample_map.tsv )
shopt -u nullglob
if (( ${#maps[@]} == 0 )); then
  F "no wave maps in $MAPS_DIR — run make_wave_maps.sh first"
else
  # pre-pass: median GVCF size across all waves, for the "unusually small" flag
  szfile="$(mktemp)"
  for m in "${maps[@]}"; do
    while IFS=$'\t' read -r label dir _ || [[ -n "${label:-}" ]]; do
      [[ -z "${label// }" || "$label" == \#* ]] && continue
      g="$dir/$label.raw_vars.$CHROM.g.vcf.gz"
      [[ -f "$g" ]] && { stat -c%s "$g" 2>/dev/null || stat -f%z "$g"; } >> "$szfile"
    done < "$m"
  done
  MEDIAN=$(sort -n "$szfile" | awk '{v[NR]=$1} END{if(NR)print (NR%2)?v[(NR+1)/2]:int((v[NR/2]+v[NR/2+1])/2); else print 0}')
  SMALL_CUTOFF=$(( MEDIAN / 2 ))
  rm -f "$szfile"
  echo "  (median GVCF ${MEDIAN} b; flagging < ${SMALL_CUTOFF} b as unusually small)"

  seen_all=" "
  for m in "${maps[@]}"; do
    n=$(grep -cvE '^[[:space:]]*(#|$)' "$m")
    echo "  --- $(basename "$m")  ($n samples) ---"
    dupes=0 miss=0 noidx=0 trunc=0 badsm=0 small=0
    while IFS=$'\t' read -r label dir _ || [[ -n "${label:-}" ]]; do
      [[ -z "${label// }" || "$label" == \#* ]] && continue
      g="$dir/$label.raw_vars.$CHROM.g.vcf.gz"
      if [[ ! -f "$g" ]]; then F "  $label: GVCF not found ($g)"; ((miss++)); continue; fi
      sz=$(stat -c%s "$g" 2>/dev/null || stat -f%z "$g")
      [[ -s "$g.tbi" || -s "$g.csi" ]] || { W "  $label: no .tbi/.csi index"; ((noidx++)); }
      tailhex=$(tail -c 28 "$g" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')
      if [[ "$tailhex" != "$BGZF_EOF" ]]; then
        W "  $label: no BGZF EOF marker — TRUNCATED ($sz b)"; ((trunc++))
      elif (( MEDIAN > 0 && sz < SMALL_CUTOFF )); then
        W "  $label: valid bgzip but UNUSUALLY SMALL ($sz b vs ${MEDIAN} b median) — near-empty GVCF? check the sample's chr22 coverage"; ((small++))
      fi
      sm=$(bcftools query -l "$g" 2>/dev/null)
      if [[ $(printf '%s\n' "$sm" | grep -c .) -ne 1 ]]; then W "  $label: GVCF header does not have exactly one sample"; ((badsm++)); fi
      [[ "$sm" != "$label" && -n "$sm" ]] && echo "        note: header sample '$sm' != map label '$label' (Step 04 uses the header name)"
      case "$seen_all" in *" ${sm:-$label} "*) F "  duplicate sample across waves: ${sm:-$label}"; ((dupes++));; esac
      seen_all="$seen_all${sm:-$label} "
    done < "$m"
    (( miss+dupes == 0 )) && P "  $(basename "$m"): all GVCFs present, no cross-wave dupes"
    (( noidx == 0 ))      || echo "        ($noidx missing index)"
    (( trunc+small == 0 )) && P "  $(basename "$m"): no truncated or unusually-small GVCFs" \
                           || echo "        ($trunc truncated, $small unusually small — Step 04 warns; GENDBI_STRICT_GVCF=true hard-fails on truncation)"
  done
fi

# --- output dir + space -------------------------------------------------------
hd "Output location + free space"
case "$OUTPUT_DIR/" in
  "${HOME}"/*) F "output_dir is under \$HOME ($HOME) — FENIX home has a hard ~5 GB quota; the copy-back WILL fail after the ~1 h import. Use group storage (/mnt/data/amedina/$USER/...)." ;;
esac
if mkdir -p "$OUTPUT_DIR/genomicsdb" 2>/dev/null && \
   dd if=/dev/zero of="$OUTPUT_DIR/genomicsdb/.writetest.$$" bs=1M count=64 status=none 2>/dev/null; then
  P "output dir writable (64 MB test ok): $OUTPUT_DIR"
  rm -f "$OUTPUT_DIR/genomicsdb/.writetest.$$"
else
  rm -f "$OUTPUT_DIR/genomicsdb/.writetest.$$" 2>/dev/null || true
  F "cannot write 64 MB to $OUTPUT_DIR/genomicsdb (quota / permissions)"
fi
if [[ -e "$OUTPUT_DIR/genomicsdb/$CHROM/callset.json" ]]; then
  W "genomicsdb/$CHROM already exists — wave 1 (create) will SKIP it; delete it for a clean first run"
fi
avail_g=$(df -Pk "$OUTPUT_DIR" 2>/dev/null | awk 'NR==2{print int($4/1048576)}')
[[ -n "${avail_g:-}" ]] && { (( avail_g >= 5 )) && P "output fs has ~${avail_g} G free" || W "output fs only ~${avail_g} G free — a full-cohort per-chrom workspace can be several GB"; }
df -h "$OUTPUT_DIR" 2>/dev/null | awk 'NR==1||NR==2{print "        "$0}'
[[ -n "${scratch_root:-}" && -d "${scratch_root:-/nonexistent}" ]] && df -h "$scratch_root" 2>/dev/null | awk 'NR==2{print "        scratch: "$0}'

# --- summary ----------------------------------------------------------------
hd "Summary"
echo "  PASS=$pass  WARN=$warn  FAIL=$fail"
if (( fail > 0 )); then
  echo "  -> resolve FAILs before launching."
  exit 1
elif (( warn > 0 )); then
  echo "  -> OK to launch; review WARNs (truncated GVCF is expected for this test)."
  exit 0
else
  echo "  -> all clear."
  exit 0
fi
