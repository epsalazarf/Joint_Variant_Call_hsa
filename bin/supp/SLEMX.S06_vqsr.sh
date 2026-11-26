#!/usr/bin/env bash

# =============================================================================
# Title       : SLEmx (lupus) — Step 06 VQSR filtering launcher [FENIX]
# Description : Checks the Step 05 VCFs and the VQSR resources, then submits
#               bin/06_gatk_vqsr.sh as a chain of SLURM jobs:
#                 prep   array, one task per chromosome  ExcessHet + sites-only
#                 model  2 jobs (SNP, INDEL), after prep  ONE genome-wide model each
#                 apply  array, after both models        final <out>/slemx.chrN.vqsr.vcf.gz
#               Anything already done is skipped (and not waited for), so
#               re-running after a failure submits only what is missing.
#
#               --chroms all       the real run: model on all 22 autosomes,
#                                  output <cohort>/vqsr/
#               --chroms chr22     a TEST: model on just those chromosomes, in its
#                                  own dir <cohort>/vqsr_test-chr22/ (never mixed
#                                  with the real model). chr22 alone is enough
#                                  sites for a smoke test of the whole chain.
#               --method hardfilter  GATK generic hard filters instead (no model;
#                                  for cohorts too small for VQSR) -> <cohort>/hardfilter/
#
#               --check  : Step 05 VCFs complete + same samples, resources
#                          staged, dbSNP, free space. Nothing submitted.
#               --status : what is done per stage / chromosome, queued jobs.
#               --report : after it ends — tranches of both models and the
#                          FILTER counts per chromosome (from the apply logs).
#
# Usage       : bash SLEMX.S06_vqsr.sh --check
#               bash SLEMX.S06_vqsr.sh --chroms chr22 [--dry-run]       # test
#               bash SLEMX.S06_vqsr.sh --chroms all   [--dry-run]       # real run
#               bash SLEMX.S06_vqsr.sh --status | --report  [--chroms ...]
#                 --method vqsr|hardfilter   (default vqsr)
#                 --out DIR   cohort root (default /mnt/data/amedina/esalazarf/SLEmx_JVC)
#                             reads <out>/genotyped/slemx.chrN.vcf.gz
# =============================================================================

set -euo pipefail

# <CONFIG> ----------------------------------------------------------------
OUTPUT_ROOT="/mnt/data/amedina/esalazarf/SLEmx_JVC"
TAG="slemx"
DATA_GROUP="amedina"
AUTOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
# Resources (first guesses — no measurement yet; check sacct after the chr22
# test and after the real run). prep/apply stream the VCF: little memory.
# VariantRecalibrator holds every site of its type in the Java heap: ~18M SNPs
# genome-wide for 164 samples -> large heap. defq has no time limit.
PREP_MEM="12G";  PREP_JAVA="6G";  PREP_HOURS=12
APPLY_MEM="12G"; APPLY_JAVA="6G"; APPLY_HOURS=12
SNP_MEM="64G";   SNP_JAVA="56G";  SNP_HOURS=72
INDEL_MEM="32G"; INDEL_JAVA="26G"; INDEL_HOURS=48
HF_MEM="12G";    HF_JAVA="6G";    HF_HOURS=12
# <\CONFIG> -----------------------------------------------------------------

SELF="$(readlink -f "$0")"
HERE="$(cd "$(dirname "$SELF")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
S06="$REPO/bin/06_gatk_vqsr.sh"
CONFIG_FILE="$REPO/config/config.yaml"

MODE="submit"; CHROMS_ARG=""; METHOD="vqsr"
while [ $# -gt 0 ]; do
  case "$1" in
    --check)   MODE="check";  shift ;;
    --status)  MODE="status"; shift ;;
    --report)  MODE="report"; shift ;;
    --dry-run) MODE="dry";    shift ;;
    --chroms)  CHROMS_ARG="${2:?--chroms needs 'all' or a list, e.g. chr22}"; shift 2 ;;
    --method)  METHOD="${2:?--method needs vqsr or hardfilter}"; shift 2 ;;
    --out)     OUTPUT_ROOT="${2:?--out needs a directory}"; shift 2 ;;
    -h|--help) sed -n '4,37p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *)         echo "[X]  unknown argument: $1  (--help)"; exit 2 ;;
  esac
done
case "$METHOD" in vqsr|hardfilter) ;; *) echo "[X]  --method must be vqsr or hardfilter"; exit 2 ;; esac
OUTPUT_ROOT="${OUTPUT_ROOT%/}"
VCF_DIR="$OUTPUT_ROOT/genotyped"
LOG_DIR="$OUTPUT_ROOT/logs"

# Chromosomes + output dir: "all" = the real run; a subset = a test dir.
if [[ -z "$CHROMS_ARG" || "$CHROMS_ARG" == all ]]; then
  CHROMS="$AUTOSOMES"
else
  CHROMS="$(echo "$CHROMS_ARG" | tr ',' ' ')"
  for c in $CHROMS; do
    case " $AUTOSOMES " in *" $c "*) ;; *) echo "[X]  not an autosome: $c"; exit 2 ;; esac
  done
  # keep genome order, drop duplicates
  _sel=""; for a in $AUTOSOMES; do case " $CHROMS " in *" $a "*) _sel="${_sel:+$_sel }$a" ;; esac; done
  CHROMS="$_sel"
fi
CHROMS_LIST=( $CHROMS )
if [[ "$CHROMS" == "$AUTOSOMES" ]]; then OUT_NAME="$METHOD"
else OUT_NAME="${METHOD}_test-$(echo "$CHROMS" | tr ' ' '_')"; fi
OUT="$OUTPUT_ROOT/$OUT_NAME"
JOB_BASE="JVC-S06-${TAG}-${OUT_NAME}"
[[ "$MODE" == submit || "$MODE" == dry ]] && [[ -z "$CHROMS_ARG" ]] \
  && { echo "[X]  --chroms is required to submit: 'all' (real run) or e.g. chr22 (test)."; exit 2; }

in_vcf()    { echo "$VCF_DIR/${TAG}.$1.vcf.gz"; }
sites_vcf() { echo "$OUT/work/${TAG}.$1.sitesonly.vcf.gz"; }
exhet_vcf() { echo "$OUT/work/${TAG}.$1.exhet.vcf.gz"; }
final_vcf() { if [[ "$METHOD" == vqsr ]]; then echo "$OUT/${TAG}.$1.vqsr.vcf.gz"; else echo "$OUT/${TAG}.$1.hardfilt.vcf.gz"; fi; }
vcf_done()  { [[ -s "$1" && -s "$1.tbi" ]]; }
prep_done() { vcf_done "$(sites_vcf "$1")" && { vcf_done "$(exhet_vcf "$1")" || vcf_done "$(final_vcf "$1")"; }; }
model_done() {   # $1 snps|indels
  local m="$OUT/model/${TAG}.$1"
  vcf_done "$m.recal.vcf.gz" && [[ -s "$m.tranches" && -s "$m.model_chroms.txt" && "$(cat "$m.model_chroms.txt")" == "$CHROMS" ]]
}
queued_any() { command -v squeue >/dev/null && squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -q "^${JOB_BASE}-"; }
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
bgzf_complete() { [[ "$(tail -c 28 "$1" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')" == "$BGZF_EOF" ]]; }
tranche_table() {
  awk -F, '/^#/ {next} /^targetTruthSensitivity/ {for (i=1;i<=NF;i++) h[$i]=i; next}
           ("targetTruthSensitivity" in h) {printf "       %7s  known=%-9s novel=%-9s TiTv=%s/%s  minVQSLod=%s\n",
             $h["targetTruthSensitivity"], $h["numKnown"], $h["numNovel"], $h["knownTiTv"], $h["novelTiTv"], $h["minVQSLod"]}' "$1"
}

# --- status -------------------------------------------------------------------
if [[ "$MODE" == status ]]; then
  echo "[&]  ${TAG} Step 06 ($METHOD) status — $OUT"
  if [[ "$METHOD" == vqsr ]]; then
    for m in snps indels; do
      if model_done "$m"; then echo "  [DONE]    model $m"; else echo "  [ ---- ]  model $m"; fi
    done
  fi
  for c in $CHROMS; do
    if vcf_done "$(final_vcf "$c")"; then printf '  [DONE]    %-6s %s\n' "$c" "$(du -h "$(final_vcf "$c")" | cut -f1)"
    elif [[ "$METHOD" == vqsr ]] && prep_done "$c"; then printf '  [prep ok] %s\n' "$c"
    else printf '  [ ---- ]  %s\n' "$c"; fi
  done
  command -v squeue >/dev/null && squeue -u "$USER" -h -o '%.18i %.50j %.8T %.12M %R' | grep "$JOB_BASE" || true
  exit 0
fi

# --- report -------------------------------------------------------------------
if [[ "$MODE" == report ]]; then
  echo "[&]  ${TAG} Step 06 ($METHOD) report — $OUT"
  if [[ "$METHOD" == vqsr ]]; then
    for m in snps indels; do
      t="$OUT/model/${TAG}.$m.tranches"
      echo; echo "[i]  $m model ($(cat "$OUT/model/${TAG}.$m.model_chroms.txt" 2>/dev/null | wc -w | tr -d ' ') chromosomes) — truth sensitivity tranches:"
      if [[ -s "$t" ]]; then tranche_table "$t"; else echo "       (not built)"; fi
    done
    echo "[i]  Expect novel Ti/Tv ~2.0-2.1 for WGS SNPs at the 99.7 level; known ~2.1."
  fi
  echo; printf '  %-6s %10s %10s %6s  %s\n' chrom records PASS %PASS "top filters"
  tot=0; totp=0
  for c in $CHROMS; do
    f="$(final_vcf "$c")"
    if ! vcf_done "$f"; then printf '  %-6s %10s\n' "$c" "-"; continue; fi
    log="$(grep -l "^\[DONE\] $f\$" "$LOG_DIR"/${TAG}-S06-${OUT_NAME}-*.out 2>/dev/null | tail -1 || true)"
    if [[ -z "$log" ]]; then printf '  %-6s %10s  (no log found)\n' "$c" "?"; continue; fi
    n=$(awk '/^\[i\]  Samples:/ {print $5}' "$log" | tail -1)
    p=$(awk '/^\[i\]  FILTER:/ {f=1; next} f && $2=="PASS" {print $1; exit}' "$log")
    top=$(awk '/^\[i\]  FILTER:/ {f=1; next} f && /^ +[0-9]+  / && $2!="PASS" {printf "%s%s %s", (k++?", ":""), $2, $1} f && !/^ +[0-9]/ {f=0}' "$log" | cut -c1-70)
    printf '  %-6s %10s %10s %5s%%  %s\n' "$c" "$n" "${p:-0}" "$(awk -v a="${p:-0}" -v b="$n" 'BEGIN{printf "%.1f", b ? 100*a/b : 0}')" "$top"
    tot=$((tot + n)); totp=$((totp + ${p:-0}))
  done
  echo "  total  $(printf '%10s %10s' "$tot" "$totp")"
  exit 0
fi

# <CHECKS> ------------------------------------------------------------------
ERR=0
fail() { echo "[X]  $*"; ERR=$((ERR + 1)); }
warn() { echo "[!]  $*"; }

if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then env_type=remote; else env_type=local; fi
eval "$(
  awk -v env="$env_type" '
    $1 ~ env":" { in_env=1; next } in_env && /^[^[:space:]]/ { in_env=0 }
    in_env && /^[[:space:]]+[a-zA-Z0-9_]+:/ { gsub(":", "=", $1); sub(/^[[:space:]]+/, "", $1); gsub(/^"/, "", $2); gsub(/"$/, "", $2); print "cfg_" $1 $2 }
  ' "$CONFIG_FILE"
)"
BUNDLE="${VQSR_BUNDLE_DIR:-${cfg_ref_vqsr_bundle:-}}"
FROZEN_MAP="$OUTPUT_ROOT/${TAG}.sample_map.S04.tsv"       # the cohort as imported in Step 04
N_SAMPLES=""
[[ -s "$FROZEN_MAP" ]] && N_SAMPLES=$(grep -cvE '^[[:space:]]*(#|$)' "$FROZEN_MAP" || true)
DBSNP="${cfg_ref_vars:-}"
if ! command -v bcftools >/dev/null && command -v module >/dev/null 2>&1; then module load bcftools 2>/dev/null || true; fi
command -v bcftools >/dev/null || { echo "[X]  bcftools not found (module load bcftools)"; exit 1; }

echo "[&]  ${TAG} Step 06 check — $METHOD   ($(date))"
echo "[i]  input  : $VCF_DIR/${TAG}.<chrom>.vcf.gz"
echo "[i]  output : $OUT"
echo "[i]  samples: ${N_SAMPLES:-? (no $FROZEN_MAP — only checked for consistency)}"
echo "[i]  chroms : $(echo "$CHROMS" | wc -w | tr -d ' ') ($([[ "$CHROMS" == "$AUTOSOMES" ]] && echo "all autosomes — real run" || echo "TEST: $CHROMS"))"
[[ -f "$S06" ]] || { echo "[X]  Step 06 script not found: $S06"; exit 1; }

echo
ref_ck=""; in_kb=0; n_rec=0
for c in $CHROMS; do
  f="$(in_vcf "$c")"
  if ! vcf_done "$f"; then fail "$c: Step 05 VCF (+ .tbi) missing: $f"; continue; fi
  bgzf_complete "$f" || fail "$c: Step 05 VCF truncated (no BGZF EOF): $f"
  ck="$(bcftools query -l "$f" | cksum)"; ns=$(bcftools query -l "$f" | wc -l | tr -d ' ')
  [[ -z "$N_SAMPLES" || "$ns" == "$N_SAMPLES" ]] || fail "$c: $ns samples, S04 map has $N_SAMPLES"
  if [[ -z "$ref_ck" ]]; then ref_ck="$ck"; elif [[ "$ck" != "$ref_ck" ]]; then fail "$c: sample names/order differ from ${CHROMS_LIST[0]}"; fi
  n=$(bcftools index -n "$f" 2>/dev/null || echo 0); n_rec=$((n_rec + n))
  in_kb=$((in_kb + $(du -sk "$f" | cut -f1)))
done
(( ERR )) || echo "[OK] Step 05 VCFs: $(echo "$CHROMS" | wc -w | tr -d ' ') complete, ${ns:-?} samples each (same order), $n_rec records, $((in_kb / 1024 / 1024)) GB"

if [[ "$METHOD" == vqsr ]]; then
  miss=0
  for stem in hapmap_3.3.hg38 1000G_omni2.5.hg38 1000G_phase1.snps.high_confidence.hg38 \
              Mills_and_1000G_gold_standard.indels.hg38 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38; do
    vcf_done "$BUNDLE/$stem.canon_chr.vcf.gz" || { fail "VQSR resource missing: $BUNDLE/$stem.canon_chr.vcf.gz"; miss=1; }
  done
  (( miss )) && echo "     stage them first (login node, ~10 min): bash $REPO/bin/supp/STAGE.vqsr_bundle_hg38.sh" \
             || echo "[OK] VQSR resources: $BUNDLE"
  vcf_done "$DBSNP" && echo "[OK] dbSNP (known): $DBSNP" || fail "dbSNP (config ref_vars) missing: $DBSNP"
fi

# Space: ExcessHet copy (deleted after apply) + final = ~2x the input.
mkdir -p "$OUTPUT_ROOT" 2>/dev/null || true
free_kb=$(df -Pk "$OUTPUT_ROOT" | awk 'NR==2 {print $4}')
need_kb=$(( in_kb * 25 / 10 ))
if (( free_kb < need_kb )); then fail "free space $((free_kb / 1024 / 1024)) GB < ~$((need_kb / 1024 / 1024)) GB needed"
else echo "[OK] free space: $((free_kb / 1024 / 1024)) GB (need ~$((need_kb / 1024 / 1024)) GB)"; fi

shopt -s nullglob
tmps=( "$OUT"/*.tmp.vcf.gz "$OUT"/work/*.tmp.vcf.gz "$OUT"/model/*.tmp.* )
left=( "/scratch/groups/$(id -gn)/$USER"/vqsr_* )
shopt -u nullglob
(( ${#tmps[@]} )) && warn "${#tmps[@]} unfinished .tmp file(s) in $OUT (from killed jobs; rewritten on the next run)"
(( ${#left[@]} )) && warn "${#left[@]} leftover Step 06 scratch dir(s): rm -rf /scratch/groups/$(id -gn)/$USER/vqsr_*  (once no S06 job runs)"

echo
if (( ERR )); then echo "[X]  Check FAILED: $ERR problem(s) above. Nothing submitted."; exit 1; fi
echo "[OK] ready."
[[ "$MODE" == check ]] && { echo "[i]  --check: nothing submitted."; exit 0; }
# <\CHECKS> -----------------------------------------------------------------

# --- submit -------------------------------------------------------------------
if [[ "$MODE" != dry ]] && ! command -v sbatch >/dev/null; then echo "[X]  sbatch not found — run this on a FENIX login node."; exit 1; fi
if queued_any; then echo "[X]  ${JOB_BASE}-* jobs are already queued/running — wait, or scancel them first."; exit 1; fi
if [[ "$MODE" != dry ]]; then
  mkdir -p "$OUT/work" "$OUT/model" "$LOG_DIR"
  chgrp "$DATA_GROUP" "$OUT" "$OUT/work" "$OUT/model" 2>/dev/null || true
  chmod g+rwxs "$OUT" "$OUT/work" "$OUT/model" 2>/dev/null || true   # new files inherit the group
fi
exclude="$(awk '/^remote:/{r=1;next} /^[^[:space:]]/{r=0} r && $1=="sbatch_exclude:"{gsub(/"/,"",$2); print $2}' "$CONFIG_FILE")"
common=(--parsable --nodes=1 --ntasks=1 --cpus-per-task=2)
[[ -n "$exclude" ]] && common+=(--exclude="$exclude")
submit() { if [[ "$MODE" == dry ]]; then echo "DRY"; else sbatch "$@"; fi; }

# Everything a stage job needs; $1 java heap, then the stage + target.
stage_wrap() {   # $1 java mem  $2 stage  $3 target (may be the literal \$c)
  cat <<EOF
set -euo pipefail; umask 0027
CHRS=( $CHROMS )
c="\${CHRS[\$(( \${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}"
_scr="/scratch/groups/\$(id -gn)"
mkdir -p "\$_scr/\$USER" || { echo "<ERROR> scratch not writable on \$(hostname): \$_scr/\$USER" >&2; exit 1; }
echo "[i] $2 $3 on node \$(hostname)"
env VQSR_PREFIX='$TAG' VQSR_CHROMS='$CHROMS' VQSR_SCRATCH_ROOT="\$_scr" VQSR_JAVA_MEM='$1' \\
  bash '$S06' '$2' '$VCF_DIR' '$OUT' "$3"
EOF
}
idx_of() { local i=1 x; for x in $CHROMS; do [[ "$x" == "$1" ]] && { echo "$i"; return; }; i=$((i + 1)); done; }
ranges() {   # "1,2,3,5" -> "1-3,5"
  echo "$1" | tr ',' '\n' | awk 'NR==1 {s=p=$1; next} $1==p+1 {p=$1; next} {printf "%s%s,", s, (s==p ? "" : "-" p); s=p=$1} END {if (NR) printf "%s%s", s, (s==p ? "" : "-" p)}'
}

echo
echo "[&]  ${TAG} Step 06 — $METHOD   (jobs ${JOB_BASE}-*, logs ${LOG_DIR}/${TAG}-S06-${OUT_NAME}-*)"

if [[ "$METHOD" == hardfilter ]]; then
  todo=""; for c in $CHROMS; do vcf_done "$(final_vcf "$c")" || todo="${todo:+$todo,}$(idx_of "$c")"; done
  if [[ -z "$todo" ]]; then echo "[SKIP] all $(echo "$CHROMS" | wc -w | tr -d ' ') chromosomes already filtered"; exit 0; fi
  log="$LOG_DIR/${TAG}-S06-${OUT_NAME}-hardfilter-%A_%a.out"
  jid=$(submit "${common[@]}" --job-name="${JOB_BASE}-hardfilter" --array="$(ranges "$todo")" \
          --mem="$HF_MEM" --time="${HF_HOURS}:00:00" --output="$log" --wrap "$(stage_wrap "$HF_JAVA" hardfilter '$c')")
  echo "[*]  hardfilter : array $jid  [$(ranges "$todo")]   $HF_MEM / ${HF_HOURS}h"
else
  # prep
  todo=""; for c in $CHROMS; do prep_done "$c" || todo="${todo:+$todo,}$(idx_of "$c")"; done
  dep_prep=()
  if [[ -n "$todo" ]]; then
    log="$LOG_DIR/${TAG}-S06-${OUT_NAME}-prep-%A_%a.out"
    jid=$(submit "${common[@]}" --job-name="${JOB_BASE}-prep" --array="$(ranges "$todo")" \
            --mem="$PREP_MEM" --time="${PREP_HOURS}:00:00" --output="$log" --wrap "$(stage_wrap "$PREP_JAVA" prep '$c')")
    echo "[*]  prep        : array $jid  [$(ranges "$todo")]   $PREP_MEM / ${PREP_HOURS}h each"
    dep_prep=(--dependency="afterok:$jid" --kill-on-invalid-dep=yes)
  else echo "[SKIP] prep: done for all chromosomes"; fi
  # models
  model_ids=""
  for m in snp indel; do
    if model_done "${m}s"; then echo "[SKIP] model $m: done"; continue; fi
    if [[ "$m" == snp ]]; then mem="$SNP_MEM"; jm="$SNP_JAVA"; h="$SNP_HOURS"; else mem="$INDEL_MEM"; jm="$INDEL_JAVA"; h="$INDEL_HOURS"; fi
    if [[ -s "$OUT/model/${TAG}.${m}s.model_chroms.txt" ]]; then
      echo "[X]  $OUT/model has a ${m}s model for other chromosomes ($(cat "$OUT/model/${TAG}.${m}s.model_chroms.txt")) — delete $OUT/model/${TAG}.${m}s.* to retrain."; exit 1
    fi
    log="$LOG_DIR/${TAG}-S06-${OUT_NAME}-model-$m-%j.out"
    jid=$(submit "${common[@]}" --job-name="${JOB_BASE}-model-$m" "${dep_prep[@]+"${dep_prep[@]}"}" \
            --mem="$mem" --time="${h}:00:00" --output="$log" --wrap "$(stage_wrap "$jm" model "$m")")
    printf '[*]  %-12s: job %s   %s / %sh   (after prep)\n' "model $m" "$jid" "$mem" "$h"
    model_ids="${model_ids:+$model_ids:}$jid"
  done
  # apply
  todo=""; for c in $CHROMS; do vcf_done "$(final_vcf "$c")" || todo="${todo:+$todo,}$(idx_of "$c")"; done
  if [[ -n "$todo" ]]; then
    deps="${model_ids}"; [[ ${#dep_prep[@]} -gt 0 ]] && deps="${deps:+$deps:}${dep_prep[0]#--dependency=afterok:}"
    dep_apply=(); [[ -n "$deps" ]] && dep_apply=(--dependency="afterok:$deps" --kill-on-invalid-dep=yes)
    log="$LOG_DIR/${TAG}-S06-${OUT_NAME}-apply-%A_%a.out"
    jid=$(submit "${common[@]}" --job-name="${JOB_BASE}-apply" "${dep_apply[@]+"${dep_apply[@]}"}" --array="$(ranges "$todo")" \
            --mem="$APPLY_MEM" --time="${APPLY_HOURS}:00:00" --output="$log" --wrap "$(stage_wrap "$APPLY_JAVA" apply '$c')")
    echo "[*]  apply       : array $jid  [$(ranges "$todo")]   $APPLY_MEM / ${APPLY_HOURS}h each   (after both models)"
  else echo "[SKIP] apply: all chromosomes done"; fi
fi

echo
if [[ "$MODE" == dry ]]; then echo "[i]  --dry-run: nothing submitted."
else
  echo "[i]  Status : bash $0 --status --chroms $CHROMS_ARG${METHOD:+ --method $METHOD}"
  echo "[i]  Report : bash $0 --report --chroms $CHROMS_ARG${METHOD:+ --method $METHOD}   (after it ends)"
  echo "[i]  A job fails -> later stages are cancelled; run the same command again (only missing parts rerun)."
fi

#EOF
