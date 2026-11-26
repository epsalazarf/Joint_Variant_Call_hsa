#!/usr/bin/env bash

# =============================================================================
# Title       : GATK VQSR variant filtering [FENIX]
# Description : Filters the Step 05 joint-genotyped VCFs (one per chromosome)
#               with GATK Variant Quality Score Recalibration, as in the Broad
#               joint-genotyping workflow:
#                 prep  <chrom>   hard-filter ExcessHet > 54.69 + sites-only copy
#                 model snp|indel VariantRecalibrator, ONE model over all the
#                                 chromosomes' sites-only VCFs (a per-chromosome
#                                 model would see too few training sites)
#                 apply <chrom>   ApplyVQSR: INDEL model, then SNP model
#               Or, for small test cohorts (too few sites for VQSR):
#                 hardfilter <chrom>  GATK's generic hard filters + ExcessHet
#               Every site is kept: the FILTER column says PASS or why not
#               (soft filtering). Downstream: bcftools view -f PASS.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-11-26
# Version     : 1.0
# Usage       : 06_gatk_vqsr.sh <stage> <vcf_dir> <output_path> <chrom | snp | indel>
#             : stage       — prep | model | apply | hardfilter
#             : vcf_dir     — Step 05 output: <vcf_dir>/<prefix>.<chrom>.vcf.gz (+ .tbi)
#             : output_path — group storage, never $HOME. Layout:
#             :   work/<prefix>.<chrom>.exhet.vcf.gz       prep: ExcessHet-filtered, genotypes
#             :   work/<prefix>.<chrom>.sitesonly.vcf.gz   prep: same, no genotypes (model input)
#             :   model/<prefix>.{snps,indels}.recal.vcf.gz + .tranches + .model_chroms.txt
#             :   <prefix>.<chrom>.vqsr.vcf.gz             apply: final (+ .tbi)
#             :   <prefix>.<chrom>.hardfilt.vcf.gz         hardfilter: final (+ .tbi)
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
#               Broad JointGenotyping WDL (resources, priors, annotations, tranches)
# =============================================================================
#
# Options (environment variables):
#   VQSR_PREFIX       file prefix, input and output        (default: cohort)
#   VQSR_CHROMS       chromosomes the model trains on      (default: chr1..chr22)
#                     `apply` refuses a chromosome that was not in the model.
#   VQSR_SNP_LEVEL    SNP truth-sensitivity filter level   (default 99.7)
#   VQSR_INDEL_LEVEL  INDEL truth-sensitivity level        (default 99.0)
#   VQSR_SNP_GAUSSIANS / VQSR_INDEL_GAUSSIANS
#                     --max-gaussians tried in order until the model converges
#                     (default "6 4 2" / "4 2 1"; few sites -> fewer gaussians)
#   VQSR_JAVA_MEM     Java -Xmx                            (default 8G remote / 4G local)
#   VQSR_PLOTS        true: VariantRecalibrator also writes R plots (needs R +
#                     ggplot2 in PATH — the model fails if Rscript fails). Default false.
#   VQSR_SCRATCH_ROOT overrides config `scratch_root` for GATK temp files
#   VQSR_BUNDLE_DIR   overrides config `ref_vqsr_bundle` (training VCFs, made by
#                     bin/supp/STAGE.vqsr_bundle_hg38.sh)
#   VQSR_HOUSEKEEP    true (default): `apply` deletes work/<chrom>.exhet.* once
#                     the final VCF is verified
#
# Skip logic: each stage keeps an existing, complete output. Outputs are written
# as *.tmp.* and renamed only after the checks pass (BGZF end-of-file marker,
# sample count, record count = Step 05 record count).
# =============================================================================

set -euo pipefail

# <ARGUMENTS> ---------------------------------------------------------------

USAGE="Usage: $(basename "$0") <prep|model|apply|hardfilter> <vcf_dir> <output_path> <chrom|snp|indel>"
STAGE="${1:?$USAGE}"
VCF_DIR="${2:?$USAGE}"
OUTPUT_PATH="${3:?$USAGE}"
TARGET="${4:?$USAGE}"

# <\ARGS> -------------------------------------------------------------------

# <ENVIRONMENT> -------------------------------------------------------------

echo
echo "[$] GATK VQSR — ${STAGE} ${TARGET} [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

PREFIX="${VQSR_PREFIX:-cohort}"
MODEL_CHROMS="${VQSR_CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22}"
MODEL_CHROMS="$(echo "$MODEL_CHROMS" | tr ',' ' ' | xargs)"
SNP_LEVEL="${VQSR_SNP_LEVEL:-99.7}"
INDEL_LEVEL="${VQSR_INDEL_LEVEL:-99.0}"
SNP_GAUSSIANS="${VQSR_SNP_GAUSSIANS:-6 4 2}"
INDEL_GAUSSIANS="${VQSR_INDEL_GAUSSIANS:-4 2 1}"
PLOTS="${VQSR_PLOTS:-false}"
HOUSEKEEP="${VQSR_HOUSEKEEP:-true}"
ALL_CHROMS=" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM "

case "$STAGE" in
  prep|apply|hardfilter)
    case "$ALL_CHROMS" in *" $TARGET "*) ;; *) echo "[X]  CANCELLED. Invalid chromosome: '$TARGET'"; exit 1 ;; esac ;;
  model)
    case "$TARGET" in snp|indel) ;; *) echo "[X]  CANCELLED. model needs snp or indel (got '$TARGET')"; exit 1 ;; esac ;;
  *) echo "[X]  CANCELLED. Unknown stage '$STAGE'. $USAGE"; exit 1 ;;
esac
for c in $MODEL_CHROMS; do
  case "$ALL_CHROMS" in *" $c "*) ;; *) echo "[X]  CANCELLED. Invalid chromosome in VQSR_CHROMS: '$c'"; exit 1 ;; esac
done

[[ -d "$VCF_DIR" ]] || { echo "[X]  CANCELLED. vcf_dir not found: $VCF_DIR"; exit 1; }
VCF_DIR="$(readlink -f "$VCF_DIR")"

home_guard() {
  case "$1/" in
    "${HOME}"/*) echo "[X]  CANCELLED. Output is under \$HOME (${HOME}) — home has a hard quota on FENIX."; exit 1 ;;
  esac
}
home_guard "$OUTPUT_PATH"
mkdir -p "$OUTPUT_PATH/work" "$OUTPUT_PATH/model"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"
home_guard "$OUTPUT_PATH"
WORK="$OUTPUT_PATH/work"; MODEL="$OUTPUT_PATH/model"

CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"; MEM="${VQSR_JAVA_MEM:-8G}"
else
  env_type="local";  MEM="${VQSR_JAVA_MEM:-4G}"
fi
echo "[i]  Environment: $env_type   java-mem=${MEM}   prefix=${PREFIX}"

eval "$(
  awk -v env="$env_type" '
    BEGIN { in_env=0 }
    $1 ~ env":" { in_env=1; next }
    in_env && /^[^[:space:]]/ { in_env=0 }
    in_env && /^[[:space:]]+[a-zA-Z0-9_]+:/ {
      gsub(":", "=", $1)
      sub(/^[[:space:]]+/, "", $1)
      gsub(/^"/, "", $2); gsub(/"$/, "", $2)
      print $1 $2
    }
  ' "$CONFIG_FILE"
)"

if [[ "$env_type" == "remote" ]]; then
  echo "[i]  Loading modules..."
  module unload oracle-java 2>/dev/null || true
  module load oracle-java/25.0.2
  module load gatk
  module load bcftools
fi

[[ -n "${ref_gnm:-}" && -f "${ref_gnm}" ]] || { echo "[X]  Reference genome not found: ${ref_gnm:-<unset>}"; exit 1; }
echo "[i]  Reference  : ${ref_gnm}"

scratch_root="${VQSR_SCRATCH_ROOT:-${scratch_root:-}}"
if [[ -n "$scratch_root" && "$env_type" == "remote" ]]; then
  SCRATCH_BASE="${scratch_root}/${USER}/vqsr_${SLURM_JOB_ID:-$$}_$(date +%s)"
else
  SCRATCH_BASE="${OUTPUT_PATH}/.vqsr_work_${SLURM_JOB_ID:-$$}_$(date +%s)"
fi
mkdir -p "$SCRATCH_BASE/tmp" || { echo "[X]  Cannot create temp dir: $SCRATCH_BASE/tmp"; exit 1; }
export TMPDIR="$SCRATCH_BASE/tmp"
cleanup() { rm -rf "$SCRATCH_BASE"; }
trap cleanup EXIT

# <\ENV> --------------------------------------------------------------------

# <FUNCTIONS> ---------------------------------------------------------------

BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
bgzf_complete() { [[ "$(tail -c 28 "$1" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')" == "$BGZF_EOF" ]]; }
done_vcf() { [[ -s "$1" && -s "$1.tbi" ]]; }
hms() { echo $(( $(date +%s) - $1 )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp'; }
n_samples() { bcftools query -l "$1" | wc -l | tr -d '[:space:]'; }
# Record count from the index (bcftools-made indexes carry counts; GATK's do
# not -> empty).
n_indexed() { bcftools index -n "$1" 2>/dev/null || true; }
gatk_run() { gatk --java-options "-Xms2G -Xmx${MEM}" "$@" --tmp-dir "$TMPDIR"; }

input_vcf() {   # $1 chrom -> Step 05 VCF, checked
  local f="$VCF_DIR/${PREFIX}.$1.vcf.gz"
  done_vcf "$f" || { echo "<ERROR> Step 05 VCF (+ .tbi) not found: $f"; exit 1; }
  bgzf_complete "$f" || { echo "<ERROR> Step 05 VCF looks truncated: $f"; exit 1; }
  echo "$f"
}

# Count records and FILTER values in one full read (also proves every block
# decompresses). Prints "total<TAB>table".
filter_table() {
  bcftools query -f '%FILTER\n' "$1" | sort | uniq -c | sort -k1,1nr \
    | awk '{t += $1; s = s sprintf("       %12d  %s\n", $1, $2)} END {printf "%d\t%s", t, s}'
}

finish_vcf() {   # $1 tmp $2 final $3 expected records $4 expected samples
  local tmp="$1" out="$2" want="$3" ns="$4" tab total
  [[ -s "$tmp" && -s "$tmp.tbi" ]] || { echo "<ERROR> no VCF + index written: $tmp"; exit 1; }
  bgzf_complete "$tmp" || { echo "<ERROR> VCF looks truncated (no BGZF EOF): $tmp"; exit 1; }
  [[ "$(n_samples "$tmp")" == "$ns" ]] || { echo "<ERROR> $tmp has $(n_samples "$tmp") samples, expected $ns"; exit 1; }
  tab="$(filter_table "$tmp")"; total="${tab%%$'\t'*}"
  if [[ -n "$want" && "$total" != "$want" ]]; then
    echo "<ERROR> $tmp has $total records, input has $want"; exit 1
  fi
  mv -f "$tmp.tbi" "$out.tbi"; mv -f "$tmp" "$out"
  echo "[DONE] $out"
  echo "[i]  Samples: $ns   Records: $total   Size: $(du -h "$out" | cut -f1)"
  echo "[i]  FILTER:"; printf '%s\n' "${tab#*$'\t'}"
}

bundle_file() {   # $1 stem -> path in the VQSR bundle, checked
  local d="${VQSR_BUNDLE_DIR:-${ref_vqsr_bundle:-}}" f
  [[ -n "$d" ]] || { echo "<ERROR> No VQSR bundle: set config ref_vqsr_bundle (or VQSR_BUNDLE_DIR)"; exit 1; }
  f="$d/$1.canon_chr.vcf.gz"
  done_vcf "$f" || { echo "<ERROR> VQSR resource missing (+ .tbi): $f — run bin/supp/STAGE.vqsr_bundle_hg38.sh"; exit 1; }
  echo "$f"
}

## prep: ExcessHet hard filter (site fails HWE for excess heterozygosity —
## mostly paralog/mapping artefacts) and a sites-only copy for the model.
stage_prep() {
  local c="$TARGET" in exhet sites t; t=$(date +%s)
  in="$(input_vcf "$c")"
  exhet="$WORK/${PREFIX}.$c.exhet.vcf.gz"; sites="$WORK/${PREFIX}.$c.sitesonly.vcf.gz"
  echo "[<]  Input      : $in"
  if done_vcf "$sites" && { done_vcf "$exhet" || done_vcf "$OUTPUT_PATH/${PREFIX}.$c.vqsr.vcf.gz"; }; then
    echo "[SKIP] prep $c done: $sites"; return 0   # exhet is deleted once apply has used it
  fi
  local want ns; want="$(n_indexed "$in")"; ns="$(n_samples "$in")"

  echo; echo "<START> VariantFiltration ExcessHet ${c}  ($(date))"
  rm -f "$WORK/${PREFIX}.$c.exhet.tmp.vcf.gz"*
  set -o xtrace
  gatk_run VariantFiltration \
    --reference "$ref_gnm" --variant "$in" \
    --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet \
    --create-output-variant-index true \
    --output "$WORK/${PREFIX}.$c.exhet.tmp.vcf.gz"
  { set +o xtrace; } 2>/dev/null
  local etmp="$WORK/${PREFIX}.$c.exhet.tmp.vcf.gz"
  bgzf_complete "$etmp" && [[ -s "$etmp.tbi" ]] || { echo "<ERROR> VariantFiltration output incomplete: $etmp"; exit 1; }
  [[ "$(n_samples "$etmp")" == "$ns" ]] || { echo "<ERROR> $etmp lost samples"; exit 1; }

  echo "<START> sites-only ${c}"
  local stmp="$WORK/${PREFIX}.$c.sitesonly.tmp.vcf.gz"
  bcftools view -G --threads 2 -Oz -o "$stmp" "$etmp"
  bcftools index -t -f "$stmp"           # bcftools index = carries the record count
  bgzf_complete "$stmp" || { echo "<ERROR> sites-only VCF truncated: $stmp"; exit 1; }
  local n; n="$(n_indexed "$stmp")"
  if [[ -n "$want" && "$n" != "$want" ]]; then echo "<ERROR> sites-only has $n records, Step 05 VCF has $want"; exit 1; fi
  [[ -n "$want" ]] || echo "[!]  Step 05 index has no record count — sites-only count ($n) not cross-checked"
  mv -f "$etmp.tbi" "$exhet.tbi"; mv -f "$etmp" "$exhet"
  mv -f "$stmp.tbi" "$sites.tbi"; mv -f "$stmp" "$sites"
  local nf; nf=$(bcftools view -H -f ExcessHet "$sites" | wc -l | tr -d '[:space:]')
  echo "[DONE] $exhet"
  echo "[DONE] $sites"
  echo "[i]  Records: $n   ExcessHet-filtered: $nf   [$(hms "$t")]"
}

## model: one VariantRecalibrator over every chromosome in VQSR_CHROMS.
stage_model() {
  local mode tag t; t=$(date +%s)
  if [[ "$TARGET" == snp ]]; then mode=SNP; tag=snps; else mode=INDEL; tag=indels; fi
  local recal="$MODEL/${PREFIX}.$tag.recal.vcf.gz" tranches="$MODEL/${PREFIX}.$tag.tranches"
  local chroms_file="$MODEL/${PREFIX}.$tag.model_chroms.txt"
  if done_vcf "$recal" && [[ -s "$tranches" && -s "$chroms_file" ]]; then
    if [[ "$(cat "$chroms_file")" == "$MODEL_CHROMS" ]]; then echo "[SKIP] $mode model done: $recal"; return 0; fi
    echo "<ERROR> $mode model exists but was trained on: $(cat "$chroms_file")"
    echo "        requested VQSR_CHROMS: $MODEL_CHROMS — delete $MODEL/${PREFIX}.$tag.* to retrain"; exit 1
  fi

  local inputs=() c f n_sites=0 n
  for c in $MODEL_CHROMS; do
    f="$WORK/${PREFIX}.$c.sitesonly.vcf.gz"
    done_vcf "$f" || { echo "<ERROR> prep not done for $c: $f"; exit 1; }
    inputs+=(--variant "$f"); n="$(n_indexed "$f")"; n_sites=$((n_sites + ${n:-0}))
  done
  echo "[<]  Model input: $(echo "$MODEL_CHROMS" | wc -w | tr -d ' ') chromosome(s), $n_sites sites"

  local dbsnp="${ref_vars:-}"
  [[ -n "$dbsnp" ]] && done_vcf "$dbsnp" || { echo "<ERROR> dbSNP (config ref_vars) not found: ${dbsnp:-<unset>}"; exit 1; }
  local res=() an=() tr=() gaussians
  if [[ "$mode" == SNP ]]; then
    res=(--resource:hapmap,known=false,training=true,truth=true,prior=15 "$(bundle_file hapmap_3.3.hg38)"
         --resource:omni,known=false,training=true,truth=true,prior=12 "$(bundle_file 1000G_omni2.5.hg38)"
         --resource:1000G,known=false,training=true,truth=false,prior=10 "$(bundle_file 1000G_phase1.snps.high_confidence.hg38)"
         --resource:dbsnp,known=true,training=false,truth=false,prior=7 "$dbsnp")
    an=(-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP)
    for x in 100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0; do tr+=(-tranche "$x"); done
    gaussians="$SNP_GAUSSIANS"
  else
    res=(--resource:mills,known=false,training=true,truth=true,prior=12 "$(bundle_file Mills_and_1000G_gold_standard.indels.hg38)"
         --resource:axiomPoly,known=false,training=true,truth=false,prior=10 "$(bundle_file Axiom_Exome_Plus.genotypes.all_populations.poly.hg38)"
         --resource:dbsnp,known=true,training=false,truth=false,prior=2 "$dbsnp")
    an=(-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP)
    for x in 100.0 99.95 99.9 99.5 99.0 97.0 96.0 95.0 94.0 93.5 93.0 92.0 91.0 90.0; do tr+=(-tranche "$x"); done
    gaussians="$INDEL_GAUSSIANS"
  fi
  local plot=()
  [[ "$PLOTS" == true ]] && plot=(--rscript-file "$MODEL/${PREFIX}.$tag.plots.R")

  local rtmp="$MODEL/${PREFIX}.$tag.recal.tmp.vcf.gz" ttmp="$MODEL/${PREFIX}.$tag.tmp.tranches" mg ok=false
  for mg in $gaussians; do
    echo; echo "<START> VariantRecalibrator $mode  --max-gaussians $mg  ($(date))"
    rm -f "$rtmp" "$rtmp.tbi" "$ttmp"
    set -o xtrace
    if gatk_run VariantRecalibrator \
         --reference "$ref_gnm" "${inputs[@]}" \
         --trust-all-polymorphic -mode "$mode" --max-gaussians "$mg" \
         "${an[@]}" "${tr[@]}" "${res[@]}" "${plot[@]+"${plot[@]}"}" \
         --create-output-variant-index true \
         --output "$rtmp" --tranches-file "$ttmp"; then
      { set +o xtrace; } 2>/dev/null; ok=true; break
    fi
    { set +o xtrace; } 2>/dev/null
    echo "[!]  VariantRecalibrator $mode failed with --max-gaussians $mg (too few sites per cluster?) — trying fewer"
  done
  [[ "$ok" == true ]] || { echo "<ERROR> VariantRecalibrator $mode failed for every --max-gaussians in: $gaussians"; exit 1; }
  bgzf_complete "$rtmp" && [[ -s "$rtmp.tbi" && -s "$ttmp" ]] || { echo "<ERROR> recal/tranches incomplete: $rtmp"; exit 1; }
  mv -f "$rtmp.tbi" "$recal.tbi"; mv -f "$rtmp" "$recal"; mv -f "$ttmp" "$tranches"
  echo "$MODEL_CHROMS" > "$chroms_file"         # written last = model complete
  echo "[DONE] $recal"
  echo "[DONE] $tranches   (max-gaussians $mg)"
  echo "[i]  Tranches (truth sensitivity: known / novel sites passing, Ti/Tv known / novel, min VQSLOD):"
  awk -F, '/^#/ {next} /^targetTruthSensitivity/ {for (i=1;i<=NF;i++) h[$i]=i; next}
           ("targetTruthSensitivity" in h) {printf "       %7s  known=%-9s novel=%-9s TiTv=%s/%s  minVQSLod=%s\n",
             $h["targetTruthSensitivity"], $h["numKnown"], $h["numNovel"], $h["knownTiTv"], $h["novelTiTv"], $h["minVQSLod"]}' "$tranches" || true
  echo "[i]  [$(hms "$t")]"
}

## apply: INDEL model first, then SNP model on its output (Broad order).
stage_apply() {
  local c="$TARGET" t; t=$(date +%s)
  local exhet="$WORK/${PREFIX}.$c.exhet.vcf.gz" sites="$WORK/${PREFIX}.$c.sitesonly.vcf.gz"
  local out="$OUTPUT_PATH/${PREFIX}.$c.vqsr.vcf.gz" tmp="$OUTPUT_PATH/${PREFIX}.$c.vqsr.tmp.vcf.gz"
  if done_vcf "$out"; then echo "[SKIP] $out exists"; return 0; fi
  done_vcf "$exhet" && done_vcf "$sites" || { echo "<ERROR> prep not done for $c: $exhet"; exit 1; }
  local tag
  for tag in indels snps; do
    done_vcf "$MODEL/${PREFIX}.$tag.recal.vcf.gz" && [[ -s "$MODEL/${PREFIX}.$tag.tranches" && -s "$MODEL/${PREFIX}.$tag.model_chroms.txt" ]] \
      || { echo "<ERROR> $tag model not built yet: $MODEL/${PREFIX}.$tag.*"; exit 1; }
    case " $(cat "$MODEL/${PREFIX}.$tag.model_chroms.txt") " in *" $c "*) ;;
      *) echo "<ERROR> $c was not in the $tag model ($(cat "$MODEL/${PREFIX}.$tag.model_chroms.txt"))"; exit 1 ;; esac
  done
  local want ns; want="$(n_indexed "$sites")"; ns="$(n_samples "$exhet")"
  local mid="$WORK/${PREFIX}.$c.indel_recal.tmp.vcf.gz"
  echo "[<]  Input      : $exhet  ($want records, $ns samples)"
  echo "[i]  Levels     : INDEL ${INDEL_LEVEL}   SNP ${SNP_LEVEL}"

  echo; echo "<START> ApplyVQSR INDEL ${c}  ($(date))"
  rm -f "$mid" "$mid.tbi" "$tmp" "$tmp.tbi"
  set -o xtrace
  gatk_run ApplyVQSR --reference "$ref_gnm" --variant "$exhet" --intervals "$c" \
    --recal-file "$MODEL/${PREFIX}.indels.recal.vcf.gz" --tranches-file "$MODEL/${PREFIX}.indels.tranches" \
    --truth-sensitivity-filter-level "$INDEL_LEVEL" -mode INDEL \
    --create-output-variant-index true --output "$mid"
  { set +o xtrace; } 2>/dev/null
  echo; echo "<START> ApplyVQSR SNP ${c}  ($(date))"
  set -o xtrace
  gatk_run ApplyVQSR --reference "$ref_gnm" --variant "$mid" --intervals "$c" \
    --recal-file "$MODEL/${PREFIX}.snps.recal.vcf.gz" --tranches-file "$MODEL/${PREFIX}.snps.tranches" \
    --truth-sensitivity-filter-level "$SNP_LEVEL" -mode SNP \
    --create-output-variant-index true --output "$tmp"
  { set +o xtrace; } 2>/dev/null
  echo; echo "[i]  Verifying output..."
  finish_vcf "$tmp" "$out" "$want" "$ns"
  rm -f "$mid" "$mid.tbi"
  if [[ "$HOUSEKEEP" == true ]]; then rm -f "$exhet" "$exhet.tbi"; echo "[i]  removed $exhet (VQSR_HOUSEKEEP=false keeps it)"; fi
  echo "[i]  [$(hms "$t")]"
}

## hardfilter: GATK's generic thresholds, for cohorts too small for VQSR.
## Mixed SNP+indel sites count as indels. A missing annotation never fails a
## site (VariantFiltration default).
stage_hardfilter() {
  local c="$TARGET" t in; t=$(date +%s)
  in="$(input_vcf "$c")"
  local out="$OUTPUT_PATH/${PREFIX}.$c.hardfilt.vcf.gz" tmp="$OUTPUT_PATH/${PREFIX}.$c.hardfilt.tmp.vcf.gz"
  if done_vcf "$out"; then echo "[SKIP] $out exists"; return 0; fi
  local want ns; want="$(n_indexed "$in")"; ns="$(n_samples "$in")"
  echo "[<]  Input      : $in  (${want:-?} records, $ns samples)"
  echo; echo "<START> VariantFiltration hard filters ${c}  ($(date))"
  rm -f "$tmp" "$tmp.tbi"
  set -o xtrace
  gatk_run VariantFiltration --reference "$ref_gnm" --variant "$in" \
    --filter-expression "vc.isSNP() && QD < 2.0"             --filter-name SNP_QD2 \
    --filter-expression "vc.isSNP() && QUAL < 30.0"          --filter-name SNP_QUAL30 \
    --filter-expression "vc.isSNP() && SOR > 3.0"            --filter-name SNP_SOR3 \
    --filter-expression "vc.isSNP() && FS > 60.0"            --filter-name SNP_FS60 \
    --filter-expression "vc.isSNP() && MQ < 40.0"            --filter-name SNP_MQ40 \
    --filter-expression "vc.isSNP() && MQRankSum < -12.5"    --filter-name SNP_MQRankSum-12.5 \
    --filter-expression "vc.isSNP() && ReadPosRankSum < -8.0" --filter-name SNP_ReadPosRankSum-8 \
    --filter-expression "!vc.isSNP() && QD < 2.0"            --filter-name INDEL_QD2 \
    --filter-expression "!vc.isSNP() && QUAL < 30.0"         --filter-name INDEL_QUAL30 \
    --filter-expression "!vc.isSNP() && FS > 200.0"          --filter-name INDEL_FS200 \
    --filter-expression "!vc.isSNP() && SOR > 10.0"          --filter-name INDEL_SOR10 \
    --filter-expression "!vc.isSNP() && ReadPosRankSum < -20.0" --filter-name INDEL_ReadPosRankSum-20 \
    --filter-expression "ExcessHet > 54.69"                  --filter-name ExcessHet \
    --create-output-variant-index true --output "$tmp"
  { set +o xtrace; } 2>/dev/null
  echo; echo "[i]  Verifying output..."
  finish_vcf "$tmp" "$out" "$want" "$ns"
  echo "[i]  [$(hms "$t")]"
}

# <\FUNCTIONS> --------------------------------------------------------------

# <MAIN> --------------------------------------------------------------------

case "$STAGE" in
  prep)       stage_prep ;;
  model)      stage_model ;;
  apply)      stage_apply ;;
  hardfilter) stage_hardfilter ;;
esac

echo
echo "[$] GATK VQSR — ${STAGE} ${TARGET} completed successfully!"
echo "[&]  Total time: $(hms "$script_timestamp")"

# <\MAIN> -------------------------------------------------------------------

#EOF
