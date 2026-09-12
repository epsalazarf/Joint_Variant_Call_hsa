#!/usr/bin/env bash

# =============================================================================
# Title       : GATK VQSR / Hard-Filtering [FENIX]
# Description : Filters the joint-genotyped cohort VCF from Step 05. Two paths,
#               selectable or auto-chosen from cohort size:
#                 vqsr        — VariantRecalibrator + ApplyVQSR, SNP and INDEL
#                               modes run separately then chained (GATK4 Best
#                               Practices). Needs the GATK hg38 resource bundle.
#                 hard-filter — GATK's standard hard-filter expressions, SNP
#                               and INDEL run separately then merged. No extra
#                               resource files; recommended for small cohorts
#                               where VQSR's Gaussian mixture model won't
#                               converge (GATK's own docs cite ~30 WGS samples
#                               as a rough lower bound).
#               Neither path removes filtered records — FILTER is annotated
#               (PASS or a semicolon-joined failure reason) and it is left to
#               downstream analysis to select on it.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-09-12
# Version     : 1.0
# Usage       : 06_gatk_vqsr.sh <joint_vcf> <output_path> [mode] [cohort_name]
#             : joint_vcf   — the gathered <cohort>.joint.vcf.gz from Step 05
#             : output_path — filtered VCF + intermediates written here
#             : mode        — vqsr | hard-filter | auto (default: auto — see
#             :               resolve_mode() / VQSR_MIN_SAMPLES below)
#             : cohort_name — label for output filenames (default: derived
#             :               from joint_vcf's own name, else "cohort")
# Source      : GATK4 Best Practices — https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
#             : Hard-filtering germline short variants — https://gatk.broadinstitute.org/hc/en-us/articles/360035890471
# =============================================================================
#
# ---------------------------------------------------------------------------
#  STATUS
# ---------------------------------------------------------------------------
#
#   DRAFT — written against GATK4's published Best Practices recipes, not yet
#   run on FENIX or against a real joint VCF. Before production use:
#     - the GATK hg38 VQSR resource bundle (hapmap, omni, 1000G SNPs, Mills
#       indels) is NOT staged on FENIX yet — config/config.yaml carries
#       EDIT_THIS placeholders for ref_hapmap/ref_omni/ref_1kg_snp/ref_mills;
#       vqsr mode refuses to run until they point at real files
#     - VQSR_MIN_SAMPLES=30 (the auto-mode SNP/INDEL threshold) is a commonly
#       cited GATK rule of thumb, not a measurement against this project's
#       actual cohorts — docs/PIPELINE_STATUS.md still lists the VQSR-vs-
#       hard-filter choice as an open decision; override --mode explicitly
#       once that decision is made instead of relying on auto
#     - annotations used are the non-allele-specific set (matches Step 05,
#       which also does not request AS annotations — see that script's header)
#     - Java heap defaults below are placeholders, not measurements
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

JOINT_VCF="${1:?Usage: $(basename "$0") <joint_vcf> <output_path> [vqsr|hard-filter|auto] [cohort_name]}"
OUTPUT_PATH="${2:?Usage: $(basename "$0") <joint_vcf> <output_path> [vqsr|hard-filter|auto] [cohort_name]}"
MODE_ARG="${3:-auto}"

case "$MODE_ARG" in
  vqsr|hard-filter|auto) ;;
  *) echo "[X]  Invalid mode: '$MODE_ARG' (expected: vqsr | hard-filter | auto)"; exit 1 ;;
esac

## Derive a cohort label from the input filename (strips .joint.vcf.gz / .vcf.gz).
default_cohort() {
  local b; b="$(basename "$JOINT_VCF")"
  b="${b%.vcf.gz}"; b="${b%.joint}"
  echo "${b:-cohort}"
}
COHORT="${4:-$(default_cohort)}"

# <\ARGS> -----------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] GATK VQSR / Hard-Filtering [FENIX] >>"
echo "[&]  Started: $(date)"
script_timestamp=$(date +%s)

echo
echo "[i]  Checking inputs..."

[ -s "$JOINT_VCF" ] || { echo "[X]  CANCELLED. Joint VCF not found or empty: $JOINT_VCF"; exit 1; }
JOINT_VCF="$(readlink -f "$JOINT_VCF")"
[[ -s "${JOINT_VCF}.tbi" || -s "${JOINT_VCF}.csi" ]] || { echo "[X]  Missing index (.tbi/.csi) for: $JOINT_VCF"; exit 1; }
mkdir -p "$OUTPUT_PATH"
OUTPUT_PATH="$(readlink -f "$OUTPUT_PATH")"

echo "[<]  Joint VCF  : $JOINT_VCF"
echo "[i]  Output     : ${OUTPUT_PATH}"
echo "[i]  Cohort     : ${COHORT}"
echo "[i]  Mode arg   : ${MODE_ARG}"

# Options (env-overridable so a launcher can vary resources per run)
#   VQSR_JAVA_MEM        — Java -Xms/-Xmx                      (default 6G remote / 4G local; PLACEHOLDER, see STATUS above)
#   VQSR_VERBOSITY       — GATK --verbosity                     (default ERROR)
#   VQSR_THREADS         — bcftools threads (sample count check) (default: $SLURM_CPUS_PER_TASK, else 2)
#   VQSR_MIN_SAMPLES     — auto-mode threshold: >= this many samples -> vqsr, else hard-filter (default 30)
#   VQSR_TS_SNP          — ApplyVQSR --truth-sensitivity-filter-level, SNP mode   (default 99.7)
#   VQSR_TS_INDEL        — ApplyVQSR --truth-sensitivity-filter-level, INDEL mode (default 99.0)
#   VQSR_SNP_MAX_GAUSSIANS   — VariantRecalibrator --max-gaussians, SNP mode   (default: GATK's own default, unset)
#   VQSR_INDEL_MAX_GAUSSIANS — VariantRecalibrator --max-gaussians, INDEL mode (default 4 — GATK's own guidance for
#                              the sparser indel training set; small cohorts commonly fail to converge above this)
VERBOSITY="${VQSR_VERBOSITY:-ERROR}"
njobs="${VQSR_THREADS:-${SLURM_CPUS_PER_TASK:-2}}"
MIN_SAMPLES="${VQSR_MIN_SAMPLES:-30}"
TS_SNP="${VQSR_TS_SNP:-99.7}"
TS_INDEL="${VQSR_TS_INDEL:-99.0}"
INDEL_MAX_GAUSSIANS="${VQSR_INDEL_MAX_GAUSSIANS:-4}"
SNP_MAX_GAUSSIANS="${VQSR_SNP_MAX_GAUSSIANS:-}"

# Config file (relative to repo root)
CONFIG_FILE="$(dirname "$(readlink -f "$0")")/../config/config.yaml"

# Detect environment
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then
  env_type="remote"
  MEM="${VQSR_JAVA_MEM:-6G}"
else
  env_type="local"
  MEM="${VQSR_JAVA_MEM:-4G}"
fi
echo "[i]  Environment: $env_type"
echo "[i]  Knobs      : java-mem=${MEM}  verbosity=${VERBOSITY}  min-samples(auto)=${MIN_SAMPLES}  ts-snp=${TS_SNP}  ts-indel=${TS_INDEL}"

# Parse YAML config into Bash variables (embedded parser — no external tool)
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

# Load modules on remote (work-around due to faulty parser [ARC02])
if [[ "$env_type" == "remote" ]]; then
  echo "[i]  Loading modules..."
  module unload oracle-java 2>/dev/null || true
  module load oracle-java/25.0.2
  module load gatk
  module load bcftools
fi

# Guard: reference genome (needed by every mode)
if [ -z "${ref_gnm:-}" ]; then
  echo "[X]  Missing ref_gnm in config: $CONFIG_FILE"
  exit 1
fi
[ -f "${ref_gnm}" ] || { echo "[X]  Reference genome not found: ${ref_gnm}"; exit 1; }
echo "[i]  Reference  : ${ref_gnm}"

TMPDIR_LOCAL="${OUTPUT_PATH}/.vqsr_tmp_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR_LOCAL"
export TMPDIR="$TMPDIR_LOCAL"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

# <\ENV> --------------------------------------------------------------------

# <FUNCTIONS> ---------------------------------------------------------------

# Standard GATK4 tranche ladder for VariantRecalibrator (both modes).
TRANCHES=(100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0)

## Decide vqsr vs hard-filter when mode=auto; otherwise just report the choice.
resolve_mode() {
  if [[ "$MODE_ARG" != "auto" ]]; then
    MODE="$MODE_ARG"
    echo "[i]  Mode        : ${MODE} (explicit)"
    return 0
  fi
  local n_samples
  n_samples="$(bcftools query -l "$JOINT_VCF" | wc -l)"
  if (( n_samples >= MIN_SAMPLES )); then
    MODE="vqsr"
  else
    MODE="hard-filter"
  fi
  echo "[i]  Mode        : ${MODE} (auto: ${n_samples} samples vs VQSR_MIN_SAMPLES=${MIN_SAMPLES})"
}

## Guard: the four required VQSR resource files (ref_axiom is optional).
guard_vqsr_resources() {
  if [ -z "${ref_vars:-}" ] || [ -z "${ref_hapmap:-}" ] || [ -z "${ref_omni:-}" ] \
     || [ -z "${ref_1kg_snp:-}" ] || [ -z "${ref_mills:-}" ]; then
    echo "[X]  CANCELLED. VQSR requires ref_vars/ref_hapmap/ref_omni/ref_1kg_snp/ref_mills in config."
    echo "[i]  These are EDIT_THIS placeholders until the GATK hg38 resource bundle is staged — see config/config.yaml."
    exit 1
  fi
  local f
  for f in "$ref_vars" "$ref_hapmap" "$ref_omni" "$ref_1kg_snp" "$ref_mills"; do
    [ -f "$f" ] || { echo "[X]  CANCELLED. VQSR resource file not found: $f"; exit 1; }
  done
  if [[ -n "${ref_axiom:-}" ]]; then
    [ -f "$ref_axiom" ] || { echo "[X]  CANCELLED. ref_axiom set but not found: $ref_axiom"; exit 1; }
    echo "[i]  Indel resources: Mills + Axiom Exome Plus"
  else
    echo "[i]  Indel resources: Mills only (ref_axiom unset — skipped)"
  fi
  echo "[i]  SNP resources  : HapMap, Omni, 1000G, dbSNP"
}

_tranche_opts() { local t; for t in "${TRANCHES[@]}"; do printf -- '--truth-sensitivity-tranche %s ' "$t"; done; }

## VariantRecalibrator, SNP mode.
recal_snp() {
  local step_name="VariantRecalibrator (SNP)"
  local outdir="${OUTPUT_PATH}/vqsr_work"
  local recal="${outdir}/${COHORT}.snp.recal"
  local tranches="${outdir}/${COHORT}.snp.tranches"
  local step_timestamp=$EPOCHSECONDS

  echo; echo "[*]  $step_name"; echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ -s "$recal" && -s "$tranches" ]]; then
    echo "[i]  Already completed ($recal exists)"; return 0
  fi
  mkdir -p "$outdir"

  local gaussians_opt=()
  [[ -n "$SNP_MAX_GAUSSIANS" ]] && gaussians_opt=(--max-gaussians "$SNP_MAX_GAUSSIANS")

  set -o xtrace
  gatk --java-options "-Xms${MEM} -Xmx${MEM}" VariantRecalibrator \
    --reference "$ref_gnm" \
    --variant "$JOINT_VCF" \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$ref_hapmap" \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 "$ref_omni" \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 "$ref_1kg_snp" \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$ref_vars" \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    --mode SNP \
    $(_tranche_opts) \
    ${gaussians_opt[@]+"${gaussians_opt[@]}"} \
    --tmp-dir "$TMPDIR" \
    --verbosity "$VERBOSITY" \
    --tranches-file "$tranches" \
    --output "$recal"
  set +o xtrace

  [[ -s "$recal" && -s "$tranches" ]] || { echo "[X]  CANCELLED: $step_name produced no output"; exit 1; }
  echo "[>]  $recal"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## VariantRecalibrator, INDEL mode.
recal_indel() {
  local step_name="VariantRecalibrator (INDEL)"
  local outdir="${OUTPUT_PATH}/vqsr_work"
  local recal="${outdir}/${COHORT}.indel.recal"
  local tranches="${outdir}/${COHORT}.indel.tranches"
  local step_timestamp=$EPOCHSECONDS

  echo; echo "[*]  $step_name"; echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ -s "$recal" && -s "$tranches" ]]; then
    echo "[i]  Already completed ($recal exists)"; return 0
  fi
  mkdir -p "$outdir"

  local axiom_opt=()
  [[ -n "${ref_axiom:-}" ]] && axiom_opt=(--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 "$ref_axiom")

  set -o xtrace
  gatk --java-options "-Xms${MEM} -Xmx${MEM}" VariantRecalibrator \
    --reference "$ref_gnm" \
    --variant "$JOINT_VCF" \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 "$ref_mills" \
    ${axiom_opt[@]+"${axiom_opt[@]}"} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$ref_vars" \
    -an QD -an DP -an FS -an ReadPosRankSum -an MQRankSum -an SOR \
    --mode INDEL \
    --max-gaussians "$INDEL_MAX_GAUSSIANS" \
    $(_tranche_opts) \
    --tmp-dir "$TMPDIR" \
    --verbosity "$VERBOSITY" \
    --tranches-file "$tranches" \
    --output "$recal"
  set +o xtrace

  [[ -s "$recal" && -s "$tranches" ]] || { echo "[X]  CANCELLED: $step_name produced no output"; exit 1; }
  echo "[>]  $recal"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## ApplyVQSR, INDEL then SNP, chained onto the same output file.
apply_vqsr() {
  local step_name="ApplyVQSR (INDEL -> SNP)"
  local outdir="${OUTPUT_PATH}/vqsr_work"
  local indel_recal="${outdir}/${COHORT}.indel.recal" indel_tranches="${outdir}/${COHORT}.indel.tranches"
  local snp_recal="${outdir}/${COHORT}.snp.recal" snp_tranches="${outdir}/${COHORT}.snp.tranches"
  local indel_out="${outdir}/${COHORT}.indel_filtered.vcf.gz"
  local final="${OUTPUT_PATH}/${COHORT}.filtered.vcf.gz"
  local tmp="${OUTPUT_PATH}/.${COHORT}.filtered.tmp.vcf.gz"
  local step_timestamp=$EPOCHSECONDS

  echo; echo "[*]  $step_name"; echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ -s "$final" && -s "${final}.tbi" ]]; then
    echo "[i]  Already completed ($final exists)"; return 0
  fi

  if [[ ! -s "$indel_out" ]]; then
    set -o xtrace
    gatk --java-options "-Xms${MEM} -Xmx${MEM}" ApplyVQSR \
      --reference "$ref_gnm" \
      --variant "$JOINT_VCF" \
      --recal-file "$indel_recal" \
      --tranches-file "$indel_tranches" \
      --truth-sensitivity-filter-level "$TS_INDEL" \
      --mode INDEL \
      --create-output-variant-index \
      --tmp-dir "$TMPDIR" \
      --verbosity "$VERBOSITY" \
      --output "$indel_out"
    set +o xtrace
    [[ -s "$indel_out" ]] || { echo "[X]  CANCELLED: ApplyVQSR (INDEL) produced no output"; exit 1; }
    echo "[>]  $indel_out"
  else
    echo "[i]  ApplyVQSR (INDEL) already completed ($indel_out exists)"
  fi

  rm -f "$tmp" "${tmp}.tbi"
  set -o xtrace
  gatk --java-options "-Xms${MEM} -Xmx${MEM}" ApplyVQSR \
    --reference "$ref_gnm" \
    --variant "$indel_out" \
    --recal-file "$snp_recal" \
    --tranches-file "$snp_tranches" \
    --truth-sensitivity-filter-level "$TS_SNP" \
    --mode SNP \
    --create-output-variant-index \
    --tmp-dir "$TMPDIR" \
    --verbosity "$VERBOSITY" \
    --output "$tmp"
  set +o xtrace

  [[ -s "$tmp" ]] || { echo "[X]  CANCELLED: ApplyVQSR (SNP) produced no output"; exit 1; }
  mv "$tmp" "$final"
  mv "${tmp}.tbi" "${final}.tbi"

  echo "[>]  $final"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## SelectVariants + VariantFiltration for one variant type, hard-filter mode.
## $1 = SNP|INDEL   $2- = filter-expression/filter-name pairs
hardfilter_one_type() {
  local vtype="$1"; shift
  local step_name="Hard-filter (${vtype})"
  local outdir="${OUTPUT_PATH}/hardfilter_work"
  local selected="${outdir}/${COHORT}.${vtype,,}.vcf.gz"
  local filtered="${outdir}/${COHORT}.${vtype,,}_filtered.vcf.gz"
  local step_timestamp=$EPOCHSECONDS

  echo; echo "[*]  $step_name"; echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ -s "$filtered" && -s "${filtered}.tbi" ]]; then
    echo "[i]  Already completed ($filtered exists)"
    HARDFILTER_OUTPUTS+=("$filtered")
    return 0
  fi
  mkdir -p "$outdir"

  if [[ ! -s "$selected" ]]; then
    set -o xtrace
    gatk --java-options "-Xms${MEM} -Xmx${MEM}" SelectVariants \
      --reference "$ref_gnm" \
      --variant "$JOINT_VCF" \
      --select-type-to-include "$vtype" \
      --tmp-dir "$TMPDIR" \
      --verbosity "$VERBOSITY" \
      --output "$selected"
    set +o xtrace
    [[ -s "$selected" ]] || { echo "[X]  CANCELLED: SelectVariants (${vtype}) produced no output"; exit 1; }
  fi

  local filter_args=()
  local expr name
  while (( "$#" )); do
    expr="$1"; name="$2"; shift 2
    filter_args+=(--filter-expression "$expr" --filter-name "$name")
  done

  local tmp="${outdir}/.${COHORT}.${vtype,,}_filtered.tmp.vcf.gz"
  rm -f "$tmp" "${tmp}.tbi"
  set -o xtrace
  gatk --java-options "-Xms${MEM} -Xmx${MEM}" VariantFiltration \
    --reference "$ref_gnm" \
    --variant "$selected" \
    ${filter_args[@]+"${filter_args[@]}"} \
    --tmp-dir "$TMPDIR" \
    --verbosity "$VERBOSITY" \
    --output "$tmp"
  set +o xtrace

  [[ -s "$tmp" ]] || { echo "[X]  CANCELLED: VariantFiltration (${vtype}) produced no output"; exit 1; }
  mv "$tmp" "$filtered"
  mv "${tmp}.tbi" "${filtered}.tbi"
  HARDFILTER_OUTPUTS+=("$filtered")

  echo "[>]  $filtered"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Merge the SNP and INDEL hard-filtered VCFs into the final callset.
merge_hardfiltered() {
  local step_name="Merge (hard-filter)"
  local final="${OUTPUT_PATH}/${COHORT}.filtered.vcf.gz"
  local step_timestamp=$EPOCHSECONDS

  echo; echo "[*]  $step_name"; echo "[&]  $(date +%Y%m%d-%H%M)"

  if [[ -s "$final" && -s "${final}.tbi" ]]; then
    echo "[i]  Already completed ($final exists)"; return 0
  fi
  (( ${#HARDFILTER_OUTPUTS[@]} == 2 )) || { echo "[X]  CANCELLED: expected 2 filtered inputs (SNP+INDEL), got ${#HARDFILTER_OUTPUTS[@]}"; exit 1; }

  local tmp="${OUTPUT_PATH}/.${COHORT}.filtered.tmp.vcf.gz"
  rm -f "$tmp" "${tmp}.tbi"
  local inputs=()
  local f; for f in "${HARDFILTER_OUTPUTS[@]}"; do inputs+=(--INPUT "$f"); done

  set -o xtrace
  gatk MergeVcfs ${inputs[@]+"${inputs[@]}"} --OUTPUT "$tmp"
  set +o xtrace

  [[ -s "$tmp" ]] || { echo "[X]  CANCELLED: MergeVcfs produced no output"; exit 1; }
  mv "$tmp" "$final"
  [[ -s "${tmp}.tbi" ]] && mv "${tmp}.tbi" "${final}.tbi" || bcftools index --tbi --threads "$njobs" --force "$final"

  echo "[>]  $final"
  echo "[!]  $step_name"
  echo "[&]  Step time: $(echo $(( EPOCHSECONDS - step_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
}

## Finisher: verify final output and report.
finisher() {
  echo
  local final="${OUTPUT_PATH}/${COHORT}.filtered.vcf.gz"
  if [[ -s "$final" ]]; then
    echo "[>]  $final"
    echo "[$] GATK VQSR / Hard-Filtering [FENIX] completed successfully! (mode: ${MODE})"
    echo "[i]  FILTER breakdown, e.g.:"
    echo "[i]    bcftools view -H '${final}' | cut -f7 | sort | uniq -c | sort -rn"
    echo "[i]  PASS-only extraction for downstream analysis, e.g.:"
    echo "[i]    bcftools view -f PASS -Oz -o '${OUTPUT_PATH}/${COHORT}.pass.vcf.gz' '${final}'"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 0
  else
    echo "[X]  GATK VQSR / Hard-Filtering [FENIX] INCOMPLETE. Final output not found: ${final}"
    echo "[&]  Total time: $(echo $(( EPOCHSECONDS - script_timestamp )) | dc -e '?60~r60~r[[0]P]szn[:]ndZ2>zn[:]ndZ2>zp')"
    exit 1
  fi
}

# <\FUNCTIONS> --------------------------------------------------------------

# <MAIN> ------------------------------------------------------------------------

HARDFILTER_OUTPUTS=()

main() {
  resolve_mode

  if [[ "$MODE" == "vqsr" ]]; then
    guard_vqsr_resources
    recal_snp
    recal_indel
    apply_vqsr
  else
    hardfilter_one_type SNP \
      "QD < 2.0" "QD2" \
      "FS > 60.0" "FS60" \
      "MQ < 40.0" "MQ40" \
      "MQRankSum < -12.5" "MQRankSum-12.5" \
      "ReadPosRankSum < -8.0" "ReadPosRankSum-8" \
      "SOR > 3.0" "SOR3"
    hardfilter_one_type INDEL \
      "QD < 2.0" "QD2" \
      "FS > 200.0" "FS200" \
      "ReadPosRankSum < -20.0" "ReadPosRankSum-20" \
      "SOR > 10.0" "SOR10"
    merge_hardfiltered
  fi

  finisher
}

main "$@"

# <\MAIN> ---------------------------------------------------------------------

#EOF
