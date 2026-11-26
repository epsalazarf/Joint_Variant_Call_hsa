#!/usr/bin/env bash

# =============================================================================
# Title       : Stage the GATK hg38 VQSR training resources [FENIX]
# Description : Downloads the five VQSR resource VCFs from the Broad's public
#               hg38 bundle, checks each download against the md5 the server
#               reports, keeps the canonical chromosomes (chr1-22, X, Y, M),
#               rewrites the ##contig header lines from our reference .fai (so
#               GATK sees the same contig set as hg38.fa, like the dbSNP157
#               canon_chr file) and indexes the result.
#
#               Output names are the ones Step 06 expects inside config key
#               `ref_vqsr_bundle`:
#                 hapmap_3.3.hg38.canon_chr.vcf.gz                       SNP truth+training
#                 1000G_omni2.5.hg38.canon_chr.vcf.gz                    SNP truth+training
#                 1000G_phase1.snps.high_confidence.hg38.canon_chr.vcf.gz  SNP training
#                 Mills_and_1000G_gold_standard.indels.hg38.canon_chr.vcf.gz  INDEL truth+training
#                 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.canon_chr.vcf.gz  INDEL training
#               (dbSNP, the "known" resource, is config `ref_vars`.)
#
#               ~2 GB download (1.9 GB of it the 1000G SNPs), ~10 min. Run it on
#               a LOGIN node (compute nodes may have no internet). Safe to re-run:
#               finished files are skipped, partial downloads resume.
#
# Usage       : bash STAGE.vqsr_bundle_hg38.sh [dest_dir]
#                 dest_dir  default: config remote `ref_vqsr_bundle`
#                           (/mnt/data/amedina/esalazarf/References/gatk_bundle_hg38)
#               Env: STAGE_BASE_URL  bundle URL (default: Broad public bucket)
#                    STAGE_REF_FAI   reference .fai (default: config ref_gnm + .fai)
# Source      : https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
# =============================================================================

set -euo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
CONFIG_FILE="$HERE/../../config/config.yaml"
BASE_URL="${STAGE_BASE_URL:-https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0}"
FILES="hapmap_3.3.hg38 1000G_omni2.5.hg38 1000G_phase1.snps.high_confidence.hg38 Mills_and_1000G_gold_standard.indels.hg38 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38"
CANON="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

# Environment + config (same embedded parser as the pipeline steps)
if [[ -n "${SSH_CLIENT:-}${SSH_TTY:-}${SSH_CONNECTION:-}" ]]; then env_type="remote"; else env_type="local"; fi
eval "$(
  awk -v env="$env_type" '
    BEGIN { in_env=0 }
    $1 ~ env":" { in_env=1; next }
    in_env && /^[^[:space:]]/ { in_env=0 }
    in_env && /^[[:space:]]+[a-zA-Z0-9_]+:/ {
      gsub(":", "=", $1); sub(/^[[:space:]]+/, "", $1)
      gsub(/^"/, "", $2); gsub(/"$/, "", $2); print $1 $2
    }
  ' "$CONFIG_FILE"
)"
DEST="${1:-${ref_vqsr_bundle:-}}"
[[ -n "$DEST" ]] || { echo "[X]  No destination: pass one, or set ref_vqsr_bundle in config.yaml ($env_type)"; exit 2; }
DEST="${DEST%/}"
case "$DEST/" in "$HOME"/*) echo "[X]  $DEST is under \$HOME (hard quota on FENIX) — use group storage."; exit 1 ;; esac
FAI="${STAGE_REF_FAI:-${ref_gnm:-}.fai}"
[[ -s "$FAI" ]] || { echo "[X]  Reference index not found: $FAI  (config ref_gnm + .fai, or STAGE_REF_FAI)"; exit 1; }

if [[ "$env_type" == remote ]] && command -v module >/dev/null 2>&1; then module load bcftools; fi
for t in curl bcftools awk; do command -v "$t" >/dev/null || { echo "[X]  '$t' not found"; exit 1; }; done

mkdir -p "$DEST/raw"
echo "[&]  VQSR bundle -> $DEST   ($(date))"
echo "[i]  source     : $BASE_URL"
echo "[i]  contigs    : from $FAI"

# Contigs to keep = canonical names that exist in the reference; their header
# lines, in reference order, go to a file the reheader step splices in.
keep=""; CONTIG_HDR="$DEST/.contigs.hdr"; : > "$CONTIG_HDR"
for c in $CANON; do
  len="$(awk -v c="$c" '$1==c {print $2; exit}' "$FAI")"
  [[ -n "$len" ]] || continue
  keep="${keep:+$keep,}$c"
  echo "##contig=<ID=${c},length=${len}>" >> "$CONTIG_HDR"
done
[[ -n "$keep" ]] || { echo "[X]  none of chr1..chrM found in $FAI"; exit 1; }

# md5 the server reports (base64), "" if it reports none (e.g. a file:// test).
remote_md5() {
  curl -sfI "$1" 2>/dev/null | tr -d '\r' | awk 'tolower($1)=="x-goog-hash:" && $2 ~ /^md5=/ {sub(/^md5=/,"",$2); print $2}'
}
local_md5() {   # base64 of the binary md5, like x-goog-hash
  if command -v openssl >/dev/null; then openssl dgst -md5 -binary "$1" | base64
  else md5sum "$1" | cut -d' ' -f1 | xxd -r -p | base64; fi
}

fetch() {   # $1 file name -> raw/<name>, verified
  local url="$BASE_URL/$1" out="$DEST/raw/$1" want got
  want="$(remote_md5 "$url" || true)"
  if [[ -s "$out" && -n "$want" && "$(local_md5 "$out")" == "$want" ]]; then echo "     [SKIP] $1 (already downloaded, md5 ok)"; return 0; fi
  echo "     [<]  $1"
  curl -fsSL --retry 5 --retry-delay 10 -C - -o "$out" "$url" \
    || { rm -f "$out"; curl -fsSL --retry 5 --retry-delay 10 -o "$out" "$url"; }   # resume refused -> start over
  if [[ -n "$want" ]]; then
    got="$(local_md5 "$out")"
    [[ "$got" == "$want" ]] || { echo "<ERROR> md5 mismatch for $1 (got $got, server $want) — deleted, re-run"; rm -f "$out"; exit 1; }
    echo "          md5 ok"
  else
    echo "     [!]  server gave no md5 for $1 — not verified"
  fi
}

t0=$(date +%s); n_ok=0
for stem in $FILES; do
  out="$DEST/${stem}.canon_chr.vcf.gz"
  echo
  echo "[*]  $stem"
  if [[ -s "$out" && -s "$out.tbi" ]]; then echo "     [SKIP] $(basename "$out") exists"; n_ok=$((n_ok + 1)); continue; fi
  fetch "${stem}.vcf.gz"
  fetch "${stem}.vcf.gz.tbi"
  raw="$DEST/raw/${stem}.vcf.gz"; tmp="$DEST/${stem}.canon_chr.tmp.vcf.gz"
  # Canonical chromosomes only, then the reference's contig lines in the header.
  bcftools view -r "$keep" -Oz -o "$tmp.body.vcf.gz" "$raw"
  bcftools view -h "$tmp.body.vcf.gz" \
    | awk -v add="$CONTIG_HDR" '/^##contig=/ {next} /^#CHROM/ {while ((getline l < add) > 0) print l} {print}' > "$tmp.hdr"
  bcftools reheader -h "$tmp.hdr" -o "$tmp" "$tmp.body.vcf.gz"
  rm -f "$tmp.body.vcf.gz" "$tmp.hdr"
  bcftools index -t -f "$tmp"
  n_raw=$(bcftools index -n "$raw" 2>/dev/null || echo "?")
  n=$(bcftools index -n "$tmp")
  (( n > 0 )) || { echo "<ERROR> no records left in $tmp"; exit 1; }
  mv -f "$tmp.tbi" "$out.tbi"; mv -f "$tmp" "$out"
  echo "     [DONE] $(basename "$out")   $n records (of $n_raw), $(du -h "$out" | cut -f1)"
  n_ok=$((n_ok + 1))
done

rm -f "$CONTIG_HDR"
chgrp -R "$(stat -c %G "$DEST" 2>/dev/null || id -gn)" "$DEST" 2>/dev/null || true
chmod -R g+rX "$DEST" 2>/dev/null || true
echo
echo "[OK] $n_ok/5 resources in $DEST   [$(( $(date +%s) - t0 ))s]"
echo "[i]  raw downloads kept in $DEST/raw (re-runs skip them; delete to save ~2 GB)"
[[ "${ref_vqsr_bundle:-}" == "$DEST" ]] || echo "[!]  config.yaml ($env_type) ref_vqsr_bundle is '${ref_vqsr_bundle:-}', not $DEST — update it."

#EOF
