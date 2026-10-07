#!/usr/bin/env bash
# =============================================================================
# Title       : GVCF Pair Identity Check [FENIX]
# Description : Decides whether two sample dirs hold the SAME individual by
#               comparing their per-sample GVCF genotypes on one chromosome
#               (default chr20). Meant for duplicate IDs sequenced twice: same
#               person → merge the datasets at Step 02; different → rename.
#               Read-only.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-10-06
# Version     : 1.0
# Usage       : IDCHECK.gvcf_pair_concordance.sh <pairs.tsv> [-o out.tsv] [-c chrom] [-q min_gq]
#               pairs.tsv — TAB-separated, one comparison per line:
#                           <label> <sample_dir_A> <sample_dir_B>
#                           (blank lines and #comments ignored)
#               -o  output TSV (default: ./<pairs basename>.idcheck.<date>.tsv)
#               -c  chromosome to compare (default chr20)
#               -q  minimum GQ for a genotype to be compared (default 20)
#
# Method: SNP sites called non-ref (GQ >= q) in EITHER sample form the site
#   list. Each sample's genotype at every site is read from its GVCF: variant
#   records give the called genotype, hom-ref blocks give 0/0 at the block's GQ,
#   so "variant in A, reference in B" is a real comparison, not a missing site.
#   Per pair:
#     discord_pct — % of compared sites whose alt-allele dosage (0/1/2) differs
#     ibs0_pct    — % of compared sites where one is 0/0 and the other 1/1
#                   (opposite homozygotes: should be ~0 for the same person,
#                   even at low depth where hets get miscalled as homs)
#   Verdict (heuristic):
#     SAME      ibs0 < 0.5% and discord < 25%
#     RELATIVE? ibs0 < 2%   and discord >= 25%  (parent/child share an allele at
#               every site, so ibs0 stays ~0 while discordance is high)
#     DIFFERENT ibs0 >= 2%
#     UNCLEAR   anything else (ibs0 0.5-2% with low discordance)
#   Simulation (GQ>=20): same person 6x/6x and 6x/22x → discord ~1.4%, ibs0 0;
#   parent/child → 52%, 0.5%; unrelated → 66-70%, 6-8%. CALIBRATE on real data:
#   include a couple of known-different pairs (same batch/depth) as controls.
# =============================================================================

set -uo pipefail   # no -e: one bad pair must not abort the others

# <ARGUMENTS> -----------------------------------------------------------------

PAIRS=""
OUT_TSV=""
CHROM="chr20"
MIN_GQ=20
while [ $# -gt 0 ]; do
  case "$1" in
    -o) OUT_TSV="${2:?-o needs a file name}"; shift 2 ;;
    -c) CHROM="${2:?-c needs a chromosome}"; shift 2 ;;
    -q) MIN_GQ="${2:?-q needs a number}"; shift 2 ;;
    -h|--help) sed -n '2,38p' "$0"; exit 0 ;;
    *)  PAIRS="$1"; shift ;;
  esac
done
[ -n "$PAIRS" ] || { echo "Usage: $(basename "$0") <pairs.tsv> [-o out.tsv] [-c chrom] [-q min_gq]"; exit 1; }
[ -s "$PAIRS" ] || { echo "<ERROR> Pairs file missing or empty: $PAIRS"; exit 1; }
[ -n "$OUT_TSV" ] || OUT_TSV="$PWD/$(basename "${PAIRS%.*}").idcheck.$(date +%F).tsv"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

if ! command -v bcftools >/dev/null; then
  module load bcftools >/dev/null 2>&1 || true
  command -v bcftools >/dev/null || { echo "<ERROR> bcftools not found in PATH"; exit 1; }
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/idcheck.XXXXXX") || { echo "<ERROR> Cannot create temp dir"; exit 1; }
trap 'rm -rf "$WORK"' EXIT HUP TERM

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## GVCF holding $CHROM for a sample dir: chrom_gvcf/ shard, else canon GVCF (needs .tbi)
## Prints "<file>\t<region arg>" — region empty for a shard (whole file is the chrom)
find_gvcf() {
  local d="$1" f
  f=$(find -H "$d/chrom_gvcf" -maxdepth 1 \( -type f -o -type l \) -name "*.raw_vars.${CHROM}.g.vcf.gz" 2>/dev/null | sort | head -1)
  [ -n "$f" ] && { printf '%s\t\n' "$f"; return; }
  f=$(find -H "$d" -maxdepth 1 \( -type f -o -type l \) -name '*.raw_variants.canon_chr.g.vcf.gz' 2>/dev/null | sort | head -1)
  [ -n "$f" ] && [ -s "${f}.tbi" ] && { printf '%s\t%s\n' "$f" "$CHROM"; return; }
  return 1
}

## Stream one GVCF as "POS END ALT TGT GQ" (END = POS for non-block records)
gvcf_stream() {   # $1 = file  $2 = region ("" = whole file)
  bcftools query ${2:+-r "$2"} -f '%POS\t%INFO/END\t%ALT\t[%TGT]\t[%GQ]\n' "$1" 2>/dev/null \
    | awk -F'\t' -v OFS='\t' '{ if ($2 == ".") $2 = $1; print }'
}

## Genotype dosage of one sample at the union site list → "POS ALT DOSAGE GQ"
## Sites file: "POS REF ALT" sorted by POS. Dosage = count of ALT allele; a
## genotype carrying any third allele gets -1 (counts as discordant).
dosages() {   # $1 = file  $2 = region  $3 = sites file
  gvcf_stream "$1" "$2" | awk -F'\t' -v OFS='\t' -v sites="$3" '
    BEGIN { n = 0; while ((getline l < sites) > 0) { split(l, s, "\t"); n++; P[n] = s[1]; R[n] = s[2]; A[n] = s[3] } i = 1 }
    {
      pos = $1; end = $2
      while (i <= n && P[i] < pos) i++                      # sites before this record: no call
      if ($3 == "<NON_REF>") {                              # hom-ref block covers [pos, end]
        if ($5 == ".") next
        for (j = i; j <= n && P[j] <= end; j++) print P[j], A[j], 0, $5
        next
      }
      for (j = i; j <= n && P[j] == pos; j++) {             # variant record at a site
        m = split($4, g, /[\/|]/); if (m != 2 || g[1] == "." || $5 == ".") continue
        d = 0
        for (k = 1; k <= 2; k++) { if (g[k] == A[j]) d++; else if (g[k] != R[j]) d = -99 }
        print P[j], A[j], (d < 0 ? -1 : d), $5
      }
    }'
}

## Compare one pair → one TSV row
compare_pair() {   # $1 label  $2 dirA  $3 dirB
  local label="$1" da="$2" db="$3" fa fb ra rb
  local na nb
  IFS=$'\t' read -r fa ra <<< "$(find_gvcf "$da")"
  IFS=$'\t' read -r fb rb <<< "$(find_gvcf "$db")"
  if [ -z "${fa:-}" ] || [ -z "${fb:-}" ]; then
    printf '%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNO_GVCF\n' "$label" "$da" "$db"
    return
  fi

  # Confident SNP sites per sample → union "POS REF ALT" list
  local sa="$WORK/a.sites" sb="$WORK/b.sites" su="$WORK/union.sites"
  snp_sites "$fa" "$ra" > "$sa"; snp_sites "$fb" "$rb" > "$sb"
  na=$(wc -l < "$sa" | tr -d ' '); nb=$(wc -l < "$sb" | tr -d ' ')
  sort -k1,1n -k3,3 -u "$sa" "$sb" > "$su"

  # one genotype per site (overlapping records, e.g. inside a deletion, keep the first)
  dosages "$fa" "$ra" "$su" | awk -F'\t' '!seen[$1 FS $2]++' > "$WORK/a.dos"
  dosages "$fb" "$rb" "$su" | awk -F'\t' '!seen[$1 FS $2]++' > "$WORK/b.dos"

  join -t $'\t' -j1 <(awk -F'\t' '{print $1":"$2"\t"$3"\t"$4}' "$WORK/a.dos" | sort -t $'\t' -k1,1) \
                    <(awk -F'\t' '{print $1":"$2"\t"$3"\t"$4}' "$WORK/b.dos" | sort -t $'\t' -k1,1) \
    | awk -F'\t' -v q="$MIN_GQ" -v OFS='\t' -v label="$label" -v da="$da" -v db="$db" -v na="$na" -v nb="$nb" '
        $3 < q || $5 < q { next }
        { n++; if ($2 != $4) dis++; if (($2 == 0 && $4 == 2) || ($2 == 2 && $4 == 0)) ibs0++ }
        END {
          if (n == 0) { print label, da, db, na, nb, 0, "NA", "NA", "NO_OVERLAP"; exit }
          dp = 100 * dis / n; ip = 100 * ibs0 / n
          if      (ip < 0.5 && dp < 25) v = "SAME"
          else if (ip < 2 && dp >= 25)  v = "RELATIVE?"
          else if (ip >= 2)             v = "DIFFERENT"
          else                          v = "UNCLEAR"
          printf "%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.3f\t%s\n", label, da, db, na, nb, n, dp, ip, v
        }'
}

## Confident SNP sites of one sample → "POS REF ALT"
snp_sites() {   # $1 = file  $2 = region
  bcftools query ${2:+-r "$2"} -i "GT=\"alt\" && FMT/GQ>=${MIN_GQ}" \
      -f '%POS\t%REF\t[%TGT]\n' "$1" 2>/dev/null \
    | awk -F'\t' -v OFS='\t' '
        length($2) != 1 { next }                            # SNP records only (REF one base)
        {
          if (split($3, g, /[\/|]/) != 2) next
          alt = ""
          for (k = 1; k <= 2; k++) {
            if (g[k] !~ /^[ACGT]$/) next                    # missing / indel / symbolic
            if (g[k] != $2) { if (alt != "" && alt != g[k]) next; alt = g[k] }
          }
          if (alt != "") print $1, $2, alt
        }'
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

echo
echo "[$] GVCF Pair Identity Check >>"
echo "[i]  Pairs : $PAIRS"
echo "[i]  Chrom : $CHROM   (min GQ $MIN_GQ)"

{
  printf 'label\tdir_a\tdir_b\tsnps_a\tsnps_b\tcompared\tdiscord_pct\tibs0_pct\tverdict\n'
  while IFS=$'\t' read -r label da db _rest; do
    [ -z "${label// }" ] && continue
    case "$label" in \#*) continue ;; esac
    if [ -z "${db:-}" ]; then echo "[!]  Skipping malformed line: $label" >&2; continue; fi
    echo "[&]  $label" >&2
    compare_pair "$label" "$da" "$db"
  done < "$PAIRS"
} > "${OUT_TSV}.tmp"

mv "${OUT_TSV}.tmp" "$OUT_TSV"
echo "[>]  $OUT_TSV"
echo "[i]  Verdicts: $(cut -f9 "$OUT_TSV" | tail -n +2 | sort | uniq -c | awk '{printf "%s=%s ", $2, $1}')"

# <\MAIN> ---------------------------------------------------------------------

#EOF
