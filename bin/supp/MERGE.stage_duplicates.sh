#!/usr/bin/env bash
# =============================================================================
# Title       : Stage duplicate-ID samples for a Step 02 re-merge [FENIX]
# Description : For samples sequenced in separate experimental batches (two
#               independent datasets under one ID), gather every Step 01 mapped
#               BAM of the same sample into ONE folder (symlinks, no copies) so
#               Step 02 merges them (the launcher passes all *.sort.bam in a
#               folder to 02 as one multiplexed input) and Step 03 calls once.
#
#               Dry-run by default: prints the plan on stdout (label, source
#               BAM, link path, status) and the pre-merge checks on stderr.
#               Nothing is created without --apply.
#
#               A group is HELD (not staged) when it has fewer than 2 BAMs, a
#               BAM fails samtools quickcheck, the @RG SM tags differ (merged
#               BAM would be multi-sample), or the @SQ dictionaries differ
#               (Picard cannot merge different reference dictionaries).
#
#               Links are named <label>_<BatchNN>[_n].sort.bam so Step 02 derives
#               the sample name <label> (text before the first underscore).
#
# Usage       : bash MERGE.stage_duplicates.sh [--apply] [--dest DIR] [--pairs FILE]
#                 --dest   default .../SLEmx-b38/BAMQC/Batch25
#                 --pairs  TSV: label <TAB> dir_A <TAB> dir_B  (default:
#                          test/idcheck/SLEmx_duplicates.pairs.tsv; lines with
#                          "_vs_" in the label are controls and are skipped)
#               Save the path list:  ... > Batch25.merge_plan.tsv
# Then        : bash bin/supp/BATCH.submit_S01-S03.sh <dest>
# =============================================================================

set -uo pipefail

HERE="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PAIRS="$REPO/test/idcheck/SLEmx_duplicates.pairs.tsv"
DEST="/mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch25"
APPLY=false

while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=true; shift ;;
    --dest)  DEST="${2:?}"; shift 2 ;;
    --pairs) PAIRS="${2:?}"; shift 2 ;;
    -h|--help) sed -n '2,29p' "$0"; exit 0 ;;
    *) echo "<ERROR> unknown argument: $1" >&2; exit 1 ;;
  esac
done

[ -f "$PAIRS" ] || { echo "<ERROR> pairs file not found: $PAIRS" >&2; exit 1; }
HAVE_ST=true
command -v samtools >/dev/null 2>&1 || { echo "<ERROR> samtools not on PATH (module load samtools): SM/@SQ checks are required" >&2; exit 1; }

echo "[i]  dest : $DEST   mode: $($APPLY && echo APPLY || echo dry-run)" >&2
printf '# label\tsource_bam\tlink_path\tstatus\n'

n_ok=0; HELD=()

while IFS=$'\t' read -r label dirA dirB _; do
  case "$label" in ''|'#'*|*_vs_*) continue ;; esac

  srcs=(); dsts=(); problems=""
  for dir in "$dirA" "$dirB"; do
    tag="$(basename "$(dirname "$dir")")"; k=0; found=0
    while IFS= read -r b; do
      [ -n "$b" ] || continue
      found=1; k=$((k + 1))
      suffix=""; [ "$k" -gt 1 ] && suffix="_$k"
      srcs+=("$b"); dsts+=("$DEST/$label/${label}_${tag}${suffix}.sort.bam")
    done < <(find "$dir" -maxdepth 1 \( -name '*.sort.bam' -o -name '*.sorted.bam' \) 2>/dev/null | sort)
    [ "$found" -eq 1 ] || problems="${problems}no S01 BAM in $dir; "
  done

  [ "${#srcs[@]}" -ge 2 ] || problems="${problems}fewer than 2 BAMs to merge; "

  if $HAVE_ST && [ "${#srcs[@]}" -ge 1 ]; then
    sm_all=""; sq_all=""
    for b in "${srcs[@]}"; do
      samtools quickcheck "$b" 2>/dev/null || problems="${problems}quickcheck failed: $(basename "$b"); "
      sm="$(samtools view -H "$b" 2>/dev/null | awk -F'\t' '$1=="@RG"{for(i=2;i<=NF;i++) if($i ~ /^SM:/) print substr($i,4)}' | sort -u | paste -sd, -)"
      sq="$(samtools view -H "$b" 2>/dev/null | grep '^@SQ' | cksum | cut -d' ' -f1)"
      echo "     $label  SM=[${sm:-none}]  SQ=$sq  $(basename "$b")" >&2
      sm_all="${sm_all}${sm}\n"; sq_all="${sq_all}${sq}\n"
    done
    [ "$(printf "$sm_all" | sort -u | wc -l | tr -d ' ')" -eq 1 ] || problems="${problems}SM tags differ across BAMs; "
    case "$(printf "$sm_all" | head -1)" in *,*) problems="${problems}a BAM holds >1 sample; " ;; esac
    [ "$(printf "$sq_all" | sort -u | wc -l | tr -d ' ')" -eq 1 ] || problems="${problems}@SQ dictionaries differ; "
  fi

  if [ -n "$problems" ]; then
    status="HELD: ${problems%; }"; HELD+=("$label  ${problems%; }")
  else
    status="OK"; n_ok=$((n_ok + 1))
  fi

  i=0
  for b in "${srcs[@]:-}"; do
    [ -n "$b" ] || continue
    printf '%s\t%s\t%s\t%s\n' "$label" "$b" "${dsts[$i]}" "$status"
    if $APPLY && [ "$status" = "OK" ]; then
      mkdir -p "$DEST/$label"
      if [ -e "${dsts[$i]}" ] || [ -L "${dsts[$i]}" ]; then
        echo "     [SKIP] link exists: ${dsts[$i]}" >&2
      else
        ln -s "$b" "${dsts[$i]}"
        for idx in "${b}.bai" "${b%.bam}.bai"; do
          [ -f "$idx" ] && { ln -s "$idx" "${dsts[$i]}.bai"; break; }
        done
      fi
    fi
    i=$((i + 1))
  done
done < "$PAIRS"

echo >&2
echo "[i]  groups OK: $n_ok   HELD: ${#HELD[@]}" >&2
if [ "${#HELD[@]}" -gt 0 ]; then
  echo "[!]  HELD (fix before merging):" >&2
  printf '       %s\n' "${HELD[@]}" >&2
fi
$APPLY || echo "[i]  dry-run only. Re-run with --apply to create the links." >&2
$APPLY && echo "[i]  next: bash $REPO/bin/supp/BATCH.submit_S01-S03.sh $DEST" >&2
#EOF
