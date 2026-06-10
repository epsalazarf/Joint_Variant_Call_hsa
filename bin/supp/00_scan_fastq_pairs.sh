#!/usr/bin/env bash

# =============================================================================
# Title       : FASTQ Pair Synchronization Scanner [FENIX]
# Description : Scans a directory of FASTQ pairs for read count and name
#               synchronization issues. Outputs a TSV of flagged samples.
# Author      : Pavel Salazar-Fernandez (epsalazarf@gmail.com)
# Institution : LIIGH (UNAM-J)
# Date        : 2026-03-26
# Version     : 1.0
# Usage       : 00_scan_fastq_pairs.sh [input_dir] [output_tsv]
# =============================================================================

set -euo pipefail

# <ARGUMENTS> -----------------------------------------------------------------

INPUT_DIR="${1:-$PWD}"
OUTPUT_TSV="${2:-$PWD/fastq_scan_flagged.tsv}"

# <\ARGS> ---------------------------------------------------------------------

# <ENVIRONMENT> ---------------------------------------------------------------

echo
echo "[$] FASTQ Pair Synchronization Scanner [FENIX] >>"
echo "[&]  Started: $(date)"

INPUT_DIR=$(readlink -f "$INPUT_DIR")

# <\ENV> ----------------------------------------------------------------------

# <FUNCTIONS> -----------------------------------------------------------------

## Read count helper: returns read count from a gzipped FASTQ, or 0 on error
count_reads() {
  local file="$1"
  local lines
  lines=$(zcat "$file" 2>/dev/null | wc -l) || { echo 0; return 0; }
  echo $(( lines / 4 ))
}

## Name sample helper: extracts read name at a given read position (1-based)
sample_read_name() {
  local file="$1"
  local pos="$2"
  local lineno=$(( (pos - 1) * 4 + 1 ))
  zcat "$file" 2>/dev/null \
    | awk -v line="$lineno" 'NR==line {print $1; exit}' \
    | sed 's|/[12]$||'
}

## Main scan
run_scan() {
  echo
  echo "[*]  Scanning: $INPUT_DIR"

  # Discover R1 files
  local r1_list
  r1_list=$(find "$INPUT_DIR" -maxdepth 1 \
    \( -name "*_R1.f*q.gz" -o -name "*_1.f*q.gz" \) \
    | sort)

  if [ -z "$r1_list" ]; then
    echo "[X]  No R1 FASTQ files found in: $INPUT_DIR"
    exit 1
  fi

  local r1_count
  r1_count=$(echo "$r1_list" | wc -l | tr -d ' ')
  echo "[i]  R1 files found: $r1_count"

  local n_ok=0 n_flagged=0

  while IFS= read -r r1; do
    local r1_abs
    r1_abs=$(readlink -f "$r1")

    # Derive sample name: strip _R1 or _1 suffix variants
    local sample_name
    sample_name=$(basename "$r1" | sed -E 's/_R?1\.f[^.]*q\.gz$//')

    # Derive R2 path
    local r2
    r2=$(echo "$r1_abs" | sed -E 's/1(\.f[^.]*q\.gz)$/2\1/')

    local reason="" details=""

    # CHECK 1 — Missing R2
    if [ ! -f "$r2" ]; then
      reason="MISSING_R2"
      details="r2_path=${r2}"
    fi

    # CHECK 2 — Read count parity
    if [ -z "$reason" ]; then
      local count1 count2 read_error=""
      count1=$(count_reads "$r1_abs")
      count2=$(count_reads "$r2")

      if [ "$count1" -eq 0 ] || [ "$count2" -eq 0 ]; then
        # Distinguish a true zero-read file from a read error
        local raw1 raw2
        raw1=$(zcat "$r1_abs" 2>/dev/null | wc -l || echo "ERR")
        raw2=$(zcat "$r2"     2>/dev/null | wc -l || echo "ERR")
        if [ "$raw1" = "ERR" ] || [ "$raw2" = "ERR" ]; then
          reason="COUNT_MISMATCH"
          details="r1_reads=${count1},r2_reads=${count2},note=read_error"
        elif [ "$count1" -ne "$count2" ]; then
          reason="COUNT_MISMATCH"
          details="r1_reads=${count1},r2_reads=${count2}"
        fi
      elif [ "$count1" -ne "$count2" ]; then
        reason="COUNT_MISMATCH"
        details="r1_reads=${count1},r2_reads=${count2}"
      fi
    fi

    # CHECK 3 — Sampled name check
    if [ -z "$reason" ]; then
      local total=$count1
      local positions=()
      positions+=( 1 )
      positions+=( $(( total / 4 )) )
      positions+=( $(( total / 2 )) )
      positions+=( $(( (total * 3) / 4 )) )

      for pos in "${positions[@]}"; do
        [ "$pos" -lt 1 ] && pos=1
        local name1 name2
        name1=$(sample_read_name "$r1_abs" "$pos")
        name2=$(sample_read_name "$r2"     "$pos")
        if [ "$name1" != "$name2" ]; then
          reason="NAME_MISMATCH"
          details="pos=${pos},r1=${name1},r2=${name2}"
          break
        fi
      done
    fi

    # Report and record
    if [ -z "$reason" ]; then
      echo "[i]  OK        : $sample_name"
      n_ok=$(( n_ok + 1 ))
    else
      echo "[!]  FLAGGED   : $sample_name  ($reason)"
      printf "%s\t%s\t%s\t%s\t%s\n" \
        "$sample_name" "$r1_abs" "$r2" "$reason" "$details" \
        >> "$OUTPUT_TSV"
      n_flagged=$(( n_flagged + 1 ))
    fi

  done <<< "$r1_list"

  echo
  echo "[i]  Samples scanned : $(( n_ok + n_flagged ))"
  echo "[i]  Samples OK      : $n_ok"
  echo "[i]  Samples flagged : $n_flagged"
  if [ "$n_flagged" -gt 0 ]; then
    echo "[>]  $OUTPUT_TSV"
  fi
  echo "[$] Scan complete."
}

# <\FUNCTIONS> ----------------------------------------------------------------

# <MAIN> ----------------------------------------------------------------------

run_scan

# <\MAIN> ---------------------------------------------------------------------

#EOF
