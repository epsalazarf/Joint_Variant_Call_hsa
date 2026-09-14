#!/usr/bin/env bash
# =============================================================================
# Title       : Batch-wide S01→S03 submission [FENIX]
# Description : Loops PIPELINE.single_sample_01-03.sh over every sample
#               subdirectory of a batch folder (e.g. .../BAMQC/Batch12/).
#               Each sample resumes independently via the launcher's own
#               S01/S02/S03 skip logic — safe to re-run on a batch that's
#               partly done, or on one that's already fully complete (every
#               sample just gets skipped).
# Usage       : bash bin/supp/BATCH.submit_S01-S03.sh <batch_dir>
#                 e.g. bash bin/supp/BATCH.submit_S01-S03.sh \
#                        /mnt/data/amedina/esalazarf/SLEmx-b38/BAMQC/Batch12
# Note        : "reprocess" here means resume-to-completion, not force-redo.
#               To force a specific step to redo for a sample, remove its
#               output file(s) first (see INSTRUCTIONS.md), then run this.
# =============================================================================

set -uo pipefail   # no -e: one sample's failure must not abort the rest

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
LAUNCHER="$SCRIPT_DIR/../PIPELINE.single_sample_01-03.sh"
BATCH_DIR="${1:?Usage: $(basename "$0") <batch_dir>}"

[ -f "$LAUNCHER" ] || { echo "<ERROR> Launcher not found: $LAUNCHER"; exit 1; }
[ -d "$BATCH_DIR" ] || { echo "<ERROR> Batch directory not found: $BATCH_DIR"; exit 1; }
BATCH_DIR="$(cd "$BATCH_DIR" && pwd)"   # normalize, strips a trailing slash

SAMPLES=()
for d in "$BATCH_DIR"/*/; do
  [ -d "$d" ] && SAMPLES+=("${d%/}")
done
[ ${#SAMPLES[@]} -gt 0 ] || { echo "<ERROR> No sample subdirectories in $BATCH_DIR"; exit 1; }

echo "[i] $(basename "$BATCH_DIR"): ${#SAMPLES[@]} sample dir(s)"
echo

FAILED=()
for d in "${SAMPLES[@]}"; do
  s="$(basename "$d")"
  echo "== $s =="
  bash "$LAUNCHER" "$d" || { echo "   [X] launcher exited non-zero for $s"; FAILED+=("$s"); }
  echo
done

echo "=============================================================="
echo "[i] Looped ${#SAMPLES[@]} sample(s) in $(basename "$BATCH_DIR")."
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "[!] Launcher reported a problem for:"
  printf '      %s\n' "${FAILED[@]}"
fi
echo "[i] Monitor: squeue -u \$USER"
#EOF
