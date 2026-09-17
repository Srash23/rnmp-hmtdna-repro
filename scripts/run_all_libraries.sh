#!/usr/bin/env bash
# Run the FS310 code path across every ribose-seq library in the sample sheet.
#
#   bash scripts/run_all_libraries.sh
#
# The per-library code path is IDENTICAL to the pilot; only the run accession and
# library id change. Library-specific inputs (fragmentation enzymes, published
# totals) are looked up from config/samples_xu2024.tsv and
# data/external/table_S2_libraries.tsv by the downstream analysis, not hardcoded.
#
# Barcode/UMI pattern are deliberately NOT passed: the deposited reads are already
# demultiplexed and UMI-extracted for this BioProject (see docs/dataset_notes.md).
#
# Results are appended to results/tables/per_library_rnmp.csv after EACH library,
# so an interrupted sweep still leaves a usable table. Re-running skips libraries
# whose coordinates already exist and replaces their summary row.
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"
mkdir -p logs results/tables
FAILED=logs/failed_libraries.txt
: > "$FAILED"

# awk, not mapfile: macOS ships bash 3.2, which has no mapfile builtin
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
             $h["role"]=="ribose_seq"{print $h["run"]"\t"$h["library"]}' \
    config/samples_xu2024.tsv > logs/.sweep_queue
N=$(wc -l < logs/.sweep_queue | tr -d ' ')
echo "=== $N ribose-seq libraries ==="

i=0
while IFS=$'\t' read -r RUN LIB; do
  i=$((i+1))
  COORD=results/ribosemap/results/${LIB}/coordinate0/${LIB}.bed

  printf '[%2d/%2d] %-6s %-12s ' "$i" "$N" "$LIB" "$RUN"
  if [[ -s $COORD ]]; then
    echo "(coordinates present, skipping pipeline)"
  else
    echo "running"
    if ! bash scripts/run_library.sh "$RUN" "$LIB" >> "logs/${LIB}.log" 2>&1; then
      echo "        FAILED -- see logs/${LIB}.log"
      echo "$LIB $RUN" >> "$FAILED"; continue
    fi
  fi
  python scripts/summarize_library.py "$LIB" || { echo "$LIB $RUN (summarize)" >> "$FAILED"; }
  # keep the workspace small: raw+trimmed FASTQ are re-fetchable, coordinates are not
  rm -f data/raw/*.fastq.gz data/trimmed/*.fq.gz
done < logs/.sweep_queue

echo
echo "=== done: $(wc -l < results/tables/per_library_rnmp.csv | awk '{print $1-1}') libraries summarized ==="
[[ -s $FAILED ]] && { echo "failures:"; cat "$FAILED"; } || echo "no failures"
