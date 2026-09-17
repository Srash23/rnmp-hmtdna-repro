#!/usr/bin/env bash
# Download SRA runs listed in a sample sheet and convert them to FASTQ.
#
#   bash scripts/01_download_sra.sh config/samples_xu2024.tsv           # everything
#   bash scripts/01_download_sra.sh config/samples_xu2024.tsv ribose_seq
#   bash scripts/01_download_sra.sh config/samples_berglund2017.tsv
#
# Arg 2 optionally filters on the `role` column (Xu sheet) so you can pull the
# 32 ribose-seq libraries without the 17 WGS control runs.
#
# Sizes: PRJNA941970 ~747 MB (49 runs). PRJNA354396 ~20 GB (44 runs) -- the
# Berglund extension is ~27x larger, so stage it by cell_type if disk is tight.
#
# Layout is read per run from the sheet; --split-files is applied for PAIRED runs
# only, so downstream globbing (_1/_2 vs bare) stays predictable.
set -euo pipefail

SHEET=${1:?usage: 01_download_sra.sh <sample_sheet.tsv> [role_filter]}
ROLE_FILTER=${2:-}
OUT=data/raw
THREADS=${THREADS:-4}
mkdir -p "$OUT" logs

# resolve column indices by name so the script survives column reordering
hdr=$(head -1 "$SHEET")
col() { echo "$hdr" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
C_RUN=$(col run); C_LAYOUT=$(col layout)
C_ROLE=$(col role 2>/dev/null || true)

tail -n +2 "$SHEET" | while IFS=$'\t' read -r -a f; do
  run=${f[$((C_RUN-1))]}
  layout=${f[$((C_LAYOUT-1))]}
  if [[ -n "$ROLE_FILTER" && -n "${C_ROLE:-}" ]]; then
    [[ "${f[$((C_ROLE-1))]}" == "$ROLE_FILTER" ]] || continue
  fi

  if compgen -G "$OUT/${run}*.fastq.gz" > /dev/null; then
    echo "[skip] $run already present"; continue
  fi

  echo "[get ] $run ($layout)"
  # ENA first: sra-tools' SDL name resolver is unreachable on restricted
  # networks, while ENA serves the same runs as gzipped FASTQ over HTTPS with a
  # published md5. Falls back to prefetch where that works.
  meta=$(curl -fsSL --max-time 120 \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv" \
    | tail -n +2 || true)
  if [[ -n "$meta" ]]; then
    n=0
    for u in $(echo "$meta" | cut -f2 | tr ';' '\n'); do
      n=$((n+1)); m=$(echo "$meta" | cut -f3 | tr ';' '\n' | sed -n "${n}p")
      f="$OUT/$(basename "$u")"
      curl -fsSL --max-time 1800 "https://${u#https://}" -o "$f"
      got=$(md5 -q "$f" 2>/dev/null || md5sum "$f" | cut -d' ' -f1)
      [[ "$got" == "$m" ]] || { echo "  md5 mismatch for $f" >&2; exit 1; }
    done
    echo "[ok  ] $run from ENA ($n file(s), md5 verified)"
    continue
  fi
  echo "[warn] $run not in ENA; falling back to prefetch"
  prefetch --max-size u --output-directory "$OUT/.sra" "$run" >> logs/prefetch.log 2>&1

  if [[ "$layout" == "PAIRED" ]]; then
    fasterq-dump --split-files --threads "$THREADS" --outdir "$OUT" "$OUT/.sra/$run" \
      >> logs/fasterq.log 2>&1
  else
    fasterq-dump --concatenate-reads --threads "$THREADS" --outdir "$OUT" "$OUT/.sra/$run" \
      >> logs/fasterq.log 2>&1
  fi

  # fasterq-dump writes uncompressed; compress immediately (20 GB -> ~6 GB for Berglund)
  for fq in "$OUT/$run".fastq "$OUT/$run"_?.fastq; do
    [[ -f "$fq" ]] && gzip -f "$fq"
  done
  rm -rf "$OUT/.sra/$run"
done

printf '\n--- downloaded ---\n'
ls -lh "$OUT" | tail -n +2 | awk '{print $9, $5}'
