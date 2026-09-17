#!/usr/bin/env bash
# End-to-end rNMP pipeline for ONE ribose-seq library.
#
#   bash scripts/run_library.sh <RUN_ACCESSION> <LIBRARY_ID>
#   bash scripts/run_library.sh SRR23726637 FS310
#
# Stages: fetch -> trim -> Ribose-Map Alignment -> Ribose-Map Coordinate.
# Downstream filtering, strand split and figures are done by the notebook via
# src/rnmp/pipeline.py, so this script stops at rNMP coordinates.
#
# Requires: conda env rnmp-pipeline, and scripts/00_setup_reference.sh already run.
#
# WHY ENA AND NOT prefetch: sra-tools resolves runs through NCBI's SDL name
# service, which is not reachable from a restricted network. ENA mirrors the same
# runs as gzipped FASTQ over plain HTTPS and publishes an md5 per file, so this
# script fetches from ENA and verifies the checksum. `prefetch` is kept as a
# fallback for environments where it does work.
set -euo pipefail

RUN=${1:?usage: run_library.sh <RUN_ACCESSION> <LIBRARY_ID>}
LIB=${2:?usage: run_library.sh <RUN_ACCESSION> <LIBRARY_ID>}
THREADS=${THREADS:-4}
REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"

REF=data/reference
ADAPTER=AGTTGCGACACGGATCTATCA      # Methods; = 5' end of Adapter.S (Suppl. Table S1)
mkdir -p data/raw data/trimmed results/ribosemap logs config

[[ -f $REF/chrM.1.bt2 ]] || { echo "ERROR: run scripts/00_setup_reference.sh first" >&2; exit 1; }

# ------------------------------------------------------------------ 1. fetch
# ENA names submitted files inconsistently across this BioProject: some runs
# expose a single <RUN>.fastq.gz even when SRA calls them PAIRED, others expose
# <RUN>_1/_2. Take the FIRST listed file, which is read 1 in both cases --
# Ribose-Map's ribose-seq coordinate path uses read 1 only.
# ENA's filereport ALWAYS prepends run_accession, whatever `fields` asks for,
# so the requested fields start at column 2. Request run_accession explicitly so
# the column positions are stated rather than implied.
META=$(curl -fsSL --max-time 120 \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${RUN}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv" \
  | tail -n +2)
URL=$(echo "$META" | cut -f2 | tr ';' '\n' | head -1)
MD5=$(echo "$META" | cut -f3 | tr ';' '\n' | head -1)
[[ "$URL" == *.fastq.gz ]] || { echo "ERROR: ENA gave no FASTQ URL for $RUN (got '$URL')" >&2; exit 1; }
[[ -n "$URL" ]] || { echo "ERROR: ENA has no FASTQ for $RUN" >&2; exit 1; }
FQ=data/raw/$(basename "$URL")
BASE=$(basename "$FQ" .fastq.gz)
if [[ ! -f $FQ ]]; then
  echo "[1/4] fetching $RUN from ENA ($(basename "$URL"))"
  curl -fsSL --max-time 3600 "https://${URL#https://}" -o "$FQ"
  GOT=$(md5 -q "$FQ" 2>/dev/null || md5sum "$FQ" | cut -d' ' -f1)
  [[ "$GOT" == "$MD5" ]] || { echo "ERROR: md5 mismatch ($GOT != $MD5)" >&2; exit 1; }
  echo "       md5 verified: $GOT"
fi

# ------------------------------------------------------------------ 2. trim
# Paper Methods: trim_galore -a <adapter> -q 15 --length 62
TRIM=data/trimmed/${BASE}_trimmed.fq.gz
if [[ ! -f $TRIM ]]; then
  echo "[2/4] trimming"
  trim_galore -a "$ADAPTER" -q 15 --length 62 "$FQ" -o data/trimmed/ >> logs/trim.log 2>&1
fi

# ------------------------------------------------- 3. Ribose-Map config + Alignment
# NOTE ON UMIs: reads deposited in SRA/ENA for this BioProject have already been
# adapter-demultiplexed and UMI-extracted by the authors -- the 11-nt
# NNNNNNXXXNN prefix is absent (no 3-mer dominates read positions 7-9, and 99.8%
# of RAW reads align to chrM). `pattern`/`barcode` are therefore deliberately
# left OUT of the config so Ribose-Map's Alignment module skips umi_tools
# extract, demultiplexing and dedup. Setting them against this data would
# produce an empty FASTQ. See docs/dataset_notes.md.
CFG=config/${LIB}.config
cat > "$CFG" <<CFGEOF
sample='${LIB}'
technique='ribose-seq'
repository='${REPO}/results/ribosemap'
basename='${REPO}/${REF}/chrM'
fasta='${REPO}/${REF}/chrM.fa'
read1='${REPO}/${TRIM}'
quality='0'
threads='${THREADS}'
seed='123456789'
mismatches='0'
percentile='0.99'
units='chrM'
CFGEOF

echo "[3/4] Ribose-Map Alignment"
ribose-map/modules/ribosemap alignment "$CFG"

# ------------------------------------------------------------- 4. Coordinate
# The "chromosomes"/divide-by-zero messages come from subset.sh and refFreqs.sh
# running their all-chromosomes pass on a chrM-only reference. Harmless.
echo "[4/4] Ribose-Map Coordinate"
ribose-map/modules/ribosemap coordinate "$CFG" 2>&1 | grep -v 'chromosomes\|Divide by zero' || true

OUT=results/ribosemap/results/${LIB}/coordinate0/${LIB}.bed
printf '\n--- %s: %s rNMP coordinates ---\n' "$LIB" "$(wc -l < "$OUT")"
awk '{c[$6]++} END {printf "  light (+): %d\n  heavy (-): %d\n", c["+"], c["-"]}' "$OUT"
echo "  -> $OUT"
