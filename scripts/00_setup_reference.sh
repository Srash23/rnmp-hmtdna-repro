#!/usr/bin/env bash
# Build every reference artefact the pipeline needs, and pin Ribose-Map.
#
#   bash scripts/00_setup_reference.sh
#
# Produces, under data/reference/:
#   chrM.fa                 GRCh38 chrM (Ensembl MT, renamed to chrM)
#   chrM.fa.fai
#   chrM.chrom.sizes        required by Ribose-Map coordinate.sh (emRiboSeq/HydEn paths)
#   chrM.{1..4,rev.1,rev.2}.bt2
#   chrM_rot<OFFSET>.fa + index   rotated copy, used to recover rNMPs spanning position 1
#   chrM_genes.bed          mtDNA gene models (CDS / rRNA / tRNA) from the Ensembl GTF
#   ribose-map/             cloned toolkit, commit recorded in ribose-map.commit
set -euo pipefail

REF=data/reference
ENS_RELEASE=108          # paper Methods: "Ensembl genome browser 108"
OFFSET=8284              # chrM length 16569 / 2, puts the 1/16569 junction mid-sequence
mkdir -p "$REF"

# ---------------------------------------------------------------- chrM FASTA
if [[ ! -f $REF/chrM.fa ]]; then
  curl -fsSL "https://ftp.ensembl.org/pub/release-${ENS_RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz" \
    | gunzip -c | sed '1s/.*/>chrM/' > "$REF/chrM.fa"
fi
samtools faidx "$REF/chrM.fa"
cut -f1,2 "$REF/chrM.fa.fai" > "$REF/chrM.chrom.sizes"
LEN=$(cut -f2 "$REF/chrM.chrom.sizes")
[[ "$LEN" == "16569" ]] || { echo "ERROR: chrM length $LEN != 16569 (wrong assembly?)" >&2; exit 1; }

# ------------------------------------------------- rotated chrM for circular recovery
# Ribose-Map processes chrM as linear, so reads spanning position 1 are unmapped and
# their rNMPs are lost. Align a second time against a rotated copy and merge
# (scripts/04_circularize.py maps rotated coordinates back to canonical ones).
if [[ ! -f $REF/chrM_rot${OFFSET}.fa ]]; then
  { echo ">chrM"
    { samtools faidx "$REF/chrM.fa" "chrM:$((OFFSET+1))-${LEN}"
      samtools faidx "$REF/chrM.fa" "chrM:1-${OFFSET}"; } | grep -v '^>' | tr -d '\n'
    echo
  } > "$REF/chrM_rot${OFFSET}.fa"
fi
samtools faidx "$REF/chrM_rot${OFFSET}.fa"
cut -f1,2 "$REF/chrM_rot${OFFSET}.fa.fai" > "$REF/chrM_rot${OFFSET}.chrom.sizes"
echo "$OFFSET" > "$REF/rotation_offset.txt"

# ---------------------------------------------------------------- bowtie2 indices
for FA in "$REF/chrM.fa" "$REF/chrM_rot${OFFSET}.fa"; do
  BASE="${FA%.fa}"
  [[ -f ${BASE}.1.bt2 ]] || bowtie2-build --quiet "$FA" "$BASE"
done

# ---------------------------------------------------------------- gene annotation
# Regenerate when absent OR empty: a 0-byte file left by a failed fetch would
# otherwise satisfy `-f` forever and silently mask the failure.
if [[ ! -s $REF/chrM_genes.bed ]]; then
  curl -fsSL "https://ftp.ensembl.org/pub/release-${ENS_RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENS_RELEASE}.chr.gtf.gz" \
    | gunzip -c | awk -F'\t' -v OFS='\t' '$1=="MT" && ($3=="CDS" || $3=="gene")' \
    | awk -F'\t' -v OFS='\t' '{
        name="."; bio=".";
        if (match($9, /gene_name "[^"]+"/))      { name=substr($9, RSTART+11, RLENGTH-12) }
        if (match($9, /gene_biotype "[^"]+"/))   { bio =substr($9, RSTART+14, RLENGTH-15) }
        print "chrM", $4-1, $5, name, $3"|"bio, $7
      }' | sort -k1,1 -k2,2n -u > "$REF/chrM_genes.bed"
  if [[ ! -s $REF/chrM_genes.bed ]]; then
    echo "ERROR: gene annotation is empty -- the Ensembl GTF fetch produced no MT" >&2
    echo "       records. Check network access to ftp.ensembl.org. Nothing in the" >&2
    echo "       current analysis consumes this file, so you may continue without" >&2
    echo "       it, but do not treat the empty file as valid annotation." >&2
    rm -f "$REF/chrM_genes.bed"
  fi
fi

# ---------------------------------------------------------------- Ribose-Map
if [[ ! -d ribose-map ]]; then
  git clone --quiet https://github.com/agombolay/ribose-map.git
fi
git -C ribose-map rev-parse HEAD > "$REF/ribose-map.commit"
chmod +x ribose-map/modules/* || true

printf '\n--- reference ready ---\n'
ls -1 "$REF"
printf 'chrM length: %s\nrotation offset: %s\nribose-map commit: %s\n' \
  "$LEN" "$OFFSET" "$(cat "$REF/ribose-map.commit")"
printf 'gene records: %s\n' "$(wc -l < "$REF/chrM_genes.bed")"
