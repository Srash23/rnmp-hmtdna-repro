"""Single-library rNMP pipeline: post-alignment coordinate calling and QC.

The read-level work (download, trim, align, Ribose-Map Coordinate) lives in
`scripts/run_library.sh`, which shells out to the real tools. This module holds
the parts that are ours rather than Ribose-Map's:

  * `call_rnmps`      -- parse Ribose-Map's alignment SAM into rNMP coordinates
                         and apply the reference-mismatch filter
  * `restriction_mask`-- positions inside a restriction-enzyme recognition site
  * `strand_summary`  -- light/heavy split of an rNMP set
  * `composition`     -- rNMP base composition from Ribose-Map's per-base BEDs

Coordinate convention (Ribose-Map `modules/coordinate.sh`, ribose-seq branch):
a read on '-' yields an rNMP at (end-1, end) on '+'; a read on '+' yields one at
(start, start+1) on '-'. Both are the read's 5'-most aligned base, reported on
the opposite strand. Because SAM SEQ is stored in reference orientation, that
base is SEQ[-1] for a reverse-strand read and SEQ[0] for a forward-strand read.
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

__all__ = ["load_chrm", "call_rnmps", "restriction_mask", "strand_summary",
           "composition", "RESTRICTION_SETS", "IUPAC"]

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "N": "ACGT"}
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "N": "N"}

# Suppl. Table S2. RE4 is used by liver libraries FS261-FS264 and is NOT
# described in the paper's Methods, which list only RE1-RE3.
RESTRICTION_SETS = {
    "RE1": {"HpyCH4V": "TGCA", "Hpy166II": "GTNNAC", "Eco53KI": "GAGCTC",
            "RsaI": "GTAC", "StuI": "AGGCCT"},
    "RE2": {"AleI": "CACNNNNGTG", "AluI": "AGCT", "PvuII": "CAGCTG",
            "DraI": "TTTAAA", "HaeIII": "GGCC", "SspI": "AATATT"},
    "RE3": {"CviKI-1": "RGCY", "MlyI": "GAGTC", "MscI": "TGGCCA", "MslI": "CAYNNNNRTG"},
    "RE4": {"AluI": "AGCT", "Hpy166II": "GTNNAC", "PmlI": "CACGTG", "StuI": "AGGCCT",
            "Cac8I": "GCNNGC", "HpyCH4V": "TGCA", "PmeI": "CACGTC"},
}


def load_chrm(fasta: str | Path) -> str:
    seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">")).upper()
    if len(seq) != 16569:
        raise ValueError(f"chrM length {len(seq)} != 16569 -- wrong assembly?")
    return seq


def _open_text(path):
    """Open a path as text, transparently handling .gz -- and falling back to a
    .gz sibling when the plain file is absent.

    The repository ships the pilot library's alignment gzipped (78 kB vs 1.6 MB)
    so that Part A of the notebook runs on a fresh clone; the full SAM set is
    11 GB and is never committed.
    """
    import gzip
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    if not path.exists() and path.with_suffix(path.suffix + ".gz").exists():
        return gzip.open(path.with_suffix(path.suffix + ".gz"), "rt")
    return open(path)


def call_rnmps(sam_path: str | Path, chrm: str) -> pd.DataFrame:
    """rNMP coordinates from an alignment SAM, with a reference-mismatch flag.

    Returns one row per aligned read: `pos` (1-based), `strand` (the rNMP's
    strand, opposite the read's), `read_base`, `ref_base`, `mismatch`.
    Paper Methods: "rNMPs that mismatch the nucleotide at their alignment
    locations in the reference genome were removed."
    """
    rows = []
    for line in _open_text(sam_path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        flag, pos, cigar, seq = int(f[1]), int(f[3]), f[5], f[9]
        if flag & 4:
            continue
        alen = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
                   if op in "MDN=X")
        if flag & 16:
            gpos, read_base, strand = pos + alen - 1, seq[-1], "+"
        else:
            gpos, read_base, strand = pos, seq[0], "-"
        ref_base = chrm[gpos - 1]
        rows.append((gpos, strand, read_base, ref_base, read_base != ref_base))
    df = pd.DataFrame(rows, columns=["pos", "strand", "read_base", "ref_base", "mismatch"])
    df.attrs["n_aligned"] = len(df)
    return df


def restriction_mask(chrm: str, enzymes: dict[str, str]) -> set[int]:
    """1-based chrM positions lying inside any recognition site, either strand.

    NOTE: this masks the whole recognition site. The paper's own filter comes
    from the authors' Zenodo deposit (10.5281/zenodo.8121711) and is narrower --
    on FS310 the whole-site mask removes 355 rNMPs where the published count
    implies ~33. Use this as an upper bound, and see docs/dataset_notes.md.
    """
    masked: set[int] = set()
    for motif in enzymes.values():
        for m in (motif, "".join(_COMP[c] for c in reversed(motif))):
            pat = "".join(f"[{IUPAC[c]}]" for c in m)
            for hit in re.finditer(f"(?=({pat}))", chrm):
                masked.update(range(hit.start() + 1, hit.start() + 1 + len(m)))
    return masked


def strand_summary(coords: pd.DataFrame, label: str = "") -> dict:
    """Light (+) / heavy (-) split. Figure 1A of the paper."""
    light = int((coords.strand == "+").sum())
    heavy = int((coords.strand == "-").sum())
    total = light + heavy
    return dict(stage=label, light=light, heavy=heavy, total=total,
                light_pct=100 * light / total, heavy_pct=100 * heavy / total)


def composition(coord_dir: str | Path, sample: str, unit: str = "chrM") -> pd.Series:
    """rNMP base composition from Ribose-Map's per-base coordinate BEDs."""
    counts = {}
    for base in "ACGU":
        p = Path(coord_dir) / f"{sample}-{unit}.{base}.bed"
        if p.exists() or p.with_suffix(".bed.gz").exists():
            counts[base] = sum(1 for _ in _open_text(p))
        else:
            counts[base] = 0
    s = pd.Series(counts, dtype=float)
    return 100 * s / s.sum()
