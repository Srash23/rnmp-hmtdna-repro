"""rNMP composition and dinucleotide-pattern quantities (paper Figures 5 and 6A).

Thin assembly layer: the formulas live in `rnmp.metrics` -- this module only maps
a per-position count table onto them. It deliberately re-uses
`composition_normalized_frequency`, `dinucleotide_normalized_frequency` and
`strand_bias_contribution` rather than re-deriving them; those carry the Methods
transcription and are pinned by tests/test_metrics.py. Kept separate from
plotting so the tables can be validated without importing matplotlib.

The one quantity added here is the strand-aware background: the existing metric
functions take background dicts as arguments, and nothing previously computed
them from the reference.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import metrics as M

CHRM_LEN = 16569
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
NR_PATTERNS = [n + r for r in "ACGT" for n in "ACGT"]
_U = lambda b: "U" if b == "T" else b


def rnmp_base(chrm: str, pos: int, strand: str) -> str:
    """rNMP identity at a 1-based position on `strand` ('+' = light)."""
    b = chrm[pos - 1]
    return b if strand == "+" else COMP.get(b, "N")


def nr_pattern(chrm: str, pos: int, strand: str) -> str:
    """The NR dinucleotide: rNMP plus its 5' dNMP neighbour, on the rNMP's own
    strand. On the light strand the 5' neighbour is at pos-1; on the heavy
    strand, which runs antiparallel, it is the complement of pos+1. Circular.
    """
    if strand == "+":
        return chrm[(pos - 2) % CHRM_LEN] + chrm[pos - 1]
    return COMP.get(chrm[pos % CHRM_LEN], "N") + COMP.get(chrm[pos - 1], "N")


def dimer_background(seq: str) -> dict:
    """Circular 2-mer counts, keyed N+R (rNMP second). Non-ACGT 2-mers dropped."""
    bg: dict = {}
    n = len(seq)
    for i in range(n):
        p = seq[i] + seq[(i + 1) % n]
        if all(c in "ACGT" for c in p):
            bg[p] = bg.get(p, 0) + 1
    return bg


def heavy_strand(chrm: str) -> str:
    """The heavy strand read 5'->3', i.e. the reverse complement."""
    return "".join(COMP.get(b, "N") for b in reversed(chrm))


def nr_background(chrm: str) -> dict:
    """NR dinucleotide background over BOTH strands, each on its own sequence."""
    bg = dimer_background(chrm)
    for k, v in dimer_background(heavy_strand(chrm)).items():
        bg[k] = bg.get(k, 0) + v
    return bg


def per_library_composition(counts: pd.DataFrame, chrm: str, cat_of: dict) -> pd.DataFrame:
    """Normalized rNMP frequency per library for both / light / heavy strands.

    Returns long-form rows with a `panel` column in {'both','light','heavy'}.
    """
    bl, bh = M.background_counts(chrm, "+"), M.background_counts(chrm, "-")
    b_both = {b: bl[b] + bh[b] for b in "ACGT"}
    out = []
    for lib, d in counts.groupby("library"):
        base = np.array([rnmp_base(chrm, p, s) for p, s in zip(d.pos, d.strand)])
        w, st = d["count"].to_numpy(), d.strand.to_numpy()
        lc = {b: int(w[(base == b) & (st == "+")].sum()) for b in "ACGT"}
        hc = {b: int(w[(base == b) & (st == "-")].sum()) for b in "ACGT"}
        for tag, cts, bg in [("both", {b: lc[b] + hc[b] for b in "ACGT"}, b_both),
                             ("light", lc, bl), ("heavy", hc, bh)]:
            nf = M.composition_normalized_frequency({_U(b): cts[b] for b in "ACGT"}, bg)
            out.append(dict(library=lib, cat=cat_of[lib], panel=tag,
                            **{f"r{b}MP": float(nf[b]) for b in "ACGU"}))
    return pd.DataFrame(out)


def per_library_contribution(counts: pd.DataFrame, chrm: str, cat_of: dict) -> pd.DataFrame:
    """Each rNMP type's contribution to a library's strand bias (Figure 5C)."""
    out = []
    for lib, d in counts.groupby("library"):
        base = np.array([rnmp_base(chrm, p, s) for p, s in zip(d.pos, d.strand)])
        w, st = d["count"].to_numpy(), d.strand.to_numpy()
        lc = {b: int(w[(base == b) & (st == "+")].sum()) for b in "ACGT"}
        hc = {b: int(w[(base == b) & (st == "-")].sum()) for b in "ACGT"}
        bl, bh = M.background_counts(chrm, "+"), M.background_counts(chrm, "-")
        df = M.strand_bias_contribution(
            {"+": {_U(b): lc[b] for b in "ACGT"}, "-": {_U(b): hc[b] for b in "ACGT"}},
            {"+": bl, "-": bh})
        out.append(dict(library=lib, cat=cat_of[lib],
                        preferred=df.attrs["preferred_strand"],
                        bias_pp=100 * df.attrs["strand_pct_difference"],
                        **{f"r{b}MP": 100 * float(df.contribution[b]) for b in "ACGU"}))
    return pd.DataFrame(out)


def per_library_nr(counts: pd.DataFrame, chrm: str, cat_of: dict) -> pd.DataFrame:
    """Normalized NR dinucleotide frequency per library (Figure 6A).

    Normalisation is within the four patterns sharing an rNMP, per the Methods.
    """
    bg = nr_background(chrm)
    out = []
    for lib, d in counts.groupby("library"):
        pats = [nr_pattern(chrm, p, s) for p, s in zip(d.pos, d.strand)]
        cts: dict = {}
        for p, c in zip(pats, d["count"].to_numpy()):
            if all(ch in "ACGT" for ch in p):
                cts[p] = cts.get(p, 0) + int(c)
        nf = M.dinucleotide_normalized_frequency(
            {k: cts.get(k, 0) for k in NR_PATTERNS}, bg, rnmp_position="second")
        out.append(dict(library=lib, cat=cat_of[lib],
                        **{k: float(nf[k]) for k in NR_PATTERNS}))
    return pd.DataFrame(out)
