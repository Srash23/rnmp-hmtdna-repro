"""Restriction-enzyme background subtraction (paper Methods, "RE-recognition sites").

DERIVED WORK -- GNU GPL v3.0
============================
The position-selection logic in `build_masks` is a direct port of `re_build()` in
`reBuild.py` from ArtificialRiboseDetection by Penghao Xu (Storici lab),
https://github.com/xph9876/ArtificialRiboseDetection, deposited as
doi:10.5281/zenodo.8121711 and released under GNU GPL v3.0. Xu, Yang et al. 2024
name that deposit as the source of "the scripts used to remove rNMPs at
restriction enzyme locations", and their Methods say the protocol follows
Balachander et al., Nat Commun 2020;11:2447 (their ref. 21).

Because this module is a derivative of GPL-v3 code, THIS FILE is GPL v3.0.
`data/external/res_all.list` is copied verbatim from the same repository.
The port was verified against the authors' `re_build()` run directly on GRCh38
chrM: identical position sets for all 18 enzymes used in the study (see
tests/test_re_sites.py).

What the filter actually does
----------------------------
The rule is much narrower than masking whole recognition sites -- which is worth
stating because masking whole sites removes far too much (it undershoots the
published rNMP counts, while no masking overshoots them). For each match of a
recognition pattern, the authors take the *midpoint* of the matched site and mask
exactly two positions, split by strand into two classes:

    mid = (match_start + match_end) // 2        # 0-based, = base 3' of a blunt cut

    noda : (mid-1, '+') and (mid, '-')   -- removed unconditionally
    da   : (mid-1, '-') and (mid, '+')   -- removed ONLY if the read carries a dA tail

`noda` positions are the 3'-terminal base of a blunt-cut fragment end on each
strand: an undigested fragment end that survives T5 exonuclease and is then
miscalled as an rNMP. `da` positions are one base further in -- where the called
rNMP lands when the dA-tailing step added a non-templated dA to that same end --
so they are only artefacts if the read actually shows the tail. The confirmation
pattern is `'T' + recognition_pattern[cut_position:]`, matched against the start
of the raw read.

Two faithfulness notes, both deliberate:

1. The `cut-after` column of `res_all.list` is used ONLY to build the dA
   confirmation pattern. The masked positions come from the recognition-site
   midpoint regardless. For the 17 blunt cutters used here midpoint == cut site,
   but for MlyI (GAGTCNNNNN, cuts after base 10) it does not -- MlyI cuts outside
   its recognition sequence. This is reproduced as written rather than
   "corrected", because the published numbers were produced by this code.
2. The authors match against raw FASTQ reads and skip `d = 11` leading bases
   (3-nt barcode + 8-nt UMI). The reads deposited for PRJNA941970 are already
   demultiplexed and UMI-extracted, so `dA_offset` defaults to 0 here. See
   docs/dataset_notes.md S2.
"""
from __future__ import annotations

import re
from pathlib import Path

__all__ = ["load_enzymes", "build_masks", "da_patterns", "IUPAC_REGEX"]

# exactly the substitution table in checkRes.py
IUPAC_REGEX = {"W": "[AT]", "S": "[CG]", "M": "[AC]", "K": "[GT]", "R": "[AG]",
               "Y": "[CT]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
               "V": "[ACG]", "N": "[ACGT]"}


def _to_regex(iupac: str) -> str:
    out = iupac
    for k, v in IUPAC_REGEX.items():
        out = out.replace(k, v)
    return out


def load_enzymes(path="data/external/res_all.list") -> dict[str, tuple[str, int]]:
    """name.lower() -> (IUPAC recognition pattern, cut-after position)."""
    table = {}
    for line in Path(path).read_text().splitlines():
        w = line.split("\t")
        if len(w) != 3:
            continue
        table[w[0].strip().lower()] = (w[1].strip(), int(w[2]))
    return table


def canonical(name: str) -> str:
    """Enzyme-name key tolerant of the punctuation/case drift between S2 and res_all.list."""
    return name.lower().replace("-", "").replace("_", "").replace(" ", "")


def build_masks(seq: str, enzymes, table):
    """Port of re_build(). Returns (noda, da), each enzyme -> set of (bed_start, strand).

    `bed_start` is 0-based, matching BED column 2 -- i.e. 1-based position minus 1.
    """
    noda = {e: set() for e in enzymes}
    da = {e: set() for e in enzymes}
    for e in enzymes:
        iupac = table[e][0]
        pat, L = _to_regex(iupac), len(iupac)
        # lookahead gives the overlapped matching that the authors get from
        # regex.finditer(..., overlapped=True); every pattern is fixed length
        for m in re.finditer(f"(?=({pat}))", seq):
            s = m.start()
            mid = (s + s + L) // 2
            noda[e].add((mid - 1, "+"))
            noda[e].add((mid, "-"))
            da[e].add((mid - 1, "-"))
            da[e].add((mid, "+"))
    return noda, da


def da_patterns(enzymes, table) -> dict[str, re.Pattern]:
    """enzyme -> compiled 'T' + pattern[cut:] confirmation regex (checkRes.py)."""
    out = {}
    for e in enzymes:
        iupac, cut = table[e]
        out[e] = re.compile(_to_regex("T" + iupac[cut:]))
    return out
