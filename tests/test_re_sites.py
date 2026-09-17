"""Port-fidelity tests for the RE-site subtraction.

src/rnmp/re_sites.py is a port of re_build() from the authors' GPL-v3
ArtificialRiboseDetection (doi:10.5281/zenodo.8121711). These tests pin the
properties that made the port verifiable, so a later edit cannot silently change
the convention that produced the published numbers.
"""
import re

import pandas as pd
import pytest

from rnmp.re_sites import build_masks, canonical, da_patterns, load_enzymes

LIST = "data/external/res_all.list"
CHRM_FA = "data/reference/chrM.fa"


@pytest.fixture(scope="module")
def table():
    return load_enzymes(LIST)


def test_enzyme_list_parses(table):
    assert len(table) == 41
    assert table["alui"] == ("AGCT", 2)
    assert table["cviki-1"] == ("RGCY", 2)
    assert table["mlyi"] == ("GAGTCNNNNN", 10)


def test_mask_is_two_positions_per_site_per_class(table):
    """AGCTAGCT has two AluI sites; each contributes one position per strand."""
    noda, da = build_masks("AGCTAGCT", ["alui"], table)
    assert sorted(noda["alui"]) == [(1, "+"), (2, "-"), (5, "+"), (6, "-")]
    assert sorted(da["alui"]) == [(1, "-"), (2, "+"), (5, "-"), (6, "+")]


def test_noda_and_da_are_disjoint_by_strand(table):
    """Same two coordinates, opposite strand assignment -- that split IS the rule."""
    noda, da = build_masks("AGCTAGCT", ["alui"], table)
    assert not (noda["alui"] & da["alui"])
    assert {p for p, _ in noda["alui"]} == {p for p, _ in da["alui"]}


def test_overlapping_sites_are_all_found(table):
    """The authors use regex overlapped=True; the lookahead port must match that."""
    noda, _ = build_masks("RGCY".replace("R", "A").replace("Y", "C") * 1 + "GCGCGC", ["cviki-1"], table)
    seq = "AGCCGCGCGC"
    noda2, _ = build_masks(seq, ["cviki-1"], table)
    naive = len(re.findall("[AG]GC[CT]", seq))          # non-overlapping
    overlapped = len(re.findall("(?=[AG]GC[CT])", seq))
    assert overlapped >= naive
    assert len(noda2["cviki-1"]) == 2 * overlapped


def test_da_pattern_is_T_plus_downstream_half(table):
    p = da_patterns(["alui", "rsai", "msci"], table)
    assert p["alui"].pattern == "TCT"          # AGCT, cut 2 -> T + CT
    assert p["rsai"].pattern == "TAC"          # GTAC, cut 2 -> T + AC
    assert p["msci"].pattern == "TCCA"         # TGGCCA, cut 3 -> T + CCA


def test_mlyi_midpoint_is_not_its_cut_site(table):
    """Faithfulness guard: MlyI cuts OUTSIDE its recognition sequence, and the
    authors' code masks the recognition midpoint regardless. Reproduced as
    written -- this test exists so nobody "fixes" it and changes the numbers."""
    noda, _ = build_masks("GAGTCAAAAA", ["mlyi"], table)
    assert sorted(noda["mlyi"]) == [(4, "+"), (5, "-")]      # midpoint, not base 10
    assert da_patterns(["mlyi"], table)["mlyi"].pattern == "T"


def test_canonical_bridges_S2_and_res_all_naming():
    assert canonical("CviKI-1") == canonical("cviki-1") == "cviki1"
    assert canonical("Eco53KI") == canonical("Eco53kI")


def test_mask_never_covers_a_published_hotspot():
    """Required, not lucky: the authors ran this filter BEFORE calling hotspots,
    so a mask that hit one of the 48 published hotspots would be wrong."""
    import os
    if not os.path.exists(CHRM_FA):
        pytest.skip("chrM.fa absent; run scripts/00_setup_reference.sh")
    chrm = "".join(l.strip() for l in open(CHRM_FA) if not l.startswith(">")).upper()
    tbl = load_enzymes(LIST)
    by_canon = {canonical(k): k for k in tbl}
    s2 = pd.read_csv("data/external/table_S2_libraries.tsv", sep="\t")
    s4 = pd.read_csv("data/external/table_S4_common_hotspots_ef.tsv", sep="\t")
    keys = set(zip(s4.pos - 1, s4.strand))
    for _, row in s2.iterrows():
        enz = [by_canon[canonical(e)] for e in str(row.enzymes).split(";")
               if e and canonical(e) in by_canon]
        if not enz:
            continue
        noda, da = build_masks(chrm, enz, tbl)
        masked = set().union(*noda.values()) | set().union(*da.values())
        assert not (masked & keys), f"{row.library}: mask covers a published hotspot"
