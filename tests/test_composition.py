"""Tests for the strand-aware background and the NR dinucleotide orientation.

The formulas themselves (`composition_normalized_frequency`,
`dinucleotide_normalized_frequency`, `strand_bias_contribution`) are pinned by
tests/test_metrics.py; this file covers what `rnmp.composition` adds on top --
the background computed from the reference, and the mapping from positions to
dinucleotide patterns.

The orientation test matters most: getting the 5'-neighbour side wrong produces
a complete, plausible-looking heatmap, so it can only be caught by pinning it
against the paper's own column ordering.
"""
import numpy as np
import pandas as pd
import pytest

from rnmp import composition as CMP
from rnmp import metrics as M

CHRM_FA = "data/reference/chrM.fa"
S5_NR = "data/external/table_S5_nr_pvalues.tsv"


@pytest.fixture(scope="module")
def chrm():
    import os
    if not os.path.exists(CHRM_FA):
        pytest.skip("chrM.fa absent; run scripts/00_setup_reference.sh")
    return "".join(l.strip() for l in open(CHRM_FA) if not l.startswith(">")).upper()


# ---------------------------------------------------------------- background --
def test_heavy_background_is_the_complement_of_the_light(chrm):
    bl, bh = M.background_counts(chrm, "+"), M.background_counts(chrm, "-")
    assert bh == {"A": bl["T"], "C": bl["G"], "G": bl["C"], "T": bl["A"]}


def test_backgrounds_are_one_short_of_the_genome_length(chrm):
    """GRCh38 chrM carries a single N at 3107 -- the rCRS placeholder kept so
    that GRCh38 chrM numbering matches rCRS. It is excluded from every
    background, so the four counts sum to 16568, not 16569."""
    assert chrm[3106] == "N"
    assert sum(M.background_counts(chrm, "+").values()) == 16568


def test_strand_specific_background_is_not_optional(chrm):
    """chrM is C-rich on the light strand, so one shared background would
    manufacture a composition difference out of the reference alone."""
    bl = M.background_counts(chrm, "+")
    assert bl["C"] > 2 * bl["G"]


def test_background_rejects_a_bad_strand(chrm):
    with pytest.raises(ValueError):
        M.background_counts(chrm, "light")


# ---------------------------------------- formulas, via the existing metrics --
def test_composition_uses_the_strand_background_it_is_given():
    """Equal counts on unequal backgrounds must not give equal frequencies."""
    bg = {"A": 100, "C": 200, "G": 50, "T": 150}
    nf = M.composition_normalized_frequency({"A": 10, "C": 10, "G": 10, "U": 10}, bg)
    assert pytest.approx(nf.sum()) == 1.0
    assert nf["G"] > nf["C"]


def test_composition_matches_U_to_reference_T():
    bg = {"A": 100, "C": 100, "G": 100, "T": 400}
    nf = M.composition_normalized_frequency({"A": 10, "C": 10, "G": 10, "U": 40}, bg)
    assert nf["U"] == pytest.approx(nf["A"])


def test_nr_normalization_is_within_the_rnmp_group():
    keys = [n + r for r in "AC" for n in "ACGT"]
    nf = M.dinucleotide_normalized_frequency(
        {k: (5 if k[-1] == "A" else 9) for k in keys},
        {k: 1000 for k in keys}, rnmp_position="second")
    for r in "AC":
        assert pytest.approx(nf[[k for k in keys if k[-1] == r]].sum()) == 1.0


def test_contribution_sums_to_the_bias_level(chrm):
    bl, bh = M.background_counts(chrm, "+"), M.background_counts(chrm, "-")
    df = M.strand_bias_contribution(
        {"+": {"A": 900, "C": 300, "G": 200, "U": 50},
         "-": {"A": 100, "C": 200, "G": 150, "U": 20}}, {"+": bl, "-": bh})
    assert df.attrs["preferred_strand"] == "+"
    assert df.contribution.sum() == pytest.approx(df.attrs["strand_pct_difference"])


def test_types_opposing_the_bias_contribute_zero(chrm):
    """The Methods build the contribution percentage from the positive
    differences only, so a type opposing the bias must contribute zero while
    remaining visible in `difference`."""
    bl, bh = M.background_counts(chrm, "+"), M.background_counts(chrm, "-")
    df = M.strand_bias_contribution(
        {"+": {"A": 5000, "C": 1, "G": 1, "U": 1},
         "-": {"A": 10, "C": 900, "G": 1, "U": 1}}, {"+": bl, "-": bh})
    assert df.difference["C"] < 0
    assert df.contribution["C"] == 0.0
    assert df.contribution["A"] == df.contribution.max()


def test_a_balanced_library_has_zero_contribution(chrm):
    bl, bh = M.background_counts(chrm, "+"), M.background_counts(chrm, "-")
    same = {"A": 100, "C": 100, "G": 100, "U": 100}
    df = M.strand_bias_contribution({"+": same, "-": dict(same)}, {"+": bl, "-": bh})
    assert df.attrs["strand_pct_difference"] == pytest.approx(0.0)
    assert df.contribution.sum() == pytest.approx(0.0)


# --------------------------------------------------------- NR orientation ---
def test_nr_pattern_puts_the_rnmp_second_on_the_light_strand(chrm):
    p = CMP.nr_pattern(chrm, 5, "+")
    assert p[1] == chrm[4]          # rNMP is the base AT the position
    assert p[0] == chrm[3]          # neighbour is 5', i.e. pos-1


def test_nr_pattern_uses_the_other_side_on_the_heavy_strand(chrm):
    """The heavy strand runs antiparallel, so its 5' neighbour is the
    complement of pos+1, not pos-1."""
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    p = CMP.nr_pattern(chrm, 5, "-")
    assert p[1] == comp[chrm[4]]
    assert p[0] == comp[chrm[5]]


def test_nr_pattern_wraps_circularly(chrm):
    assert CMP.nr_pattern(chrm, 1, "+") == chrm[-1] + chrm[0]


def test_table_s5_column_keys_confirm_the_rnmp_is_second():
    """Table S5's own ordering (AA, CA, GA, TA, AC, ...) varies the FIRST
    character fastest within a block sharing the second -- so the second
    character is the rNMP. This is the external evidence for the orientation."""
    import os
    if not os.path.exists(S5_NR):
        pytest.skip("Table S5 not exported")
    cols = [c for c in pd.read_csv(S5_NR, sep="\t").columns if c != "Celltype"]
    assert cols[:4] == ["AA", "CA", "GA", "TA"]
    assert [c[1] for c in cols] == list("AAAACCCCGGGGUUUU")


def test_nr_background_covers_both_strands(chrm):
    bg = CMP.nr_background(chrm)
    assert len(bg) == 16
    # two strands, circular, minus the 2-mers touching the N on each strand
    assert sum(bg.values()) == 2 * (len(chrm) - 2)


def test_heavy_strand_is_the_reverse_complement():
    assert CMP.heavy_strand("AAGG") == "CCTT"


# ------------------------------------------------- shipped-pilot readability --
def test_pilot_alignment_is_readable_gzipped():
    """The repo ships FS310's alignment gzipped so Part A of the notebook runs on
    a fresh clone. This test is the guard that a clone is actually sufficient:
    it reads the shipped file and checks it reproduces the published count."""
    import os
    import sys
    sys.path.insert(0, "src")
    from rnmp import pipeline as P
    sam = "results/ribosemap/results/FS310/alignment/aligned.sam.gz"
    fa = "data/reference/chrM.fa"
    if not (os.path.exists(sam) and os.path.exists(fa)):
        pytest.skip("pilot alignment or reference absent")
    coords = P.call_rnmps(sam, P.load_chrm(fa))
    assert len(coords) == 4866
    assert int((~coords.mismatch).sum()) == 4798      # Table S2 publishes 4,765


def test_pilot_composition_reads_gzipped_beds():
    import os
    import sys
    sys.path.insert(0, "src")
    from rnmp import pipeline as P
    d = "results/ribosemap/results/FS310/coordinate0"
    if not os.path.isdir(d):
        pytest.skip("pilot coordinates absent")
    comp = P.composition(d, "FS310")
    assert comp.sum() == pytest.approx(100.0)
    assert comp.idxmax() == "G"                        # Table S2: rGMP-dominant
