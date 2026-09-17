"""Unit tests for the transcribed PPB / EF formulas.

These are deliberately hand-computable: the point is to catch a formula
transcription error, not to test numpy. Run with `pytest -q` in rnmp-analysis.
"""
import numpy as np
import pytest

from rnmp.metrics import (
    ppb, enrichment_factor, circular_moving_average,
    bin_counts, bin_enrichment, common_rez, composition_normalized_frequency,
    dinucleotide_normalized_frequency, strand_bias_contribution,
)

CHRM = 16569


def test_ppb_hand_computed():
    # 50 rNMPs in a 1000-nt region, 2000 rNMPs total -> 50 / (1000 * 2000)
    assert ppb(50, 1000, 2000) == pytest.approx(2.5e-5)


def test_ef_hand_computed():
    # (50 / 2000) * (16569 / 1000)
    assert enrichment_factor(50, 1000, 2000, CHRM) == pytest.approx(0.025 * 16.569)


def test_ef_is_ppb_times_genome_length():
    """The algebraic identity between the paper's two formulas."""
    rng = np.random.default_rng(0)
    r = rng.integers(0, 500, 200)
    ln = rng.integers(50, 3000, 200)
    tot = 12345
    np.testing.assert_allclose(
        enrichment_factor(r, ln, tot, CHRM), ppb(r, ln, tot) * CHRM, rtol=1e-12)


def test_ef_equals_one_for_uniform_region():
    """A region whose rNMP share equals its length share has EF exactly 1."""
    assert enrichment_factor(200, 200, CHRM, CHRM) == pytest.approx(1.0)


def test_degenerate_inputs_are_nan_not_raise():
    assert np.isnan(ppb(5, 0, 100))
    assert np.isnan(enrichment_factor(5, 100, 0, CHRM))


def test_bin_counts_matches_paper_bin_number():
    b = bin_counts([1, 16569], ["+", "-"], CHRM, bin_size=200)
    assert len(b) == 166, "paper: 'divided into 200-nt bins (N = 166)'"
    assert (b.groupby("strand").size() == 83).all()
    # last bin is truncated, not 200 nt
    assert b[b.bin_index == 82].length.unique().tolist() == [16569 - 82 * 200]
    assert b.groupby("strand").length.sum().unique().tolist() == [CHRM]


def test_bin_counts_assigns_edges_correctly():
    b = bin_counts([1, 200, 201, 16569], ["+"] * 4, CHRM)
    plus = b[b.strand == "+"].set_index("bin_index")["count"]
    assert plus[0] == 2 and plus[1] == 1 and plus[82] == 1


def _uniform_bins(n=20000, seed=1):
    rng = np.random.default_rng(seed)
    return bin_counts(rng.integers(1, CHRM + 1, n),
                      rng.choice(["+", "-"], n), CHRM)


def test_bin_enrichment_genome_normalisation_baseline_is_half():
    """Literal R_total,L (both strands pooled) puts a uniform library at EF=0.5.

    166 bins of 200 nt cover 2 x 16569 nt, so pooling the strand totals halves
    the baseline. This is the consequence flagged in bin_enrichment's docstring;
    the test pins it so the behaviour cannot drift silently.
    """
    e = bin_enrichment(_uniform_bins(), CHRM, normalize_by="genome")
    assert np.average(e.ef, weights=e.length) == pytest.approx(0.5)
    assert e.ef.median() == pytest.approx(0.5, abs=0.05)


def test_bin_enrichment_strand_normalisation_baseline_is_one():
    e = bin_enrichment(_uniform_bins(), CHRM, normalize_by="strand")
    for strand, grp in e.groupby("strand"):
        assert np.average(grp.ef, weights=grp.length) == pytest.approx(1.0)
    assert e.ef.median() == pytest.approx(1.0, abs=0.1)


def test_bin_enrichment_modes_differ_under_strand_bias():
    """The two normalisations must disagree on a strand-biased library --
    otherwise the flag would be cosmetic and the ambiguity would not matter."""
    pos = np.arange(1, CHRM + 1)
    strand = np.where(np.arange(CHRM) % 10 < 8, "+", "-")   # 80/20 light bias
    b = bin_counts(pos, strand, CHRM)
    g = bin_enrichment(b, CHRM, normalize_by="genome")
    s = bin_enrichment(b, CHRM, normalize_by="strand")
    assert g.enriched.sum() != s.enriched.sum()


def test_bin_enrichment_explicit_total_overrides():
    b = _uniform_bins()
    e = bin_enrichment(b, CHRM, total_count=1)
    assert e.total_count.unique().tolist() == [1.0]


def test_common_rez_any_rule_requires_every_category_and_80pct():
    """Bin 0 is enriched in every library; bin 1 misses one category entirely;
    bin 2 is enriched in 3/5 libraries only."""
    libs, cats = {}, {}
    plan = {                      # library -> (category, bins enriched)
        "L1": ("catA", {0, 1, 2}), "L2": ("catA", {0, 1, 2}),
        "L3": ("catB", {0, 1, 2}), "L4": ("catB", {0}),
        "L5": ("catC", {0}),
    }
    for lib, (cat, hot) in plan.items():
        b = bin_counts([], [], CHRM)
        b["ef"] = [2.0 if i in hot else 0.1 for i in b.bin_index]
        libs[lib], cats[lib] = b, cat
    res = common_rez(libs, cats, category_rule="any")
    plus = res.xs("+", level="strand")
    assert bool(plus.loc[0, "is_common_rez"]) is True
    assert bool(plus.loc[1, "is_common_rez"]) is False   # catC never enriched there
    assert bool(plus.loc[2, "is_common_rez"]) is False   # only 3/5 libraries
    assert plus.loc[2, "frac_libraries"] == pytest.approx(0.6)


def test_common_rez_all_rule_is_stricter():
    """Under category_rule='all', a bin enriched in 4/5 libraries fails even when
    every category has at least one enriched library.

    This is why "any" must be the paper's rule: requiring every library of every
    category would force frac_libraries == 1.0, making the separate >80%
    criterion vacuous, yet Table S3 lists REZs at 26/32 libraries.
    """
    libs, cats = {}, {}
    plan = {"L1": ("catA", {0}), "L2": ("catA", {0}), "L3": ("catB", {0}),
            "L4": ("catB", set()), "L5": ("catC", {0})}
    for lib, (cat, hot) in plan.items():
        b = bin_counts([], [], CHRM)
        b["ef"] = [2.0 if i in hot else 0.1 for i in b.bin_index]
        libs[lib], cats[lib] = b, cat
    any_ = common_rez(libs, cats, category_rule="any").xs("+", level="strand")
    all_ = common_rez(libs, cats, category_rule="all").xs("+", level="strand")
    assert any_.loc[0, "frac_libraries"] == pytest.approx(0.8)
    assert bool(any_.loc[0, "is_common_rez"]) is True
    assert bool(all_.loc[0, "is_common_rez"]) is False


def test_common_rez_rejects_bad_category_rule():
    b = bin_counts([], [], CHRM); b["ef"] = 2.0
    with pytest.raises(ValueError):
        common_rez({"L1": b}, {"L1": "catA"}, category_rule="most")


def test_composition_sums_to_one():
    s = composition_normalized_frequency(
        {"A": 100, "C": 400, "G": 50, "U": 150},
        {"A": 5124, "C": 5181, "G": 2169, "T": 4094})   # GRCh38 chrM-like background
    assert s.sum() == pytest.approx(1.0)
    assert s.idxmax() == "C", "paper reports a strong rCMP bias"


def test_dinucleotide_rows_sum_to_one_per_rnmp():
    counts = {f"{n}{r}": 10 + i for i, (n, r) in
              enumerate((n, r) for r in "ACGU" for n in "ACGT")}
    bg = {k: 100 for k in counts}
    s = dinucleotide_normalized_frequency(counts, bg, rnmp_position="second")
    for r in "ACGU":
        assert sum(v for k, v in s.items() if k[1] == r) == pytest.approx(1.0)


def test_strand_bias_contribution_shape_and_normalisation():
    df = strand_bias_contribution(
        {"+": {"A": 100, "C": 500, "G": 60, "U": 140},
         "-": {"A": 80,  "C": 120, "G": 90, "U": 110}},
        {"+": {"A": 5124, "C": 5181, "G": 2169, "T": 4094},
         "-": {"A": 4094, "C": 2169, "G": 5181, "T": 5124}})
    assert df.attrs["preferred_strand"] == "+"
    assert df.contribution_pct.sum() == pytest.approx(1.0)
    assert df.contribution.sum() == pytest.approx(df.attrs["strand_pct_difference"])
    assert (df.contribution >= 0).all()
