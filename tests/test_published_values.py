"""Regression tests against the paper's own published numbers.

These are the tests that make this a *reproduction* rather than a reimplementation:
each one checks our code against a value Xu, Yang et al. actually printed, using
the supplementary tables extracted to data/external/. They require no sequencing
data and run in under a second, so they gate every later pipeline change.

Source: gkad1204_supplemental_files.zip, obtained from NLM's PMC Open Access S3
bucket (s3://pmc-oa-opendata/PMC10853789.1/), Tables S1-S4 and S7.
"""
import numpy as np
import pandas as pd
import pytest

from rnmp.metrics import enrichment_factor, ppb, bin_counts, composition_normalized_frequency

CHRM = 16569
EXT = "data/external"


@pytest.fixture(scope="module")
def s2():
    return pd.read_csv(f"{EXT}/table_S2_libraries.tsv", sep="\t").set_index("library")


@pytest.fixture(scope="module")
def s3():
    return pd.read_csv(f"{EXT}/table_S3_common_rez.tsv", sep="\t")


@pytest.fixture(scope="module")
def s4():
    return pd.read_csv(f"{EXT}/table_S4_common_hotspots_ef.tsv", sep="\t")


# --------------------------------------------------------------- EF normalisation
def test_published_single_nucleotide_ef_inverts_to_integer_counts(s2, s4):
    """THE decisive test for the EF normalisation ambiguity.

    At a single nucleotide Length(G) = 1, so the paper's formula reduces to
        EF = (R_pos / R_total,L) * 16569
    Inverting it must give back an integer rNMP count. It does -- for all 48
    published hotspot positions x 32 libraries -- when R_total,L is the library's
    WHOLE-mtDNA rNMP count from Table S2 (both strands pooled).

    This is what settles `bin_enrichment(normalize_by=...)`: the literal
    "genome" reading is the paper's convention. A per-strand R_total,L would
    make these residuals nonzero.
    """
    libs = [c for c in s4.columns if c.startswith("FS")]
    assert len(libs) == 32

    residuals, n = [], 0
    for lib in libs:
        total = int(s2.loc[lib, "n_rnmp"])
        ef = pd.to_numeric(s4[lib], errors="coerce")
        implied = ef * total / CHRM
        ok = implied.notna() & (ef > 0)
        residuals.append(np.abs(implied[ok] - implied[ok].round()))
        n += int(ok.sum())

    resid = pd.concat(residuals)
    assert n == 1536, f"expected 48 positions x 32 libraries, got {n}"
    assert resid.max() < 1e-6, f"max residual {resid.max():.3e} -- pooled R_total,L is wrong"


def test_our_ef_reproduces_published_values(s2, s4):
    """Round-trip: recover the published EF from the integer counts it implies."""
    libs = [c for c in s4.columns if c.startswith("FS")]
    for lib in libs[:6]:
        total = int(s2.loc[lib, "n_rnmp"])
        published = pd.to_numeric(s4[lib], errors="coerce")
        counts = (published * total / CHRM).round()
        ours = enrichment_factor(counts, 1, total, CHRM)
        np.testing.assert_allclose(ours, published.to_numpy(), rtol=1e-9)


def test_ef_and_ppb_agree_on_published_data(s2, s4):
    lib = "FS185"
    total = int(s2.loc[lib, "n_rnmp"])
    counts = (pd.to_numeric(s4[lib]) * total / CHRM).round()
    np.testing.assert_allclose(
        enrichment_factor(counts, 1, total, CHRM), ppb(counts, 1, total) * CHRM, rtol=1e-12)


# ------------------------------------------------------------------ bin geometry
def test_our_bin_geometry_matches_published_rez_coordinates(s3):
    """Table S3's REZ coordinates must fall exactly on our bin edges.

    This independently confirms 0-based edges at multiples of 200 and a
    truncated final bin -- including REZ8L/REZ5H at 16400-16569 (169 nt, index 82).
    """
    bins = bin_counts([], [], CHRM, bin_size=200).set_index(["strand", "bin_index"])
    for _, r in s3.iterrows():
        b = bins.loc[(r["strand"], r["bin_index"])]
        assert b["start"] == r["start"], r["name"]
        assert b["end"] == r["end"], r["name"]
        assert b["length"] == r["length"], r["name"]


def test_published_rez_counts_and_thresholds(s3):
    assert (s3.strand == "+").sum() == 8, "paper: eight common REZs on the light strand"
    assert (s3.strand == "-").sum() == 5, "paper: five common REZs on the heavy strand"
    assert (s3.ef_mean > 1.0).all(), "every common REZ must clear the EF > 1 threshold"
    assert (s3.frac_libraries >= 0.80).all(), "paper: enriched in at least 80% of libraries"
    assert s3.length.isin([200, 169]).all()


# ----------------------------------------------------------------- composition
def test_published_composition_percentages_are_consistent(s2):
    pct = s2[["rAMP_pct", "rCMP_pct", "rGMP_pct", "rUMP_pct"]]
    np.testing.assert_allclose(pct.sum(axis=1).to_numpy(), 100.0, atol=1e-6)


def test_composition_helper_reproduces_a_published_library(s2):
    """hESC-H9 FS197 is rCMP-dominant in Table S2; our helper must agree on rank
    order once the chrM dNMP background is divided out.

    Uses the published percentages as rNMP counts scaled by the library total,
    so this checks the normalisation arithmetic rather than the background model.
    """
    lib = "FS197"
    total = int(s2.loc[lib, "n_rnmp"])
    counts = {b: s2.loc[lib, f"r{b if b != 'U' else 'U'}MP_pct"] / 100 * total for b in "ACGU"}
    # uniform background isolates the normalisation step
    s = composition_normalized_frequency(counts, {"A": 1, "C": 1, "G": 1, "T": 1})
    assert s.sum() == pytest.approx(1.0)
    assert s.idxmax() == "C"
    assert s["C"] == pytest.approx(s2.loc[lib, "rCMP_pct"] / 100, abs=1e-9)


# ----------------------------------------------------------------- adapters
def test_adapter_tags_and_library_barcodes_are_consistent():
    """Every Table S2 barcode must be an Adapter.L tag or its reverse complement.

    Table S2 is not internally consistent about which form it reports: the liver
    libraries (FS257-FS264) give the reverse complement while the rest give the
    tag as synthesised. That is why scripts/03_ribosemap.sh verifies the barcode
    against the FASTQ instead of trusting this column -- see docs/dataset_notes.md.
    """
    s1 = pd.read_csv(f"{EXT}/table_S1_adapters.tsv", sep="\t")
    samples = pd.read_csv(f"{EXT}/table_S2_libraries.tsv", sep="\t")
    allowed = set(s1.tag) | set(s1.barcode_revcomp)
    assert len(s1) == 8
    unknown = set(samples.barcode) - allowed
    assert not unknown, f"barcodes absent from Table S1: {unknown}"
    assert (samples.barcode.str.len() == 3).all()


def test_max_counts_cover_every_library(s2):
    s7 = pd.read_csv(f"{EXT}/table_S7_max_counts.tsv", sep="\t")
    assert len(s7) == 32
    assert set(s7.library) == set(s2.index)
    assert (s7[["max_count_light", "max_count_heavy"]] > 0).all().all()


# --------------------------------------------------------------------------
# Figure 2: conventions confirmed against Table S4's own annotation columns.
# These need only the reference and the supplementary table, no sequencing data.
# --------------------------------------------------------------------------
CHRM_FA = "data/reference/chrM.fa"
_COMP = str.maketrans("ACGT", "TGCA")


def _chrm():
    import os
    if not os.path.exists(CHRM_FA):
        pytest.skip("chrM.fa absent; run scripts/00_setup_reference.sh")
    return "".join(l.strip() for l in open(CHRM_FA) if not l.startswith(">")).upper()


def test_rnmp_base_matches_table_S4():
    """Our strand convention must reproduce Table S4's `rnmp` column exactly.

    A '+' (light-strand) rNMP is the reference base; a '-' rNMP is its
    complement. Getting this backwards would flip 48/48 of these.
    """
    chrm = _chrm()
    s4 = pd.read_csv(f"{EXT}/table_S4_common_hotspots_ef.tsv", sep="\t")
    ours = [chrm[p-1] if st == "+" else chrm[p-1].translate(_COMP)
            for p, st in zip(s4.pos, s4.strand)]
    assert list(ours) == list(s4.rnmp)


def test_7nt_context_matches_table_S4():
    """And Table S4's `Pattern(NNNRNNN)` column, which additionally pins orientation.

    A reverse-complement error leaves the central base right but reverses the
    7-mer, so this catches what the base test alone cannot.
    """
    chrm = _chrm()
    n = len(chrm)
    s4 = pd.read_csv(f"{EXT}/table_S4_common_hotspots_ef.tsv", sep="\t")
    ours = []
    for p, st in zip(s4.pos, s4.strand):
        seg = "".join(chrm[i % n] for i in range(p-4, p+3))
        ours.append(seg if st == "+" else seg.translate(_COMP)[::-1])
    assert ours == list(s4.context_7nt)


def test_single_nucleotide_ef_is_invertible_to_integer_counts(s2, s4):
    """The pooled-convention proof, restated as a test over all 1,536 values.

    At one nucleotide Length(G) = 1, so EF = (R_pos / R_total,L) * 16569. With
    R_total,L pooled across strands, every published EF must invert to an integer.
    """
    resid = []
    for lib in s2.index:
        implied = s4[lib].astype(float) * s2.n_rnmp[lib] / CHRM
        resid.extend((implied - implied.round()).abs().tolist())
    assert len(resid) == 48 * 32
    assert max(resid) < 1e-6, f"max residual {max(resid):.3g} -- pooled convention violated"
