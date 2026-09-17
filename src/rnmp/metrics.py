"""rNMP-embedment frequency metrics.

Transcribed from Xu, Yang et al., Nucleic Acids Res. 2024;52(3):1207-1225,
Methods sections "Calculate the normalized frequencies for rNMP embedment" and
"Strand bias contribution of each type of rNMPs".

The two published formulas, verbatim from the paper:

    PPB_{G,L} = R_{G,L} / ( Length(G) * R_total,L )

    EF_{G,L}  = ( R_{G,L} / R_total,L ) * ( Length(Genome) / Length(G) )

where R_{G,L} is the rNMP count in genomic region G in library L and R_total,L
is the total rNMP count in mtDNA in library L.

Note the exact algebraic relationship, which `tests/test_metrics.py` asserts and
which is the cheapest way to catch a transcription error in either function:

    EF_{G,L} == PPB_{G,L} * Length(Genome)

So EF is PPB expressed in units of "relative to the genome-wide average", and a
bin with EF > 1 is a bin whose PPB exceeds the genome mean PPB. This is what
makes the paper's REZ threshold (EF > 1) a threshold on relative enrichment
rather than on an absolute rate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "ppb", "enrichment_factor", "circular_moving_average",
    "background_counts",
    "bin_counts", "bin_enrichment", "common_rez", "composition_normalized_frequency",
    "dinucleotide_normalized_frequency", "strand_bias_contribution",
]


# ----------------------------------------------------------------- core formulas
def ppb(region_count, region_length, total_count):
    """rNMP-embedment probability per base (PPB).

    Paper: "The PPB represents the rNMP-embedment probability at each nucleotide,
    which is comparable among the different genomic regions and different rNMP
    libraries."

        PPB = R_{G,L} / ( Length(G) * R_total,L )

    Parameters accept scalars or array-likes and broadcast elementwise.
    Returns NaN where `region_length` or `total_count` is zero rather than
    raising, so that empty libraries/regions propagate instead of aborting a
    whole-cohort run.
    """
    R = np.asarray(region_count, dtype=float)
    L = np.asarray(region_length, dtype=float)
    T = np.asarray(total_count, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = R / (L * T)
    return np.where((L > 0) & (T > 0), out, np.nan)


def enrichment_factor(region_count, region_length, total_count, genome_length):
    """rNMP enrichment factor (EF).

    Paper: "The enrichment factor is similar to PPB but measures the rNMP
    embedment frequency in a region."

        EF = ( R_{G,L} / R_total,L ) * ( Length(Genome) / Length(G) )

    `genome_length` is the length of the reference used for normalisation --
    16569 for GRCh38 chrM. Pass the *single-strand* length even when binning
    both strands: the paper divides both strands into 200-nt bins over the same
    16569-nt coordinate system (166 bins = 83 per strand x 2), so each bin's
    length is 200 (or 169 for the last, truncated bin) and the genome length
    stays 16569.
    """
    R = np.asarray(region_count, dtype=float)
    L = np.asarray(region_length, dtype=float)
    T = np.asarray(total_count, dtype=float)
    G = float(genome_length)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = (R / T) * (G / L)
    return np.where((L > 0) & (T > 0), out, np.nan)


def circular_moving_average(values, window):
    """Moving average over a circular coordinate system.

    Paper: "The moving average of PPB (window size = 51 nt) is used to measure
    the rNMP embedment frequency in the control region."

    The control region straddles position 1/16569, so a linear moving average
    would truncate exactly the region being measured. `window` must be odd so
    the average is centred on each position.
    """
    v = np.asarray(values, dtype=float)
    if window % 2 == 0:
        raise ValueError(f"window must be odd to stay centred; got {window}")
    if window > v.size:
        raise ValueError(f"window {window} exceeds series length {v.size}")
    half = window // 2
    padded = np.concatenate([v[-half:], v, v[:half]])
    kernel = np.ones(window) / window
    return np.convolve(padded, kernel, mode="valid")


# ------------------------------------------------------------------- REZ binning
def bin_counts(positions, strands, genome_length, bin_size=200):
    """Bin single-nucleotide rNMP positions into per-strand bins.

    Returns a DataFrame with one row per bin: `strand`, `bin_index`, `start`,
    `end` (0-based half-open), `length`, `count`.

    Paper: "Both strands of hmtDNA were divided into 200-nt bins (N = 166)."
    ceil(16569 / 200) = 83 bins per strand, 166 total; the last bin of each
    strand is 169 nt, not 200, so `length` is carried per bin and used by
    `bin_enrichment` rather than assuming `bin_size`.
    """
    positions = np.asarray(positions, dtype=np.int64)   # 1-based
    strands = np.asarray(strands)
    if positions.size and (positions.min() < 1 or positions.max() > genome_length):
        raise ValueError("positions must be 1-based within [1, genome_length]")

    n_bins = int(np.ceil(genome_length / bin_size))
    edges = np.arange(n_bins) * bin_size
    lengths = np.minimum(edges + bin_size, genome_length) - edges

    rows = []
    for strand in ("+", "-"):
        sel = positions[strands == strand]
        idx = (sel - 1) // bin_size
        counts = np.bincount(idx, minlength=n_bins)
        rows.append(pd.DataFrame({
            "strand": strand, "bin_index": np.arange(n_bins),
            "start": edges, "end": edges + lengths,
            "length": lengths, "count": counts,
        }))
    return pd.concat(rows, ignore_index=True)


def bin_enrichment(binned, genome_length, normalize_by="genome",
                   total_count=None, ef_threshold=1.0):
    """Add PPB / EF / `enriched` columns to the output of `bin_counts`.

    Paper: "For each rNMP library, the rNMP enrichment factors (EF) were
    calculated in each bin. All bins with EF > 1 represent enriched zones."

    RESOLVED: use the default, normalize_by="genome". R_total,L is the library's
    WHOLE-mtDNA rNMP count, both strands pooled -- the literal reading of the
    printed formula.

    The Methods were ambiguous on this point (they define R_total,L as "the total
    rNMP count in mtDNA in library L" but then bin *both strands* into 166 bins
    over one 16569-nt coordinate system, and those statements imply different
    EF=1 baselines). It was settled against the paper's own published numbers,
    not by argument -- see tests/test_published_values.py and
    docs/dataset_notes.md S3:

        At a single nucleotide Length(G) = 1, so EF = (R_pos / R_total,L) * 16569.
        Supplementary Table S4 prints that EF for 48 hotspot positions x 32
        libraries, and Table S2 prints each library's total rNMP count. Inverting
        every one of those 1536 published values with the pooled total recovers
        an exact integer rNMP count (max residual 2.8e-14, i.e. float noise). A
        per-strand R_total,L does not.

    Consequences to keep in mind when reading REZ output:

      normalize_by="genome"  (verified, default)
          A uniform library gives every bin EF = 0.5, not 1.0, because 166 bins
          of 200 nt cover 2 x 16569 nt of sequence. The paper's "EF > 1"
          threshold therefore means "> 2x the genome-average per-nucleotide
          rate" -- a strict criterion, and the reason only 13 of 166 bins
          qualify as common REZs.

      normalize_by="strand"  (retained for sensitivity analysis only)
          R_total,L = rNMPs on that bin's own strand, so a uniform library sits
          at EF = 1.0. Do NOT use this to reproduce the paper; it is here to
          quantify how much of the REZ set depends on the choice.

    `total_count` overrides both modes with an explicit R_total,L -- use it when
    a library's total must come from pre-filter counts.
    """
    if normalize_by not in ("genome", "strand"):
        raise ValueError("normalize_by must be 'genome' or 'strand'")
    out = binned.copy()

    if total_count is not None:
        T = pd.Series(float(total_count), index=out.index)
    elif normalize_by == "genome":
        T = pd.Series(float(out["count"].sum()), index=out.index)
    else:
        T = out.groupby("strand")["count"].transform("sum").astype(float)

    out["total_count"] = T
    out["ppb"] = ppb(out["count"], out["length"], T)
    out["ef"] = enrichment_factor(out["count"], out["length"], T, genome_length)
    out["enriched"] = out["ef"] > ef_threshold
    return out


def common_rez(per_library_bins, library_categories, ef_threshold=1.0,
               library_fraction=0.80, category_rule="any"):
    """Identify the paper's "common REZs" across a cohort of libraries.

    Paper: "The common REZs are the bins that are enriched in all cell categories
    and at least in 80% of the libraries." Figure 1C caption restates this as
    "REZs in all cell categories and genotypes, and more than 80% of the
    libraries studied".

    Both conditions must hold, and they are not redundant: the first protects
    against a zone driven by one over-represented cell type (HEK293T contributes
    11 of the 32 libraries), the second against a zone that is marginal
    everywhere.

    `category_rule` resolves what "enriched in a cell category" means, which the
    paper does not spell out:

        "any"  (default, the paper's convention) -- a category counts if AT LEAST
               ONE of its libraries has EF > 1.
        "all"  -- a category counts only if EVERY one of its libraries does.

    "any" is the paper's rule, and this is forced rather than chosen: under "all",
    a common REZ would need every library enriched, i.e. 100% of libraries, which
    would make the separate ">80% of libraries" criterion vacuous. Supplementary
    Table S3 settles it directly -- its common REZs include zones enriched in only
    26 of 32 libraries (81.25%), so not every library can be required. Running
    `category_rule="all"` against this cohort yields zero common REZs.

    Parameters
    ----------
    per_library_bins : dict[str, pandas.DataFrame]
        library id -> output of `bin_enrichment` for that library.
    library_categories : dict[str, str]
        library id -> cell category (the paper's 10 categories).

    Returns a DataFrame indexed by (strand, bin_index) with `n_libraries`,
    `frac_libraries`, `n_categories_enriched` and `is_common_rez`.
    """
    if category_rule not in ("any", "all"):
        raise ValueError("category_rule must be 'any' or 'all'")
    frames = []
    for lib, df in per_library_bins.items():
        d = df[["strand", "bin_index", "ef"]].copy()
        d["library"] = lib
        d["category"] = library_categories[lib]
        d["enriched"] = d["ef"] > ef_threshold
        frames.append(d)
    if not frames:
        raise ValueError("no libraries supplied")
    long = pd.concat(frames, ignore_index=True)

    n_lib = long["library"].nunique()
    categories = set(long["category"])

    g = long.groupby(["strand", "bin_index"])
    res = pd.DataFrame({
        "n_libraries": g["enriched"].sum(),
        "frac_libraries": g["enriched"].mean(),
    })
    # what "enriched in a cell category" means -- see the docstring
    agg = "any" if category_rule == "any" else "all"
    by_cat = long.groupby(["strand", "bin_index", "category"])["enriched"].agg(agg)
    res["n_categories_enriched"] = by_cat.groupby(level=[0, 1]).sum()
    res["is_common_rez"] = (
        (res["n_categories_enriched"] == len(categories))
        & (res["frac_libraries"] >= library_fraction)
    )
    res.attrs["n_libraries"] = n_lib
    res.attrs["n_categories"] = len(categories)
    return res


# ------------------------------------------------- composition / pattern frequency
def composition_normalized_frequency(rnmp_counts, background_counts):
    """Normalised rNMP composition (paper Fig 5 left panels).

    Paper: "the count of rNMPs ... are divided by the corresponding background
    dNMP count in the reference genome to get the frequency. For composition,
    the frequencies are further divided by the sum of frequencies to get the
    normalized frequencies."

    `rnmp_counts` / `background_counts` are dicts keyed by 'A','C','G','U'
    (rNMP) and 'A','C','G','T' (reference dNMP). Keys are matched with U == T.
    Returns a Series summing to 1.0.
    """
    key = {"A": "A", "C": "C", "G": "G", "U": "T", "T": "T"}
    freq = {}
    for base, count in rnmp_counts.items():
        bg = background_counts[key[base]]
        freq[base] = count / bg if bg else np.nan
    s = pd.Series(freq, dtype=float)
    return s / s.sum()


def dinucleotide_normalized_frequency(dinuc_counts, background_dinuc_counts,
                                      rnmp_position="second"):
    """Normalised dinucleotide pattern frequency (paper Fig 5 heatmaps).

    Paper: "For dinucleotide and trinucleotide patterns, the frequencies are
    divided by the sum of frequencies that share the same rNMP to get the
    normalized frequencies."

    So each *rNMP identity* is normalised independently and the four
    neighbour-base values for a given rNMP sum to 1.0 -- which is why the
    paper's null expectation for a dinucleotide is 0.25 (and 0.0625 for
    trinucleotides, where 16 neighbour contexts share one rNMP).

    `dinuc_counts` is keyed by the two-character pattern as written in the
    paper's figures. `rnmp_position` says which character of that key is the
    rNMP. Getting it wrong transposes the heatmap without changing its row sums,
    so it will not trip an assertion -- the values below are taken from the
    paper's own sheet layouts, not guessed:

        rnmp_position="second"  -> the paper's "NR" patterns (5' neighbour).
            Supplementary Table S5 sheet `human_nr` is keyed
            AA, CA, GA, TA, AC, CC, GC, TC, ... so the SECOND character is the
            rNMP and the first is the 5'-adjacent reference base.

        rnmp_position="first"   -> the paper's "RN" patterns (3' neighbour).
            Sheet `human_rn` is keyed AA, AC, AG, AT, CA, CC, CG, CT, ...
            so the FIRST character is the rNMP.

    Trinucleotide ("NNR") patterns use the same per-rNMP normalisation over 16
    neighbour contexts, so their null expectation is 1/16 = 0.0625, as stated in
    the Methods. Supplementary Table S5's footnote says 0.125; that is a typo --
    sheet `human_nnr` carries all 64 trinucleotides, i.e. 16 per rNMP. See
    config/params.yaml.

    NOTE: the paper computes these with RibosePreferenceAnalysis (ref. 33).
    Cross-check one library against that package before publishing numbers.
    """
    if rnmp_position not in ("first", "second"):
        raise ValueError("rnmp_position must be 'first' or 'second'")
    i = 0 if rnmp_position == "first" else 1

    freq = {}
    for pat, count in dinuc_counts.items():
        bg = background_dinuc_counts[pat]
        freq[pat] = count / bg if bg else np.nan
    s = pd.Series(freq, dtype=float)
    groups = pd.Series([p[i] for p in s.index], index=s.index)
    return s / s.groupby(groups).transform("sum")


# ------------------------------------------------------------ strand-bias decomposition
def strand_bias_contribution(rnmp_counts_by_strand, background_counts_by_strand):
    """Per-rNMP-type contribution to the observed strand bias (paper Fig 1 panels).

    Implements Methods, "Strand bias contribution of each type of rNMPs", in the
    four steps the paper describes:

      1. embedment rate per strand = rNMP count / background dNMP count
      2. strand-incorporation difference = rate(preferred strand) - rate(other),
         where "preferred" is the strand carrying more rNMPs overall in this library
      3. contribution % = each positive difference / sum of positive differences
      4. final contribution = contribution % * (difference in rNMP embedment
         percentage between the two strands), normalising for bias magnitude

    Inputs are nested dicts: {'+': {'A': n, ...}, '-': {...}}, where '+' is the
    light strand and '-' the heavy strand (see config/params.yaml: `strands`).

    Returns a DataFrame indexed by rNMP with columns `rate_light`, `rate_heavy`,
    `difference`, `contribution_pct`, `contribution`. Types whose difference is
    <= 0 (i.e. that oppose the library's overall bias) get contribution 0, per
    "the contribution percentage is calculated using all the positive differences".
    """
    key = {"A": "A", "C": "C", "G": "G", "U": "T", "T": "T"}
    light, heavy = rnmp_counts_by_strand["+"], rnmp_counts_by_strand["-"]
    bg_l, bg_h = background_counts_by_strand["+"], background_counts_by_strand["-"]

    tot_l, tot_h = sum(light.values()), sum(heavy.values())
    total = tot_l + tot_h
    if total == 0:
        raise ValueError("library contains no rNMPs")
    preferred = "+" if tot_l >= tot_h else "-"

    # step 4's scale factor: difference in rNMP embedment percentage between strands
    pct_diff = abs(tot_l - tot_h) / total

    df = pd.DataFrame(index=sorted(light))
    df["rate_light"] = [light[b] / bg_l[key[b]] if bg_l[key[b]] else np.nan for b in df.index]
    df["rate_heavy"] = [heavy[b] / bg_h[key[b]] if bg_h[key[b]] else np.nan for b in df.index]
    df["difference"] = (df["rate_light"] - df["rate_heavy"]) * (1 if preferred == "+" else -1)

    positive = df["difference"].clip(lower=0)
    denom = positive.sum()
    df["contribution_pct"] = positive / denom if denom else np.nan
    df["contribution"] = df["contribution_pct"] * pct_diff
    df.attrs["preferred_strand"] = preferred
    df.attrs["strand_pct_difference"] = pct_diff
    return df


# --------------------------------------------------------------------------- #
# Strand-aware background base counts
# --------------------------------------------------------------------------- #

def background_counts(chrm: str, strand: str) -> dict:
    """Background dNMP counts for one strand of a double-stranded reference.

    The light strand ('+') is the reference sequence as given; the heavy strand
    ('-') is its complement, so a heavy-strand rNMP of type X sits opposite a
    reference base complement(X). Counting the complement rather than
    reverse-complementing is equivalent here (composition is order-free) and
    avoids an unnecessary string reversal.
    """
    comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
    if strand == "+":
        return {b: chrm.count(b) for b in "ACGT"}
    if strand == "-":
        return {b: chrm.count(comp[b]) for b in "ACGT"}
    raise ValueError(f"strand must be '+' or '-', got {strand!r}")
