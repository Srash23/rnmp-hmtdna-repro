#!/usr/bin/env python
"""Apply the authors' restriction-enzyme background subtraction to every library.

    python scripts/filter_re_sites.py

Reads each library's Ribose-Map alignment, calls rNMPs, applies the
reference-mismatch filter and then the RE-site subtraction implemented in
src/rnmp/re_sites.py (a GPL-v3 port of the authors' own code -- see that module).

Writes:
  results/coordinates/per_position_counts_refiltered.csv.gz
  results/tables/re_filter_audit.csv        per-library before/after vs published

The `noda` class is removed positionally. The `da` class needs the read sequence,
so those rows are resolved individually -- there are few of them, and the bulk of
each library is handled vectorised.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rnmp import pipeline as P                                    # noqa: E402
from rnmp import re_sites as R                                    # noqa: E402

DA_OFFSET = 0        # deposited reads are already barcode/UMI-stripped; authors use 11
_RC = str.maketrans("ACGTN", "TGCAN")
CIG = __import__("re").compile(r"(\d+)([MIDNSHP=X])")


def rnmps_from_sam(sam, chrm):
    """Per aligned read: 1-based rNMP position, strand, read base, original read seq."""
    rows = []
    with open(sam) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 10)
            flag, pos, cigar, seq = int(f[1]), int(f[3]), f[5], f[9]
            if flag & 4:
                continue
            rev = bool(flag & 16)
            alen = sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X")
            if rev:
                g, rb, st = pos + alen - 1, seq[-1], "+"
                original = seq.translate(_RC)[::-1]
            else:
                g, rb, st = pos, seq[0], "-"
                original = seq
            rows.append((g, st, rb, original))
    d = pd.DataFrame(rows, columns=["pos", "strand", "read_base", "read"])
    d["ref_base"] = chrm_take(chrm, d.pos.to_numpy())
    return d


def chrm_take(chrm, pos1):
    arr = np.frombuffer(chrm.encode(), dtype="S1")
    return np.char.decode(arr[pos1 - 1])


def main():
    chrm = P.load_chrm(ROOT / "data/reference/chrM.fa")
    table = R.load_enzymes(ROOT / "data/external/res_all.list")
    by_canon = {R.canonical(k): k for k in table}
    samples = pd.read_csv(ROOT / "config/samples_xu2024.tsv", sep="\t", keep_default_na=False)
    s2 = pd.read_csv(ROOT / "data/external/table_S2_libraries.tsv", sep="\t").set_index("library")
    libs = samples.loc[samples.role == "ribose_seq", "library"].tolist()

    frames, audit = [], []
    for i, lib in enumerate(sorted(libs), 1):
        enz = [by_canon[R.canonical(e)] for e in str(s2.enzymes[lib]).split(";")
               if e and R.canonical(e) in by_canon]
        d = rnmps_from_sam(ROOT / f"results/ribosemap/results/{lib}/alignment/aligned.sam", chrm)
        n_aligned = len(d)
        d = d[d.read_base == d.ref_base]                       # reference-mismatch filter
        n_mm = n_aligned - len(d)
        d["bed"] = d.pos - 1

        n_noda = n_da_cand = n_da_conf = 0
        if enz:
            noda, da = R.build_masks(chrm, enz, table)
            pats = R.da_patterns(enz, table)
            noda_all = set().union(*noda.values())
            key = list(zip(d.bed, d.strand))
            is_noda = np.fromiter((k in noda_all for k in key), bool, len(d))
            n_noda = int(is_noda.sum())
            d = d[~is_noda]
            # dA class: positional candidate, then confirmed on the read itself
            key = list(zip(d.bed, d.strand))
            cand = np.fromiter((any(k in da[e] for e in enz) for k in key), bool, len(d))
            n_da_cand = int(cand.sum())
            if n_da_cand:
                sub = d[cand]
                conf = []
                for k, rd in zip(zip(sub.bed, sub.strand), sub.read):
                    hit = any(k in da[e] and pats[e].match(rd[DA_OFFSET:]) for e in enz)
                    conf.append(hit)
                conf = np.array(conf, bool)
                n_da_conf = int(conf.sum())
                drop = sub.index[conf]
                d = d.drop(index=drop)

        g = (d.groupby(["pos", "strand"]).size().rename("count").reset_index())
        g.insert(0, "library", lib)
        frames.append(g)

        pub = int(s2.n_rnmp[lib])
        before = n_aligned - n_mm
        audit.append(dict(library=lib, enzymes=";".join(enz), n_enzymes=len(enz),
                          n_aligned=n_aligned, n_mismatch=n_mm,
                          n_before=before, n_noda_removed=n_noda,
                          n_da_candidates=n_da_cand, n_da_removed=n_da_conf,
                          n_after=int(g["count"].sum()), n_published=pub,
                          delta_pct_before=100 * (before - pub) / pub,
                          delta_pct_after=100 * (int(g["count"].sum()) - pub) / pub))
        a = audit[-1]
        print(f"[{i:2d}/32] {lib}  before {a['n_before']:>8,} "
              f"-noda {n_noda:>6,} -dA {n_da_conf:>5,} -> {a['n_after']:>8,} "
              f"| pub {pub:>8,}  {a['delta_pct_before']:+6.2f}% -> {a['delta_pct_after']:+6.2f}%",
              flush=True)

    out = pd.concat(frames, ignore_index=True)
    (ROOT / "results/coordinates").mkdir(parents=True, exist_ok=True)
    out.to_csv(ROOT / "results/coordinates/per_position_counts_refiltered.csv.gz", index=False)
    ad = pd.DataFrame(audit)
    ad.round(4).to_csv(ROOT / "results/tables/re_filter_audit.csv", index=False)
    tot, pub = int(out["count"].sum()), int(ad.n_published.sum())
    print(f"\ncohort: {tot:,} rNMPs vs {pub:,} published ({100*(tot-pub)/pub:+.2f}%)")
    print(f"was {ad.n_before.sum():,} ({100*(ad.n_before.sum()-pub)/pub:+.2f}%)")
    print(f"libraries within 3% of own published: {int((ad.delta_pct_after.abs()<3).sum())}/32 "
          f"(was {int((ad.delta_pct_before.abs()<3).sum())}/32)")


if __name__ == "__main__":
    main()
