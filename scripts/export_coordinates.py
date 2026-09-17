#!/usr/bin/env python
"""Export mismatch-filtered rNMP counts per (library, position, strand).

    python scripts/export_coordinates.py

Writes results/coordinates/per_position_counts.csv.gz -- the input to every
downstream analysis (REZ, hotspots, per-nucleotide maps). Collapsing to counts
keeps it small (<= 2 x 16569 rows per library) while losing nothing those
analyses need.

Vectorised rather than line-by-line: the cohort is ~23M aligned reads.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnmp import pipeline as P                                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/coordinates/per_position_counts.csv.gz"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def aligned_len(cigar: str) -> int:
    if cigar.endswith("M") and cigar[:-1].isdigit():          # fast path
        return int(cigar[:-1])
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X")


def counts_for(lib: str, chrm_arr: np.ndarray) -> pd.DataFrame:
    sam = ROOT / f"results/ribosemap/results/{lib}/alignment/aligned.sam"
    df = pd.read_csv(sam, sep="\t", comment="@", header=None, usecols=[1, 3, 5, 9],
                     names=["flag", "pos", "cigar", "seq"], dtype={"cigar": str, "seq": str},
                     quoting=3, on_bad_lines="skip")
    df = df[(df.flag & 4) == 0]
    rev = (df.flag & 16).to_numpy() > 0
    alen = np.array([aligned_len(c) for c in df.cigar.to_numpy()])
    gpos = np.where(rev, df.pos.to_numpy() + alen - 1, df.pos.to_numpy())
    seq = df.seq.to_numpy()
    read_base = np.where(rev, [s[-1] for s in seq], [s[0] for s in seq])
    ok = (gpos >= 1) & (gpos <= chrm_arr.size)
    gpos, read_base, rev = gpos[ok], read_base[ok], rev[ok]
    keep = read_base == chrm_arr[gpos - 1]                     # reference-mismatch filter
    gpos, rev = gpos[keep], rev[keep]
    strand = np.where(rev, "+", "-")                           # rNMP strand is opposite the read
    out = (pd.DataFrame({"pos": gpos, "strand": strand})
             .value_counts().rename("count").reset_index())
    out.insert(0, "library", lib)
    return out


def main():
    chrm_arr = np.array(list(P.load_chrm(ROOT / "data/reference/chrM.fa")))
    samples = pd.read_csv(ROOT / "config/samples_xu2024.tsv", sep="\t", keep_default_na=False)
    libs = samples.loc[samples.role == "ribose_seq", "library"].tolist()
    frames = []
    for i, lib in enumerate(libs, 1):
        d = counts_for(lib, chrm_arr)
        frames.append(d)
        print(f"[{i:2d}/{len(libs)}] {lib}: {int(d['count'].sum()):>9,} rNMPs "
              f"at {len(d):>6,} (pos,strand) sites", flush=True)
    all_ = pd.concat(frames, ignore_index=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    all_.to_csv(OUT, index=False)
    print(f"\n{OUT}: {len(all_):,} rows, {int(all_['count'].sum()):,} rNMPs total")


if __name__ == "__main__":
    main()
