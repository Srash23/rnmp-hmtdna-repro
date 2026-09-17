#!/usr/bin/env python
"""Circular recovery of rNMPs on reads spanning the chrM 1/16569 junction.

    python scripts/04_circularize.py

Ribose-Map treats chrM as a linear contig, so a read crossing position 1 cannot
align end-to-end and is lost -- together with its rNMP. The paper's Methods
handle this by aligning a second time against a rotated copy of chrM and merging.

WHY THIS DOES NOT RE-FETCH ANY FASTQ
------------------------------------
Bowtie 2 runs in end-to-end mode here, so a junction-spanning read fails to align
and is written to `aligned.sam` with flag 4 rather than discarded. Every read this
step needs is therefore already on disk in the existing alignments. We extract the
unmapped records, align just those against `chrM_rot8284`, and merge back.

WHAT IS KEPT
------------
Only alignments that actually cross the junction. An unmapped read can fail
linear alignment for unrelated reasons (adapter remnant, contamination, too many
mismatches), and such a read aligning somewhere in the middle of the rotated
reference would be a *new* rNMP call that the linear pass rejected -- not a
recovered one. Admitting those would inflate counts rather than repair a known
loss, so the junction-spanning test is the whole point of the filter:

    rotated 1-based interval [p, p+alen-1] crosses the junction
    iff  p <= 8285 <= 8286 <= p+alen-1

because canonical 16569 sits at rotated 8285 and canonical 1 at rotated 8286.

Recovered rNMPs then pass the same two filters as the main path: the
reference-mismatch filter, and the restriction-site subtraction of
`src/rnmp/re_sites.py`.

Writes:
  results/coordinates/per_position_counts_circular.csv.gz   merged counts
  results/tables/circular_recovery_audit.csv                per-library audit
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from rnmp import pipeline as P                                    # noqa: E402
from rnmp import re_sites as R                                    # noqa: E402

OFFSET, CHRM_LEN = 8284, 16569
J_LO, J_HI = CHRM_LEN - OFFSET, CHRM_LEN - OFFSET + 1             # 8285, 8286
SEED, MISMATCHES = "123456789", "0"
_RC = str.maketrans("ACGTN", "TGCAN")
CIG = __import__("re").compile(r"(\d+)([MIDNSHP=X])")


def to_canonical(rot_pos: int) -> int:
    """Rotated 1-based position -> canonical 1-based position."""
    return ((rot_pos - 1 + OFFSET) % CHRM_LEN) + 1


def extract_unmapped(sam: Path, fq: Path) -> int:
    n = 0
    with open(sam) as fh, open(fq, "w") as out:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.rstrip("\n").split("\t")
            if not (int(f[1]) & 4):
                continue
            out.write(f"@{f[0]}\n{f[9]}\n+\n{f[10]}\n")
            n += 1
    return n


def recover(rot_sam: Path, chrm: str):
    """Junction-spanning alignments -> canonical rNMP rows."""
    rows = []
    n_aligned = n_spanning = 0
    for line in open(rot_sam):
        if line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t")
        flag, pos, cigar, seq = int(f[1]), int(f[3]), f[5], f[9]
        if flag & 4:
            continue
        n_aligned += 1
        alen = sum(int(n) for n, op in CIG.findall(cigar) if op in "MDN=X")
        if not (pos <= J_LO and pos + alen - 1 >= J_HI):
            continue                                   # not a junction read
        n_spanning += 1
        if flag & 16:
            rot_g, read_base, strand = pos + alen - 1, seq[-1], "+"
            orig = seq.translate(_RC)[::-1]
        else:
            rot_g, read_base, strand = pos, seq[0], "-"
            orig = seq
        g = to_canonical(rot_g)
        rows.append((g, strand, read_base, chrm[g - 1], orig))
    df = pd.DataFrame(rows, columns=["pos", "strand", "read_base", "ref_base", "read"])
    return df, n_aligned, n_spanning


def main():
    chrm = P.load_chrm(ROOT / "data/reference/chrM.fa")
    rot_idx = ROOT / "data/reference/chrM_rot8284"
    if not (rot_idx.with_suffix(".1.bt2")).exists():
        sys.exit("ERROR: rotated index missing; run scripts/00_setup_reference.sh")

    tbl = R.load_enzymes(ROOT / "data/external/res_all.list")
    by_canon = {R.canonical(k): k for k in tbl}
    s2 = pd.read_csv(ROOT / "data/external/table_S2_libraries.tsv", sep="\t").set_index("library")
    base = pd.read_csv(ROOT / "results/coordinates/per_position_counts_refiltered.csv.gz")
    libs = sorted(base.library.unique())

    tmp = ROOT / "results/coordinates/.circ"
    tmp.mkdir(parents=True, exist_ok=True)
    audit, extra = [], []

    for i, lib in enumerate(libs, 1):
        sam = ROOT / f"results/ribosemap/results/{lib}/alignment/aligned.sam"
        fq, rsam = tmp / f"{lib}.unmapped.fq", tmp / f"{lib}.rot.sam"
        n_un = extract_unmapped(sam, fq)
        if n_un:
            subprocess.run(["bowtie2", "--threads", "4", "--seed", SEED, "-N", MISMATCHES,
                            "-x", str(rot_idx), "-U", str(fq), "-S", str(rsam)],
                           check=True, capture_output=True)
            df, n_al, n_span = recover(rsam, chrm)
        else:
            df, n_al, n_span = pd.DataFrame(columns=["pos","strand","read_base","ref_base","read"]), 0, 0

        n_mm = int(df.read_base.ne(df.ref_base).sum()) if len(df) else 0
        df = df[df.read_base == df.ref_base] if len(df) else df

        n_re = 0
        if len(df):
            enz = [by_canon[R.canonical(e)] for e in str(s2.enzymes.get(lib, "")).split(";")
                   if e and R.canonical(e) in by_canon]
            if enz:
                noda, da = R.build_masks(chrm, enz, tbl)
                pats = R.da_patterns(enz, tbl)
                noda_all = set().union(*noda.values())
                keys = list(zip(df.pos - 1, df.strand))
                keep = [k not in noda_all for k in keys]
                df = df[keep]
                keys = list(zip(df.pos - 1, df.strand))
                conf = [any(k in da[e] and pats[e].match(rd) for e in enz)
                        for k, rd in zip(keys, df.read)]
                n_re = int(sum(conf)) + (len(keep) - sum(keep))
                df = df[[not c for c in conf]]

        if len(df):
            g = df.groupby(["pos", "strand"]).size().rename("count").reset_index()
            g["library"] = lib
            extra.append(g[["library", "pos", "strand", "count"]])
        audit.append(dict(library=lib, n_unmapped=n_un, n_realigned=n_al,
                          n_junction_spanning=n_span, n_mismatch_removed=n_mm,
                          n_re_removed=n_re, n_recovered=int(len(df))))
        print(f"[{i:2d}/{len(libs)}] {lib}  unmapped {n_un:>6,} -> realigned {n_al:>6,} "
              f"-> spanning {n_span:>4,} -> recovered {len(df):>4,}", flush=True)

    ad = pd.DataFrame(audit)
    ad.to_csv(ROOT / "results/tables/circular_recovery_audit.csv", index=False)
    merged = pd.concat([base] + extra, ignore_index=True) if extra else base
    merged = merged.groupby(["library", "pos", "strand"], as_index=False)["count"].sum()
    merged.to_csv(ROOT / "results/coordinates/per_position_counts_circular.csv.gz", index=False)

    print(f"\nunmapped reads across cohort : {ad.n_unmapped.sum():,}")
    print(f"junction-spanning alignments  : {ad.n_junction_spanning.sum():,}")
    print(f"rNMPs recovered               : {ad.n_recovered.sum():,}")
    print(f"cohort total {int(base['count'].sum()):,} -> {int(merged['count'].sum()):,}")


if __name__ == "__main__":
    main()
