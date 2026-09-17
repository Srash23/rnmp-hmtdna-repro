#!/usr/bin/env python
"""Append one library's rNMP totals, strand split and composition to the running
per-library table. Idempotent: re-running for an existing library replaces its row,
so an interrupted sweep can simply be restarted.

    python scripts/summarize_library.py FS310
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rnmp import pipeline as P                                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/tables/per_library_rnmp.csv"
COLS = ["library", "run", "sample", "cell_category", "genotype", "fragmentation",
        "n_aligned", "n_mismatch", "n_rnmp", "n_rnmp_published", "delta_pct",
        "light", "heavy", "light_pct", "heavy_pct",
        "A_pct", "C_pct", "G_pct", "U_pct",
        "A_pct_pub", "C_pct_pub", "G_pct_pub", "U_pct_pub", "comp_max_dev_pp"]


def summarize(lib: str) -> dict:
    samples = pd.read_csv(ROOT / "config/samples_xu2024.tsv", sep="\t", keep_default_na=False)
    row = samples[samples.library == lib].iloc[0]
    # cell_category in the sheet is derived from Table S2, NOT from the SRA
    # SampleName: SRA mislabels FS303 as HEK293T-WT when the paper (Tables S2
    # and S4, and the N values in Fig 1A) assigns it to hESC-H9.
    s2 = pd.read_csv(ROOT / "data/external/table_S2_libraries.tsv", sep="\t").set_index("library").loc[lib]

    chrm = P.load_chrm(ROOT / "data/reference/chrM.fa")
    coords = P.call_rnmps(ROOT / f"results/ribosemap/results/{lib}/alignment/aligned.sam", chrm)
    kept = coords[~coords.mismatch]
    strand = P.strand_summary(kept, "after mismatch filter")
    comp = P.composition(ROOT / f"results/ribosemap/results/{lib}/coordinate0", lib)
    pub = {b: float(s2[f"r{b}MP_pct"]) for b in "ACGU"}
    n_pub = int(s2["n_rnmp"])

    return dict(
        library=lib, run=row.run, sample=row["sample"], cell_category=row.cell_category,
        genotype=s2["genotype"].strip(), fragmentation=row.fragmentation,
        n_aligned=len(coords), n_mismatch=int(coords.mismatch.sum()), n_rnmp=len(kept),
        n_rnmp_published=n_pub, delta_pct=100 * (len(kept) - n_pub) / n_pub,
        light=strand["light"], heavy=strand["heavy"],
        light_pct=strand["light_pct"], heavy_pct=strand["heavy_pct"],
        **{f"{b}_pct": float(comp[b]) for b in "ACGU"},
        **{f"{b}_pct_pub": pub[b] for b in "ACGU"},
        comp_max_dev_pp=max(abs(float(comp[b]) - pub[b]) for b in "ACGU"),
    )


if __name__ == "__main__":
    rec = summarize(sys.argv[1])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(OUT) if OUT.exists() else pd.DataFrame(columns=COLS)
    df = df[df.library != rec["library"]]
    df = pd.concat([df, pd.DataFrame([rec])], ignore_index=True)[COLS]
    df.sort_values(["cell_category", "library"]).to_csv(OUT, index=False)
    print(f"  {rec['library']}: {rec['n_rnmp']:,} rNMPs "
          f"({rec['delta_pct']:+.1f}% vs published) | "
          f"light {rec['light_pct']:.1f}% | comp dev {rec['comp_max_dev_pp']:.2f} pp")
