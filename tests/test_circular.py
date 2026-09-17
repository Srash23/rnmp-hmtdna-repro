"""Tests for circular recovery of junction-spanning rNMPs (scripts/04_circularize.py).

These pin the coordinate arithmetic and the junction test. Both are easy to get
subtly wrong (off-by-one on the rotation, or a junction test that admits reads
which never cross it), and either error would silently change Figure 1C.
"""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("circ", ROOT / "scripts/04_circularize.py")
circ = importlib.util.module_from_spec(spec)
spec.loader.exec_module(circ)

CHRM_LEN = 16569


def test_rotation_offset_is_half_the_genome():
    assert circ.OFFSET == 8284
    assert circ.CHRM_LEN == CHRM_LEN


def test_junction_sits_between_the_two_flagged_rotated_positions():
    """canonical 16569 is at rotated 8285, canonical 1 at rotated 8286."""
    assert circ.to_canonical(circ.J_LO) == CHRM_LEN
    assert circ.to_canonical(circ.J_HI) == 1


def test_coordinate_map_is_a_bijection():
    seen = {circ.to_canonical(i) for i in range(1, CHRM_LEN + 1)}
    assert seen == set(range(1, CHRM_LEN + 1))


def test_rotated_reference_matches_the_declared_offset():
    fa = ROOT / "data/reference/chrM_rot8284.fa"
    ref = ROOT / "data/reference/chrM.fa"
    if not (fa.exists() and ref.exists()):
        pytest.skip("reference absent; run scripts/00_setup_reference.sh")
    rd = lambda p: "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()
    chrm, rot = rd(ref), rd(fa)
    assert len(rot) == CHRM_LEN
    assert rot == chrm[circ.OFFSET:] + chrm[:circ.OFFSET]
    # the map must recover the same base from either coordinate system
    for r in (1, 4242, circ.J_LO, circ.J_HI, CHRM_LEN):
        assert rot[r - 1] == chrm[circ.to_canonical(r) - 1]


@pytest.mark.parametrize("pos,alen,crosses", [
    (8285, 1, False),      # ends exactly at canonical 16569
    (8286, 1, False),      # starts exactly at canonical 1
    (8285, 2, True),       # spans 16569 -> 1
    (8200, 100, True),     # comfortably across
    (8187, 100, True),     # latest start a 100-nt read can have and still cross
    (8186, 100, False),    # ends at 8285 -- reaches the junction but does not cross it
    (8186, 101, True),     # one base longer, now crosses
    (8187, 99, False),     # also ends at 8285
    (1, 100, False),       # nowhere near
    (16000, 150, False),
])
def test_junction_span_test(pos, alen, crosses):
    assert (pos <= circ.J_LO and pos + alen - 1 >= circ.J_HI) is crosses


def test_recovered_positions_are_within_one_read_length_of_the_junction():
    """Structural guarantee: a junction-crossing read can only place an rNMP
    within its own length of the junction, so nothing far away may change."""
    import pandas as pd
    p = ROOT / "results/coordinates/per_position_counts_circular.csv.gz"
    q = ROOT / "results/coordinates/per_position_counts_refiltered.csv.gz"
    if not (p.exists() and q.exists()):
        pytest.skip("coordinate files absent; run the pipeline first")
    a = pd.read_csv(p); b = pd.read_csv(q)
    m = a.merge(b, on=["library", "pos", "strand"], how="left", suffixes=("", "_pre"))
    chg = m[m.count_pre.isna() | (m["count"] != m.count_pre)]
    d = (chg.pos - 1).abs().combine((CHRM_LEN - (chg.pos - 1).abs()), min)
    assert d.max() <= 300, f"a count changed {d.max()} nt from the junction"


def test_missing_hotspots_are_not_reachable_by_circular_recovery():
    """(+,7018) and (-,14923) sit >1,600 nt from the junction, so this step
    cannot recover them -- documented in docs/dataset_notes.md S13 and S16."""
    for pos in (7018, 14923):
        d = min(abs(pos - 1), CHRM_LEN - abs(pos - 1))
        assert d > 300, f"position {pos} is only {d} nt from the junction"
