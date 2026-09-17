"""Figure style for this project's plots.

Kept in the repo rather than relying on an interactive helper, so the notebook
runs standalone from a fresh clone. Values are the ones used for every figure in
figures/: a light typographic grid, open frames (no top/right spines), and a
single neutral grey for secondary annotation.
"""
from __future__ import annotations

import matplotlib as mpl

GREY = "#5A6169"          # secondary text: captions, reference-line labels
LIGHT_STRAND = "#C0392B"  # paper's Fig 1A encoding
HEAVY_STRAND = "#2471A3"
BASE_COLOURS = {"A": "#2E9B57", "C": "#2471A3", "G": "#E8A33D", "T": "#C0392B", "U": "#C0392B"}

__all__ = ["apply_figure_style", "set_frame", "GREY",
           "LIGHT_STRAND", "HEAVY_STRAND", "BASE_COLOURS"]


def apply_figure_style(sizes: tuple[float, float, float] = (9, 8, 7)) -> None:
    """Set rcParams. `sizes` is (title, axis-label, tick-label) in points."""
    title, label, tick = sizes
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": label,
        "axes.titlesize": title,
        "axes.labelsize": label,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "legend.fontsize": tick,
        "axes.titlelocation": "left",
        "axes.titlepad": 8.0,
        "axes.labelcolor": "#1A1D20",
        "text.color": "#1A1D20",
        "axes.edgecolor": "#33383D",
        "axes.linewidth": 0.8,
        "xtick.color": "#33383D",
        "ytick.color": "#33383D",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def set_frame(ax, style: str = "open") -> None:
    """`open` keeps left+bottom spines; `none` removes all four."""
    keep = {"open": {"left", "bottom"}, "none": set()}
    if style not in keep:
        raise ValueError("style must be 'open' or 'none'")
    for name, sp in ax.spines.items():
        sp.set_visible(name in keep[style])
