"""Publication figure style for the RTENSO2026 7a-7d notebooks.

Figures are drawn at their printed size (journal column widths in inches) with 7-9 pt text, so
they need no rescaling in a manuscript. ``save`` writes a vector PDF (fonts embedded as TrueType)
and a 600-dpi PNG through ``rtenso_case.save_figure``.
"""

from __future__ import annotations

import string

import matplotlib as mpl
import matplotlib.pyplot as plt

from . import rtenso_case as rc

WIDTH = {"single": 3.5, "onehalf": 5.5, "double": 7.2}   # inches (AGU/AMS column widths)
FONT = 8
SAVE_EXTS = ("pdf", "png")
SAVE_DPI = 600


def use_style(font=FONT):
    mpl.rcParams.update({
        "font.size": font, "axes.titlesize": font, "axes.labelsize": font,
        "xtick.labelsize": font - 1, "ytick.labelsize": font - 1, "legend.fontsize": font - 1,
        "figure.titlesize": font + 1, "axes.titlelocation": "left", "axes.titlepad": 3,
        "axes.linewidth": 0.6, "lines.linewidth": 1.2, "lines.markersize": 3.5, "patch.linewidth": 0.5,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 3, "ytick.major.size": 3, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.linewidth": 0.4, "grid.alpha": 0.4,
        "legend.frameon": False, "legend.handlelength": 1.8, "legend.columnspacing": 1.2, "legend.borderaxespad": 0.3,
        "figure.dpi": 150, "savefig.dpi": SAVE_DPI, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular",
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.h_pad": 0.03, "figure.constrained_layout.w_pad": 0.03,
    })


def panel_label(ax, k, text=None, **kwargs):
    """Bold (a), (b), ... in the top-left corner above the axes, optionally followed by ``text``."""
    label = f"({string.ascii_lowercase[k]})"
    ax.set_title(label if text is None else f"{label} {text}", loc="left", **kwargs)


def zero_line(ax, axis="y", **kwargs):
    style = {"color": "0.5", "lw": 0.6, "zorder": 0} | kwargs
    (ax.axhline if axis == "y" else ax.axvline)(0.0, **style)


def legend_below(fig, handles=None, labels=None, ncol=4, ax=None, **kwargs):
    """Figure legend centred below the axes."""
    if handles is None:
        handles, labels = (ax or fig.axes[0]).get_legend_handles_labels()
    return fig.legend(handles, labels, loc="outside lower center", ncol=ncol, **kwargs)


def save(fig, notebook, *parts, exts=SAVE_EXTS, dpi=SAVE_DPI):
    """Save as fig_<parts>.<ext> under OUTPUT_ROOT/<notebook>/figures/ and print the paths."""
    for ext in exts:
        print(rc.save_figure(fig, notebook, *parts, dpi=dpi, ext=ext))


__all__ = ["WIDTH", "FONT", "use_style", "panel_label", "zero_line", "legend_below", "save"]
