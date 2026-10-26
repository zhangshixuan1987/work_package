"""Plots of the NSIDCv6 spin-up sea-ice files against the JRA55-FOSIRL archive.

``05_reanalysis_seaice_check.ipynb``: the sea-ice files in
``initial_conditions/Reanalysis_ICs/YYYYMMDD/`` (provided by D. Comeau) are
compared with the BruteForce (JRA55-FOSIRL) restart of the same date where it
exists (2015; ``archive_states``, reusing the key-field cache of notebooks
01-03), through hemispheric totals (``hemispheric_totals``, ``plot_totals``) and native-mesh
``mosaic`` maps (``plot_mpassi_rows``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import netCDF4 as nc
import numpy as np
import pandas as pd

import clim_restart as cr
import ic_eval as ice

QUANTITIES = ("ice concentration", "ice thickness", "snow on ice")
EXTENT_THRESHOLD = 0.15      # concentration above which a cell counts toward the ice extent
HEMISPHERES = {"NH": [-180, 180, 45, 90], "SH": [-180, 180, -90, -45]}

# Fixed source order and colors (categorical slots 1-2; markers as second encoding).
SOURCES = ("NSIDCv6 spin-up", "JRA55-FOSIRL")
SOURCE_STYLE = {
    "NSIDCv6 spin-up": {"color": "#2a78d6", "marker": "o"},
    "JRA55-FOSIRL": {"color": "#1baf7a", "marker": "^"},
}


def ice_file(root: Path, pattern: str, date: str) -> Path:
    """Sea-ice file of ``date`` (YYYY-MM-DD); ``pattern`` takes the YYYYMMDD tag as ``{tag}``."""
    return Path(root) / pattern.format(tag=date.replace("-", ""))


def read_mesh(path: Path) -> dict[str, np.ndarray]:
    """Cell latitude (degrees) and area (m2)."""
    with nc.Dataset(path, "r") as ds:
        return {"lat": np.degrees(cr._read(ds, "latCell")), "area": cr._read(ds, "areaCell").astype(np.float64)}


def read_state(path: Path) -> dict[str, np.ndarray]:
    """``{quantity: values[nCells]}``: category sums of concentration, ice and snow thickness."""
    fields = cr.read_key_fields("mpassi", Path(path))
    return {q: fields[f"NH {q}"][0] for q in QUANTITIES}    # NH and SH entries hold the same values


def archive_states(src_root: Path, prefix: str, mmdd: str, years: Sequence[int],
                   cache_dir: Path | None = None, log=print) -> dict[int, dict[str, np.ndarray]]:
    """Archive states ``{year: {quantity: values}}`` on ``mmdd``.

    Taken from ``<cache_dir>/mpassi_<MM-DD>.npz`` (written by notebooks 01-03)
    where it holds the year, otherwise read from the restarts.
    """
    out: dict[int, dict[str, np.ndarray]] = {}
    cache = Path(cache_dir) / f"mpassi_{mmdd}.npz" if cache_dir else None
    if cache is not None and cache.is_file():
        with np.load(cache) as z:
            names = list(z["names"])
            cached = [int(y) for y in z["years"]]
            arrays = {q: z[f"f{names.index(f'NH {q}')}"] for q in QUANTITIES}
            for y in years:
                if y in cached:
                    out[y] = {q: arrays[q][cached.index(y)].astype(np.float64) for q in QUANTITIES}
    for y in years:
        if y not in out:
            path = cr.source_file(src_root, prefix, "mpassi", y, mmdd)
            if path.is_file():
                out[y] = read_state(path)
                log(f"  archive {y}-{mmdd}: read {path.name}")
    return out


def hemispheric_totals(state: Mapping[str, np.ndarray], mesh: Mapping[str, np.ndarray]) -> dict[str, float]:
    """NH/SH ice area and extent (1e6 km2) and ice and snow volume (1e3 km3) of one state."""
    out = {}
    km2, km3 = 1e6 * 1e6, 1e9 * 1e3
    for hemi in HEMISPHERES:
        w = np.where(mesh["lat"] > 0 if hemi == "NH" else mesh["lat"] < 0, mesh["area"], 0.0)
        c = state["ice concentration"]
        out[f"{hemi} area"] = float(np.nansum(c * w)) / km2
        out[f"{hemi} extent"] = float(np.nansum(np.where(c >= EXTENT_THRESHOLD, w, 0.0))) / km2
        out[f"{hemi} ice volume"] = float(np.nansum(state["ice thickness"] * w)) / km3
        out[f"{hemi} snow volume"] = float(np.nansum(state["snow on ice"] * w)) / km3
    return out


def plot_totals(table: pd.DataFrame,
                quantities: Sequence[str] = ("area", "extent", "ice volume", "snow volume"),
                fontsize: float = ice.BASE_FONTSIZE, panel_size: tuple[float, float] | None = None):
    """Hemispheric totals against date: one row per hemisphere, one column per quantity.

    ``table`` has columns ``date``, ``source`` and the ``hemispheric_totals`` keys.
    """
    import matplotlib.pyplot as plt

    f = ice.font_sizes(fontsize)
    pw, ph = ice._panel(panel_size, (3.6, 2.8), fontsize)
    dates = list(dict.fromkeys(table["date"]))
    x = np.arange(len(dates))
    fig, axes = plt.subplots(2, len(quantities), figsize=(pw * len(quantities), ph * 2),
                             constrained_layout=True, squeeze=False)
    for r, hemi in enumerate(HEMISPHERES):
        for c, q in enumerate(quantities):
            ax, key = axes[r, c], f"{hemi} {q}"
            for src in SOURCES:
                sub = table[table["source"] == src].set_index("date").reindex(dates)
                if sub[key].isna().all():
                    continue
                st = SOURCE_STYLE[src]
                ax.plot(x, sub[key], color=st["color"], marker=st["marker"], markersize=6, linewidth=2, label=src)
            unit = "10$^6$ km$^2$" if q in ("area", "extent") else "10$^3$ km$^3$"
            ax.set_title(f"{hemi} {q} ({unit})", fontsize=f["title"])
            ax.set_xticks(x, [d[2:] for d in dates], rotation=45, ha="right", fontsize=f["tick"])
            ax.tick_params(axis="y", labelsize=f["tick"])
            ax.set_ylim(bottom=0)
            ax.grid(axis="y", color="0.9", linewidth=0.6)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels), fontsize=f["label"], frameon=False)
    fig.suptitle("MPAS-Seaice hemispheric totals by initialization date", fontsize=f["suptitle"],
                 fontweight="bold")
    return fig


def plot_mpassi_rows(rows: Sequence[tuple[str, Sequence[tuple | None]]], descriptors: Mapping[str, object],
                     title: str, fontsize: float = ice.BASE_FONTSIZE,
                     panel_size: tuple[float, float] | None = None):
    """Native-mesh polar maps: ``rows`` of ``(hemisphere, panels)``.

    Each panel is ``(title, values, cmap, lo, hi, colorbar label)`` or None (left
    blank); a symmetric ``lo, hi`` is drawn as a diverging scale. Columns share a
    scale down the rows, with one colorbar at the bottom of each column.
    """
    import cartopy.crs as ccrs
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    import mosaic

    ice.ensure_proj_data()
    theta = np.linspace(0, 2 * np.pi, 100)
    circle = mpath.Path(np.vstack([np.sin(theta), np.cos(theta)]).T * 0.5 + 0.5)
    f = ice.font_sizes(fontsize)
    pw, ph = ice._panel(panel_size, (3.6, 3.6), fontsize)
    ncol = max(len(p) for _, p in rows)
    fig = plt.figure(figsize=(pw * ncol, ph * len(rows) + 0.6), constrained_layout=True)
    colls: dict[int, tuple[object, str, bool]] = {}
    col_axes: dict[int, list] = {}      # every axes of a column, blank ones too, so its colorbar sits at the bottom
    for r, (hemi, panels) in enumerate(rows):
        desc = descriptors[hemi]
        for c, panel in enumerate(panels):
            ax = fig.add_subplot(len(rows), ncol, r * ncol + c + 1, projection=desc.projection)
            col_axes.setdefault(c, []).append(ax)
            if panel is None:
                ax.axis("off")
                continue
            ptitle, arr, cmap, lo, hi, label = panel
            ax.set_facecolor(ice.MISSING["mpassi"])
            coll = mosaic.polypcolor(ax, desc, np.ma.masked_invalid(arr), **ice.discrete_cmap(cmap, lo, hi, arr),
                                     antialiaseds=False)
            mosaic.coastlines(ax, desc, color="0.25", linewidth=0.5)
            ax.set_extent(HEMISPHERES[hemi], crs=ccrs.PlateCarree())
            ax.set_boundary(circle, transform=ax.transAxes)
            ax.set_title(ptitle, fontsize=f["title"])
            colls[c] = (coll, label, bool(np.isclose(lo, -hi)))
    for c, (coll, label, diverging) in colls.items():
        ice._cbar(fig, coll, fontsize, ax=col_axes[c], orientation="horizontal", fraction=0.045, pad=0.03,
                  aspect=24, shrink=0.9, skip_ticks=diverging, label=label)
    fig.suptitle(title, fontsize=f["suptitle"], fontweight="bold")
    return fig
