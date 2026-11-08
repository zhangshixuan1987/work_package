"""North-polar map helpers for the polar-vortex notebooks (10 and 11).

The maps are north-polar stereographic (North America at the bottom), cut to a circle at
``lat_min``. Shading is an anomaly with fixed levels; contours are a total field with fixed
absolute levels, so a change between panels or frames is a change of the field, never of the
plot settings. ``fill_map`` redraws only the field artists of a panel, so the same axes serve
every animation frame.
"""

from __future__ import annotations

import subprocess

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
from cartopy.util import add_cyclic_point
from matplotlib import animation
from matplotlib.colors import BoundaryNorm

from esp_lab.utils import colormap_utils as mycolors

PROJ = ccrs.NorthPolarStereo(central_longitude=-90)   # North America at the bottom, Eurasia at the top
_THETA = np.linspace(0, 2 * np.pi, 120)
CIRCLE = mpath.Path(np.vstack([np.sin(_THETA), np.cos(_THETA)]).T * 0.5 + [0.5, 0.5])


def shading(style):
    """(levels, cmap, norm) of an anomaly style {"range", "ci", "cmap"}; fixed for every panel and frame."""
    levels = np.arange(-style["range"], style["range"] + style["ci"] / 2, style["ci"])
    cmap = mycolors.blue2red_cmap(levels) if style.get("cmap", "blue2red") == "blue2red" else plt.get_cmap(style["cmap"], len(levels) - 1)
    return levels, cmap, BoundaryNorm(levels, ncolors=cmap.N)


def setup_map(ax, title, *, lat_min, color="k", title_size=10.5, badge_size=9.5, label_size=10.5, land="0.88",
              coast_lw=0.5, dashed_frame=False):
    """Polar map axes cut to a circle at ``lat_min``; returns the panel state used by ``fill_map``.

    The extent is set in projected coordinates: ``set_extent([-180, 180, ...])`` collapses the y range
    of a polar-stereographic map in cartopy 0.26.
    """
    x, y = PROJ.transform_point(0.0, lat_min, ccrs.PlateCarree())
    r = float(np.hypot(x, y))
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_boundary(CIRCLE, transform=ax.transAxes)
    ax.add_feature(cfeature.LAND, facecolor=land, zorder=0.5)
    ax.coastlines(lw=coast_lw, color="0.35", zorder=3)
    ax.gridlines(crs=ccrs.PlateCarree(), lw=0.3, color="0.5", alpha=0.5, xlocs=range(-180, 181, 45), ylocs=[40, 60, 80])
    if dashed_frame:
        ax.spines["geo"].set_linestyle("--")
        ax.spines["geo"].set_edgecolor("0.4")
    ax.set_title(title, loc="left", fontsize=title_size, color=color, fontweight="bold")
    return {"artists": [],
            "lead": ax.set_title("", loc="right", fontsize=title_size, color=color),
            "empty": ax.text(0.5, 0.5, "", transform=ax.transAxes, ha="center", va="center", fontsize=label_size,
                             color="0.35", zorder=6, bbox={"facecolor": "white", "edgecolor": "0.6", "boxstyle": "round,pad=0.4"}),
            "badge": ax.text(0.5, -0.04, "", transform=ax.transAxes, ha="center", va="top", fontsize=badge_size)}


def fill_map(ax, state, anom=None, total=None, *, style=None, total_levels=None, empty="", badge="", badge_color="k",
             lead="", total_color="k", total_lw=0.8):
    """Shade ``anom`` (lat, lon) with ``style`` and contour ``total`` at ``total_levels``; a grey panel with
    ``empty`` when both are None."""
    for artist in state["artists"]:
        artist.remove()
    state["artists"] = []
    state["badge"].set_text(badge)
    state["badge"].set_color(badge_color)
    state["lead"].set_text(lead)
    if anom is None and total is None:
        ax.set_facecolor("0.93")
        state["empty"].set_text(empty)
        state["empty"].set_visible(True)
        return
    ax.set_facecolor("white")
    state["empty"].set_visible(False)
    if anom is not None and style is not None:
        levels, cmap, norm = shading(style)
        values, lon = add_cyclic_point(anom.values, coord=anom["lon"].values)
        state["artists"].append(ax.contourf(lon, anom["lat"].values, values, levels=levels, cmap=cmap, norm=norm,
                                            extend="both", transform=ccrs.PlateCarree(), zorder=1))
    if total is not None:
        values, lon = add_cyclic_point(total.values, coord=total["lon"].values)
        state["artists"].append(ax.contour(lon, total["lat"].values, values, levels=total_levels, colors=total_color,
                                           linewidths=total_lw, transform=ccrs.PlateCarree(), zorder=2))


def u60_badge(value, prefix="U60"):
    """Keyword arguments (badge, badge_color) of ``fill_map`` for a U(60N, 10 hPa) value; red when easterly."""
    if value is None or not np.isfinite(value):
        return {"badge": "", "badge_color": "k"}
    return {"badge": f"{prefix}(10 hPa) {value:+.0f} m/s" + ("  easterly" if value < 0 else ""),
            "badge_color": "tab:red" if value < 0 else "k"}


def video_writer(fps):
    """(extension, writer factory) for the best video codec the available ffmpeg has, or (None, None):
    H.264 MP4 with libx264 (e.g. from imageio-ffmpeg), else VP9 WebM."""
    try:
        import imageio_ffmpeg   # bundled ffmpeg with H.264
        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        pass
    if not animation.writers.is_available("ffmpeg"):
        return None, None
    encoders = subprocess.run([plt.rcParams["animation.ffmpeg_path"], "-hide_banner", "-encoders"],
                              capture_output=True, text=True).stdout
    even = ["-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p"]
    if "libx264" in encoders:
        return "mp4", lambda: animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=even + ["-crf", "20"])
    if "libvpx-vp9" in encoders:
        return "webm", lambda: animation.FFMpegWriter(fps=fps, codec="libvpx-vp9", extra_args=even + ["-crf", "30", "-b:v", "0"])
    return None, None
