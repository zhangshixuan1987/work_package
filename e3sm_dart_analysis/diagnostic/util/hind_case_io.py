"""Run plumbing shared by the hind_lead_eval case notebooks (05, 06, 07): products cached on disk with
the settings they were built from, and figure saving.

A product is read from its file when the file exists, carries the same settings stamp and no recompute
is forced; otherwise it is built and written (to a temporary file first, so an interrupted run never
leaves a truncated product behind).
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Callable

import pandas as pd
import xarray as xr


def _replace(tmp: Path, path: Path) -> None:
    os.replace(tmp, path)
    print(f"[WRITE] {path}")


def cached_dataset(path: Path, build: Callable[[], xr.Dataset], stamp: str, *, force: bool = False) -> xr.Dataset:
    """``build()`` cached in the NetCDF file ``path`` (stamp kept in the ``stamp`` attribute)."""
    path = Path(path)
    if path.exists() and not force:
        ds = xr.load_dataset(path)
        if ds.attrs.get("stamp") == stamp:
            print(f"[READ] {path}")
            return ds
        print(f"[STALE] {path}: built with other settings, recomputing")
    ds = build()
    ds.attrs["stamp"] = stamp
    tmp = path.with_name(path.name + ".tmp")
    ds.to_netcdf(tmp)
    _replace(tmp, path)
    return ds


def cached_table(path: Path, build: Callable[[], pd.DataFrame], stamp: str, *, force: bool = False) -> pd.DataFrame:
    """``build()`` (a table with a default index) cached in the CSV file ``path``; the stamp is kept
    next to it in ``<name>.stamp``. The table is returned as read back from the file, so a fresh and a
    cached run give the same dtypes."""
    path = Path(path)
    stamp_file = path.with_name(path.name + ".stamp")
    if path.exists() and stamp_file.exists() and not force:
        if stamp_file.read_text() == stamp:
            print(f"[READ] {path}")
            return pd.read_csv(path)
        print(f"[STALE] {path}: built with other settings, recomputing")
    write_table(build(), path)
    stamp_file.write_text(stamp)
    return pd.read_csv(path)


def write_dataset(ds: xr.Dataset, path: Path) -> None:
    """Write ``ds`` (a figure's data) to ``path`` through a temporary file."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    ds.to_netcdf(tmp)
    _replace(tmp, path)


def write_table(table: pd.DataFrame, path: Path, *, index: bool = False) -> None:
    """Write ``table`` to the CSV file ``path`` through a temporary file."""
    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    table.to_csv(tmp, index=index)
    _replace(tmp, path)


def figure_saver(figure_dir: Path, *, dpi: int = 200, show: bool = True) -> Callable:
    """``save(fig, name)``: write ``figure_dir/name``; close the figure unless ``show``."""
    def save(fig, name: str) -> Path:
        import matplotlib.pyplot as plt

        out = Path(figure_dir) / name
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        print(f"[FIGURE] {out}")
        if not show:
            plt.close(fig)
        return out
    return save


__all__ = ["cached_dataset", "cached_table", "write_dataset", "write_table", "figure_saver"]
