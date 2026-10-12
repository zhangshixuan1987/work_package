"""Time labels for post-processed EAM time series."""
from __future__ import annotations


def label_by_interval_start(ds, field):
    """Label interval means by their start and drop zero-width records.

    Forecast ts files stamp each mean at its interval end (the Jan 1 daily mean
    is stamped Jan 2 00Z) and begin with a zero-width snapshot of the initial
    state. Series without positive-width records (e.g. analysis-run daily
    files, which hold one 09Z snapshot per day) are returned unchanged.
    """
    name = ds.time.attrs.get("bounds") or next((n for n in ("time_bnds", "time_bounds") if n in ds), None)
    if name is None:
        return field
    bounds = ds[name].load()
    start = bounds.isel({bounds.dims[-1]: 0})
    end = bounds.isel({bounds.dims[-1]: 1})
    keep = (end > start).values
    if not keep.any():
        return field
    return field.isel(time=keep).assign_coords(time=start.values[keep])
