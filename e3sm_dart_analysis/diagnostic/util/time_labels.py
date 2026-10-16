"""Time labels for post-processed EAM/ELM time series."""
from __future__ import annotations

import numpy as np


def label_by_interval_start(ds, field):
    """Label interval means by their start and drop partial leading records.

    Forecast ts files stamp each mean at its interval end (the Jan 1 daily mean
    is stamped Jan 2 00Z), so matching the raw stamp with a reference compares
    each model day with the next reference day. EAM series begin with a
    zero-width snapshot of the initial state and ELM daily series with a
    one-step sliver (e.g. 23:30-00Z); records shorter than half the median
    record width are dropped. Series not stamped at their interval ends (e.g.
    observations stamped at the interval start or midpoint, whose floor is
    already the right day) are returned unchanged. Instantaneous fields (``time: point``, e.g. the
    00Z snapshots in the forecast daily files of PSL, TREFHT, T850) are already
    stamped at their instant and are returned unchanged, as are series without
    positive-width records (e.g. analysis-run daily files, which hold one 09Z
    snapshot per day).
    """
    if "time: point" in str(field.attrs.get("cell_methods", "")):
        return field
    name = ds.time.attrs.get("bounds") or next((n for n in ("time_bnds", "time_bounds") if n in ds), None)
    if name is None:
        return field
    bounds = ds[name].load()
    start = bounds.isel({bounds.dims[-1]: 0})
    end = bounds.isel({bounds.dims[-1]: 1})
    # timedelta64 for standard calendars, datetime.timedelta objects for cftime
    width = np.asarray((end - start).values).astype("timedelta64[ms]").astype(float)
    if not (width > 0).any():
        return field
    keep = width >= 0.5 * np.median(width[width > 0])
    if not (field.time.values[keep] == end.values[keep]).all():
        return field  # not stamped at interval end (e.g. GPCP, centred stamps)
    _check_daily_windows(field, start.values[keep], width[keep])
    return field.isel(time=keep).assign_coords(time=start.values[keep])


def _check_daily_windows(field, start, width_ms):
    """Refuse daily series whose means do not cover 00Z-00Z days.

    Hindcast ts/daily files built from 6-hourly output before the post_lib.sh
    daily_member fix (2026-10-12; e.g. LHFLX, SHFLX, FSNS) average 18Z-18Z and
    carry 03Z-09Z bounds; relabelling them would hide a 6-h offset.
    """
    if len(start) < 2:
        return
    spacing = np.median(np.diff(start.astype("datetime64[ms]").astype(float)))
    day_ms = 86400e3
    if abs(spacing - day_ms) > 3600e3:
        return  # not a daily series
    starts = start.astype("datetime64[ms]")
    misaligned = (np.abs(width_ms - day_ms) > 1e3) | (starts != starts.astype("datetime64[D]"))
    if misaligned.any():
        first = int(np.flatnonzero(misaligned)[0])
        raise ValueError(
            f"{field.name}: daily records do not cover 00Z-00Z days (record {first} starts "
            f"{starts[first]} and spans {width_ms[first] / 3600e3:g} h); this daily file was "
            "built from 6-hourly output with the old post_lib.sh daily_member. Use ts/6hourly "
            "or regenerate ts/daily.")
