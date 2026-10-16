#!/usr/bin/env python
"""6-hourly and daily time series for one member of a continuous (free-run) ensemble.

A free run writes every history tape at the 6-h interval, four records per
daily file. Each variable carries its own cell_methods:
  time: mean / maximum / minimum   record stamped at the interval end, with
                                   bounds; the file named for day d spans 18Z
                                   of d-1 to 18Z of d, and a hybrid start adds
                                   a one-timestep record at the start time
  time: point                      instantaneous state at the stamped time
Records are therefore selected by time, not file names: interval variables keep
the full-length intervals whose lower bound lies in [start, end), point
variables keep the states at [start, end). Every 6-h step must be present.

Products, in the same layout as the cycled runs (cycled_split.py) and the
hindcasts, under <archive>/post/<atm|lnd>/180x360_aave/:
  ts/6hourly     interval records labelled at the interval midpoint (with
                 bounds); point records at their valid time
  ts/daily_mean  day D 00Z-24Z: mean of the four 6-h means; for point
                 variables the mean of the 00/06/12/18Z states
  ts/daily_max   max of the four 6-h maxima;  ts/daily_min  min of the minima
  ts/daily_00z   point variables at 00Z of D
  monthly/       mean of ts/daily_mean over the days of each month in range,
                 VAR.ENxx.YYYY-MM.nc
ts files are VAR.ENxx.YYYY.nc; daily values are labelled 12Z with bounds
[D, D+1) (daily_00z: 00Z, bounds [D, D]). A free run has no analysis, so there
is no ts/6hourly_analysis. The first run on a member moves the products of the
earlier version (ts/daily, monthly, clim) to legacy/ and leaves a LAYOUT marker.

Every variable on the component's tapes with a time and horizontal dimension is
processed. A name on more than one tape is taken from the first tape listed.
Native fields are reduced in time first, then regridded with ncremap.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

COMPONENTS = {
    # component: (history directory, tapes in priority order, horizontal dims)
    "atm": ("atm", ("eam.h2", "eam.h1", "eam.h0"), {"ncol"}),
    "lnd": ("lnd", ("elm.h1", "elm.h0"), {"lat", "lon"}),
}
KINDS = {"time: mean": "mean", "time: maximum": "max", "time: minimum": "min", "time: point": "point"}
DAILY_CELL_METHODS = {
    "mean": "time: mean (00Z to 24Z of the day, from four 6-h windows)",
    "max": "time: maximum (00Z to 24Z of the day, from four 6-h windows)",
    "min": "time: minimum (00Z to 24Z of the day, from four 6-h windows)",
    "point": "time: mean (average of the states at 00, 06, 12 and 18Z)",
}
FILL_VALUE = np.float32(1.0e36)


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def tape_files(hist, case_name, tape, start, end):
    # The 18-24Z interval of the last day is the first record of the end-date file.
    first, last = start.strftime("%Y-%m-%d"), end.strftime("%Y-%m-%d")
    files = [path for path in sorted(hist.glob(f"{case_name}.{tape}.*.nc"))
             if first <= path.name.split(f".{tape}.")[1][:10] <= last]
    if not files:
        fail(f"no {case_name}.{tape} files from {first} to {last} in {hist}")
    return files


def classify(dataset, horizontal, taken):
    """Map each new variable on this tape to its kind, by cell_methods."""
    kinds = {}
    for name, var in dataset.data_vars.items():
        if name in taken or var.dims[:1] != ("time",) or not horizontal <= set(var.dims):
            continue
        kind = KINDS.get(var.attrs.get("cell_methods", "").strip())
        if kind is None:
            print(f"WARNING: {name}: unsupported cell_methods {var.attrs.get('cell_methods')!r}; skipped", file=sys.stderr)
            continue
        kinds[name] = kind
    return kinds


def select_records(dataset, start, end, interval, point):
    """Indices and 6-hourly labels of the records in [start, end), plus bounds."""
    expected = pd.date_range(start, end, freq=interval, inclusive="left")
    if point:
        times = pd.DatetimeIndex(dataset["time"].values)
        valid = (times >= start) & (times < end)
        keep = pd.Series(np.flatnonzero(valid), index=times[valid])
    else:
        bounds = dataset[dataset["time"].attrs.get("bounds", "time_bnds")].load()
        lower = pd.DatetimeIndex(bounds.isel({bounds.dims[-1]: 0}).values)
        upper = pd.DatetimeIndex(bounds.isel({bounds.dims[-1]: 1}).values)
        # Full-length intervals only (drops the hybrid-start timestep record).
        valid = (np.abs((upper - lower) - interval) <= pd.Timedelta("1min")) & (lower >= start) & (lower < end)
        keep = pd.Series(np.flatnonzero(valid), index=lower[valid])
    # A restart that rewrites a record keeps the last copy.
    keep = keep[~keep.index.duplicated(keep="last")].sort_index()
    absent = expected.difference(keep.index)
    if len(absent):
        listed = ", ".join(f"{t:%Y-%m-%d %H}Z" for t in absent[:8])
        fail(f"{len(absent)} of {len(expected)} 6-h {'states' if point else 'intervals'} missing (starting {listed})")
    if point:
        return keep.to_numpy(), expected, None
    return keep.to_numpy(), expected + interval / 2, np.stack([expected, expected + interval], axis=1)


def with_bounds(fields, labels, bounds):
    fields = fields.assign_coords(time=labels.to_numpy())
    fields["time"].attrs.update(long_name="time")
    if bounds is not None:
        fields["time_bnds"] = xr.DataArray(bounds, dims=("time", "nbnd"))
        fields["time"].attrs["bounds"] = "time_bnds"
    return fields


def daily_series(fields, kinds, per_day, days):
    """Reduce consecutive groups of ``per_day`` 6-h records to days."""
    reduced = {}
    for name in kinds:
        grouped = fields[name].coarsen(time=per_day)
        kind = kinds[name]
        reduced[name] = (grouped.max if kind == "max" else grouped.min if kind == "min" else grouped.mean)(keep_attrs=True)
        reduced[name].attrs["cell_methods"] = DAILY_CELL_METHODS[kind]
    daily = xr.Dataset(reduced).assign_coords(time=(days + pd.Timedelta("12h")).to_numpy())
    return with_bounds(daily, days + pd.Timedelta("12h"), np.stack([days, days + pd.Timedelta("1D")], axis=1))


def varying_mask(subset, names):
    """Variables whose missing values change between records. ncremap -P elm
    applies the first record's mask to every record, so the fill value leaks
    into cells that become missing later (seen with SNORDSL); such variables
    are regridded one record at a time."""
    out = []
    for name in names:
        var = subset[name]
        if "time" in var.dims and var.sizes["time"] > 1:
            mask = np.isnan(var.values)
            if (mask != mask[:1]).any():
                out.append(name)
    return out


def regrid(subset, names, encoding, component, map_file, workdir, native, regridded):
    """Regrid ``names`` of ``subset``; returns {name: regridded file}."""
    others = [v for v in subset.variables if v not in names and v != "time"]   # bounds, landfrac
    varying = varying_mask(subset, names) if component == "lnd" else []
    steady = [name for name in names if name not in varying]

    def remap(ds, src, dst):
        # time must be the record dimension: NCO tools average only along it
        ds.to_netcdf(src, encoding={k: v for k, v in encoding.items() if k in ds.variables}, unlimited_dims=["time"])
        command = ["ncremap", "-m", str(map_file), "-i", str(src), "-o", str(dst)]
        if component == "lnd":
            command[1:1] = ["-P", "elm"]
        subprocess.run(command, check=True, cwd=workdir)
        src.unlink()

    out = {}
    if steady:
        remap(subset[steady + others], native, regridded)
        out.update(dict.fromkeys(steady, regridded))
    if varying:
        parts = []
        for i in range(subset.sizes["time"]):
            part = regridded.with_name(f"{regridded.stem}_r{i:04d}.nc")
            remap(subset[varying + others].isel(time=[i]), native.with_name(f"{native.stem}_r{i:04d}.nc"), part)
            parts.append(part)
        joined = regridded.with_name(f"{regridded.stem}_varying.nc")
        subprocess.run(["ncrcat", "-O"] + [str(p) for p in parts] + [str(joined)], check=True)
        for part in parts:
            part.unlink()
        out.update(dict.fromkeys(varying, joined))
    return out


def regrid_period(dataset, names, period, select, member, component, map_file, ts_dest, workdir, tag, extras, overwrite):
    """Regrid the records ``select`` of ``names`` and split them into per-variable
    files <ts_dest>/<VAR>.<member>.<period>.nc; returns the target paths."""
    targets = {name: ts_dest / f"{name}.{member}.{period}.nc" for name in names}
    ts_dest.mkdir(parents=True, exist_ok=True)
    if not overwrite and all(path.is_file() for path in targets.values()):
        print(f"Reusing {ts_dest.name} {tag} {period}: all {len(targets)} files exist")
        return list(targets.values())
    native = workdir / f"native_{tag}_{period}.nc"
    regridded = workdir / f"regridded_{tag}_{period}.nc"
    subset = dataset[names + [v for v in ("time_bnds",) if v in dataset]].isel(time=np.flatnonzero(select))
    for name in extras:
        subset[name] = extras[name]
    units = f"days since {period[:4]}-01-01 00:00:00"
    encoding = {"time": {"units": units, "calendar": "gregorian", "dtype": "float64", "_FillValue": None}}
    if "time_bnds" in subset:
        encoding["time_bnds"] = {"units": units, "calendar": "gregorian", "dtype": "float64", "_FillValue": None}
    for name in names + list(extras):
        encoding[name] = {"_FillValue": FILL_VALUE, "dtype": "float32"}
    sources = regrid(subset, names, encoding, component, map_file, workdir, native, regridded)
    for name, target in targets.items():
        partial = target.with_name(target.name + ".tmp")
        subprocess.run(["ncks", "-O", "-4", "-L", "1", "-v", name, str(sources[name]), str(partial)], check=True)
        partial.replace(target)
    for path in set(sources.values()):
        path.unlink()
    print(f"Wrote {len(targets)} {ts_dest.name} series ({tag}) for {member} {period}")
    return list(targets.values())


def at_00z(fields, names, per_day, days):
    """Point variables at 00Z of each day (bounds [D, D])."""
    out = fields[names].isel(time=slice(0, None, per_day))
    out = with_bounds(out.drop_vars([v for v in ("time_bnds",) if v in out]), days, np.stack([days, days], axis=1))
    for name in names:
        out[name].attrs["cell_methods"] = "time: point (00Z of the day)"
    return out


def monthly_means(daily, names, days):
    """Mean of the daily values over the days of each month in the range."""
    months = days.to_period("M")
    labels, bounds, values = [], [], {name: [] for name in names}
    for month in months.unique():
        sel = np.flatnonzero(months == month)
        first, last = days[sel[0]], days[sel[-1]] + pd.Timedelta("1D")
        labels.append(first + (last - first) / 2)
        bounds.append([first, last])
        for name in names:
            values[name].append(daily[name].isel(time=sel).mean("time", keep_attrs=True))
    out = xr.Dataset({name: xr.concat(values[name], dim="time") for name in names})
    for name in names:
        out[name].attrs["cell_methods"] = "time: mean (mean of the daily means of the days in range)"
    return with_bounds(out, pd.DatetimeIndex(labels), np.array(bounds, dtype="datetime64[ns]")), \
        np.array([m.strftime("%Y-%m") for m in months.unique()])


def retire_legacy(root):
    """Move the earlier free-run products aside once, so no directory mixes layouts."""
    marker = root / "LAYOUT"
    if marker.is_file():
        return
    moved = []
    for sub in ("ts/daily", "monthly", "clim"):
        src, dst = root / sub, root / "legacy" / sub
        if src.exists():
            if dst.exists():
                fail(f"cannot move {src} to {dst}: target exists")
            dst.parent.mkdir(parents=True, exist_ok=True)
            src.rename(dst)
            moved.append(sub)
    root.mkdir(parents=True, exist_ok=True)
    marker.write_text("layout=free_split\n" + "".join(f"moved_to_legacy={sub}\n" for sub in moved))
    if moved:
        print(f"Moved earlier products to {root / 'legacy'}: {' '.join(moved)}")


def batches(names, size):
    for index in range(0, len(names), size):
        yield names[index:index + size]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--component", choices=sorted(COMPONENTS), required=True)
    parser.add_argument("--archive", type=Path, required=True, help="member archive directory (<case>/ENxx/archive)")
    parser.add_argument("--case-name", required=True, help="history file prefix, <case>.ENxx")
    parser.add_argument("--member", required=True, help="member label used in output names, ENxx")
    parser.add_argument("--start", required=True, help="first day, YYYY-MM-DD")
    parser.add_argument("--end", required=True, help="day after the last day, YYYY-MM-DD")
    parser.add_argument("--map-file", type=Path, required=True)
    parser.add_argument("--interval-hours", type=int, default=6)
    parser.add_argument("--batch-size", type=int, default=20, help="variables per regridding call")
    parser.add_argument("--manifest", type=Path, required=True, help="list of every file written, for the driver")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    start, end = pd.Timestamp(args.start), pd.Timestamp(args.end)
    interval = pd.Timedelta(hours=args.interval_hours)
    if start.normalize() != start or end.normalize() != end or end <= start:
        fail("--start and --end must be dates with start < end")
    if pd.Timedelta("1D") % interval:
        fail("--interval-hours must divide 24")
    if not args.map_file.is_file():
        fail(f"map file not found: {args.map_file}")
    per_day = pd.Timedelta("1D") // interval
    days = pd.date_range(start, end, freq="1D", inclusive="left")

    hist_dir, tapes, horizontal = COMPONENTS[args.component]
    hist = args.archive / hist_dir / "hist"
    root = args.archive / "post" / hist_dir / "180x360_aave"
    retire_legacy(root)
    dest = {name: root / "ts" / name for name in ("6hourly", "daily_mean", "daily_max", "daily_min", "daily_00z")}
    dest["monthly"] = root / "monthly"
    years = np.array([f"{t.year:04d}" for t in days])
    written = []
    workdir = Path(tempfile.mkdtemp(prefix=f"tmp_free_hfreq_{args.component}_{args.member}_", dir=args.archive / "post"))

    taken = set()
    try:
        for tape in tapes:
            dataset = xr.open_mfdataset(
                tape_files(hist, args.case_name, tape, start, end), combine="nested", concat_dim="time",
                data_vars="minimal", coords="minimal", compat="override", chunks={"time": 1},
            )
            kinds = classify(dataset, horizontal, taken)
            if not kinds:
                print(f"{tape}: no new variables")
                continue
            # ncremap -P elm weights by the land fraction, so it travels with every file.
            extras = {"landfrac": dataset["landfrac"].load()} if args.component == "lnd" else {}
            for point in (False, True):
                names = sorted(name for name, kind in kinds.items() if (kind == "point") == point)
                if not names:
                    continue
                indices, labels, bounds = select_records(dataset, start, end, interval, point)
                group = "point" if point else "interval"
                print(f"{tape}: {len(names)} {group} variable(s), {len(indices)} records")
                for number, chunk in enumerate(batches(names, args.batch_size)):
                    fields = with_bounds(dataset[chunk].isel(time=indices), labels, bounds).load()
                    daily = daily_series(fields, {name: kinds[name] for name in chunk}, per_day, days)
                    tag = f"{tape}_{group}{number}"
                    common = (args.member, args.component, args.map_file)
                    six_years = np.array([f"{t.year:04d}" for t in labels])
                    for year in sorted(set(six_years)):
                        written += regrid_period(fields, chunk, year, six_years == year, *common,
                                                 dest["6hourly"], workdir, tag, extras, args.overwrite)
                    # daily values: means (and point-variable means) to daily_mean, extremes to daily_max/min
                    by_dest = {"daily_mean": [n for n in chunk if kinds[n] in ("mean", "point")],
                               "daily_max": [n for n in chunk if kinds[n] == "max"],
                               "daily_min": [n for n in chunk if kinds[n] == "min"]}
                    for product, product_names in by_dest.items():
                        for year in sorted(set(years)) if product_names else ():
                            written += regrid_period(daily, product_names, year, years == year, *common,
                                                     dest[product], workdir, tag, extras, args.overwrite)
                    if point:
                        at00 = at_00z(fields, chunk, per_day, days)
                        for year in sorted(set(years)):
                            written += regrid_period(at00, chunk, year, years == year, *common,
                                                     dest["daily_00z"], workdir, tag, extras, args.overwrite)
                    if by_dest["daily_mean"]:
                        monthly, months = monthly_means(daily, by_dest["daily_mean"], days)
                        for month in months:
                            written += regrid_period(monthly, by_dest["daily_mean"], month, months == month, *common,
                                                     dest["monthly"], workdir, tag, extras, args.overwrite)
            taken.update(kinds)
            dataset.close()
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    partial = args.manifest.with_name(args.manifest.name + ".tmp")
    partial.write_text("".join(f"{path}\n" for path in written))
    partial.replace(args.manifest)
    counts = pd.Series([path.parent.name for path in written]).value_counts().sort_index()
    for product, count in counts.items():
        print(f"  {product:<12} {count:6d} file(s)")
    print(f"Done: {len(taken)} variables for {args.member} ({args.component}); manifest {args.manifest}")


if __name__ == "__main__":
    main()
