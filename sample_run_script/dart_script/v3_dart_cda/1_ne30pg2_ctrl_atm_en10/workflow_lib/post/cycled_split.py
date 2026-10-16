#!/usr/bin/env python
"""Analysis, forecast, daily and monthly series for one member of a DA-cycled ensemble.

Every DA cycle is a separate 6-h model run started from the DART analysis, and
each history tape writes one file per cycle, named for the cycle start:
  record 0   nstep 0, zero-width time bounds. Point fields hold the analysis
             state. Mean/max/min fields hold a single-timestep value right
             after the analysis is inserted (spin-up), not a time statistic.
  record 1   only on tapes written every 6 h (nhtfrq=-6): bounds span the
             cycle. Mean/max/min fields are statistics over the 6-h forecast;
             point fields are the background state at the cycle end.
A daily (nhtfrq=-24) tape never reaches its output time in a 6-h run, so it
holds record 0 only; a monthly tape has no per-cycle files and is skipped.
Records are selected by their time bounds, not their position in the file.

Products, under <archive>/post/<atm|lnd>/180x360_aave/:
  ts/6hourly           forecast: 6-h means/max/min labelled at the window
                       middle (with bounds); point fields (background) at the
                       window end
  ts/6hourly_analysis  analysis state at every cycle start: record-0 point
                       fields, plus the state variables in STATE that a tape
                       labels "time: mean" (their nstep-0 value is the state)
  ts/daily_mean        day D, 00Z-24Z: mean of the four 6-h forecast means;
                       fields without a forecast mean use the mean of the
                       00/06/12/18Z analysis states (stated in cell_methods)
  ts/daily_max/min     max/min of the four 6-h forecast maxima/minima
  ts/daily_00z         analysis state at 00Z of D
  monthly/             mean of ts/daily_mean over the days of each month in
                       the range, VAR.ENxx.YYYY-MM.nc
ts files are VAR.ENxx.YYYY.nc. Daily and monthly values are labelled at their
middle with bounds [D, D+1) (daily_00z: 00Z, bounds [D, D]). Record-0 values
of mean/max/min fields are never used. Every cycle in [start, end) must be
present, and every per-cycle tape must hold the 6-h records unless it is
listed in --analysis-only (archives of runs set up with nhtfrq=-24 or with
mfilt=1, where the 6-h records were never kept). The first run on a member moves the outputs of the earlier cycled
post-processing (ts/6hourly, ts/6hourly_da, ts/daily, monthly, clim) to
legacy/ and leaves a LAYOUT marker.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr

COMPONENTS = {
    # component: (history directory, tapes in priority order, horizontal dims)
    "atm": ("atm", ("eam.h2", "eam.h1", "eam.h0"), {"ncol"}),
    "lnd": ("lnd", ("elm.h1", "elm.h0"), {"lat", "lon"}),
}
# State variables that tapes label "time: mean"; at nstep 0 their value is the
# instantaneous model state, so record 0 is kept for them as analysis.
STATE = {
    "atm": re.compile(r"^((U|V|T|Q|Z|OMEGA)\d{3,4}|PS|PSL)$"),
    "lnd": re.compile(r"^(TSOI|TSOI_10CM|TSOI_ICE|TG|H2OSOI|SOILLIQ|SOILICE|SOILLIQ_ICE|SOILICE_ICE|"
                      r"SOILWATER_10CM|H2OSNO|SNOW_DEPTH|FSNO|H2OSFC|TWS|ZWT|TLAI|HC|HCSOI)$"),
}
KINDS = {"time: mean": "mean", "time: maximum": "max", "time: minimum": "min", "time: point": "point"}
CELL = {"mean": "time: mean", "max": "time: maximum", "min": "time: minimum"}
LEGACY = ("ts/6hourly", "ts/6hourly_da", "ts/daily", "monthly", "clim")
FILL_VALUE = np.float32(1.0e36)
BATCH_BYTES = 2 * 1024 ** 3
TOL = pd.Timedelta("1min")


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def stamp(t):
    return f"{t:%Y-%m-%d}-{t.hour * 3600 + t.minute * 60 + t.second:05d}"


def cycle_files(hist, case_name, tape, stamps):
    """Per-cycle files of a tape, or None if it has none (e.g. a monthly tape)."""
    paths = [hist / f"{case_name}.{tape}.{stamp(t)}.nc" for t in stamps]
    present = [path.is_file() for path in paths]
    if not any(present):
        return None
    if not all(present):
        missing = [path.name for path, ok in zip(paths, present) if not ok]
        fail(f"{tape}: {len(missing)} of {len(paths)} cycle files missing (first: {missing[0]})")
    return paths


def bounds_name(dataset):
    for name in (dataset["time"].attrs.get("bounds"), dataset["time"].encoding.get("bounds"), "time_bnds", "time_bounds"):
        if name and name in dataset:
            return name
    fail("history has no time bounds variable")


def select(dataset, stamps, interval, tape):
    """Record indices of the analysis states and of the forecast windows, one
    per cycle start; the forecast indices are None when the tape writes no 6-h
    records. The nstep-0 record ends at the cycle start and has zero width
    (EAM) or one model timestep (ELM); a forecast record spans the cycle."""
    bounds = dataset[bounds_name(dataset)].load().values
    if not np.issubdtype(bounds.dtype, np.datetime64):
        fail(f"{tape}: time bounds were not decoded to dates")
    lower, upper = pd.DatetimeIndex(bounds[:, 0]), pd.DatetimeIndex(bounds[:, 1])
    width = upper - lower

    def by_time(mask, times):
        found = pd.Series(np.flatnonzero(mask), index=times[mask])
        # a rerun cycle rewrites its records: keep the last copy
        found = found[~found.index.duplicated(keep="last")].reindex(stamps)
        return found, found.isna().to_numpy()

    analysis, absent = by_time(np.asarray(width < interval / 2), upper)
    if absent.any():
        fail(f"{tape}: {absent.sum()} of {len(stamps)} analysis records missing (first {stamps[absent][0]:%Y-%m-%d %H}Z)")
    forecast, absent = by_time(np.asarray(abs(width - interval) <= TOL), lower)
    if absent.all():
        return analysis.to_numpy(int), None
    if absent.any():
        fail(f"{tape}: {absent.sum()} of {len(stamps)} 6-h forecast records missing (first {stamps[absent][0]:%Y-%m-%d %H}Z)")
    return analysis.to_numpy(int), forecast.to_numpy(int)


def classify(dataset, horizontal):
    """Map each history field on the tape to its kind, by cell_methods."""
    kinds = {}
    for name, var in dataset.data_vars.items():
        if var.dims[:1] != ("time",) or not horizontal <= set(var.dims):
            continue
        kind = KINDS.get(var.attrs.get("cell_methods", "").strip())
        if kind is None:
            print(f"WARNING: {name}: unsupported cell_methods {var.attrs.get('cell_methods')!r}; skipped", file=sys.stderr)
            continue
        kinds[name] = kind
    return kinds


def batches(dataset, names, nrec):
    """Groups of variables of at most BATCH_BYTES each once loaded."""
    group, size = [], 0
    for name in names:
        var = dataset[name]
        nbytes = 4 * nrec * int(np.prod([var.sizes[d] for d in var.dims if d != "time"]))
        if group and size + nbytes > BATCH_BYTES:
            yield group
            group, size = [], 0
        group.append(name)
        size += nbytes
    if group:
        yield group


def relabel(fields, labels, bounds, cell_methods):
    """Give fields new time labels (and bounds) and per-variable cell_methods."""
    fields = fields.drop_vars([v for v in ("time_bnds", "time_bounds") if v in fields])
    fields = fields.assign_coords(time=pd.DatetimeIndex(labels).to_numpy())
    fields["time"].attrs = {"long_name": "time"}
    if bounds is not None:
        fields["time_bnds"] = xr.DataArray(np.asarray(bounds, dtype="datetime64[ns]"), dims=("time", "nbnd"))
        fields["time"].attrs["bounds"] = "time_bnds"
    for name, method in cell_methods.items():
        fields[name].attrs["cell_methods"] = method
    return fields


def reduce_days(fields, names, how, per_day):
    """Consecutive groups of per_day records reduced to one value per day."""
    out = {}
    for name in names:
        grouped = fields[name].drop_vars("time").coarsen(time=per_day, boundary="exact")
        out[name] = getattr(grouped, how[name])(keep_attrs=True)
    return xr.Dataset(out)


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
    return relabel(out, labels, bounds, {}), np.array([m.strftime("%Y-%m") for m in months.unique()])


def write_periods(fields, names, dest, periods, ctx, tag):
    """Regrid fields and write one file per variable and period
    (dest/VAR.ENxx.<period>.nc); returns the target paths."""
    dest.mkdir(parents=True, exist_ok=True)
    written = []
    for period in sorted(set(periods)):
        targets = {name: dest / f"{name}.{ctx.member}.{period}.nc" for name in names}
        written += targets.values()
        if not ctx.overwrite and all(path.is_file() for path in targets.values()):
            continue
        keep = names + [v for v in ("time_bnds",) if v in fields]
        subset = fields[keep].isel(time=np.flatnonzero(periods == period))
        for name in ctx.extras:
            subset[name] = ctx.extras[name]
        for name in subset.variables:
            subset[name].encoding = {}
        subset.attrs = {"source_case": ctx.case_name, "product": f"cycled_split.py {dest.relative_to(ctx.root)}"}
        units = f"days since {period[:4]}-01-01 00:00:00"
        encoding = {"time": {"units": units, "calendar": ctx.calendar, "dtype": "float64", "_FillValue": None}}
        if "time_bnds" in subset:
            encoding["time_bnds"] = {"units": units, "calendar": ctx.calendar, "dtype": "float64", "_FillValue": None}
        for name in names + list(ctx.extras):
            encoding[name] = {"_FillValue": FILL_VALUE, "dtype": "float32"}
        native = ctx.workdir / f"native_{tag}.nc"
        regridded = ctx.workdir / f"regridded_{tag}.nc"
        # time must be the record dimension for NCO tools downstream
        subset.to_netcdf(native, encoding=encoding, unlimited_dims=["time"])
        command = ["ncremap", "-m", str(ctx.map_file), "-i", str(native), "-o", str(regridded)]
        if ctx.component == "lnd":
            command[1:1] = ["-P", "elm"]
        subprocess.run(command, check=True, cwd=ctx.workdir, stdout=subprocess.DEVNULL)
        for name, target in targets.items():
            partial = target.with_name(target.name + ".tmp")
            subprocess.run(["ncks", "-O", "-4", "-L", "1", "-v", name, str(regridded), str(partial)], check=True)
            partial.replace(target)
        native.unlink()
        regridded.unlink()
    return written


def retire_legacy(root):
    """Move the earlier cycled products aside once, so no directory mixes layouts."""
    marker = root / "LAYOUT"
    if marker.is_file():
        return
    moved = []
    for sub in LEGACY:
        src, dst = root / sub, root / "legacy" / sub
        if src.exists():
            if dst.exists():
                fail(f"cannot move {src} to {dst}: target exists")
            dst.parent.mkdir(parents=True, exist_ok=True)
            src.rename(dst)
            moved.append(sub)
    root.mkdir(parents=True, exist_ok=True)
    marker.write_text("layout=cycled_split\n" + "".join(f"moved_to_legacy={sub}\n" for sub in moved))
    if moved:
        print(f"Moved earlier products to {root / 'legacy'}: {' '.join(moved)}")


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
    parser.add_argument("--manifest", type=Path, required=True, help="list of every file written, for the driver")
    parser.add_argument("--analysis-only", default="",
                        help="comma-separated tapes allowed to have no 6-h records, e.g. eam.h1,elm.h1")
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
    stamps = pd.date_range(start, end, freq=interval, inclusive="left")
    days = pd.date_range(start, end, freq="1D", inclusive="left")
    years = np.array([f"{t.year:04d}" for t in stamps])
    day_years = np.array([f"{t.year:04d}" for t in days])

    analysis_only = {tape for tape in args.analysis_only.replace(" ", ",").split(",") if tape}
    hist_dir, tapes, horizontal = COMPONENTS[args.component]
    state = STATE[args.component]
    hist = args.archive / hist_dir / "hist"
    root = args.archive / "post" / hist_dir / "180x360_aave"
    retire_legacy(root)
    dest = {name: root / "ts" / name for name in
            ("6hourly", "6hourly_analysis", "daily_mean", "daily_max", "daily_min", "daily_00z")}
    dest["monthly"] = root / "monthly"
    workdir = Path(tempfile.mkdtemp(prefix=f"tmp_cycled_{args.component}_{args.member}_", dir=args.archive / "post"))
    ctx = SimpleNamespace(member=args.member, case_name=args.case_name, component=args.component,
                          map_file=args.map_file, workdir=workdir, root=root, overwrite=args.overwrite,
                          calendar="gregorian", extras={})

    # Pass 1: open every per-cycle tape and decide where each field comes from.
    plans, analysis_of, forecast_of = [], {}, {}
    for tape in tapes:
        files = cycle_files(hist, args.case_name, tape, stamps)
        if files is None:
            print(f"{tape}: no per-cycle files; skipped")
            continue
        dataset = xr.open_mfdataset(files, combine="nested", concat_dim="time", data_vars="minimal",
                                    coords="minimal", compat="override", chunks={"time": 1})
        ctx.calendar = dataset["time"].encoding.get("calendar", ctx.calendar)
        kinds = classify(dataset, horizontal)
        an_idx, fc_idx = select(dataset, stamps, interval, tape)
        unused = []
        for name, kind in kinds.items():
            if kind == "point" or state.match(name):
                analysis_of.setdefault(name, (len(plans), kind))
            if fc_idx is not None:
                forecast_of.setdefault(name, (len(plans), kind))
            elif kind != "point" and not state.match(name):
                unused.append(name)
        if fc_idx is None:
            if tape not in analysis_only:
                fail(f"{tape}: no record spans a {args.interval_hours}-h cycle, so its means, maxima and minima "
                     f"are lost (tape written with nhtfrq other than -{args.interval_hours}, or mfilt=1 and "
                     f"overwritten). If this archive is known to lack them, list {tape} in ANALYSIS_ONLY_TAPES "
                     f"of 7_run_post_hist.sh")
            print(f"{tape}: no 6-h forecast records (record 0 only; allowed by --analysis-only); analysis states only")
        if unused:
            print(f"{tape}: {len(unused)} single-timestep field(s) not used: {' '.join(sorted(unused))}")
        if args.component == "lnd" and not ctx.extras:
            # ncremap -P elm weights by the land fraction, so it travels with every file
            ctx.extras = {"landfrac": dataset["landfrac"].load()}
        plans.append(SimpleNamespace(tape=tape, dataset=dataset, an_idx=an_idx, fc_idx=fc_idx))
    if not plans:
        fail(f"no per-cycle history files for {args.case_name} from {args.start} in {hist}")
    has_forecast_mean = {name for name, (_, kind) in forecast_of.items() if kind == "mean"}

    # Pass 2: build and write the products, tape by tape, in memory-sized batches.
    written = []
    try:
        for number, plan in enumerate(plans):
            ds, n = plan.dataset, len(stamps)
            interval_names = sorted(k for k, (p, kind) in forecast_of.items() if p == number and kind != "point")
            point_names = sorted(k for k, (p, kind) in forecast_of.items() if p == number and kind == "point")
            analysis_names = sorted(k for k, (p, _) in analysis_of.items() if p == number)
            kind = {name: forecast_of[name][1] for name in interval_names}
            print(f"{plan.tape}: forecast {len(interval_names)} statistic + {len(point_names)} background, "
                  f"analysis {len(analysis_names)} field(s)")

            for b, chunk in enumerate(batches(ds, interval_names, n)):
                tag = f"{plan.tape}_fc{b}"
                fields = relabel(ds[chunk].isel(time=plan.fc_idx).load(), stamps + interval / 2,
                                 np.stack([stamps, stamps + interval], axis=1), {})
                written += write_periods(fields, chunk, dest["6hourly"], years, ctx, tag)
                daily = reduce_days(fields, chunk, kind, per_day)
                daily = relabel(daily, days + pd.Timedelta("12h"), np.stack([days, days + pd.Timedelta("1D")], axis=1),
                                {name: f"{CELL[kind[name]]} (00Z to 24Z of the day, from four 6-h forecast windows)"
                                 for name in chunk})
                for stat in ("mean", "max", "min"):
                    names = [name for name in chunk if kind[name] == stat]
                    if names:
                        written += write_periods(daily, names, dest[f"daily_{stat}"], day_years, ctx, tag)
                means = [name for name in chunk if kind[name] == "mean"]
                if means:
                    monthly, months = monthly_means(daily, means, days)
                    written += write_periods(monthly, means, dest["monthly"], months, ctx, tag)

            for b, chunk in enumerate(batches(ds, point_names, n)):
                fields = relabel(ds[chunk].isel(time=plan.fc_idx).load(), stamps + interval, None,
                                 {name: "time: point (background: end of the 6-h forecast)" for name in chunk})
                written += write_periods(fields, chunk, dest["6hourly"], years, ctx, f"{plan.tape}_bg{b}")

            for b, chunk in enumerate(batches(ds, analysis_names, n)):
                tag = f"{plan.tape}_an{b}"
                fields = relabel(ds[chunk].isel(time=plan.an_idx).load(), stamps, None,
                                 {name: "time: point (analysis: model state at the cycle start)" for name in chunk})
                written += write_periods(fields, chunk, dest["6hourly_analysis"], years, ctx, tag)
                at00 = relabel(fields.isel(time=slice(0, None, per_day)), days, np.stack([days, days], axis=1),
                               {name: "time: point (analysis at 00Z of the day)" for name in chunk})
                written += write_periods(at00, chunk, dest["daily_00z"], day_years, ctx, tag)
                means = [name for name in chunk if name not in has_forecast_mean]
                if means:
                    daily = reduce_days(fields, means, dict.fromkeys(means, "mean"), per_day)
                    daily = relabel(daily, days + pd.Timedelta("12h"), np.stack([days, days + pd.Timedelta("1D")], axis=1),
                                    {name: "time: mean (average of the analysis states at 00, 06, 12 and 18Z)"
                                     for name in means})
                    written += write_periods(daily, means, dest["daily_mean"], day_years, ctx, tag)
                    monthly, months = monthly_means(daily, means, days)
                    written += write_periods(monthly, means, dest["monthly"], months, ctx, tag)
            ds.close()
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    partial = args.manifest.with_name(args.manifest.name + ".tmp")
    partial.write_text("".join(f"{path}\n" for path in written))
    partial.replace(args.manifest)
    counts = pd.Series([path.parent.name for path in written]).value_counts().sort_index()
    for product, count in counts.items():
        print(f"  {product:<18} {count:6d} file(s)")
    print(f"Done: {len(written)} files for {args.member} ({args.component}); manifest {args.manifest}")


if __name__ == "__main__":
    main()
