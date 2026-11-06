"""Complete an ERA5-derived EAM initial file with the fields it lacks (``complete_era5_initial_states.ipynb``).

The ERA5 initial files (``v3.LR.ne30pg2.ERA5.eam.i.<date>.nc``) hold only T, U, V, Q, PS, CLDLIQ and
CLDICE on the GLL dynamics columns. A full E3SMv3 initial file also holds ozone, the Linoz tracers,
chemistry, MAM aerosols, hydrometeor numbers and surface fields; EAM starts any missing field from
zero. The radiation reads its gases from tracers (``rad_climate = 'A:H2OLNZ:H2O', 'A:O3:O3',
'A:N2OLNZ:N2O', 'A:CH4LNZ:CH4'``), so a run from the ERA5 file starts radiatively dry, without ozone,
N2O or CH4.

The completed file takes the native layout of a DART member file (dynamics fields on ``ncol_d``,
physics-grid surface fields on ``ncol``), every state field from the mean of the given DART members,
and then the ERA5 values of the ERA5 fields. ``H2OLNZ`` (radiative water vapour) follows the ERA5
moisture with the model's own relative difference from ``Q`` (multiplicative, so it stays positive
where ERA5 is drier than the DART mean):

    H2OLNZ = Q_ERA5 * (mean H2OLNZ_DART / mean Q_DART), the ratio bounded to ``ratio_bounds``.
"""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

ERA5_VARS = ("T", "U", "V", "Q", "PS", "CLDLIQ", "CLDICE")
GRAVITY = 9.80616
DU_KG_M2 = 2.1415e-5  # one Dobson unit of ozone in kg m-2


def dart_member_file(dart_root: Path, case: str, member: str, date_tag: str) -> Path:
    """``<dart_root>/<case>/<member>/archive/rest/<date>/<case>.<member>.eam.i.<date>.nc``."""
    return Path(dart_root) / case / member / "archive" / "rest" / date_tag / f"{case}.{member}.eam.i.{date_tag}.nc"


def _is_state(da: xr.DataArray) -> bool:
    """A field averaged over members: floating point, with time and a column dimension."""
    return da.dtype.kind == "f" and "time" in da.dims and bool({"ncol", "ncol_d"} & set(da.dims))


def member_mean(files: list[Path]) -> xr.Dataset:
    """The first file with every state field replaced by the mean over ``files``.

    One member file is open at a time (holding all 20 open exhausts memory through their caches).
    """
    with xr.open_dataset(files[0], decode_times=False) as first:
        out = first.load()
    state = [v for v in out.data_vars if _is_state(out[v])]
    acc = {v: out[v].values.astype("float64", copy=True) for v in state}
    for f in files[1:]:
        with xr.open_dataset(f, decode_times=False) as ds:
            for v in state:
                acc[v] += ds[v].values
    for v in state:
        out[v].values[...] = acc.pop(v) / len(files)
    return out


def column_ozone_du(ds: xr.Dataset, o3: str = "O3") -> xr.DataArray:
    """Total ozone column (DU) per dynamics column from the O3 mass mixing ratio."""
    ps = ds["PS"].isel(time=0)
    pi = ds["hyai"] * 1.0e5 + ds["hybi"] * ps            # (ilev, ncol_d) interface pressure, Pa
    dp = pi.diff("ilev").rename(ilev="lev").assign_coords(lev=ds["lev"])
    return (ds[o3].isel(time=0) * dp).sum("lev") / GRAVITY / DU_KG_M2


def complete_era5_ic(era5_file: Path, dart_files: list[Path], out_file: Path, *, era5_vars=ERA5_VARS,
                     ratio_bounds=(0.2, 5.0), overwrite: bool = False) -> dict:
    """Write the completed initial file and return its checks (see ``validate``)."""
    out_file = Path(out_file)
    if out_file.exists() and not overwrite:
        raise FileExistsError(f"{out_file} exists; pass overwrite=True to replace it")
    era5 = xr.open_dataset(era5_file, decode_times=False).load()
    if "ncol" not in era5.dims or "ncol_d" in era5.dims:
        raise ValueError(f"{era5_file}: expected the GLL columns on 'ncol'")
    out = member_mean(dart_files)
    if era5.sizes["ncol"] != out.sizes["ncol_d"]:
        raise ValueError(f"ERA5 has {era5.sizes['ncol']} columns, the DART files {out.sizes['ncol_d']}")
    if np.abs(era5["lat"].values - out["lat_d"].values).max() > 1e-8 or \
            np.abs((era5["lon"].values - out["lon_d"].values + 180) % 360 - 180).max() > 1e-8:
        raise ValueError("ERA5 and DART columns are not in the same order")
    for c in ("hyam", "hybm", "hyai", "hybi"):
        if not np.array_equal(era5[c].values, out[c].values):
            raise ValueError(f"{c} differs between ERA5 and DART files")

    q_mean, h_mean = out["Q"].values.copy(), out["H2OLNZ"].values.copy()
    for v in era5_vars:
        new = era5[v].rename(ncol="ncol_d").transpose(*out[v].dims).values
        out[v].values[...] = new
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(q_mean > 0, h_mean / q_mean, 1.0)
    lo, hi = ratio_bounds
    n_bounded = int(((ratio < lo) | (ratio > hi)).sum())
    out["H2OLNZ"].values[...] = out["Q"].values * np.clip(ratio, lo, hi)

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    out.attrs.update(
        title="ERA5 initial state completed with DART ensemble-mean fields",
        completed_era5_source=str(era5_file),
        completed_era5_vars=",".join(era5_vars),
        completed_dart_members=",".join(str(f) for f in dart_files),
        completed_h2olnz="Q_ERA5 * (mean H2OLNZ_DART / mean Q_DART), ratio bounded to [%g, %g]" % ratio_bounds,
        history=f"{stamp} complete_era5_ic.py; " + out.attrs.get("history", ""),
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)
    encoding = {v: {"zlib": True, "complevel": 1} for v in out.data_vars if out[v].dtype.kind == "f" and out[v].ndim > 1}
    tmp = out_file.with_suffix(".tmp.nc")
    out.to_netcdf(tmp, format="NETCDF4", encoding=encoding)
    tmp.rename(out_file)
    return validate(out_file, era5_file, dart_files[0], era5_vars=era5_vars) | {"h2olnz_ratio_bounded_points": n_bounded}


def validate(out_file: Path, era5_file: Path, dart_file: Path, *, era5_vars=ERA5_VARS) -> dict:
    """Checks of a completed file: ERA5 fields exact, nothing missing or non-finite, ozone column."""
    out = xr.open_dataset(out_file, decode_times=False)
    era5 = xr.open_dataset(era5_file, decode_times=False)
    ref = xr.open_dataset(dart_file, decode_times=False)
    try:
        missing = sorted(set(ref.data_vars) - set(out.data_vars))
        era5_maxdiff = {v: float(np.abs(out[v].values - era5[v].rename(ncol="ncol_d")
                                        .transpose(*out[v].dims).values).max()) for v in era5_vars}
        nonfinite = [v for v in out.data_vars if out[v].dtype.kind == "f" and not np.isfinite(out[v].values).all()]
        toz = column_ozone_du(out)
        hq = (out["H2OLNZ"] / out["Q"]).isel(time=0)
        low = hq.sel(lev=slice(500, None))
        return dict(file=str(out_file), missing_vs_dart=missing, era5_max_abs_diff=era5_maxdiff,
                    nonfinite_vars=nonfinite, toz_du_mean=float(toz.mean()), toz_du_min=float(toz.min()),
                    toz_du_max=float(toz.max()), h2olnz_over_q_below_500hPa_median=float(low.median()),
                    n_vars=len(out.data_vars))
    finally:
        out.close(); era5.close(); ref.close()


__all__ = ["ERA5_VARS", "dart_member_file", "member_mean", "column_ozone_du", "complete_era5_ic", "validate"]
