"""
Surface density flux and deep convection in the North Atlantic deep-water formation regions, per decade, from the
remapped (0.5°) decadal MPAS-Ocean climatologies ``post/time_series/ocn_2d/<case>_<ANN|03>_YYYY..._climo.nc``.

Surface density flux (kg m-2 s-1, positive = the surface water gets denser), with TEOS-10 alpha, beta at the
surface temperature and salinity of the climatology:

    F_T = -(alpha / c_p) * Q_net                          thermal part
    F_S = -beta * S * F_fw / (1 - S/1000)  +  beta * 1000 * F_salt      haline part
    Q_net  = SW + LW_down + LW_up + latent + sensible + sea-ice heat     (W m-2, positive into the ocean)
    F_fw   = rain + snow + evaporation + river + ice runoff + sea-ice fresh water (kg m-2 s-1, positive into the ocean)
    F_salt = sea-ice salt flux (kg m-2 s-1, positive into the ocean)

Regional values are area-weighted means over the ocean cells of each box. March mixed-layer depth
(``dThreshMLD``) from the March climatologies measures deep convection (mean, maximum, and the area deeper than
``conv_depth``). MPAS's own ``surfaceBuoyancyForcing`` is kept as a check (opposite sign to the density flux).
"""
import os

import numpy as np
import xarray as xr
import gsw
import matplotlib.pyplot as plt

CP = 3996.0          # J kg-1 K-1 (MPAS-Analysis)
R_EARTH = 6371229.0  # m

# name -> (lat_min, lat_max, lon_min, lon_max), lon in -180..180
REGIONS = {
    "Labrador Sea":            (53.0, 65.0, -64.0, -45.0),
    "Irminger Sea":            (57.0, 66.0, -45.0, -25.0),
    "GIN Seas":                (65.0, 80.0, -20.0, 15.0),
    "Subpolar North Atlantic": (45.0, 65.0, -65.0, 0.0),
}
P = "timeMonthly_avg_"
HEAT = ["shortWaveHeatFlux", "longWaveHeatFluxDown", "longWaveHeatFluxUp", "latentHeatFlux", "sensibleHeatFlux",
        "seaIceHeatFlux"]
FRESH = ["rainFlux", "snowFlux", "evaporationFlux", "riverRunoffFlux", "iceRunoffFlux", "seaIceFreshWaterFlux"]
BOX = (45.0, 80.0, -70.0, 20.0)   # read window covering all regions


def climo_file(case_dir, case, period, y0, y1):
    """<case>_ANN_YYYY01_YYYY12_climo.nc (period "ANN") or <case>_MM_YYYYMM_YYYYMM_climo.nc (period "03" ...)."""
    if period == "ANN":
        return os.path.join(case_dir, f"{case}_ANN_{y0:04d}01_{y1:04d}12_climo.nc")
    return os.path.join(case_dir, f"{case}_{period}_{y0:04d}{period}_{y1:04d}{period}_climo.nc")


def _box(ds):
    la0, la1, lo0, lo1 = BOX
    return ds.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))


def _regional(field, area, regions):
    out = {}
    for name, (la0, la1, lo0, lo1) in regions.items():
        f = field.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
        a = area.sel(lat=slice(la0, la1), lon=slice(lo0, lo1)).where(f.notnull())
        out[name] = float((f * a).sum() / a.sum())
    return out


def density_flux_fields(path):
    """Surface density-flux fields (lat, lon) of one annual climatology over BOX."""
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = _box(ds).isel(Time=0)
        pt = ds[P + "activeTracers_temperature"].isel(nVertLevels=0).load()
        sp = ds[P + "activeTracers_salinity"].isel(nVertLevels=0).load()
        q = sum(ds[P + v].load() for v in HEAT)
        fw = sum(ds[P + v].load() for v in FRESH)
        fs = ds[P + "seaIceSalinityFlux"].load()
        b = ds[P + "surfaceBuoyancyForcing"].load()
        area = ds["area"].load() * R_EARTH ** 2
    lat2 = pt["lat"].broadcast_like(pt).values
    lon2 = pt["lon"].broadcast_like(pt).values
    sa_ = gsw.SA_from_SP(sp.values, 0.0, lon2, lat2)
    ct_ = gsw.CT_from_pt(sa_, pt.values)
    like = lambda a: xr.DataArray(a, dims=pt.dims, coords=pt.coords)
    alpha, beta = like(gsw.alpha(sa_, ct_, 0.0)), like(gsw.beta(sa_, ct_, 0.0))
    f_t = -(alpha / CP) * q
    f_s = -beta * sp * fw / (1.0 - sp / 1000.0) + beta * 1000.0 * fs
    sigma0 = like(gsw.sigma0(sa_, ct_))
    return xr.Dataset(dict(Q_net=q, F_fw=fw * 86400.0 * 365.0, F_T=f_t, F_S=f_s, F_rho=f_t + f_s, SST=pt, SSS=sp,
                           sigma0=sigma0, buoyancy=b)), area


def mld_fields(path):
    """March mixed-layer depth (lat, lon) over BOX."""
    with xr.open_dataset(path, decode_times=False) as ds:
        ds = _box(ds).isel(Time=0)
        return ds[P + "dThreshMLD"].load(), ds["area"].load() * R_EARTH ** 2


def regional_decades(case_dir, case, years, regions=REGIONS, conv_depth=1000.0, mld_month="03", decade=10,
                     verbose=True):
    """Dataset (region, decade) of regional means for every decade of ``years`` that has its climatologies."""
    rows = []
    for y0 in range(years[0], years[1] + 1, decade):
        y1 = y0 + decade - 1
        fa, fm = climo_file(case_dir, case, "ANN", y0, y1), climo_file(case_dir, case, mld_month, y0, y1)
        if not os.path.isfile(fa):
            if verbose:
                print(f"[skip] {os.path.basename(fa)} missing")
            continue
        flds, area = density_flux_fields(fa)
        rec = {v: _regional(flds[v], area, regions) for v in flds.data_vars}
        if os.path.isfile(fm):
            mld, am = mld_fields(fm)
            rec["MLD_mean"] = _regional(mld, am, regions)
            rec["MLD_max"], rec["conv_area"] = {}, {}
            for name, (la0, la1, lo0, lo1) in regions.items():
                m = mld.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
                a = am.sel(lat=slice(la0, la1), lon=slice(lo0, lo1))
                rec["MLD_max"][name] = float(m.max())
                rec["conv_area"][name] = float(a.where(m > conv_depth).sum()) / 1e9      # 10^3 km2
        rows.append((0.5 * (y0 + y1), rec))
        if verbose:
            print(f"  {case} {y0:04d}-{y1:04d} done")
    names = list(regions)
    keys = list(rows[0][1]) if rows else []
    data = {k: (("decade", "region"), np.array([[r[k].get(n, np.nan) for n in names] for _, r in rows]))
            for k in keys}
    out = xr.Dataset(data, coords={"decade": [d for d, _ in rows], "region": names})
    units = dict(Q_net="W m-2", F_fw="mm yr-1 (kg m-2 yr-1)", F_T="kg m-2 s-1", F_S="kg m-2 s-1", F_rho="kg m-2 s-1",
                 SST="degC", SSS="psu", sigma0="kg m-3", buoyancy="m2 s-3", MLD_mean="m", MLD_max="m", conv_area="10^3 km2")
    for k in out.data_vars:
        out[k].attrs["units"] = units.get(k, "")
    return out


def plot_region_panels(results, columns, *, segments=None, markers=(), fontz=11.0, figsize=None, save=None):
    """
    Rows = regions, columns = quantities; one line per run against the decade centre.

    results: {label: dict(ds=regional_decades() Dataset, color=...)}
    columns: [(title, [(variable, scale, linestyle), ...], unit)] — several variables in one column share the axis
    segments: optional {label: [(y0, y1), ...]} shaded per run
    """
    regions = list(next(iter(results.values()))["ds"]["region"].values)
    nr, nc = len(regions), len(columns)
    figsize = figsize or (4.6 * nc, 2.5 * nr)
    fig, axes = plt.subplots(nr, nc, figsize=figsize, sharex=True, squeeze=False, layout="constrained")
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    for i, reg in enumerate(regions):
        for j, (title, vars_, unit) in enumerate(columns):
            ax = axes[i, j]
            for lab, r in results.items():
                for v, scale, ls in vars_:
                    ax.plot(r["ds"]["decade"], r["ds"][v].sel(region=reg) * scale, ls, color=r["color"], lw=1.4, ms=3,
                            label=lab if (v, ls) == (vars_[0][0], vars_[0][2]) else None)
                for y0, y1 in (segments or {}).get(lab, []):
                    ax.axvspan(y0, y1, color=r["color"], alpha=0.07, lw=0)
            for yr in markers:
                ax.axvline(yr, color="0.4", ls=":", lw=1)
            ax.set_title(f"({next(letters)}) {reg}: {title}", fontsize=fontz, loc="left")
            ax.set_ylabel(unit, fontsize=0.9 * fontz)
            ax.tick_params(labelsize=0.85 * fontz)
            ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel("Model year (decade centre)", fontsize=0.9 * fontz)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="outside lower center", ncol=len(l), fontsize=fontz, frameon=False)
    if save:
        fig.savefig(save, dpi=300, bbox_inches="tight")
    return fig, axes
