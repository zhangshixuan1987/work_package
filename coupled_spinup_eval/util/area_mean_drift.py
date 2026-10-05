"""
Decadal drift / recoupling-shock metrics from EAM area-mean time series.

Input: the zppy-style ``post/time_series/atm/{VAR}_{YYYY}01_{YYYY}12.nc`` files, each holding
monthly area means ``VAR(time, rgn)`` over rgn = Global, Northern Hemisphere, Southern Hemisphere.
Because these are area means (no lat/lon), spatial PCC/RMSE cannot be formed; the metrics here are
the area-mean analogues, all measured against a fixed piControl climatology:

  * seasonal bias        (window mean - piControl mean) / |piControl mean| * 100        [%]
  * seasonal bias / σ    (window mean - piControl mean) / σ_dec                          [σ]
  * difference from base (window mean - base-exp window mean) / (√2 σ_dec)               [σ]
  * annual cycle         correlation and RMSE of the 12-month window climatology vs the piControl
                         climatology; RMSE also as % change from the piControl decadal RMSE level

σ_dec is the standard deviation of the (linearly detrended) non-overlapping window means of the
piControl run, i.e. the internal decadal variability for that variable, region and season.

Seasonal means use the December of the same year (ncclimo ``sdd``), so a window never reaches into
the year before it (e.g. the FOSI phase before a recoupling).
"""

import glob
import os
import re

import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Polygon, Rectangle

DPM = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], float)   # noleap
CUM_DAYS = np.concatenate([[0], np.cumsum(DPM)])
SEASON_MONTHS = {"DJF": [12, 1, 2], "MAM": [3, 4, 5], "JJA": [6, 7, 8], "SON": [9, 10, 11],
                 "ANN": list(range(1, 13))}
REGIONS = ["Global", "NH", "SH"]       # order of rgn in the files

# unit conversion applied on read: {var: (scale, new_units)}
UNIT_SCALE = {
    "PRECC":  (8.64e7, "mm/day"), "PRECL":  (8.64e7, "mm/day"),
    "PRECSC": (8.64e7, "mm/day"), "PRECSL": (8.64e7, "mm/day"),
    "QFLX":   (86400.0, "mm/day"),
    "PS":     (0.01, "hPa"), "PSL": (0.01, "hPa"),
}
# derived variables: {name: (function of dict of DataArrays, units, inputs)}
DERIVED = {
    "PRECT":  (lambda d: d["PRECC"] + d["PRECL"], "mm/day", ("PRECC", "PRECL")),
    "RESTOM": (lambda d: d["FSNT"] - d["FLNT"],   "W/m2",   ("FSNT", "FLNT")),
    # latent heat flux of evaporation, L_v * QFLX (QFLX already in mm/day = kg m-2 day-1)
    "LHFLX":  (lambda d: d["QFLX"] * (2.501e6 / 86400.0), "W/m2", ("QFLX",)),
}

_FNAME = re.compile(r"^(?P<var>.+)_(?P<y0>\d{4})01_(?P<y1>\d{4})12\.nc$")


# ======================================================================
# reading
# ======================================================================
_UNITS = re.compile(r"days since\s+(\d+)-(\d+)-(\d+)")


def _month_year(time_bnds, units):
    """(year, month) of each record from noleap time bounds in 'days since Y-M-D' units."""
    m = _UNITS.match(units.strip())
    if m is None:
        raise ValueError(f"unsupported time units: {units!r}")
    y, mo, d = (int(g) for g in m.groups())
    offset = (y - 1) * 365 + CUM_DAYS[mo - 1] + (d - 1)          # days since 0001-01-01
    mid = np.asarray(time_bnds, float).mean(axis=1) + offset
    year = (mid // 365).astype(int) + 1
    doy = mid % 365
    month = np.searchsorted(CUM_DAYS[1:], doy, side="right") + 1
    return year, month


def read_area_means(ts_dir, variables, years, derived=("PRECT", "RESTOM"), verbose=True):
    """
    Read monthly area means into one DataArray ``(variable, year, month, region)``.

    ``years`` = (first, last), inclusive. Missing months stay NaN (and are reported).
    ``derived`` names from DERIVED are added after the unit conversion.
    """
    y0, y1 = int(years[0]), int(years[1])
    yrs = np.arange(y0, y1 + 1)
    need = list(dict.fromkeys(list(variables) + [i for d in derived for i in DERIVED[d][2]]))
    out, units = {}, {}
    for v in need:
        arr = np.full((yrs.size, 12, len(REGIONS)), np.nan)
        files = []
        for f in glob.glob(os.path.join(ts_dir, f"{v}_*.nc")):
            m = _FNAME.match(os.path.basename(f))
            if m and m["var"] == v and int(m["y1"]) >= y0 and int(m["y0"]) <= y1:
                files.append(f)
        for f in sorted(files):
            with xr.open_dataset(f, decode_times=False) as ds:
                yr, mo = _month_year(ds["time_bnds"].values, ds["time"].attrs["units"])
                vals = ds[v].values
                units.setdefault(v, ds[v].attrs.get("units", ""))
            keep = (yr >= y0) & (yr <= y1)
            arr[yr[keep] - y0, mo[keep] - 1, :] = vals[keep]
        scale, new_units = UNIT_SCALE.get(v, (1.0, None))
        out[v] = arr * scale
        if new_units:
            units[v] = new_units
        nmiss = int(np.isnan(arr[..., 0]).sum())
        if nmiss and verbose:
            print(f"[WARNING] {v}: {nmiss} of {arr.shape[0] * 12} months missing in {ts_dir}")
    for d in derived:
        fn, u, _ = DERIVED[d]
        out[d] = fn(out)
        units[d] = u
    names = list(variables) + list(derived)
    da = xr.DataArray(
        np.stack([out[v] for v in names]),
        dims=("variable", "year", "month", "region"),
        coords={"variable": names, "year": yrs, "month": np.arange(1, 13), "region": REGIONS},
        name="area_mean",
    )
    da.attrs["units"] = "; ".join(f"{v}={units.get(v, '')}" for v in names)
    da = da.assign_coords(units=("variable", [units.get(v, "") for v in names]))
    return da


def load_or_read(cache_nc, ts_dir, variables, years, derived=("PRECT", "RESTOM"), force=False):
    """read_area_means() with a netCDF cache (rebuilt if variables/years differ or force=True)."""
    names = list(variables) + list(derived)
    if os.path.isfile(cache_nc) and not force:
        da = xr.open_dataarray(cache_nc).load()
        if (list(da["variable"].values) == names
                and int(da["year"][0]) == years[0] and int(da["year"][-1]) == years[1]):
            return da
    da = read_area_means(ts_dir, variables, years, derived=derived)
    os.makedirs(os.path.dirname(cache_nc), exist_ok=True)
    da.to_netcdf(cache_nc)
    return da


# ======================================================================
# metrics
# ======================================================================
def window_starts(years, length):
    """Start years of consecutive non-overlapping windows inside years=(first, last)."""
    return list(range(int(years[0]), int(years[1]) - length + 2, length))


def window_label(start, length):
    return f"{start:04d}-{start + length - 1:04d}"


def monthly_climo(da, start, length):
    """12-month climatology ``(variable, month, region)`` of years start..start+length-1."""
    return da.sel(year=slice(start, start + length - 1)).mean("year", skipna=False)


def seasonal_from_climo(clim, seasons):
    """Days-weighted seasonal means ``(variable, season, region)`` from a monthly climatology."""
    out = []
    for s in seasons:
        w = xr.DataArray(DPM[np.array(SEASON_MONTHS[s]) - 1], dims="month",
                         coords={"month": SEASON_MONTHS[s]})
        sub = clim.sel(month=SEASON_MONTHS[s])
        out.append((sub * w).sum("month", skipna=False) / w.sum())
    return xr.concat(out, dim=xr.DataArray(list(seasons), dims="season", name="season"))


def _ac_corr_rmse(clim, ref):
    """Annual-cycle (12-month) centered correlation and days-weighted RMSE along 'month'."""
    w = xr.DataArray(DPM, dims="month", coords={"month": np.arange(1, 13)})
    a = clim - clim.mean("month")
    b = ref - ref.mean("month")
    corr = (a * b).sum("month") / np.sqrt((a ** 2).sum("month") * (b ** 2).sum("month"))
    rmse = np.sqrt(((clim - ref) ** 2 * w).sum("month", skipna=False) / w.sum())
    return corr, rmse


def _detrended_std(x, dim):
    """Std (ddof=1) along dim after removing a least-squares linear trend."""
    n = x.sizes[dim]
    t = xr.DataArray(np.arange(n, dtype=float), dims=dim, coords={dim: x[dim]})
    t = t - t.mean()
    slope = (x - x.mean(dim)).dot(t, dims=dim) / (t ** 2).sum()
    resid = x - x.mean(dim) - slope * t
    return np.sqrt((resid ** 2).sum(dim) / (n - 1))


class AreaMeanDrift:
    """
    Window metrics for several experiments against one piControl reference.

    monthly    : {exp: DataArray(variable, year, month, region)} from read_area_means()
    reference  : DataArray(variable, year, month, region) of the piControl
    years      : (first, last) analysed in every experiment
    window     : window length (years)
    base_exp   : experiment the others are differenced against (e.g. the continuously coupled run)
    """

    def __init__(self, monthly, reference, years, window=10, seasons=("DJF", "MAM", "JJA", "SON", "ANN"),
                 base_exp=None):
        self.monthly = monthly
        self.reference = reference
        self.years = tuple(years)
        self.window = int(window)
        self.seasons = list(seasons)
        self.base_exp = base_exp
        self.starts = window_starts(self.years, self.window)
        self.labels = [window_label(s, self.window) for s in self.starts]
        self.ds = None

    def _reference_stats(self):
        ref = self.reference
        ref_clim = ref.mean("year", skipna=False)                       # (variable, month, region)
        ref_seas = seasonal_from_climo(ref_clim, self.seasons)
        # piControl non-overlapping windows -> internal decadal variability and RMSE noise level
        y0, y1 = int(ref["year"][0]), int(ref["year"][-1])
        pi_starts = window_starts((y0, y1), self.window)
        pi_seas, pi_rmse = [], []
        for s in pi_starts:
            c = monthly_climo(ref, s, self.window)
            pi_seas.append(seasonal_from_climo(c, self.seasons))
            pi_rmse.append(_ac_corr_rmse(c, ref_clim)[1])
        pi_seas = xr.concat(pi_seas, dim="pi_window")
        sigma = _detrended_std(pi_seas, "pi_window")
        rmse_pi = xr.concat(pi_rmse, dim="pi_window").mean("pi_window")
        return ref_clim, ref_seas, sigma, rmse_pi, len(pi_starts)

    def compute(self):
        ref_clim, ref_seas, sigma, rmse_pi, n_pi = self._reference_stats()
        exps = list(self.monthly)
        seas, corr, rmse = [], [], []
        for e in exps:
            s_e, c_e, r_e = [], [], []
            for s in self.starts:
                c = monthly_climo(self.monthly[e], s, self.window)
                s_e.append(seasonal_from_climo(c, self.seasons))
                cc, rr = _ac_corr_rmse(c, ref_clim)
                c_e.append(cc)
                r_e.append(rr)
            seas.append(xr.concat(s_e, dim="window"))
            corr.append(xr.concat(c_e, dim="window"))
            rmse.append(xr.concat(r_e, dim="window"))
        exp_dim = xr.DataArray(exps, dims="exp", name="exp")
        win = {"window": self.labels}
        mean = xr.concat(seas, dim=exp_dim).assign_coords(win)
        ac_corr = xr.concat(corr, dim=exp_dim).assign_coords(win)
        ac_rmse = xr.concat(rmse, dim=exp_dim).assign_coords(win)

        ds = xr.Dataset({
            "mean": mean,
            "ref_mean": ref_seas,
            "sigma_dec": sigma,
            "bias_pct": (mean - ref_seas) / np.abs(ref_seas) * 100.0,
            "bias_z": (mean - ref_seas) / sigma,
            "ac_corr": ac_corr,
            "ac_rmse": ac_rmse,
            "ac_rmse_pi": rmse_pi,
            "ac_rmse_ratio_pct": (ac_rmse - rmse_pi) / rmse_pi * 100.0,
        })
        if self.base_exp is not None:
            ds["diff_base_z"] = (mean - mean.sel(exp=self.base_exp)) / (np.sqrt(2.0) * sigma)
            ds.attrs["base_exp"] = self.base_exp
        ds.attrs.update({
            "window_years": self.window, "years": f"{self.years[0]}-{self.years[1]}",
            "n_picontrol_windows": n_pi,
            "seasons": "DJF uses December of the same year (sdd)",
        })
        self.ds = ds
        return ds

    def annual_means(self):
        """Days-weighted annual means {exp: DataArray(variable, year, region)}."""
        w = xr.DataArray(DPM, dims="month", coords={"month": np.arange(1, 13)})
        return {e: (m * w).sum("month", skipna=False) / w.sum() for e, m in self.monthly.items()}


# ======================================================================
# plotting
# ======================================================================
def _discrete_cmap(edges, cmap_name, whiten_zero=0):
    """Discrete colormap on bin edges; whiten_zero=n paints the n bins on each side of 0 white."""
    edges = np.asarray(edges, float)
    n = edges.size - 1
    cols = np.asarray(mpl.colormaps.get_cmap(cmap_name)(np.linspace(0, 1, n + 2)))
    bins = cols[1:-1].copy()
    if whiten_zero and edges[0] < 0.0 < edges[-1]:
        k0 = int(np.searchsorted(edges, 0.0, side="right") - 1)
        for k in range(k0 - int(whiten_zero), k0 + int(whiten_zero)):
            if 0 <= k < n:
                bins[k] = (1.0, 1.0, 1.0, 1.0)
    cmap = ListedColormap(bins)
    cmap.set_under(cols[0])
    cmap.set_over(cols[-1])
    cmap.set_bad((0.85, 0.85, 0.85, 1.0))
    return cmap, BoundaryNorm(edges, n, clip=False)


_TRI = {   # cell (j=x, i=y) split into four triangles, y axis pointing down
    "DJF": lambda j, i: [(j, i), (j + 1, i), (j + 0.5, i + 0.5)],
    "MAM": lambda j, i: [(j + 1, i), (j + 1, i + 1), (j + 0.5, i + 0.5)],
    "JJA": lambda j, i: [(j, i + 1), (j + 1, i + 1), (j + 0.5, i + 0.5)],
    "SON": lambda j, i: [(j, i), (j, i + 1), (j + 0.5, i + 0.5)],
}


def _draw_season_box(ax, fontz):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="0.25", linewidth=1.1))
    tri = {"DJF": [(0, 1), (1, 1), (0.5, 0.5)], "MAM": [(1, 0), (1, 1), (0.5, 0.5)],
           "JJA": [(0, 0), (1, 0), (0.5, 0.5)], "SON": [(0, 0), (0, 1), (0.5, 0.5)]}
    pos = {"DJF": (0.5, 0.8), "MAM": (0.78, 0.5), "JJA": (0.5, 0.2), "SON": (0.22, 0.5)}
    for s in tri:
        ax.add_patch(Polygon(tri[s], closed=True, facecolor="white", edgecolor="0.4", lw=1.0))
        ax.text(*pos[s], s, ha="center", va="center", fontsize=fontz)


def plot_window_heatmaps(
    da, *, exps, explabels, regions, region_titles, variables, fig_path,
    bins, cmap="RdBu_r", whiten_zero=1, cblabel="", title="",
    quad_seasons=True, fontz=16, x_step=1, figsize=None,
):
    """
    Grid of heatmaps: rows = regions, columns = experiments; x = windows, y = variables.

    da : DataArray with dims (exp, window, variable, region[, season]).
         quad_seasons=True draws DJF/MAM/JJA/SON as triangles; otherwise one square per cell
         (da must then have no season dim, or season must be selected beforehand).
    Missing values (e.g. a variable left out of a metric) are drawn grey.
    """
    windows = [str(w) for w in da["window"].values]
    V, W = len(variables), len(windows)
    nr, nc = len(regions), len(exps)
    if figsize is None:
        figsize = (4.0 + 0.55 * W * nc, 1.5 + 0.33 * V * nr)
    fig, axes = plt.subplots(nr, nc, figsize=figsize, squeeze=False, sharex=True, sharey=True)
    cm, norm = _discrete_cmap(bins, cmap, whiten_zero=whiten_zero)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cm)

    for r, reg in enumerate(regions):
        for c, e in enumerate(exps):
            ax = axes[r, c]
            sub = da.sel(exp=e, region=reg).sel(variable=variables)
            for j in range(W):
                for i in range(V):
                    if quad_seasons:
                        for s, tri in _TRI.items():
                            val = float(sub.sel(season=s).isel(window=j, variable=i))
                            fc = sm.to_rgba(val) if np.isfinite(val) else (0.85, 0.85, 0.85, 1.0)
                            ax.add_patch(Polygon(tri(j, i), closed=True, facecolor=fc,
                                                 edgecolor="grey", linewidth=0.05))
                    else:
                        val = float(sub.isel(window=j, variable=i))
                        fc = sm.to_rgba(val) if np.isfinite(val) else (0.85, 0.85, 0.85, 1.0)
                        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=fc, edgecolor="grey",
                                               linewidth=0.3))
            ax.set_xlim(0, W); ax.set_ylim(V, 0)
            ax.set_yticks(np.arange(V) + 0.5)
            ax.set_yticklabels(variables, fontsize=fontz * 0.8)
            ax.set_xticks(np.arange(0, W, x_step) + 0.5)
            ax.set_xticklabels(windows[::x_step], rotation=90, fontsize=fontz * 0.8)
            ax.tick_params(length=0)
            if r == 0:
                ax.set_title(explabels[c], fontsize=fontz * 1.05, pad=8)
            if c == 0:
                ax.set_ylabel(region_titles.get(reg, reg), fontsize=fontz)

    fig.subplots_adjust(left=0.08, right=0.86, bottom=0.10, top=0.93 if title else 0.96,
                        wspace=0.06, hspace=0.06)
    top, bot = axes[0, -1].get_position(), axes[-1, -1].get_position()
    h = (top.y1 - bot.y0) * 0.6
    y = bot.y0 + 0.5 * ((top.y1 - bot.y0) - h)
    cax = fig.add_axes([top.x1 + 0.015, y, 0.014, h])
    cb = fig.colorbar(sm, cax=cax, boundaries=bins, extend="both")
    cb.set_ticks(bins)
    cb.set_ticklabels([f"{v:g}" for v in bins])
    cb.ax.tick_params(labelsize=fontz * 0.8)
    cb.set_label(cblabel, fontsize=fontz)
    if quad_seasons:
        bax = fig.add_axes([top.x1 + 0.055, y + h + 0.01, 0.07, 0.07])
        _draw_season_box(bax, fontz * 0.7)
    if title:
        fig.suptitle(title, fontsize=fontz * 1.15)
    fig.savefig(fig_path, bbox_inches="tight")
    return fig


def plot_annual_timeseries(
    annual, *, exps, explabels, colors, variables, region, fig_path, reference=None,
    sigma_band=2.0, mark_year=None, ncols=3, fontz=13, figsize=None, var_titles=None,
    panel_labels=False, detrend_band=False,
):
    """
    Annual-mean area-mean time series, one panel per variable.

    annual    : {exp: DataArray(variable, year, region)} (AreaMeanDrift.annual_means())
    reference : DataArray(variable, year, region) of the piControl annual means; its mean ±
                sigma_band × interannual std (ddof=1) is shaded.
    detrend_band : remove a linear trend from the piControl annual means before taking the std.
    mark_year : vertical dashed line (e.g. the recoupling year).
    panel_labels : prefix titles with (a), (b), ... in reading order.
    """
    n = len(variables)
    nrows = int(np.ceil(n / ncols))
    if figsize is None:
        figsize = (5.2 * ncols, 3.0 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False, sharex=True)
    for k, v in enumerate(variables):
        ax = axes.flat[k]
        if reference is not None:
            ref = reference.sel(variable=v, region=region)
            m = float(ref.mean())
            sd = float(_detrended_std(ref, "year") if detrend_band else ref.std("year", ddof=1))
            ax.axhspan(m - sigma_band * sd, m + sigma_band * sd, color="0.88", zorder=0,
                       label=f"piControl mean ±{sigma_band:g}σ")
            ax.axhline(m, color="0.5", lw=0.8, zorder=0)
        for e, lab, col in zip(exps, explabels, colors):
            x = annual[e].sel(variable=v, region=region)
            ax.plot(x["year"], x, color=col, lw=1.2, label=lab)
        if mark_year is not None:
            ax.axvline(mark_year, color="k", ls="--", lw=0.9)
        units = str(annual[exps[0]]["units"].sel(variable=v).values) if "units" in annual[exps[0]].coords else ""
        label = f"({chr(ord('a') + k)}) " if panel_labels else ""
        ax.set_title(label + (var_titles or {}).get(v, v) + (f" [{units}]" if units else ""), fontsize=fontz)
        ax.tick_params(labelsize=fontz * 0.85)
        ax.ticklabel_format(axis="y", useOffset=False)
        ax.grid(alpha=0.3)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("Model year", fontsize=fontz)
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=fontz, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(fig_path, bbox_inches="tight")
    return fig
