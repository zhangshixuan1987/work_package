"""2-D climatology maps from MPAS-Analysis output: surface fields (SST, SSS) and the
AMOC streamfunction (latitude-depth), with reference + difference panels.

Options / separate methods:
  Surface2DClimoPlotter.plot_map_compare(ref_cmap_n, ref_darken, diff_cmap_n)
      SST: 200, 1.0, 200 (defaults)   SSS: 256, 0.8, 256
  AMOCClimoPlotter.plot_amoc_2d_compare       contour panels
  AMOCClimoPlotter.plot_amoc_2d_compare_mesh  pcolormesh panels
"""
from typing import Optional, Tuple, List
import os, re, glob, warnings, string
import numpy as np
import xarray as xr
import cftime
import json
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.ticker import (
    FixedLocator, FixedFormatter, AutoMinorLocator,
    MultipleLocator, NullFormatter
)
from matplotlib.lines import Line2D
import statsmodels.api as sm
import cmocean
from scipy import stats, signal
import matplotlib as mpl
from matplotlib.colors import BoundaryNorm, TwoSlopeNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cm as cm
import warnings
import string
import os, re, glob
from matplotlib.ticker import (
    FixedLocator, FixedFormatter, AutoMinorLocator,
    MultipleLocator, NullFormatter, FuncFormatter
)
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import TwoSlopeNorm


# ------------------- significance tests (shared by the plotters) -------------------
def control_ttest_p(d: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Two-sided p-values of a difference d between two window means, with the variance of
    one window mean taken from the control samples c (chunk axis 0): t = d / (s_c sqrt 2),
    df = n_chunks - 1.
    """
    s = np.nanstd(c, axis=0, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = d / (s * np.sqrt(2.0))
    p = 2.0 * stats.t.sf(np.abs(t), c.shape[0] - 1)
    p[~np.isfinite(t)] = np.nan
    return p


def control_resample_p(d: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Resampling p-values: |d| against the |differences| of all pairs of control samples,
    p = (1 + #{|null| >= |d|}) / (1 + n_pairs). Resolution 1/(1 + n_pairs).
    """
    i, j = np.triu_indices(c.shape[0], 1)
    null = np.abs(c[i] - c[j])
    with np.errstate(invalid="ignore"):
        p = (1.0 + np.sum(null >= np.abs(d)[None], axis=0)) / (1.0 + i.size)
    p[~np.isfinite(d) | ~np.all(np.isfinite(null), axis=0)] = np.nan
    return p


def fdr_significant(p: np.ndarray, alpha_fdr: float) -> np.ndarray:
    """Benjamini-Hochberg false-discovery-rate control over all finite p-values (Wilks 2016)."""
    ok = np.isfinite(p)
    pv = np.sort(p[ok])
    if pv.size == 0:
        return np.zeros_like(p, dtype=bool)
    below = pv <= alpha_fdr * np.arange(1, pv.size + 1) / pv.size
    p_crit = pv[below].max() if below.any() else -1.0
    return ok & (p <= p_crit)


def climo_window_error(clim_nc: str, ts_nc: str, years: Tuple[int, int]) -> float:
    """
    Check that an MPAS-Analysis climatology really averages `years`: max |mocAtlantic of the MOC
    climatology `clim_nc` (mocStreamfunction_years*.nc) - days-weighted mean of the monthly MOC
    time series `ts_nc` over `years`| in Sv. ~0 for a correct climatology; a climatology that
    also averaged other months (e.g. year 1) is off by several 0.1 Sv. All variables of an
    MPAS-Analysis climatology come from the same averaging step, so this checks SST/SSS too.
    """
    dpm = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=float)
    with xr.open_dataset(clim_nc) as dc, xr.open_dataset(ts_nc, decode_times=False) as dt:
        clim = dc["mocAtlantic"].values
        sel = (dt["year"].values >= years[0]) & (dt["year"].values <= years[1])
        if sel.sum() != 12 * (years[1] - years[0] + 1):
            raise ValueError(f"{ts_nc} does not hold all months of {years}")
        w = dpm[dt["month"].values[sel] - 1]
        ts = np.tensordot(w, dt["mocAtlantic"].values[sel][:, :clim.shape[0]], axes=(0, 0)) / w.sum()
    return float(np.nanmax(np.abs(clim - ts)))


def stipple(ax, x, y, mask, *, grid=(72, 36), size=4.0, color="k", transform=None, zorder=1.5):
    """
    Mark mask==True with dots on a regular grid of grid[0] x grid[1] positions spanning the data
    range (each dot takes the mask value of the nearest data point), so every panel gets the same
    dot density whatever its data grid (0.5° maps, 1° x 65-level MOC). Dots are plain markers
    (size in points^2): unlike hatch patterns they scale with the figure when it is included at a
    different size (e.g. by LaTeX), and unlike filled mask contours they cannot fill the globe.
    x, y: 1-D coordinates of mask's columns / rows (any monotonic order).
    """
    x = np.asarray(x, float); y = np.asarray(y, float); mask = np.asarray(mask, bool)
    def nearest(coord, n):
        lo, hi = np.nanmin(coord), np.nanmax(coord)
        pos = lo + (np.arange(n) + 0.5) * (hi - lo) / n
        return pos, np.abs(coord[None, :] - pos[:, None]).argmin(axis=1)
    xs, ix = nearest(x, int(grid[0]))
    ys, iy = nearest(y, int(grid[1]))
    m = mask[np.ix_(iy, ix)]
    if not m.any():
        return None
    X, Y = np.meshgrid(xs, ys)
    kw = {} if transform is None else {"transform": transform}
    return ax.scatter(X[m], Y[m], s=size, c=color, marker="o", linewidths=0, zorder=zorder, **kw)


def fit_cbar_to_panel(ax, cax, height_frac: float = 1.0, pad: float = 0.006):
    """
    Place colorbar axes `cax` right of panel `ax` (gap `pad`, figure fraction) with height_frac of the
    panel's drawn height. Uses the position after the panel's aspect is applied (maps, set_box_aspect),
    so the bar matches what is drawn rather than the gridspec slot.
    """
    ax.apply_aspect()
    pb, cb = ax.get_position(), cax.get_position()
    h = pb.height * float(height_frac)
    cax.set_position([pb.x1 + pad, pb.y0 + 0.5 * (pb.height - h), cb.width, h])


def arrange_panel_grid(
    fig,
    rows,
    *,
    width: float = 24.0,
    aspect: float = 0.5,
    left: float = 1.2,
    right: float = 1.6,
    top: float = 0.6,
    bottom: float = 0.9,
    col_gaps=(1.6, 0.35),
    row_gap: float = 1.25,
    cbar_width: float = 0.2,
    cbar_pad: float = 0.12,
    cbar_frac: float = 0.95,
):
    """
    Put rows of panels on one grid, in inches: every panel has the same width and shape
    (height = aspect * width; 0.5 = global lat-lon maps), columns line up across rows, and each
    visible colorbar sits cbar_pad right of its panel with cbar_frac of the panel height.
    rows: list of (axes, caxes) per row (plotter.last_axes after a plot call), same length each.
    col_gaps: space after each column but the last (room for that column's colorbar and labels);
    right: space after the last column. The figure size is set to fit (width fixed).
    """
    ncol = len(rows[0][0])
    gaps = list(col_gaps) + [right]
    if len(gaps) != ncol:
        raise ValueError(f"col_gaps needs {ncol - 1} values for {ncol} columns")
    w = (width - left - sum(gaps)) / ncol
    h = aspect * w
    height = top + len(rows) * h + (len(rows) - 1) * row_gap + bottom
    fig.set_size_inches(width, height)
    for r, (axes, caxes) in enumerate(rows):
        y0 = height - top - (r + 1) * h - r * row_gap
        x0 = left
        for k, (ax, cax) in enumerate(zip(axes, caxes)):
            ax.set_position([x0 / width, y0 / height, w / width, h / height])
            hc = cbar_frac * h
            cax.set_position([(x0 + w + cbar_pad) / width, (y0 + 0.5 * (h - hc)) / height,
                              cbar_width / width, hc / height])
            x0 += w + gaps[k]
    return fig


class Surface2DClimoPlotter:
    """
    2-D map climatology/statistics + plotting for fields like SST/SSS.

    - Works with regular (lat, lon) grids directly.
    - If MPAS-native (nCells, latCell/lonCell), bins to a regular lat-lon grid via histogram2d.
    - Computes monthly climatology, annual mean, seasonal amplitude, and month-of-max.
    - Can produce Reference | (Tgt-Ref) | (Tgt2-Ref) comparison maps with discrete colorbars.
    """

    def __init__(
        self,
        root_dir: str,
        sub_dir: str,
        file_key: str = "mpasTimeSeriesOcean",
        var_aliases: Optional[List[str]] = None,
        scale: float = 1.0,
        weighted: bool = False,                 # days-in-month weights
        engine: Optional[str] = "netcdf4",
        chunks: Optional[dict] = None,          # e.g., {"Time": 120}
        use_mfdataset: bool = True,
        combine: str = "nested",
        open_parallel: bool = True,
        time_origin: Optional[str] = None,      # fallback if needed
        calendar: str = "noleap",
        dtype: str = "float32",
        target_res: Tuple[float, float] = (1.0, 1.0),   # deg lat/lon for MPAS binning
        area_weighted: bool = True,                     # use cellArea if present
        use_cartopy: bool = True,                       # requires cartopy
    ):
        self.root_dir = root_dir
        self.sub_dir = sub_dir
        self.file_key = file_key
        self.scale = float(scale)
        self.weighted = bool(weighted)
        self.engine = engine
        self.use_mfdataset = bool(use_mfdataset)
        self.combine = combine
        self.open_parallel = bool(open_parallel)
        self.time_origin = time_origin
        self.calendar = calendar
        self.dtype = dtype
        self.chunks = chunks if chunks is not None else {"Time": 120}
        self.target_res = target_res
        self.area_weighted = bool(area_weighted)
        self.use_cartopy = bool(use_cartopy)

        # common surface aliases (customize as needed)
        self.var_aliases = var_aliases or [
            # SST
            "timeMonthly_avg_activeTracers_surfaceTemperature",
            "timeMonthly_avg_surfaceTemperature",
            "SST", "tos", "sea_surface_temperature",
            # SSS
            "timeMonthly_avg_activeTracers_surfaceSalinity",
            "timeMonthly_avg_surfaceSalinity",
            "SSS", "sos", "sea_surface_salinity",
        ]

        self._meta: List[dict] = []
        self._stats: List[dict] = []
        self._control: Optional[xr.DataArray] = None   # control-run samples (add_control_climo)

    # -------------------- helpers --------------------
    def _find_files(self, exp: str, year_token: Optional[str] = None) -> List[str]:
        """
        Find files for this experiment. Behavior:
          - If file_key ends with '.nc' => treat as exact filename.
          - If year_token is missing/empty/unparseable => look for 'file_key.nc'.
          - Else, search yearly-sliced pattern 'file_key_YYYY-YYYY.nc', with fallback to 'file_key.nc'.
        """
        base = os.path.join(self.root_dir, exp, self.sub_dir or "")

        # 1) Exact filename provided
        if self.file_key.lower().endswith(".nc"):
            hits = glob.glob(os.path.join(base, self.file_key))
            hits.sort()
            return hits

        # 2) No/empty/unparseable year_token -> just use 'file_key.nc'
        if not year_token:
            hits = glob.glob(os.path.join(base, f"{self.file_key}.nc"))
            hits.sort()
            return hits

        try:
            a, b = [int(t.strip()) for t in str(year_token).split("-")]
            if a > b:
                a, b = b, a
        except Exception:
            # Fallback if token can't be parsed
            hits = glob.glob(os.path.join(base, f"{self.file_key}.nc"))
            hits.sort()
            return hits

        # 3) Search yearly-sliced files; fallback to single file
        hits: List[str] = []
        for yy in range(a, b + 1):
            pat = os.path.join(base, f"{self.file_key}_{yy:04d}-{yy:04d}.nc")
            hits.extend(glob.glob(pat))

        if not hits:
            hits = glob.glob(os.path.join(base, f"{self.file_key}.nc"))

        hits.sort()
        return hits

    def _guess_var(self, ds: xr.Dataset) -> Optional[str]:
        for nm in self.var_aliases:
            if nm in ds.data_vars:
                return nm
        # sometimes variables are tucked in variables but not as data_vars
        for nm in self.var_aliases:
            if nm in ds.variables:
                return nm
        return None

    def _decode_time_inplace(self, ds: xr.Dataset) -> xr.Dataset:
        # Try common MPAS time decode paths
        if "Time" in ds.dims and "time" not in ds.coords:
            # CF units on Time
            if "Time" in ds and "units" in ds["Time"].attrs:
                try:
                    coder = xr.coders.CFDatetimeCoder(use_cftime=True)
                    t = xr.coding.times.decode_cf_datetime(
                        ds["Time"].values,
                        units=ds["Time"].attrs.get("units"),
                        calendar=ds["Time"].attrs.get("calendar", self.calendar),
                        time_coder=coder,   # <- replaces use_cftime=True
                    )
                    ds = ds.assign_coords(time=("Time", t)).swap_dims({"Time": "time"})
                    return ds
                except Exception:
                    pass
    
            # MPAS xtime chars
            for nm in ("xtime_startMonthly", "xtime_endMonthly", "startTime", "endTime"):
                if nm in ds:
                    arr = ds[nm].values
                    try:
                        s = ["".join(r.astype(str)).strip().replace("\x00", "") for r in arr]
                    except Exception:
                        continue
                    tt = []
                    ok = True
                    for txt in s:
                        txt = txt.replace("T", " ").replace("_", " ")
                        try:
                            y, m, d = int(txt[0:4]), int(txt[5:7]), int(txt[8:10])
                            tt.append(cftime.DatetimeNoLeap(y, m, d))
                        except Exception:
                            ok = False
                            break
                    if ok and len(tt) == ds.dims["Time"]:
                        ds = ds.assign_coords(time=("Time", np.asarray(tt, dtype=object))).swap_dims({"Time": "time"})
                        return ds
    
        if "time" in ds.dims and "time" not in ds.coords:
            ds = ds.assign_coords(time=("time", np.arange(ds.sizes["time"], dtype=int)))
        return ds

    def _open_lazy(self, files: List[str]) -> xr.Dataset:
        time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
        xr_opts = dict(engine=self.engine, chunks=self.chunks, decode_times=time_coder)
    
        if (len(files) == 1) or (not self.use_mfdataset):
            ds = xr.open_dataset(files[0], **xr_opts)
            ds = self._decode_time_inplace(ds)
            return ds
    
        try:
            ds = xr.open_mfdataset(
                files,
                combine=self.combine,
                concat_dim=("time" if self.combine == "nested" else None),
                data_vars="minimal",
                coords="minimal",
                compat="override",
                join="override",
                combine_attrs="drop",
                parallel=self.open_parallel,
                **xr_opts,
            )
            ds = self._decode_time_inplace(ds)
            return ds
    
        except Exception:
            parts = []
            try:
                for fp in files:
                    dso = xr.open_dataset(fp, **xr_opts)
                    dso = self._decode_time_inplace(dso)
                    parts.append(dso)
    
                ds = xr.concat(parts, dim="time", join="override")
    
                # Detach from on-disk file handles before closing
                ds.load()
                return ds
    
            finally:
                for dso in parts:
                    try:
                        dso.close()
                    except Exception:
                        pass

    def _is_regular_latlon(self, ds: xr.Dataset, var: str) -> bool:
        da = ds[var]
        dims = list(da.dims)
        has_lat = any("lat" == d.lower() or d.lower().endswith("lat") for d in dims)
        has_lon = any("lon" == d.lower() or d.lower().endswith("lon") for d in dims)
        return has_lat and has_lon and ("nCells" not in dims)

    def _to_regular_latlon(
        self, ds: xr.Dataset, var: str, res: Tuple[float, float], area_weighted: bool = True
    ) -> xr.DataArray:
        """
        If (lat, lon) already present, return data.
        If MPAS-native (nCells), bin-average to regular lat-lon grid at given resolution.
        """
        da = ds[var] * self.scale

        # Case 1: Already lat-lon
        if self._is_regular_latlon(ds, var):
            # normalize coordinate names to 'lat', 'lon'
            lat_name = next(d for d in da.dims if "lat" in d.lower())
            lon_name = next(d for d in da.dims if "lon" in d.lower())
            if lat_name != "lat":
                da = da.rename({lat_name: "lat"})
            if lon_name != "lon":
                da = da.rename({lon_name: "lon"})
            return da

        # Case 2: MPAS native (nCells, latCell/lonCell)
        if ("nCells" in da.dims) and (("latCell" in ds) and ("lonCell" in ds)):
            latc = np.degrees(ds["latCell"].values) if np.abs(ds["latCell"]).max() <= np.pi else ds["latCell"].values
            lonc = np.degrees(ds["lonCell"].values) if np.abs(ds["lonCell"]).max() <= np.pi else ds["lonCell"].values
            # wrap lon to [-180, 180]
            lonc = ((lonc + 180.0) % 360.0) - 180.0

            dlat, dlon = float(res[0]), float(res[1])
            lat_edges = np.arange(-90.0, 90.0 + dlat, dlat)
            lon_edges = np.arange(-180.0, 180.0 + dlon, dlon)
            lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
            lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])

            # optional area weights
            if area_weighted and ("areaCell" in ds):
                wcell = ds["areaCell"].values
            else:
                wcell = np.ones_like(lonc, dtype=float)

            # We bin each time slice separately to avoid massive memory
            time_dim = "time" if "time" in da.dims else da.dims[0]
            ntime = da.sizes[time_dim]
            out = np.full((ntime, lat_centers.size, lon_centers.size), np.nan, dtype=self.dtype)

            # Pull values to host as needed (dask-friendly iteration)
            for t in range(ntime):
                vals = (da.isel({time_dim: t}).values).astype(float)  # (nCells,)
                good = np.isfinite(vals) & np.isfinite(lonc) & np.isfinite(latc)
                if not np.any(good):
                    continue
                v = vals[good]
                wx = wcell[good]
                y = latc[good]
                x = lonc[good]

                # weighted average in bins: sum(w*v)/sum(w)
                num, _, _ = np.histogram2d(y, x, bins=[lat_edges, lon_edges], weights=wx * v)
                den, _, _ = np.histogram2d(y, x, bins=[lat_edges, lon_edges], weights=wx)
                with np.errstate(invalid="ignore", divide="ignore"):
                    out[t, :, :] = num / den

            out_da = xr.DataArray(
                out,
                dims=("time", "lat", "lon"),
                coords=dict(
                    time=da["time"].values if "time" in da.coords else np.arange(ntime),
                    lat=lat_centers, lon=lon_centers
                ),
                name=var
            )
            return out_da

        raise ValueError("Could not detect regular lat–lon or MPAS-native (nCells) layout.")

    def _cmap_with_over_under(self, base_cmap, n: int = 256, darken_factor: float = 1.0):
        def _darken(color, f=0.8):
            r, g, b, a = mpl.colors.to_rgba(color)
            return (r*f, g*f, b*f, a)
        cmap = mpl.colormaps.get_cmap(base_cmap) if isinstance(base_cmap, str) else base_cmap
        colors1 = cmap(np.linspace(0.05, 0.95, n+21))
        new = mpl.colors.LinearSegmentedColormap.from_list(f"{cmap.name}_mod", colors1[10:n+10], N=n)
        new.set_under(_darken(colors1[0], darken_factor))
        new.set_over(_darken(colors1[-1], darken_factor))
        new.set_bad("1.0")
        return new

    def _symmetric_levels(self, vmin, vmax, *, n_levels=None, step=None):
        a = float(max(abs(vmin), abs(vmax)))
        if step is not None:
            a = step * np.ceil(a / step)
            levels = np.arange(-a, a + 0.5*step, step, dtype=float)
        else:
            n = int(n_levels) if n_levels else 13
            if n % 2 == 0:
                n += 1
            levels = np.linspace(-a, a, n, dtype=float)
        zidx = np.argmin(abs(levels))
        levels[zidx] = 0.0
        return levels
    
    def _with_cyclic(self, lon, Z):
        Zc, lonc = add_cyclic_point(Z, coord=np.asarray(lon))
        return lonc, Zc
    
    def _bin_centers(self, levels):
        levels = np.asarray(levels, float)
        return 0.5*(levels[:-1] + levels[1:])
    
    # -------------------- public API --------------------
    def add_experiment(self, exp: str, var_prefer: Optional[str] = None,
                       label: Optional[str] = None, year_token: str = "0001-0030",
                       sub_dir: Optional[str] = None, clim_start: Optional[str] = None,
                       clim_end: Optional[str] = None):
        if sub_dir is not None:
            old = self.sub_dir
            self.sub_dir = sub_dir

        files = self._find_files(exp, year_token)
        if not files:
            raise FileNotFoundError(f"No files under {os.path.join(self.root_dir, exp, self.sub_dir)}")

        ds = self._open_lazy(files)
        var = var_prefer or self._guess_var(ds)
        if var is None:
            raise KeyError(f"None of var_aliases found in dataset. Tried: {self.var_aliases}")

        da_ll = self._to_regular_latlon(ds, var, res=self.target_res, area_weighted=self.area_weighted).astype(self.dtype)

        # NEW: if there is no time dimension, assume this is already an annual climatology (e.g., mpaso_ANN_*_climo.nc)
        if "time" not in da_ll.dims:
            stats = dict(
                monthly_climatology=None,
                annual_mean=da_ll,                 # <-- direct field (lat, lon)
                seasonal_amplitude=None,
                month_of_max=None,
                _var_name=var
            )
        else:
            # (existing monthly/annual climatology code stays the same)
            has_dt = "time" in da_ll.coords
            if clim_start or clim_end:
                if has_dt:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        da_ll = da_ll.sel(time=slice(clim_start, clim_end))
                else:
                    a, b = [int(t) for t in year_token.split("-")]
                    ys = int(str(clim_start)[:4]) if clim_start else a
                    ye = int(str(clim_end)[:4])   if clim_end   else b
                    m0 = max(0, (ys - a) * 12)
                    m1 = min(da_ll.sizes["time"], (ye - a + 1) * 12)
                    da_ll = da_ll.isel(time=slice(m0, m1))

            if "time" in da_ll.dims:
                if has_dt:
                    if self.weighted:
                        dpm = da_ll["time"].dt.days_in_month.astype("float")
                        w = dpm / dpm.groupby("time.month").sum()
                        w_al = xr.DataArray(w, dims=["time"]).broadcast_like(da_ll)
                        clim12 = (da_ll * w_al).groupby("time.month").sum("time", skipna=True)
                        dpm_m = dpm.groupby("time.month").mean().astype("float")
                        w_m = xr.DataArray(dpm_m / dpm_m.sum(), dims=["month"]).broadcast_like(clim12)
                        annual_mean = (clim12 * w_m).sum("month", skipna=True)
                    else:
                        clim12 = da_ll.groupby("time.month").mean("time", skipna=True)
                        annual_mean = clim12.mean("month", skipna=True)
                else:
                    n = da_ll.sizes["time"]
                    months = xr.DataArray((np.arange(n) % 12) + 1, dims=["time"], name="month")
                    clim12 = da_ll.groupby(months).mean("time", skipna=True)
                    clim12 = clim12.rename({"group": "month"}) if "group" in clim12.dims else clim12
                    annual_mean = clim12.mean("month", skipna=True)

                seasonal_amplitude = clim12.max("month") - clim12.min("month")
                month_of_max = (clim12.argmax("month") + 1).astype("int16")
                stats = dict(
                    monthly_climatology=clim12,
                    annual_mean=annual_mean,
                    seasonal_amplitude=seasonal_amplitude,
                    month_of_max=month_of_max,
                    _var_name=var
                )

        self._meta.append({
            "exp": exp,
            "label": label or exp,
            "year_token": year_token,
            "clim_start": clim_start,
            "clim_end": clim_end,
        })
        self._stats.append(stats)

        if sub_dir is not None:
            self.sub_dir = old

    # -------------------- plotting --------------------
    def add_experiment_mean(
        self,
        exp: str,
        files: List[str],
        var: str,
        label: Optional[str] = None,
        level: int = 0,
        weights: Optional[List[float]] = None,
    ):
        """
        Add an experiment as the (weighted) mean of several lat-lon annual climatologies, e.g.
        five decadal climatologies for a 50-yr window. files: paths under <root_dir>/<exp>;
        var: variable in them (level `level` of nVertLevels is taken if present);
        weights: e.g. window lengths (default: equal).
        """
        fields = []
        for name in files:
            path = os.path.join(self.root_dir, exp, name)
            if not os.path.isfile(path):
                raise FileNotFoundError(path)
            with xr.open_dataset(path, engine=self.engine) as ds:
                da = ds[var]
                if "nVertLevels" in da.dims:
                    da = da.isel(nVertLevels=level)
                fields.append(da.squeeze(drop=True).load())
        w = np.ones(len(fields)) if weights is None else np.asarray(weights, float)
        annual_mean = (sum(wk * f for wk, f in zip(w, fields)) / w.sum()).astype(self.dtype)
        self._meta.append({"exp": exp, "label": label or exp, "year_token": None,
                           "clim_start": None, "clim_end": None})
        self._stats.append(dict(monthly_climatology=None, annual_mean=annual_mean,
                                seasonal_amplitude=None, month_of_max=None, _var_name=var))

    def add_control_climo(
        self,
        exp: str,
        files: List[str],
        sub_dir: str,
        var: str,
    ):
        """
        Load non-overlapping annual climatologies of a control run (e.g. the ten 50-yr windows
        of piControl) as samples of internal variability for plot_map_compare(sig_test=True).
        files: file names under <root_dir>/<exp>/<sub_dir> (sub_dir may contain a glob);
        var: variable in those files, on the same lat-lon grid as the experiments.
        """
        chunks = []
        for name in files:
            hits = sorted(glob.glob(os.path.join(self.root_dir, exp, sub_dir, name)))
            if not hits:
                raise FileNotFoundError(f"{name} under {exp}/{sub_dir}")
            with xr.open_dataset(hits[0], engine=self.engine) as ds:
                chunks.append((ds[var].squeeze(drop=True) * self.scale).astype(self.dtype).load())
        self._control = xr.concat(chunks, dim="chunk")   # (chunk, lat, lon)

    def _resolve_field(self, idx: int, field: str, month: Optional[int]):
        stat = self._stats[idx]
        if field not in stat:
            raise KeyError(f"Field '{field}' not found. Options: {list(stat.keys())}")
        da = stat[field]
        if field == "monthly_climatology":
            if month is None:
                raise ValueError("Provide month=1..12 for 'monthly_climatology'.")
            da = da.sel(month=int(month))
        return da

    def _axes(self, ncols, figsize, wspace, cbar_width=0.018, fig=None, slot=None):
        """slot: a SubplotSpec of `fig` to draw the row into (multi-row figures); else a new figure."""
        if fig is None:
            fig = plt.figure(figsize=figsize)
        # layout: (panel|cbar) × ncols
        width_ratios = []
        for _ in range(ncols):
            width_ratios += [1.0, cbar_width]

        grid = slot.subgridspec if slot is not None else fig.add_gridspec
        gs = grid(1, ncols*2, wspace=wspace, width_ratios=width_ratios)

        axes, caxes = [], []
        for i in range(ncols):
            ax = fig.add_subplot(gs[0, 2*i], projection=ccrs.PlateCarree())
            ax.set_global()
            ax.coastlines(linewidth=0.6, zorder=3)
            ax.add_feature(cfeature.LAND, zorder=2, edgecolor="none", facecolor="0.9")
            axes.append(ax)
            caxes.append(fig.add_subplot(gs[0, 2*i + 1]))  # dedicated cbar slot for this panel
        return fig, axes, caxes

    def _add_geophysical_lines(
        self, ax, *, gridlines=False, xlocs=None, ylocs=None,
        equator=True, tropics=True, polar_circles=False,
        prime_meridian=False, date_line=False,
    ):
        if xlocs is None: xlocs = np.arange(-180, 181, 60)
        if ylocs is None: ylocs = np.arange(-90,  91, 30)
        pc = ccrs.PlateCarree()

        if gridlines:
            ax.gridlines(crs=pc, xlocs=xlocs, ylocs=ylocs,
                         linewidth=0.4, color="0.6", alpha=0.5,
                         linestyle="--", draw_labels=False)

        lons = np.linspace(-180, 180, 721)
        def _plot_lat(lat, **kw): ax.plot(lons, np.full_like(lons, lat), transform=pc, **kw)
        def _plot_lon(lon, **kw): ax.plot(np.full_like(lons, lon), np.linspace(-90, 90, lons.size), transform=pc, **kw)

        if equator:       _plot_lat(0.0,      color="k", lw=0.8, ls="-",  alpha=0.8,  zorder=4)
        if tropics:
            _plot_lat(+23.436, color="k", lw=0.6, ls=":",  alpha=0.85, zorder=4)
            _plot_lat(-23.436, color="k", lw=0.6, ls=":",  alpha=0.85, zorder=4)
        if polar_circles:
            _plot_lat(+66.563, color="k", lw=0.6, ls="--", alpha=0.8,  zorder=4)
            _plot_lat(-66.563, color="k", lw=0.6, ls="--", alpha=0.8,  zorder=4)
        if prime_meridian: _plot_lon(0.0,     color="k", lw=0.6, ls="--", alpha=0.7,  zorder=4)
        if date_line:      _plot_lon(180.0,   color="k", lw=0.6, ls="--", alpha=0.7,  zorder=4)
            
    def _shrink_cbar_axes(self, fig, cax, height_frac=0.9, align="center"):
        """
        Return a new colorbar axes with the same x-position as `cax` but a shorter height.
        height_frac in (0,1]; align is one of {"center","top","bottom"}.
        """
        height_frac = max(0.05, min(1.0, float(height_frac)))
        bb = cax.get_position()
        new_h = bb.height * height_frac
        if align == "top":
            y0 = bb.y0 + (bb.height - new_h)
        elif align == "bottom":
            y0 = bb.y0
        else:  # center
            y0 = bb.y0 + 0.5*(bb.height - new_h)
        new_cax = fig.add_axes([bb.x0, y0, bb.width, new_h])
        cax.remove()
        return new_cax
    
    def _font_sizes(self, fontz: float):
        """
        Scale all typography from one knob: `fontz` (interpreted as the base tick size).
        Tuned so things look balanced for most figures.
        """
        f = float(fontz)
        return {
            "tick":        f,        # major tick labels
            "axis_label":  f*1.0,   # "Longitude", "Latitude"
            "title":       f*1.1,   # panel titles
            "cbar_tick":   f*0.95,   # colorbar tick labels
            "cbar_label":  f*0.95,   # colorbar title
        }
    
    def _dock_cbar(self, fig, ax_panel, cax, *, pad=0.006, side="right"):
        """
        Move `cax` so it hugs `ax_panel` with a small gap.
        `pad` is in figure fraction (0–1). side = "right"|"left".
        """
        pb = ax_panel.get_position()      # panel bbox in figure coords
        cb = cax.get_position()           # current cbar bbox
        new_x0 = (pb.x1 + pad) if side == "right" else (pb.x0 - pad - cb.width)
        cax.set_position([new_x0, cb.y0, cb.width, cb.height])
        return cax
    
    def plot_map_compare(
        self,
        var_label: str,
        ref_idx: int,
        tgt_idx,
        fig_idx: int=0,
        field: str = "annual_mean",
        month: Optional[int] = None,
        cmap: str = "viridis",
        diff_cmap: str = "bwr",
        vmin=None, vmax=None,
        diff_vmin=None, diff_vmax=None,
        clevs: Optional[List[float]] = None,
        n_levels: Optional[int] = 13,
        diff_clevs: Optional[List[float]] = None,
        diff_n_levels: Optional[int] = 13,
        extend: str = "both",
        diff_extend: str = "both",
        figsize=(13.2, 6.0),
        wspace: float = 0.25,
        cbar_width: float = 0.03,     # thinner/thicker bars
        cbar_height: float = 0.55,     # 0–1 fraction of panel height
        cbar_align: str = "center",    # "center" | "top" | "bottom"
        cbar_pad: float = 0.006,   # ← distance between panel and colorbar (figure fraction)
        cbar_tick_every: int = 0,   # optional: thin out ticks on both CBs
        fontz: float = 12.0,
        cbar_label: str = "",
        diff_label: str = "Δ",
        save: Optional[str] = None,
        tick_every: Optional[int] = None,
        diff_tick_every: Optional[int] = None,
        # geophysical lines
        gridlines: bool = False,
        grid_xlocs=None,
        grid_ylocs=None,
        equator: bool = False,
        tropics: bool = False,
        polar_circles: bool = False,
        prime_meridian: bool = False,
        date_line: bool = False,
        # colormap resolution / end-color darkening (sst_2d defaults; sss_2d used 256, 0.8, 256)
        ref_cmap_n: int = 200,
        ref_darken: float = 1.0,
        diff_cmap_n: int = 200,
        # draw into a row of an existing figure (fig + slot = a SubplotSpec of it); figsize is then ignored
        fig=None,
        slot=None,
        # layout: "shared" = one colorbar for all difference panels (on the last one);
        # cbar_fit_panel = colorbar height as a fraction of the drawn panel height (None: cbar_height)
        diff_cbar: str = "each",
        cbar_fit_panel: Optional[float] = None,
        # significance of the differences, from control-run samples (add_control_climo)
        sig_test: bool = False,
        sig_method: str = "control",    # "control" | "control_resample"
        sig_alpha: float = 0.05,        # test level (or the FDR level when sig_fdr=True)
        sig_fdr: bool = False,          # Benjamini-Hochberg control of the false discovery rate
        sig_hatch: str = "..",          # unused (stippling is drawn by stipple()); kept for old calls
        sig_hatch_color: str = "k",
        sig_min_diff: Optional[float] = None,   # stipple only where also |Tgt−Ref| >= this
        sig_grid=(72, 36),              # stipple dots per panel (lon x lat), see stipple()
        sig_dot_size: float = 4.0,      # dot size (points^2, at this figure's size)
    ):
        """
        sig_test=True stipples the difference panels where Tgt−Ref is significant. The variance
        of a climatology-window mean comes from the control climatologies of add_control_climo:
          "control"          t-test, t = d / (s_c sqrt 2), df = n_chunks − 1 (control_ttest_p)
          "control_resample" |d| ranked against all pairwise control differences (control_resample_p)
        Only field="annual_mean".
        """
        FS = self._font_sizes(fontz)
        
        tgt_list = list(tgt_idx) if isinstance(tgt_idx, (list, tuple)) else [int(tgt_idx)]
        ncols = 1 + len(tgt_list)

        pc = ccrs.PlateCarree()
        
        fig, axes, caxes = self._axes(ncols, figsize, wspace, cbar_width=cbar_width, fig=fig, slot=slot)
        
        # move them closer to their panels
        for k in range(ncols):
            caxes[k] = self._shrink_cbar_axes(fig, caxes[k], height_frac=cbar_height, align=cbar_align)
            self._dock_cbar(fig, axes[k], caxes[k], pad=cbar_pad, side="right")
            if cbar_fit_panel is not None:
                fit_cbar_to_panel(axes[k], caxes[k], cbar_fit_panel, pad=cbar_pad)

        letters = string.ascii_lowercase
        panel_prefix = [f"({letters[i+fig_idx]}) " for i in range(ncols)]

        # --- reference ---
        da_ref = self._resolve_field(ref_idx, field, month)
        Zref = da_ref.values
        lat  = da_ref["lat"].values
        lon  = da_ref["lon"].values

        if (vmin is None) or (vmax is None):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rmax = float(np.nanmax(Zref)); rmin = float(np.nanmin(Zref))
            vmin = rmin if vmin is None else vmin
            vmax = rmax if vmax is None else vmax

        ref_levels = np.asarray(clevs, float) if clevs is not None else np.linspace(vmin, vmax, int(n_levels))
        cmap_ref   = self._cmap_with_over_under(cmap, n=ref_cmap_n, darken_factor=ref_darken)
        norm_ref   = BoundaryNorm(ref_levels, cmap_ref.N, clip=False)

        lonr, Zrefp = self._with_cyclic(lon, Zref)
        
        ax0 = axes[0]
        im0 = ax0.contourf(lonr, lat, Zrefp, levels=ref_levels, cmap=cmap_ref,
                           norm=norm_ref, extend=extend, transform=ccrs.PlateCarree(), zorder=1)
        ax0.coastlines(linewidth=0.6, zorder=3)
        ax0.add_feature(cfeature.LAND, zorder=2, edgecolor="none", facecolor="0.9")
        self._add_geophysical_lines(ax0, gridlines=gridlines, equator=equator, tropics=tropics)
        
        # tidy ticks/labels
        ax0.set_xticks(np.arange(-180, 181, 60), crs=pc)
        ax0.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
        ax0.set_yticks(np.arange(-90, 91, 30), crs=pc)
        ax0.yaxis.set_major_formatter(LatitudeFormatter())
        ax0.tick_params(axis="both", which="major", labelsize=FS["tick"])

        # ---- reference colorbar (replace the old fig.colorbar call) ----
        cb0 = fig.colorbar(im0, cax=caxes[0], boundaries=ref_levels, spacing="uniform", extend=extend)
        cb0.set_label(cbar_label, fontsize=FS["cbar_label"])
        cb0.ax.tick_params(labelsize=FS["cbar_tick"])
        if cbar_tick_every:                  # label every n-th level (matplotlib's own ticks otherwise)
            cb0.set_ticks(ref_levels[::int(cbar_tick_every)])
        cb0.outline.set_linewidth(0.6)
        
        lab_ref = self._meta[ref_idx].get("label", f"Exp {ref_idx}")
        title0 = f"{panel_prefix[0]}{var_label} ({lab_ref})" if field != "monthly_climatology" \
                 else f"{panel_prefix[0]}{var_label} {month:02d} ({lab_ref})"
        ax0.set_title(title0, fontsize=FS["title"],loc="left")
        
        # --- diffs ---
        diffs = []
        for ti in tgt_list:
            da_tgt = self._resolve_field(ti, field, month).interp(lat=da_ref["lat"], lon=da_ref["lon"])
            diffs.append((ti, da_tgt - da_ref))

        if (diff_vmin is None) or (diff_vmax is None):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                dmax = float(np.nanmax(np.abs(np.concatenate([np.ravel(D.values) for _, D in diffs]))))
                if not np.isfinite(dmax) or dmax == 0: dmax = 1.0
            diff_vmin = -dmax if diff_vmin is None else diff_vmin
            diff_vmax =  dmax if diff_vmax is None else diff_vmax

        # significance masks (True = significant), on the ref grid
        sig_masks = {}
        if sig_test:
            if field != "annual_mean":
                raise ValueError("sig_test supports field='annual_mean' only.")
            if self._control is None:
                raise ValueError("sig_test needs add_control_climo() first.")
            if sig_method not in ("control", "control_resample"):
                raise ValueError(f"Unknown sig_method '{sig_method}'.")
            ctrl = self._control.interp(lat=da_ref["lat"], lon=da_ref["lon"]).values.astype(float)
            test = control_ttest_p if sig_method == "control" else control_resample_p
            for ti, D in diffs:
                p = test(D.values.astype(float), ctrl)
                sig_masks[ti] = (fdr_significant(p, sig_alpha) if sig_fdr
                                 else np.isfinite(p) & (p < sig_alpha))
                if sig_min_diff is not None:
                    sig_masks[ti] &= np.abs(D.values) >= sig_min_diff

        diff_levels = np.asarray(diff_clevs, float) if diff_clevs is not None \
                      else self._symmetric_levels(diff_vmin, diff_vmax, n_levels=diff_n_levels)
        cmap_d = self._cmap_with_over_under(diff_cmap, n=diff_cmap_n, darken_factor=0.8)
        norm_d = BoundaryNorm(diff_levels, cmap_d.N, clip=False)

        for j, (ti, D) in enumerate(diffs, start=1):
            ax = axes[j]
            Zd = D.values
            lond, Zdp  = self._with_cyclic(lon, Zd)

            im = ax.contourf(lond, lat, Zdp, levels=diff_levels, cmap=cmap_d, norm=norm_d,
                             extend=diff_extend, transform=ccrs.PlateCarree(), zorder=1)
            if ti in sig_masks:
                stipple(ax, lon, lat, sig_masks[ti], grid=sig_grid, size=sig_dot_size,
                        color=sig_hatch_color, transform=ccrs.PlateCarree())
            ax.coastlines(linewidth=0.6, zorder=3)
            ax.add_feature(cfeature.LAND, zorder=2, edgecolor="none", facecolor="0.9")
            self._add_geophysical_lines(ax, gridlines=gridlines, equator=equator, tropics=tropics)
            # tidy ticks/labels
            ax.set_xticks(np.arange(-180, 181, 60), crs=pc)
            ax.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))
            ax.set_yticks(np.arange(-80, 81, 20), crs=pc)
            ax.yaxis.set_major_formatter(LatitudeFormatter())
            ax.tick_params(axis="both", which="major", labelsize=FS["tick"])

            lab_tgt = self._meta[ti].get("label", f"Exp {ti}")
            titlej = f"{panel_prefix[j]}{var_label} ({lab_tgt} − {lab_ref})"
            
            ax.set_title(titlej, fontsize=FS["title"], loc="left")
            # ---- reference colorbar (replace the old fig.colorbar call) ----
            cb = fig.colorbar(im, cax=caxes[j], boundaries=diff_levels, spacing="uniform", extend=diff_extend)
            cb.set_label(diff_label, fontsize=FS["cbar_label"])
            cb.ax.tick_params(labelsize=FS["cbar_tick"])
            every = diff_tick_every if diff_tick_every is not None else cbar_tick_every
            if every:                        # label every n-th level
                cb.set_ticks(diff_levels[::int(every)])
            cb.outline.set_linewidth(0.6)

        if diff_cbar == "shared":            # keep only the last difference colorbar
            for k in range(1, ncols - 1):
                caxes[k].set_visible(False)
        self.last_axes = (axes, caxes)       # for arrange_panel_grid

        for k, ax in enumerate(axes):
            if k > 0:
                ax.set_yticks([])            # no lat labels on middle/right
                ax.yaxis.set_ticklabels([])
        for ax in axes:
            ax.set_xlabel("Longitude", fontsize=FS["axis_label"])
        for k in range(1, len(axes)):
            axes[k].set_ylabel("")  # hide duplicate y-labels

        if save:
            fig.savefig(save, dpi=300, bbox_inches="tight")
        return fig, axes

    def save_surface2d_cache(
        self,
        out_nc: str,
        *,
        fields=("annual_mean",),
        compress: bool = True,
        complevel: int = 4,
        float_dtype: str = "float32",
    ) -> str:
        """
        Save plotter._stats fields to a NetCDF cache.

        Output schema
        -------------
        coords:
          - exp   (exp)
          - label (exp)
          - lat   (lat)
          - lon   (lon)

        data_vars (if present):
          - annual_mean(exp, lat, lon)
          - seasonal_amplitude(exp, lat, lon)
          - month_of_max(exp, lat, lon)  [int16 if present]
          - monthly_climatology(exp, month, lat, lon)  [if present]

        Notes
        -----
        For ANN climo files (no time dim), annual_mean is already (lat, lon).
        """

        if (not self._meta) or (not self._stats):
            raise ValueError("No experiments added yet (self._meta/self._stats are empty).")

        exps = [m.get("exp", f"exp{i}") for i, m in enumerate(self._meta)]
        labels = [m.get("label", e) for m, e in zip(self._meta, exps)]

        # --- Find reference lat/lon (and month if needed)
        lat = lon = None
        for st in self._stats:
            for f in fields:
                da = st.get(f, None)
                if isinstance(da, xr.DataArray) and ("lat" in da.coords) and ("lon" in da.coords):
                    lat = da["lat"].values
                    lon = da["lon"].values
                    break
            if lat is not None:
                break
        if lat is None or lon is None:
            raise ValueError("Could not find lat/lon in stored stats for requested fields.")

        ds = xr.Dataset(coords=dict(
            exp=("exp", np.asarray(exps, dtype=object)),
            label=("exp", np.asarray(labels, dtype=object)),
            lat=("lat", lat),
            lon=("lon", lon),
        ))

        # --- Helper to set encoding
        def _enc(name, dtype=None):
            e = {}
            if compress:
                e.update(dict(zlib=True, complevel=int(complevel)))
            if dtype is not None:
                e["dtype"] = dtype
            e["_FillValue"] = np.nan if (dtype is None or "float" in str(dtype)) else None
            return e

        encoding = {}

        for f in fields:
            # monthly_climatology has (month, lat, lon)
            if f == "monthly_climatology":
                # Find a month coordinate
                month_vals = None
                for st in self._stats:
                    da = st.get("monthly_climatology", None)
                    if isinstance(da, xr.DataArray) and ("month" in da.coords):
                        month_vals = da["month"].values
                        break
                if month_vals is None:
                    # nothing to write
                    continue

                arr = np.full((len(exps), len(month_vals), len(lat), len(lon)), np.nan, dtype=float_dtype)
                ok_any = False
                for i, st in enumerate(self._stats):
                    da = st.get("monthly_climatology", None)
                    if not isinstance(da, xr.DataArray):
                        continue
                    da2 = da.transpose("month", "lat", "lon")
                    arr[i, :, :, :] = da2.values.astype(float_dtype)
                    ok_any = True
                if ok_any:
                    ds = ds.assign_coords(month=("month", month_vals))
                    ds[f] = (("exp", "month", "lat", "lon"), arr)
                    encoding[f] = _enc(f, dtype=float_dtype)
                continue

            # 2D fields: (lat, lon)
            arr = np.full((len(exps), len(lat), len(lon)), np.nan, dtype=float_dtype)
            ok_any = False
            for i, st in enumerate(self._stats):
                da = st.get(f, None)
                if not isinstance(da, xr.DataArray):
                    continue
                da2 = da.transpose("lat", "lon")
                arr[i, :, :] = da2.values.astype(float_dtype)
                ok_any = True

            if ok_any:
                # Use int16 for month_of_max if present
                if f == "month_of_max":
                    ds[f] = (("exp", "lat", "lon"), arr.astype("int16"))
                    encoding[f] = _enc(f, dtype="int16")
                else:
                    ds[f] = (("exp", "lat", "lon"), arr)
                    encoding[f] = _enc(f, dtype=float_dtype)

        # --- Provenance
        ds.attrs["surface2d_meta_json"] = json.dumps(dict(
            root_dir=self.root_dir,
            file_key=self.file_key,
            scale=self.scale,
            target_res=self.target_res,
            area_weighted=self.area_weighted,
            calendar=self.calendar,
            time_origin=self.time_origin,
            per_exp_meta=self._meta,
        ), default=str)

        ds.to_netcdf(out_nc, encoding=encoding)
        return out_nc


class AMOCClimoPlotter:
    """
    AMOC-only (time, depth, lat) climatology/statistics + 2-D plotting.

    - Opens lazily with dask (xarray.open_mfdataset).
    - Normalizes the time dim to 'time' and supplies a datetime coord if possible.
    - If datetime is unavailable, uses a numeric month label (1..12) for climatology.
    """

    def __init__(
        self,
        root_dir: str,
        sub_dir: str,
        file_key: str = "mocTimeSeries",
        var_mpas: Optional[str] = "mocAtlantic",   # internal variable name in files
        scale: float = 1.0,
        weighted: bool = False,                    # days-in-month weights (requires datetime)
        engine: Optional[str] = "netcdf4",
        chunks: Optional[dict] = None,             # e.g. {"Time": 240} or {"time": 240}
        use_mfdataset: bool = True,
        time_origin: Optional[str] = "0001-01-01 00:00:00",
        calendar: str = "noleap",
        compute_interannual_std: bool = False,
        dtype: str = "float32",
        combine: str = "nested",                   # strict concat along time avoids merge errors
        open_parallel: bool = True,
    ):
        self.root_dir = root_dir
        self.sub_dir = sub_dir
        self.file_key = file_key
        self.var_mpas = var_mpas or "mocAtlantic"
        self.scale = float(scale)
        self.weighted = bool(weighted)
        self.engine = engine
        self.use_mfdataset = bool(use_mfdataset)
        self.time_origin = time_origin
        self.calendar = calendar
        self.compute_interannual_std = bool(compute_interannual_std)
        self.dtype = dtype
        self.combine = combine
        self.open_parallel = bool(open_parallel)
        self.chunks = chunks if chunks is not None else {"Time": 240}

        self._meta: List[dict] = []
        self._stats: List[dict] = []
        self._control: Optional[xr.DataArray] = None   # control-run samples (add_control_climo)
        
    def _lat_formatter(self):
        """Format latitude ticks as 30°S, 0°, 45°N."""
        def _fmt(x, pos):
            if np.isnan(x):
                return ""
            x = float(x)
            if x == 0:
                return "0°"
            hemi = "N" if x > 0 else "S"
            return f"{abs(x):g}°{hemi}"
        return FuncFormatter(_fmt)
    
    def _cmap_with_over_under(self, base_cmap, n: int = 256, darken_factor: float = 0.8):
        def _darken(color, f=0.8):
            r, g, b, a = mpl.colors.to_rgba(color)
            return (r*f, g*f, b*f, a)

        cmap = mpl.colormaps.get_cmap(base_cmap) if isinstance(base_cmap, str) else base_cmap
        colors = cmap(np.linspace(0.05,0.95, n+21))
        new = mpl.colors.LinearSegmentedColormap.from_list(f"{cmap.name}_mod", colors[10:n+10], N=n)
        new.set_under(_darken(colors[0], darken_factor))
        new.set_over(_darken(colors[-1], darken_factor))
        new.set_bad("0.9")
        return new
    
    def _symmetric_levels(self, vmin, vmax, *, n_levels=None, step=None):
        """
        Build symmetric boundaries about 0 that include 0 exactly.
        Provide either n_levels (bins) or step (bin width).
        """
        import numpy as np
        a = float(max(abs(vmin), abs(vmax)))
        if step is not None:
            a = step * np.ceil(a / step)
            levels = np.arange(-a, a + 0.5*step, step, dtype=float)
        else:
            # n_levels boundaries -> n_levels-1 bins; keep odd so 0 is a boundary
            n = int(n_levels) if n_levels else 13
            if n % 2 == 0:
                n += 1
            levels = np.linspace(-a, a, n, dtype=float)
        # ensure exact 0
        zidx = np.argmin(abs(levels))
        levels[zidx] = 0.0
        return levels

    # --------------------------- time helpers ---------------------------
    def _construct_monthly_time(self, n: int):
        """Build a cftime monthly index starting at time_origin."""
        if self.time_origin is None:
            return None
    
        s = str(self.time_origin)
        try:
            y = int(s[0:4]); m = int(s[5:7]); d = int(s[8:10])
            hh = int(s[11:13]) if len(s) >= 13 else 0
            mm = int(s[14:16]) if len(s) >= 16 else 0
            ss = int(s[17:19]) if len(s) >= 19 else 0
        except Exception:
            return None
    
        n = int(n)
        if n <= 0:
            return None
    
        t0 = cftime.DatetimeNoLeap(y, m, d, hh, mm, ss)
        return xr.date_range(
            start=t0,
            periods=n,
            freq="MS",
            use_cftime=True,
            calendar="noleap",
        )

    # ----------------------- file discovery & open -----------------------
    def _find_moc_timeseries_files(self, exp: str, year_token: str) -> List[str]:
        a, b = [int(t.strip()) for t in year_token.split("-")]
        if a > b:
            a, b = b, a
        base = os.path.join(self.root_dir, exp, self.sub_dir if self.sub_dir else "")
        hits = []
        for yy in range(a, b + 1):
            pat = os.path.join(base, f"{self.file_key}_{yy:04d}-{yy:04d}.nc")
            hits.extend(glob.glob(pat))
        if not hits:
            hits = glob.glob(os.path.join(base, f"{self.file_key}.nc"))
        hits.sort()
        return hits

    def _preprocess(self, desired_name: str, ds: xr.Dataset) -> xr.Dataset:
        # --- find an AMOC-like field, even if it's in coords ---
        cand_names = [
            desired_name, "mocAtlantic", "mocAtlantic26", "mocStreamfunction",
            "moc_atlantic", "MOCAtlantic", "MOC", "moc"
        ]
        actual = None
        for cand in cand_names:
            if cand in ds.data_vars:
                actual = cand
                break
        if actual is None:
            # some files park it as a coordinate-like variable
            for cand in cand_names:
                if cand in ds.variables:
                    actual = cand
                    break
        if actual is None:
            raise KeyError(
                f"AMOC variable not found. Available data_vars: {list(ds.data_vars)}; "
                f"variables: {list(ds.variables)}"
            )
        if actual != desired_name:
            ds = ds.rename({actual: desired_name})

        # --- time decoding branch stays the same from your version down to "keep only AMOC" ---
        # try to build a true 'time' coord from CF units on "Time"
        if "Time" in ds.dims and "time" not in ds.coords:
            if "Time" in ds and "units" in ds["Time"].attrs:
                try:
                    t = xr.coding.times.decode_cf_datetime(
                        ds["Time"].values,
                        units=ds["Time"].attrs.get("units"),
                        calendar=ds["Time"].attrs.get("calendar", "noleap"),
                        use_cftime=True,
                    )
                    ds = ds.assign_coords(time=("Time", t)).swap_dims({"Time": "time"})
                except Exception:
                    pass

        # if still no 'time', try MPAS character times
        if "Time" in ds.dims and "time" not in ds.coords:
            for nm in ("xtime_startMonthly", "xtime_endMonthly", "startTime", "endTime"):
                if nm in ds:
                    a = ds[nm].values  # (Time, StrLen)
                    s = ["".join(r.astype(str)).strip().strip("\x00").replace("\x00", "") for r in a]
                    tt = []
                    ok = True
                    for txt in s:
                        txt = txt.replace("T", " ").replace("_", " ")
                        try:
                            y, m, d = int(txt[0:4]), int(txt[5:7]), int(txt[8:10])
                            tt.append(cftime.DatetimeNoLeap(y, m, d))
                        except Exception:
                            ok = False
                            break
                    if ok and len(tt) == ds.dims["Time"]:
                        ds = ds.assign_coords(time=("Time", np.asarray(tt, dtype=object))).swap_dims({"Time": "time"})
                        break

        if "time" not in ds.dims and "Time" in ds.dims:
            ds = ds.rename({"Time": "time"})
        if "time" in ds.dims and "time" not in ds.coords:
            ds = ds.assign_coords(time=("time", np.arange(ds.sizes["time"], dtype=int)))

        drop_vars = [v for v in ("xtime_startMonthly", "xtime_endMonthly", "startTime", "endTime") if v in ds.variables]
        if drop_vars:
            ds = ds.drop_vars(drop_vars)

        # keep only AMOC field + coords
        keep = {desired_name} | set(ds.coords)
        ds = ds[list(keep & (set(ds.data_vars) | set(ds.coords)))]
        return ds

    def _infer_start_year_from_files(self, files: List[str]) -> Optional[int]:
        # match "..._YYYY-YYYY.nc"
        for fp in files:
            m = re.search(r"_(\d{4})-\1\.nc$", os.path.basename(fp))
            if m:
                return int(m.group(1))
            # match "..._YYYYMM-YYYYMM.nc"
            m = re.search(r"_(\d{6})-(\d{6})\.nc$", os.path.basename(fp))
            if m:
                yyyymm = int(m.group(1))
                return yyyymm // 100
        return None
    
    def _time_has_dt(self, da: xr.DataArray) -> bool:
        """Return True if da has a datetime-like 'time' coordinate with .dt access."""
        if "time" not in da.coords:
            return False
        try:
            _ = da["time"].dt.month
            return True
        except Exception:
            return False

    def _dedupe_time(self, da):
        # keep first occurrence for any duplicate time stamps
        try:
            idx = da.get_index("time")
            if not idx.is_unique:
                keep = ~idx.duplicated()
                da = da.isel(time=keep)
        except Exception:
            pass
        return da

    def _open_amoc_lazy(self, files: List[str], year_token: str) -> xr.DataArray:
        """
        Open lazily and concatenate strictly along time.
        Always guarantees a monthly `cftime` datetime coordinate named 'time'.
        """
        xr_opts = dict(engine=self.engine, chunks=self.chunks,
                       decode_times=True, use_cftime=True)

        def _pp(ds):
            return self._preprocess(self.var_mpas, ds)

        # --- open files ---
        if (len(files) == 1) or (not self.use_mfdataset):
            ds = xr.open_dataset(files[0], **xr_opts)
            ds = _pp(ds)
            da = ds[self.var_mpas]
        else:
            try:
                ds = xr.open_mfdataset(
                    files,
                    preprocess=_pp,
                    combine=self.combine,           # "nested" recommended
                    concat_dim=("time" if self.combine == "nested" else None),
                    data_vars="minimal",
                    coords="minimal",
                    compat="override",
                    join="override",
                    combine_attrs="drop",
                    parallel=self.open_parallel,
                    **xr_opts,
                )
                da = ds[self.var_mpas]
            except Exception:
                # fallback manual concat (ensure each file is closed)
                datasets = []
                try:
                    for fp in files:
                        dso = xr.open_dataset(fp, **xr_opts)
                        dso = _pp(dso)
                        datasets.append(dso)

                    parts = []
                    for dso in datasets:
                        da_part = dso[self.var_mpas]
                        if "time" not in da_part.dims:
                            for d in da_part.dims:
                                if d.lower().startswith("time"):
                                    da_part = da_part.rename({d: "time"})
                                    break
                        parts.append(da_part)

                    # Concatenate while handles are still valid
                    da = xr.concat(parts, dim="time", join="override")

                finally:
                    # Proactively close all open datasets
                    for dso in datasets:
                        try:
                            dso.close()
                        except Exception:
                            pass
        # --- normalize dim name ---
        if "time" not in da.dims:
            for d in da.dims:
                if d.lower().startswith("time"):
                    da = da.rename({d: "time"})
                    break

        # --- reconstruct guaranteed monthly datetime coordinate ---
        y0 = self._infer_start_year_from_files(files) or int(year_token.split("-")[0])
        ntime = int(da.sizes["time"])
        
        tt = xr.date_range(
            start=cftime.DatetimeNoLeap(y0, 1, 1),
            periods=ntime,
            freq="MS",
            use_cftime=True,
            calendar="noleap",
        )
        
        da = da.assign_coords(time=("time", tt))
        
        # --- slice to requested years ---
        a, b = [int(t) for t in year_token.split("-")]
        mask = (tt.year >= a) & (tt.year <= b)
        if bool(mask.any()):
            da = da.isel(time=mask)

        return da

    # ----------------------- public: add & compute -----------------------
    def add_experiment(
        self,
        exp: str,
        label: Optional[str] = None,
        year_token: str = "0001-0350",
        sub_dir: Optional[str] = None,
        clim_start: Optional[str] = None,   # e.g., "0321-01-01"
        clim_end: Optional[str] = None,     # e.g., "0350-12-31"
    ):
        """Open AMOC lazily, compute monthly climatology + stats."""
        if sub_dir is not None:
            old = self.sub_dir
            self.sub_dir = sub_dir

        files = self._find_moc_timeseries_files(exp, year_token)
        if not files:
            raise FileNotFoundError(f"No AMOC files under {os.path.join(self.root_dir, exp, self.sub_dir)}")

        da = self._open_amoc_lazy(files, year_token) * self.scale  # (time, depth, lat) dask-backed

        # Ensure datetime coord if possible; else we’ll fall back to numeric months
        has_dt = self._time_has_dt(da)

        # Climatology window
        if clim_start or clim_end:
            if has_dt:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    da = da.sel(time=slice(clim_start, clim_end))
            else:
                # numeric fallback based on year_token
                a, b = [int(t) for t in year_token.split("-")]
                ys = int(str(clim_start)[:4]) if clim_start else a
                ye = int(str(clim_end)[:4])   if clim_end   else b
                m0 = max(0, (ys - a) * 12)
                m1 = min(da.sizes["time"], (ye - a + 1) * 12)
                da = da.isel(time=slice(m0, m1))

        # ----- Monthly climatology -----
        if has_dt:
            if self.weighted:
                # days-in-month weights
                dpm = da["time"].dt.days_in_month.astype("float")
                w = dpm / dpm.groupby("time.month").sum()
                w_al = xr.DataArray(w, dims=["time"]).broadcast_like(da)
                clim12 = (da * w_al).groupby("time.month").sum("time", skipna=True)
                # annual mean via mean month-length weights
                dpm_m = dpm.groupby("time.month").mean().astype("float")
                w_m = dpm_m / dpm_m.sum()
                w_m = xr.DataArray(w_m, dims=["month"]).broadcast_like(clim12)
                annual_mean = (clim12 * w_m).sum("month", skipna=True)
            else:
                clim12 = da.groupby("time.month").mean("time", skipna=True)
                annual_mean = clim12.mean("month", skipna=True)
        else:
            # No datetime: group by repeating month labels (1..12)
            n = da.sizes["time"]
            months = xr.DataArray((np.arange(n) % 12) + 1, dims=["time"], name="month")
            clim12 = da.groupby(months).mean("time", skipna=True)
            clim12 = clim12.rename({"group": "month"}) if "group" in clim12.dims else clim12
            annual_mean = clim12.mean("month", skipna=True)
            if self.weighted:
                warnings.warn("weighted=True requested but time is not datetime-like; using unweighted means.")

        seasonal_amplitude = clim12.max("month") - clim12.min("month")
        month_of_max = (clim12.argmax("month") + 1).astype("int16")

        if self.compute_interannual_std:
            if has_dt:
                interannual_std = da.groupby("time.month").std("time", skipna=True)
            else:
                n = da.sizes["time"]
                months = xr.DataArray((np.arange(n) % 12) + 1, dims=["time"], name="month")
                interannual_std = da.groupby(months).std("time", skipna=True).rename({"group": "month"})
        else:
            interannual_std = None

        # ----- Annual-mean series over the climatology window (samples for significance tests) -----
        if has_dt:
            annual_series = da.groupby("time.year").mean("time", skipna=True)
        else:
            n = da.sizes["time"]
            years = xr.DataArray(np.arange(n) // 12, dims=["time"], name="year")
            annual_series = da.groupby(years).mean("time", skipna=True)

        stats = {
            "monthly_climatology": clim12,                # (month, depth, lat)
            "annual_mean": annual_mean,                   # (depth, lat)
            "seasonal_amplitude": seasonal_amplitude,     # (depth, lat)
            "month_of_max": month_of_max,                 # (depth, lat) int16
            "annual_series": annual_series,               # (year, depth, lat)
        }
        if interannual_std is not None:
            stats["monthly_interannual_std"] = interannual_std

        self._meta.append({
            "exp": exp,
            "label": label or exp,
            "year_token": year_token,
            "clim_start": clim_start,
            "clim_end": clim_end,
        })
        self._stats.append(stats)

        if sub_dir is not None:
            self.sub_dir = old  # restore

    def add_experiment_climo(
        self,
        exp: str,
        label: Optional[str] = None,
        sub_dir: Optional[str] = None,
        years: str = "0301-0350",
        region: str = "Global",
    ):
        """
        Add an experiment from the MPAS-Analysis MOC climatology
        <root_dir>/<exp>/<sub_dir>/mocStreamfunction_years<years>.nc, which holds the annual-mean
        streamfunction of every basin (region = "Global", "Atlantic", "AtlanticMed", "IndoPacific";
        the time series hold only the Atlantic). Only field="annual_mean" is available, on
        dims (depth [m], lat).
        """
        path = os.path.join(self.root_dir, exp, sub_dir if sub_dir is not None else self.sub_dir,
                            f"mocStreamfunction_years{years}.nc")
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        with xr.open_dataset(path, engine=self.engine) as ds:
            da = xr.DataArray(
                ds[f"moc{region}"].values.astype(self.dtype) * self.scale,
                dims=("depth", "lat"),
                coords={"depth": ds["depth"].values, "lat": ds[f"lat{region}"].values},
                name=f"moc{region}",
            )
        y0, y1 = years.split("-")
        self._meta.append({
            "exp": exp,
            "label": label or exp,
            "year_token": years,
            "clim_start": f"{y0}-01-01",
            "clim_end": f"{y1}-12-31",
        })
        self._stats.append({"annual_mean": da})

    def add_control_climo(
        self,
        exp: str,
        windows: List[str],
        sub_dir: str = "post/analysis/mpas_analysis/*/clim/mpas/avg",
        region: str = "Global",
        level_aligned: bool = False,
    ):
        """
        Load non-overlapping climatologies of a control run (e.g. the ten 50-yr windows of
        piControl) as samples of internal variability for sig_method="control"/"control_resample".
        windows: e.g. ["0001-0050", "0051-0100", ...]; sub_dir may contain a glob.
        level_aligned: the experiments come from the MOC time series (add_experiment), whose depth
        levels are the layer interfaces (0, 10, ... m) while the climatology files label the same
        levels by layer midpoints (5, 15, ... m); use the experiments' depths for the climatology
        levels by index so that the samples are not shifted by half a layer.
        """
        chunks = []
        for years in windows:
            hits = sorted(glob.glob(os.path.join(self.root_dir, exp, sub_dir,
                                                 f"mocStreamfunction_years{years}.nc")))
            if not hits:
                raise FileNotFoundError(f"mocStreamfunction_years{years}.nc under {exp}/{sub_dir}")
            with xr.open_dataset(hits[0], engine=self.engine) as ds:
                chunks.append(xr.DataArray(
                    ds[f"moc{region}"].values.astype(self.dtype) * self.scale,
                    dims=("depth", "lat"),
                    coords={"depth": ds["depth"].values, "lat": ds[f"lat{region}"].values},
                ))
        ctrl = xr.concat(chunks, dim="chunk")   # (chunk, depth, lat)
        if level_aligned:
            ref = self._stats[0]["annual_mean"]
            dname = next(d for d in ref.dims if "depth" in d.lower() or d.lower().startswith("z"))
            zref = ref[dname].values
            if zref.size < ctrl.sizes["depth"]:
                raise ValueError("experiments have fewer depth levels than the control climatology")
            ctrl = ctrl.assign_coords(depth=zref[:ctrl.sizes["depth"]])
        self._control = ctrl

    # ------------------------- significance test -------------------------
    _control_ttest = staticmethod(control_ttest_p)
    _control_resample_p = staticmethod(control_resample_p)
    _fdr_significant = staticmethod(fdr_significant)

    @staticmethod
    def _lag1_autocorr(x: np.ndarray) -> np.ndarray:
        """Lag-1 autocorrelation along axis 0 (NaN where the series is constant)."""
        a = x - np.nanmean(x, axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.nansum(a[1:] * a[:-1], axis=0) / np.nansum(a * a, axis=0)

    def _welch_ttest_neff(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        Two-sided Welch t-test of mean(b) - mean(a) at every point, samples along axis 0.
        Each sample size is reduced for serial correlation, n_eff = n (1 - r1) / (1 + r1)
        with r1 the lag-1 autocorrelation clipped to [0, 0.95] (Bretherton et al. 1999).
        Returns p-values (NaN where both samples have zero variance).
        """
        def _parts(x):
            n = np.sum(np.isfinite(x), axis=0).astype(float)
            r1 = np.clip(np.nan_to_num(self._lag1_autocorr(x)), 0.0, 0.95)
            neff = np.maximum(n * (1.0 - r1) / (1.0 + r1), 2.0)
            return np.nanmean(x, axis=0), np.nanvar(x, axis=0, ddof=1), neff

        ma, va, na = _parts(a)
        mb, vb, nb = _parts(b)
        with np.errstate(invalid="ignore", divide="ignore"):
            sa, sb = va / na, vb / nb
            t = (mb - ma) / np.sqrt(sa + sb)
            df = (sa + sb) ** 2 / (sa ** 2 / (na - 1.0) + sb ** 2 / (nb - 1.0))
        p = 2.0 * stats.t.sf(np.abs(t), df)
        p[~np.isfinite(t)] = np.nan
        return p

    # ----------------------------- plotting -----------------------------
    def _resolve_field(self, idx, field, month, depth_coord, lat_coord):
        stat = self._stats[idx]
        if field not in stat:
            raise KeyError(f"Field '{field}' not found. Options: {list(stat.keys())}")
        da = stat[field]
        if field == "monthly_climatology":
            if month is None:
                raise ValueError("Provide month=1..12 for 'monthly_climatology'.")
            da = da.sel(month=int(month))
        # resolve coordinate names
        dc = depth_coord if depth_coord in da.dims else next((d for d in da.dims if "depth" in d.lower() or d.lower().startswith("z")), da.dims[0])
        lc = lat_coord   if lat_coord   in da.dims else next((d for d in da.dims if "lat" in d.lower()), da.dims[-1])
        if list(da.dims) != [dc, lc]:
            da = da.transpose(dc, lc)
        return da, dc, lc
    
    def _shrink_cbar_axes(self, fig, cax, height_frac=0.9, width_frac=1.0, 
                          valign="center", halign="center"):
        """
        Return a new colorbar axes with adjusted width and height.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
        cax : Axes
            Original colorbar axes.
        height_frac : float, optional
            Fraction of original height (0–1). Default 0.9.
        width_frac : float, optional
            Fraction of original width (0–1). Default 1.0.
        valign : {"center","top","bottom"}
            Vertical alignment after shrinking.
        halign : {"center","left","right"}
            Horizontal alignment after shrinking.
        """
        # clip fractions to reasonable range
        height_frac = max(0.05, min(1.0, float(height_frac)))
        width_frac  = max(0.05, min(1.0, float(width_frac)))

        bb = cax.get_position()

        # new dimensions
        new_h = bb.height * height_frac
        new_w = bb.width  * width_frac

        # vertical alignment
        if valign == "top":
            y0 = bb.y0 + (bb.height - new_h)
        elif valign == "bottom":
            y0 = bb.y0
        else:  # center
            y0 = bb.y0 + 0.5*(bb.height - new_h)

        # horizontal alignment
        if halign == "right":
            x0 = bb.x0 + (bb.width - new_w)
        elif halign == "left":
            x0 = bb.x0
        else:  # center
            x0 = bb.x0 + 0.5*(bb.width - new_w)

        # make new axes and delete old
        new_cax = fig.add_axes([x0, y0, new_w, new_h])
        cax.remove()
        return new_cax
    
    def _dock_cbar(self, fig, ax_panel, cax, *, pad=0.006, side="right"):
        """
        Move `cax` so it hugs `ax_panel` with a small gap.
        `pad` is in figure fraction (0–1). side = "right"|"left".
        """
        pb = ax_panel.get_position()      # panel bbox in figure coords
        cb = cax.get_position()           # current cbar bbox
        new_x0 = (pb.x1 + pad) if side == "right" else (pb.x0 - pad - cb.width)
        cax.set_position([new_x0, cb.y0, cb.width, cb.height])
        return cax
    
    # --- comparison plot: Reference | (T1−Ref) | (T2−Ref) with depth in km + discrete colorbars ---
    def plot_amoc_2d_compare(
        self,
        varnam: str,
        ref_idx: int,
        tgt_idx,                          # int or list/tuple of ints
        fig_idx: int = 1,                 # start panel lettering: 1->(a), 4->(d)
        field: str = "annual_mean",
        month: Optional[int] = None,
        depth_coord: str = "depth",
        lat_coord: str = "lat",
        cmap: str = "RdBu_r",
        diff_cmap: str = "bwr",
        vmin: float = None,
        vmax: float = None,
        diff_vmin: float = None,
        diff_vmax: float = None,
        clevs: Optional[List[float]] = None,
        n_levels: Optional[int] = None,
        diff_clevs: Optional[List[float]] = None,
        diff_n_levels: Optional[int] = None,
        extend: str = "both",          # <- keep current behavior by default
        diff_extend: str = "both",     # <- keep current behavior by default
        tick_every: Optional[int] = None,
        diff_tick_every: Optional[int] = None,
        depth_units: str = "m",
        figsize: tuple = (13.2, 7.4),              # taller & a bit narrower by default
        cbar_label: str = "Sv",
        cbar_width: float = 0.55,     # thinner/thicker bars
        cbar_height: float = 0.95,    # 0–1 fraction of panel height
        cbar_align: str = "center",
        cbar_pad: float = 0.006,      # distance between panel and colorbar (figure fraction)
        diff_label: str = "ΔSv",
        yreverse: bool = True,
        save: Optional[str] = None,
        fontz: float = 12.0,
        wspace: float = 0.12,
        lat_band: Optional[tuple] = (-34.0, 90.0),
        cbar_labelpad: Optional[float] = None,   # gap between colorbar ticks and label (points)
        sig_test: bool = False,       # stipple diff panels where Tgt−Ref is significant
        sig_method: str = "welch",    # "welch" | "control" | "control_resample"
        sig_alpha: float = 0.05,      # test level (or the FDR level when sig_fdr=True)
        sig_fdr: bool = False,        # Benjamini-Hochberg control of the false discovery rate
        sig_hatch: str = "..",          # unused (stippling is drawn by stipple()); kept for old calls
        sig_hatch_color: str = "k",
        sig_min_diff: Optional[float] = None,   # stipple only where also |Tgt−Ref| >= this
        sig_grid=(72, 36),            # stipple dots per panel (lat x depth), see stipple()
        sig_dot_size: float = 4.0,    # dot size (points^2, at this figure's size)
        fig=None,                     # draw into a row of an existing figure (fig + slot = a
        slot=None,                    # SubplotSpec of it); figsize is then ignored
        diff_cbar: str = "each",      # "shared": one colorbar for all difference panels (on the last one)
        cbar_fit_panel: Optional[float] = None,   # colorbar height / drawn panel height (None: cbar_height)
        box_aspect: Optional[float] = None,       # panel height / width (e.g. 0.5 = shape of global maps)
    ):
        """
        sig_test=True stipples the difference panels where p < sig_alpha. Only field="annual_mean".
        sig_method:
          "welch"            Welch t-test on the annual means of the climatology window (serial
                             correlation reduces the sample sizes, see _welch_ttest_neff). Needs
                             experiments added with add_experiment (time series).
          "control"          t-test with the variance of a window mean taken from the control
                             climatologies of add_control_climo (see _control_ttest). For
                             experiments added with add_experiment_climo (no yearly samples).
          "control_resample" |difference| ranked against all pairwise differences of the
                             control climatologies (see _control_resample_p).
        """
        # ---- base fonts (points) scaled by fontz ----
        BASE = dict(title=1.0, label=0.95, tick=0.95, cbar=0.95)
        FS   = {k: v * float(fontz) for k, v in BASE.items()}
        lpad = {} if cbar_labelpad is None else {"labelpad": cbar_labelpad}
    
        # normalize targets
        tgt_list = list(tgt_idx) if isinstance(tgt_idx, (list, tuple)) else [int(tgt_idx)]
    
        # reference
        da_ref, dc, lc = self._resolve_field(ref_idx, field, month, depth_coord, lat_coord)
    
        # --- optional latitude band (no-op unless lat_band is not None) ---
        if lat_band is not None:
            lo, hi = float(lat_band[0]), float(lat_band[1])
            da_ref = da_ref.sel({lc: slice(lo, hi)})
    
        X = da_ref[lc].values
        Y = da_ref[dc].values
        Zref = da_ref.values
    
        # depth units
        ylabel = "Depth (m)"
        if str(depth_units).lower() == "km":
            Y = np.asarray(Y, float) / 1000.0
            ylabel = "Depth (km)"
    
        # defaults for ranges
        if vmin is None or vmax is None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                amax = np.nanmax(np.abs(Zref))
            if np.isfinite(amax) and amax > 0:
                vmin = -amax if vmin is None else vmin
                vmax =  amax if vmax is None else vmax
    
        # diffs on ref grid (already banded if lat_band applied)
        diffs = []
        for ti in tgt_list:
            da_tgt, _, _ = self._resolve_field(ti, field, month, depth_coord, lat_coord)
            da_tgt_interp = da_tgt.interp({lc: da_ref[lc], dc: da_ref[dc]})
            diffs.append((ti, da_tgt_interp - da_ref))

        # significance masks (True = significant), on the ref grid
        sig_masks = {}
        if sig_test:
            if field != "annual_mean":
                raise ValueError("sig_test supports field='annual_mean' only.")
            def _on_ref(s):
                """(sample, depth, lat) samples interpolated to the reference grid."""
                sd = next(d for d in s.dims if "depth" in d.lower() or d.lower().startswith("z"))
                sl = next(d for d in s.dims if "lat" in d.lower())
                s = s.rename({k: v for k, v in ((sd, dc), (sl, lc)) if k != v})
                s = s.transpose(s.dims[0], dc, lc)
                return s.interp({lc: da_ref[lc], dc: da_ref[dc]}).values.astype(float)

            if sig_method == "welch":
                def _series(idx):
                    if "annual_series" not in self._stats[idx]:
                        raise ValueError(f"No annual series for '{self._meta[idx]['label']}': "
                                         "sig_method='welch' needs add_experiment().")
                    return _on_ref(self._stats[idx]["annual_series"])
                s_ref = _series(ref_idx)
                pvals = {ti: self._welch_ttest_neff(s_ref, _series(ti)) for ti in tgt_list}
            elif sig_method in ("control", "control_resample"):
                if self._control is None:
                    raise ValueError(f"sig_method='{sig_method}' needs add_control_climo() first.")
                ctrl = _on_ref(self._control)
                test = self._control_ttest if sig_method == "control" else self._control_resample_p
                pvals = {ti: test(D.values.astype(float), ctrl) for ti, D in diffs}
            else:
                raise ValueError(f"Unknown sig_method '{sig_method}'.")
            for ti, p in pvals.items():
                sig_masks[ti] = (self._fdr_significant(p, sig_alpha) if sig_fdr
                                 else np.isfinite(p) & (p < sig_alpha))
            if sig_min_diff is not None:
                for ti, D in diffs:
                    sig_masks[ti] &= np.abs(D.values) >= sig_min_diff
    
        # diff range
        if diff_vmin is None or diff_vmax is None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                all_diff = np.concatenate([np.ravel(D.values) for _, D in diffs]) if diffs else np.array([0.0])
                dmax = np.nanmax(np.abs(all_diff))
            if np.isfinite(dmax) and dmax > 0:
                diff_vmin = -dmax if diff_vmin is None else diff_vmin
                diff_vmax =  dmax if diff_vmax is None else diff_vmax
    
        # --- LAYOUT (outer panels + per-panel colorbar axes) ---
        n_panels = 1 + len(tgt_list)  # 2 or 3
        if fig is None:
            fig = plt.figure(figsize=figsize)
        grid = slot.subgridspec if slot is not None else fig.add_gridspec
        outer = grid(
            1, n_panels,
            wspace=wspace,
            width_ratios=[1.0] * n_panels
        )
    
        axes, caxes = [], []
        for i in range(n_panels):
            sg = outer[0, i].subgridspec(1, 2, width_ratios=[1.0, 0.055], wspace=0.04)
            axes.append(fig.add_subplot(sg[0, 0]))
            caxes.append(fig.add_subplot(sg[0, 1]))
            if box_aspect is not None:
                axes[-1].set_box_aspect(box_aspect)
    
        # shrink + dock cbar axes
        for k in range(n_panels):
            caxes[k] = self._shrink_cbar_axes(fig, caxes[k], height_frac=cbar_height, width_frac=cbar_width)
            self._dock_cbar(fig, axes[k], caxes[k], pad=cbar_pad, side="right")
            if cbar_fit_panel is not None:
                fit_cbar_to_panel(axes[k], caxes[k], cbar_fit_panel, pad=cbar_pad)
    
        # panel lettering from fig_idx
        letters = string.ascii_lowercase
        start = max(1, int(fig_idx))  # 1-indexed
        panel_prefix = [f"({letters[i % 26]}) " for i in range(start - 1, start - 1 + n_panels)]
    
        # --- panel 1: reference ---
        ax0 = axes[0]
        cmap_ref = self._cmap_with_over_under(cmap, n=256, darken_factor=1.0)
    
        if clevs is not None:
            ref_levels = np.asarray(clevs, float)
        elif n_levels is not None and vmin is not None and vmax is not None:
            ref_levels = np.linspace(vmin, vmax, int(n_levels) + 1)
        else:
            ref_levels = np.linspace(vmin, vmax, 13)
    
        norm_ref = BoundaryNorm(ref_levels, cmap_ref.N, clip=False)
        im0 = ax0.contourf(X, Y, Zref, levels=ref_levels, cmap=cmap_ref, norm=norm_ref, extend=extend)
        cb0 = fig.colorbar(
            im0, cax=caxes[0],
            boundaries=ref_levels, norm=norm_ref,
            spacing="uniform", extend=extend
        )
        ax0.contour(X, Y, Zref, levels=ref_levels, colors="k", linewidths=0.5, alpha=0.6)
    
        if tick_every:                       # label every n-th level (matplotlib's own ticks otherwise)
            cb0.set_ticks(ref_levels[::int(tick_every)])
    
        if yreverse:
            ax0.invert_yaxis()
        ax0.tick_params(labelsize=FS["tick"])
        lab_ref = self._meta[ref_idx].get("label", f"Exp {ref_idx}")
        ax0.set_title(panel_prefix[0] + f"{varnam} ({lab_ref})", fontsize=FS["title"], loc="left")
        ax0.set_xlabel("Latitude", fontsize=FS["label"])
        ax0.xaxis.set_major_formatter(self._lat_formatter())
        ax0.set_ylabel(ylabel, fontsize=FS["label"])
        
        cb0.set_label(cbar_label, fontsize=FS["label"], **lpad)
        cb0.ax.tick_params(labelsize=FS["tick"])
        ax0.grid(True, alpha=0.2)
    
        # --- panels 2..: diffs vs ref ---
        for j, (ti, D) in enumerate(diffs, start=1):
            ax = axes[j]
            Zd = D.values
            cmap_d = self._cmap_with_over_under(diff_cmap, n=256, darken_factor=0.8)
    
            if diff_clevs is not None:
                diff_levels = np.asarray(diff_clevs, float)
            elif diff_n_levels is not None and (diff_vmin is not None and diff_vmax is not None):
                diff_levels = self._symmetric_levels(diff_vmin, diff_vmax, n_levels=int(diff_n_levels))
            else:
                diff_levels = self._symmetric_levels(diff_vmin, diff_vmax, n_levels=13)
    
            norm_d = TwoSlopeNorm(vmin=float(diff_levels.min()), vcenter=0.0, vmax=float(diff_levels.max()))
            im = ax.contourf(X, Y, Zd, levels=diff_levels, cmap=cmap_d, norm=norm_d, extend=diff_extend)
            cb = fig.colorbar(
                im, cax=caxes[j],
                boundaries=diff_levels,
                spacing="uniform", extend=diff_extend
            )
    
            ax.contour(X, Y, Zd, levels=diff_levels, colors="k", linewidths=0.5, alpha=0.6)

            if ti in sig_masks:
                stipple(ax, X, Y, sig_masks[ti], grid=sig_grid, size=sig_dot_size, color=sig_hatch_color)

            if diff_tick_every:              # label every n-th level
                cb.set_ticks(diff_levels[::int(diff_tick_every)])
    
            if yreverse:
                ax.invert_yaxis()

            ax.tick_params(labelsize=FS["tick"])
            lab_tgt = self._meta[ti].get("label", f"Exp {ti}")
            ax.set_title(panel_prefix[j] + f"{varnam} ({lab_tgt} − {lab_ref})", fontsize=FS["title"], loc="left")            
            ax.set_xlabel("Latitude", fontsize=FS["label"])
            ax.xaxis.set_major_formatter(self._lat_formatter())
            if j == 1:
                ax.set_ylabel(ylabel, fontsize=FS["label"])
                
            cb.set_label(diff_label, fontsize=FS["label"], **lpad)
            cb.ax.tick_params(labelsize=FS["tick"])
            ax.grid(True, alpha=0.2)
    
        # hide duplicate y-axis labels/ticks on panels 2+
        for k in range(1, len(axes)):
            axes[k].set_ylabel("")
            axes[k].set_yticks([])
            axes[k].yaxis.set_ticklabels([])

        if diff_cbar == "shared":            # keep only the last difference colorbar
            for k in range(1, n_panels - 1):
                caxes[k].set_visible(False)
        self.last_axes = (axes, caxes)       # for arrange_panel_grid
    
        if save:
            fig.savefig(save, dpi=300, bbox_inches="tight")
        return fig, axes
        
    def save_amoc_cache(
        self,
        out_nc: str,
        *,
        fields=("annual_mean",),          # e.g. ("annual_mean","seasonal_amplitude","month_of_max","monthly_climatology")
        compress: bool = True,
        complevel: int = 4,
        dtype: str = "float32",
    ) -> str:
        """
        Save selected fields from self._stats to a NetCDF cache.

        Expected shapes:
          - annual_mean, seasonal_amplitude: (depth, lat)
          - month_of_max: (depth, lat) integer-ish
          - monthly_climatology: (month, depth, lat)
          - monthly_interannual_std (optional): (month, depth, lat)

        Cache layout:
          coords: exp, label, depth, lat, month
          vars:
            - <field>(exp, depth, lat) or (exp, month, depth, lat)
        """
        
        if not self._meta or not self._stats:
            raise ValueError("No experiments added yet (self._meta / self._stats empty).")

        exps   = [m.get("exp", f"exp{i}") for i, m in enumerate(self._meta)]
        labels = [m.get("label", e) for m, e in zip(self._meta, exps)]

        # --- Find reference depth/lat from the first available DataArray among requested fields ---
        depth = lat = None
        depth_name = lat_name = None

        def _pick_depth_lat(da: xr.DataArray):
            nonlocal depth, lat, depth_name, lat_name
            # pick dims by name heuristics (same logic spirit as your resolver)
            dc = next((d for d in da.dims if "depth" in d.lower() or d.lower().startswith("z")), None)
            lc = next((d for d in da.dims if "lat" in d.lower()), None)
            # if monthly_climatology slice, da may have month too; that's fine
            if dc is None:
                # fallback: choose 2D dims other than month/exp
                cand = [d for d in da.dims if d not in ("month", "exp")]
                dc = cand[0] if cand else da.dims[0]
            if lc is None:
                cand = [d for d in da.dims if d not in ("month", "exp", dc)]
                lc = cand[-1] if cand else da.dims[-1]
            depth_name, lat_name = dc, lc
            depth = da[dc].values if dc in da.coords else np.arange(da.sizes[dc], dtype=float)
            lat   = da[lc].values if lc in da.coords else np.arange(da.sizes[lc], dtype=float)

        for st in self._stats:
            for f in fields:
                da = st.get(f, None)
                if isinstance(da, xr.DataArray):
                    _pick_depth_lat(da)
                    break
            if depth is not None and lat is not None:
                break

        if depth is None or lat is None:
            raise ValueError("Could not infer depth/lat coordinates from stored stats.")

        # month coord if any monthly-like field is requested / available
        month = np.arange(1, 13, dtype=np.int16)

        ds = xr.Dataset(
            coords=dict(
                exp=("exp", np.asarray(exps, dtype=object)),
                label=("exp", np.asarray(labels, dtype=object)),
                depth=("depth", np.asarray(depth)),
                lat=("lat", np.asarray(lat)),
                month=("month", month),
            )
        )

        # helper to coerce da -> desired dim order
        def _to_2d(da: xr.DataArray) -> xr.DataArray:
            da2 = da
            # rename to internal canonical dim names for cache
            rename = {}
            if depth_name in da2.dims and depth_name != "depth":
                rename[depth_name] = "depth"
            if lat_name in da2.dims and lat_name != "lat":
                rename[lat_name] = "lat"
            if rename:
                da2 = da2.rename(rename)
            # drop month if present unexpectedly
            if "month" in da2.dims:
                raise ValueError("Expected 2D (depth, lat) but got month dimension.")
            return da2.transpose("depth", "lat")

        def _to_3d_month(da: xr.DataArray) -> xr.DataArray:
            da2 = da
            rename = {}
            if depth_name in da2.dims and depth_name != "depth":
                rename[depth_name] = "depth"
            if lat_name in da2.dims and lat_name != "lat":
                rename[lat_name] = "lat"
            # month dim might be named "month" already; if not, try to detect
            if "month" not in da2.dims:
                md = next((d for d in da2.dims if d.lower() == "month"), None)
                if md and md != "month":
                    rename[md] = "month"
            if rename:
                da2 = da2.rename(rename)
            if "month" not in da2.dims:
                raise ValueError("Expected monthly field with 'month' dimension.")
            return da2.transpose("month", "depth", "lat")

        # --- write each requested field ---
        for f in fields:
            if f == "monthly_climatology" or f == "monthly_interannual_std":
                arr = np.full((len(exps), 12, len(depth), len(lat)), np.nan, dtype=np.float32)
                ok_any = False
                for i, st in enumerate(self._stats):
                    da = st.get(f, None)
                    if not isinstance(da, xr.DataArray):
                        continue
                    da3 = _to_3d_month(da)
                    # align months 1..12 if coords exist
                    if "month" in da3.coords:
                        try:
                            da3 = da3.sel(month=month)
                        except Exception:
                            pass
                    arr[i, :, :, :] = da3.values.astype(np.float32)
                    ok_any = True
                if ok_any:
                    ds[f] = (("exp", "month", "depth", "lat"), arr)

            elif f == "month_of_max":
                arr = np.full((len(exps), len(depth), len(lat)), -1, dtype=np.int16)
                ok_any = False
                for i, st in enumerate(self._stats):
                    da = st.get(f, None)
                    if not isinstance(da, xr.DataArray):
                        continue
                    da2 = _to_2d(da)
                    arr[i, :, :] = da2.values.astype(np.int16)
                    ok_any = True
                if ok_any:
                    ds[f] = (("exp", "depth", "lat"), arr)

            else:
                # treat as float 2D (depth, lat)
                arr = np.full((len(exps), len(depth), len(lat)), np.nan, dtype=np.float32)
                ok_any = False
                for i, st in enumerate(self._stats):
                    da = st.get(f, None)
                    if not isinstance(da, xr.DataArray):
                        continue
                    da2 = _to_2d(da)
                    arr[i, :, :] = da2.values.astype(np.float32)
                    ok_any = True
                if ok_any:
                    ds[f] = (("exp", "depth", "lat"), arr)

        # provenance
        ds.attrs["amoc_meta_json"] = json.dumps(
            dict(
                root_dir=self.root_dir,
                file_key=self.file_key,
                var_mpas=self.var_mpas,
                scale=self.scale,
                weighted=self.weighted,
                calendar=self.calendar,
                time_origin=self.time_origin,
                compute_interannual_std=self.compute_interannual_std,
                per_exp_meta=self._meta,
            ),
            default=str,
        )

        encoding = {}
        if compress:
            for v in ds.data_vars:
                if np.issubdtype(ds[v].dtype, np.floating):
                    encoding[v] = dict(zlib=True, complevel=int(complevel), _FillValue=np.nan)
                else:
                    encoding[v] = dict(zlib=True, complevel=int(complevel))

        ds.to_netcdf(out_nc, encoding=encoding)
        return out_nc

    # ---- from the map2d_ohc variant ----
    def plot_amoc_2d(
        self,
        exp_idx: int = 0,
        field: str = "annual_mean",
        month: Optional[int] = None,
        depth_coord: str = "depth",
        lat_coord: str = "lat",
        cmap: str = "RdBu_r",
        vmin=None, vmax=None,
        add_contours: bool = True,
        clevs: Optional[List[float]] = None,
        n_levels: Optional[int] = None,     # NEW: build discrete levels from vmin/vmax
        extend: str = "neither",            # NEW: colorbar extend for discrete
        tick_every: Optional[int] = None,   # NEW: thin discrete ticks (e.g., 2 -> every 2nd)
        depth_units: str = "m",             # NEW: "m" or "km"
        figsize=(8, 6),
        title: Optional[str] = None,
        cbar_label: str = "Sv",
        yreverse: bool = True,
        save: Optional[str] = None,
    ):
        
        da, dc, lc = self._resolve_field(exp_idx, field, month, depth_coord, lat_coord)
        Z = da.values
        X = da[lc].values
        Y = da[dc].values

        # convert depth to km if requested
        ylabel = "Depth (m)"
        if depth_units.lower() == "km":
            Y = np.asarray(Y, float) / 1000.0
            ylabel = "Depth (km)"

        # defaults for vmin/vmax if not provided
        if vmin is None or vmax is None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                amax = np.nanmax(np.abs(Z))
            if np.isfinite(amax) and amax > 0:
                vmin = -amax if vmin is None else vmin
                vmax =  amax if vmax is None else vmax

        # construct discrete levels if requested
        levels = None
        if clevs is not None:
            levels = np.asarray(clevs, float)
        elif n_levels is not None and vmin is not None and vmax is not None:
            levels = np.linspace(vmin, vmax, int(n_levels)+1)

        fig, ax = plt.subplots(figsize=figsize)

        if levels is not None:
            cmap_obj = mpl.colormaps.get_cmap(cmap)
            norm = BoundaryNorm(levels, cmap_obj.N, clip=False)
            im = ax.pcolormesh(X, Y, Z, shading="auto", cmap=cmap_obj, norm=norm)
            # discrete colorbar
            ticks = levels
            if tick_every and tick_every > 1:
                ticks = levels[::int(tick_every)]
            cb = fig.colorbar(im, ax=ax, boundaries=levels, ticks=ticks, spacing="proportional", extend=extend, pad=0.02)
        else:
            im = ax.pcolormesh(X, Y, Z, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            cb = fig.colorbar(im, ax=ax, pad=0.02)

        # optional contour lines (use same levels if provided)
        if add_contours and (levels is not None and len(levels) > 1):
            cs = ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=0.6, alpha=0.7)
            ax.clabel(cs, fmt="%g", fontsize=9, inline=True)

        ax.set_xlabel("Latitude (°)")
        ax.set_ylabel(ylabel)
        if yreverse:
            ax.invert_yaxis()

        meta = self._meta[exp_idx]
        default_title = f"{meta.get('label','Exp')} – {field.replace('_',' ').title()}"
        if field == "monthly_climatology":
            default_title += f" (Month {int(month)})"
        ax.set_title(title or default_title)

        cb.set_label(cbar_label)
        ax.grid(True, alpha=0.2)

        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")
        return fig, ax

    # ---- from the ohc_2d variant (renamed from plot_amoc_2d_compare) ----
    def plot_amoc_2d_compare_mesh(
        self,
        ref_idx: int,
        tgt_idx,                               # int or list/tuple of ints
        field: str = "annual_mean",
        month: Optional[int] = None,
        depth_coord: str = "depth",
        lat_coord: str = "lat",
        cmap: str = "RdBu_r",                  # reference cmap
        diff_cmap: str = "bwr",                # diffs cmap
        vmin=None, vmax=None,                  # reference range
        diff_vmin=None, diff_vmax=None,        # diffs range
        clevs: Optional[List[float]] = None,   # reference discrete boundaries
        n_levels: Optional[int] = None,        # build reference boundaries
        diff_clevs: Optional[List[float]] = None,  # diffs discrete boundaries
        diff_n_levels: Optional[int] = None,       # build diffs boundaries
        extend: str = "neither",               # reference colorbar extend
        diff_extend: str = "neither",          # diffs colorbar extend
        tick_every: Optional[int] = None,      # thin ticks for ref colorbar
        diff_tick_every: Optional[int] = None, # thin ticks for diff colorbar
        depth_units: str = "m",                # "m" or "km"
        figsize=(18, 6),
        cbar_label: str = "Sv",
        diff_label: str = "ΔSv",
        yreverse: bool = True,
        save: Optional[str] = None,
        fontz: float = 12.0,                    # NEW: global font size scale factor
    ):
        # ---- base fonts (points) scaled by fontz ----
        BASE = dict(title=1.0, label=0.95, tick=0.95, cbar=0.95)
        FS   = {k: v * float(fontz) for k, v in BASE.items()}
    
        # normalize targets
        tgt_list = list(tgt_idx) if isinstance(tgt_idx, (list, tuple)) else [int(tgt_idx)]

        # reference
        da_ref, dc, lc = self._resolve_field(ref_idx, field, month, depth_coord, lat_coord)
        X = da_ref[lc].values
        Y = da_ref[dc].values
        Zref = da_ref.values

        # depth units
        ylabel = "Depth (m)"
        if depth_units.lower() == "km":
            Y = np.asarray(Y, float) / 1000.0
            ylabel = "Depth (km)"

        # defaults for ranges
        if vmin is None or vmax is None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                amax = np.nanmax(np.abs(Zref))
            if np.isfinite(amax) and amax > 0:
                vmin = -amax if vmin is None else vmin
                vmax =  amax if vmax is None else vmax

        # build diffs (on ref grid)
        diffs = []
        for ti in tgt_list:
            da_tgt, _, _ = self._resolve_field(ti, field, month, depth_coord, lat_coord)
            da_tgt_interp = da_tgt.interp({lc: da_ref[lc], dc: da_ref[dc]})
            diffs.append((ti, da_tgt_interp - da_ref))

        # diff range if not provided
        if diff_vmin is None or diff_vmax is None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                all_diff = np.concatenate([np.ravel(D.values) for _, D in diffs]) if diffs else np.array([0.0])
                dmax = np.nanmax(np.abs(all_diff))
            if np.isfinite(dmax) and dmax > 0:
                diff_vmin = -dmax if diff_vmin is None else diff_vmin
                diff_vmax =  dmax if diff_vmax is None else diff_vmax

        # discrete boundaries
        ref_levels = None
        if clevs is not None:
            ref_levels = np.asarray(clevs, float)
        elif n_levels is not None and vmin is not None and vmax is not None:
            ref_levels = np.linspace(vmin, vmax, int(n_levels) + 1)

        diff_levels = None
        if diff_clevs is not None:
            diff_levels = np.asarray(diff_clevs, float)
        elif diff_n_levels is not None and diff_vmin is not None and diff_vmax is not None:
            diff_levels = np.linspace(diff_vmin, diff_vmax, int(diff_n_levels) + 1)

        # figure layout: 1 + len(targets)
        ncols = 1 + len(tgt_list)
        fig_w = max(figsize[0], 6 * ncols)
        fig, axes = plt.subplots(1, ncols, figsize=(fig_w, figsize[1]), constrained_layout=True)
        if ncols == 1:
            axes = [axes]

        # --- panel 1: reference ---
        ax0 = axes[0]
        # --- reference panel ---
        cmap_ref = self._cmap_with_over_under(cmap, n=256, darken_factor=0.8)
        if ref_levels is not None:
            norm_ref = BoundaryNorm(ref_levels, cmap_ref.N, clip=False)  # IMPORTANT
            im0 = ax0.pcolormesh(X, Y, Zref, shading="auto", cmap=cmap_ref, norm=norm_ref)
            cb0 = fig.colorbar(im0, ax=ax0, boundaries=ref_levels,
                               spacing="proportional", extend="both", pad=0.02)
        else:
            im0 = ax0.pcolormesh(X, Y, Zref, shading="auto", cmap=cmap_ref, vmin=vmin, vmax=vmax)
            cb0 = fig.colorbar(im0, ax=ax0, extend="both", pad=0.02)

        if yreverse:
            ax0.invert_yaxis()
        ax0.set_xlabel("Latitude (°)", fontsize=FS["label"])
        ax0.set_ylabel(ylabel, fontsize=FS["label"])
        ax0.tick_params(labelsize=FS["tick"])
        lab_ref = self._meta[ref_idx].get("label", f"Exp {ref_idx}")
        ttl_ref = (
            f"{lab_ref} – {field.replace('_',' ').title()}"
            if field != "monthly_climatology"
            else f"{lab_ref} – Monthly Climatology (Month {int(month)})"
        )
        ax0.set_title(ttl_ref, fontsize=FS["title"])
        cb0.set_label(cbar_label, fontsize=FS["label"])
        cb0.ax.tick_params(labelsize=FS["tick"])
        ax0.grid(True, alpha=0.2)

        # --- panels 2..: diffs vs ref ---
        for j, (ti, D) in enumerate(diffs, start=1):
            ax = axes[j]
            Zd = D.values
            cmap_d = self._cmap_with_over_under(diff_cmap, n=256, darken_factor=0.8)
            if diff_levels is not None:
                norm_d = BoundaryNorm(diff_levels, cmap_d.N, clip=False)     # IMPORTANT
                im = ax.pcolormesh(X, Y, Zd, shading="auto", cmap=cmap_d, norm=norm_d)
                cb = fig.colorbar(im, ax=ax, boundaries=diff_levels,
                                  spacing="proportional", extend="both", pad=0.02)
            else:
                im = ax.pcolormesh(X, Y, Zd, shading="auto", cmap=cmap_d, vmin=diff_vmin, vmax=diff_vmax)
                cb = fig.colorbar(im, ax=ax, extend="both", pad=0.02)

            if yreverse:
                ax.invert_yaxis()
            ax.set_xlabel("Latitude (°)", fontsize=FS["label"])
            if j == 1:
                ax.set_ylabel(ylabel, fontsize=FS["label"])
            ax.tick_params(labelsize=FS["tick"])
            lab_tgt = self._meta[ti].get("label", f"Exp {ti}")
            ax.set_title(f"{lab_tgt} − {lab_ref}", fontsize=FS["title"])
            cb.set_label(diff_label, fontsize=FS["label"])
            cb.ax.tick_params(labelsize=FS["tick"])
            ax.grid(True, alpha=0.2)

        if save:
            plt.savefig(save, dpi=300, bbox_inches="tight")
        return fig, axes
