from __future__ import annotations

import os
import string
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter

from pcmdi_mov_reader import EMOVDiagReader, ModeFileSpec
from movs_plotter import add_sig_dots, ExtrapropicalModeMapPlotter, MultimodelPCTimeSeriesPlotter


VALID_PSL_TOKENS = {"DJF", "MAM", "JJA", "SON", "monthly", "yearly"}
VALID_TS_TOKENS = {"monthly", "yearly"}


def resolve_mode_spec(
    mode_dict: Mapping[str, Mapping[str, Any]],
    mode: str,
    token: str,
    period: Tuple[int, int],
) -> tuple[Mapping[str, Any], ModeFileSpec]:
    """Validate a mode/token pair and return its PCMDI file spec."""
    if mode not in mode_dict:
        raise KeyError(f"mode={mode!r} not in mode_dict. Valid: {list(mode_dict)}")

    cfg = mode_dict[mode]
    var = cfg["var"]
    eof = cfg["eof"]
    valid_tokens = VALID_PSL_TOKENS if var == "psl" else VALID_TS_TOKENS
    if token not in valid_tokens:
        raise ValueError(f"token={token!r} invalid for var={var!r}. Valid: {sorted(valid_tokens)}")

    return cfg, ModeFileSpec(mode=mode, var=var, eof=eof, period=period, season_or_freq=token)


def format_mode_filename(
    template: str,
    *,
    mode: str,
    token: str,
    var: str,
    eof: str,
) -> str:
    return template.format(mode=mode, token=token, var=var, eof=eof)


def plot_pcmdi_mode_outputs(
    *,
    reader: EMOVDiagReader,
    map_plotter: ExtrapropicalModeMapPlotter,
    pc_plotter: Optional[MultimodelPCTimeSeriesPlotter],
    model_tags: Sequence[str],
    mode_dict: Mapping[str, Mapping[str, Any]],
    mode: str,
    token: str,
    period: Tuple[int, int],
    pattern_config: Optional[Mapping[str, Any]] = None,
    teleconnection_config: Optional[Mapping[str, Any]] = None,
    pc_config: Optional[Mapping[str, Any]] = None,
    member_labels: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Plot the standard PCMDI variability-mode products from saved PMP NetCDF files.

    The config dictionaries are intentionally notebook-friendly: they can include
    ``which`` and ``filename_template`` plus any keyword accepted by the matching
    plotter method.
    """
    if not model_tags:
        raise ValueError("model_tags is empty.")

    mode_cfg, spec = resolve_mode_spec(mode_dict, mode, token, period)
    var = mode_cfg["var"]
    eof = mode_cfg["eof"]
    lat_bnds = mode_cfg.get("lat_bnds")
    lon_bnds = mode_cfg.get("lon_bnds")
    labels = list(member_labels) if member_labels is not None else list(model_tags)

    # Fail early with explicit paths before doing any plotting.
    reader.obs_path(model_tags[0], spec)
    reader.model_path(model_tags[0], spec)

    outputs: Dict[str, Any] = {"spec": spec}

    if pattern_config is not None:
        cfg = dict(pattern_config)
        which = cfg.pop("which", "eof")
        member_dim = cfg.pop("member_dim", "member")
        filename = format_mode_filename(
            cfg.pop("filename_template", "{mode}_{token}_{var}_{eof}_MODE_PATTERN_region.pdf"),
            mode=mode,
            token=token,
            var=var,
            eof=eof,
        )
        dd = reader.build_multimodel_stack(model_tags, spec, which=which, member_dim=member_dim)
        outputs["pattern_data"] = dd
        outputs["pattern"] = map_plotter.plot_multimodel_mode_pattern_with_stats(
            mode=mode,
            token=token,
            obs_map=dd["reference"],
            model_stack=dd["hist"],
            member_labels=labels,
            member_dim=member_dim,
            filename=filename,
            region_lat_bounds=lat_bnds,
            region_lon_bounds=lon_bnds,
            **cfg,
        )

    if teleconnection_config is not None:
        cfg = dict(teleconnection_config)
        which = cfg.pop("which", "slope")
        member_dim = cfg.pop("member_dim", "member")
        filename = format_mode_filename(
            cfg.pop("filename_template", "{mode}_{token}_{var}_{eof}_TELECONNECTION_global.pdf"),
            mode=mode,
            token=token,
            var=var,
            eof=eof,
        )
        dd = reader.build_multimodel_stack(model_tags, spec, which=which, member_dim=member_dim)
        outputs["teleconnection_data"] = dd
        outputs["teleconnection"] = map_plotter.plot_multimodel_teleconnection_with_stats(
            mode=mode,
            token=token,
            obs_map=dd["reference"],
            model_stack=dd["hist"],
            member_labels=labels,
            member_dim=member_dim,
            filename=filename,
            **cfg,
        )

    if pc_config is not None:
        if pc_plotter is None:
            raise ValueError("pc_plotter is required when pc_config is provided.")
        cfg = dict(pc_config)
        which = cfg.pop("which", "pc")
        filename = format_mode_filename(
            cfg.pop("filename_template", "{mode}_{token}_{var}_{eof}_PC_timeseries.pdf"),
            mode=mode,
            token=token,
            var=var,
            eof=eof,
        )
        dd = reader.build_multimodel_stack(model_tags, spec, which=which, member_dim=cfg.get("member_dim", "member"))
        outputs["pc_data"] = dd
        outputs["pc"] = pc_plotter.plot_multimodel_pc_timeseries_with_stats(
            mode=mode,
            token=token,
            obs_pc=dd["reference"],
            model_stack=dd["hist"],
            member_labels=labels,
            filename=filename,
            **cfg,
        )

    return outputs


def dataset_var(ds: xr.Dataset, var: str, fallback: Optional[str] = None) -> xr.DataArray:
    """Return a variable from a raw-analysis result dataset with optional fallback."""
    if var in ds:
        return ds[var]
    if fallback is not None and fallback in ds:
        return ds[fallback]
    keys = list(ds.data_vars)
    msg = f"Variable {var!r} not found. Available: {keys}"
    if fallback is not None:
        msg = f"Variables {var!r}/{fallback!r} not found. Available: {keys}"
    raise KeyError(msg)


def build_raw_multimodel_stack(
    all_cases: Sequence[tuple[str, xr.Dataset]],
    *,
    obs_var: str,
    model_var: Optional[str] = None,
    obs_fallback: Optional[str] = None,
    model_fallback: Optional[str] = None,
    unit_scale: float = 1.0,
    member_dim: str = "member",
) -> Dict[str, xr.DataArray]:
    """Build the Obs + model stack expected by the shared multimodel plotters."""
    if len(all_cases) < 2:
        raise ValueError("all_cases must contain one observation and at least one model.")

    model_var = model_var or obs_var
    obs_name, obs_ds = all_cases[0]
    model_cases = list(all_cases[1:])

    obs = dataset_var(obs_ds, obs_var, obs_fallback) * unit_scale
    model_maps = []
    for _, ds in model_cases:
        m = dataset_var(ds, model_var, model_fallback) * unit_scale
        # Align model grid exactly to obs grid to prevent spurious
        # bilinear re-interpolation inside the plotter.
        if "lat" in m.coords and "lon" in m.coords:
            if not (
                m.lat.size == obs.lat.size
                and m.lon.size == obs.lon.size
                and np.array_equal(m.lat.values, obs.lat.values)
                and np.array_equal(m.lon.values, obs.lon.values)
            ):
                m = m.interp(lat=obs.lat, lon=obs.lon)
        model_maps.append(m)
    model_labels = [name for name, _ in model_cases]
    stack = xr.concat(model_maps, dim=member_dim).assign_coords({member_dim: model_labels})

    return {
        "reference": obs,
        "hist": stack,
        "reference_label": obs_name,
        "member_labels": model_labels,
    }


def build_raw_pvalue_stack(
    all_cases: Sequence[tuple[str, xr.Dataset]],
    *,
    obs_var: Optional[str],
    model_var: Optional[str],
    obs_fallback: Optional[str] = None,
    model_fallback: Optional[str] = None,
    member_dim: str = "member",
) -> Dict[str, Optional[xr.DataArray]]:
    """Build optional p-value maps for stippling; missing variables return None."""
    if obs_var is None and model_var is None:
        return {"reference": None, "hist": None}

    obs_pval = None
    if obs_var is not None:
        try:
            obs_pval = dataset_var(all_cases[0][1], obs_var, obs_fallback)
        except KeyError:
            obs_pval = None

    model_pvals = []
    model_labels = []
    if model_var is not None:
        for name, ds in all_cases[1:]:
            try:
                pv = dataset_var(ds, model_var, model_fallback)
                # Align p-value grid to obs p-value grid when available.
                if obs_pval is not None and "lat" in pv.coords and "lon" in pv.coords:
                    if not (
                        pv.lat.size == obs_pval.lat.size
                        and pv.lon.size == obs_pval.lon.size
                        and np.array_equal(pv.lat.values, obs_pval.lat.values)
                        and np.array_equal(pv.lon.values, obs_pval.lon.values)
                    ):
                        pv = pv.interp(lat=obs_pval.lat, lon=obs_pval.lon)
                model_pvals.append(pv)
                model_labels.append(name)
            except KeyError:
                model_pvals.append(None)

    if model_pvals and all(p is not None for p in model_pvals):
        stack = xr.concat(model_pvals, dim=member_dim).assign_coords({member_dim: model_labels})
    else:
        stack = None

    return {"reference": obs_pval, "hist": stack}


def plot_raw_mode_pattern(
    *,
    map_plotter: ExtrapropicalModeMapPlotter,
    all_cases: Sequence[tuple[str, xr.Dataset]],
    mode: str,
    token: str,
    filename: str,
    region_lat_bounds: Optional[Tuple[float, float]],
    region_lon_bounds: Optional[Tuple[float, float]],
    obs_var: str,
    model_var: Optional[str] = None,
    obs_fallback: Optional[str] = None,
    model_fallback: Optional[str] = None,
    pval_obs_var: Optional[str] = None,
    pval_model_var: Optional[str] = None,
    pval_obs_fallback: Optional[str] = None,
    pval_model_fallback: Optional[str] = None,
    unit_scale: float = 1.0,
    member_dim: str = "member",
    show_significance: bool = False,
    **plot_kwargs: Any,
) -> tuple[Any, Any]:
    """Plot raw-analysis map products using the same panel plotter as PCMDI files."""
    dd = build_raw_multimodel_stack(
        all_cases,
        obs_var=obs_var,
        model_var=model_var,
        obs_fallback=obs_fallback,
        model_fallback=model_fallback,
        unit_scale=unit_scale,
        member_dim=member_dim,
    )
    pvals = build_raw_pvalue_stack(
        all_cases,
        obs_var=pval_obs_var,
        model_var=pval_model_var,
        obs_fallback=pval_obs_fallback,
        model_fallback=pval_model_fallback,
        member_dim=member_dim,
    )
    return map_plotter.plot_multimodel_mode_pattern_with_stats(
        mode=mode,
        token=token,
        obs_map=dd["reference"],
        model_stack=dd["hist"],
        member_labels=dd["member_labels"],
        filename=filename,
        region_lat_bounds=region_lat_bounds,
        region_lon_bounds=region_lon_bounds,
        obs_pval_map=pvals["reference"],
        model_pval_stack=pvals["hist"],
        member_dim=member_dim,
        show_significance=show_significance and (pvals["reference"] is not None or pvals["hist"] is not None),
        **plot_kwargs,
    )


def plot_raw_teleconnection(
    *,
    map_plotter: ExtrapropicalModeMapPlotter,
    all_cases: Sequence[tuple[str, xr.Dataset]],
    mode: str,
    token: str,
    filename: str,
    obs_var: str,
    model_var: Optional[str] = None,
    obs_fallback: Optional[str] = None,
    model_fallback: Optional[str] = None,
    pval_obs_var: Optional[str] = None,
    pval_model_var: Optional[str] = None,
    pval_obs_fallback: Optional[str] = None,
    pval_model_fallback: Optional[str] = None,
    unit_scale: float = 1.0,
    member_dim: str = "member",
    show_significance: bool = False,
    **plot_kwargs: Any,
) -> tuple[Any, Any]:
    """Plot raw-analysis global teleconnection maps using the shared PCMDI plotter."""
    dd = build_raw_multimodel_stack(
        all_cases,
        obs_var=obs_var,
        model_var=model_var,
        obs_fallback=obs_fallback,
        model_fallback=model_fallback,
        unit_scale=unit_scale,
        member_dim=member_dim,
    )
    pvals = build_raw_pvalue_stack(
        all_cases,
        obs_var=pval_obs_var,
        model_var=pval_model_var,
        obs_fallback=pval_obs_fallback,
        model_fallback=pval_model_fallback,
        member_dim=member_dim,
    )
    return map_plotter.plot_multimodel_teleconnection_with_stats(
        mode=mode,
        token=token,
        obs_map=dd["reference"],
        model_stack=dd["hist"],
        member_labels=dd["member_labels"],
        filename=filename,
        obs_pval_map=pvals["reference"],
        model_pval_stack=pvals["hist"],
        member_dim=member_dim,
        show_significance=show_significance and (pvals["reference"] is not None or pvals["hist"] is not None),
        **plot_kwargs,
    )


def plot_raw_pc_timeseries(
    *,
    pc_plotter: MultimodelPCTimeSeriesPlotter,
    all_cases: Sequence[tuple[str, xr.Dataset]],
    mode: str,
    token: str,
    filename: str,
    pc_var: str = "pc_proj",
    pc_fallback: Optional[str] = None,
    member_dim: str = "member",
    **plot_kwargs: Any,
) -> tuple[Any, Any]:
    """Plot cached raw-analysis PC time series through the shared PC plotter."""
    if "dpi" in plot_kwargs and "fig_dpi" not in plot_kwargs:
        plot_kwargs["fig_dpi"] = plot_kwargs.pop("dpi")
    dd = build_raw_multimodel_stack(
        all_cases,
        obs_var=pc_var,
        model_var=pc_var,
        obs_fallback=pc_fallback,
        model_fallback=pc_fallback,
        unit_scale=1.0,
        member_dim=member_dim,
    )
    return pc_plotter.plot_multimodel_pc_timeseries_with_stats(
        mode=mode,
        token=token,
        obs_pc=dd["reference"],
        model_stack=dd["hist"],
        member_labels=dd["member_labels"],
        filename=filename,
        member_dim=member_dim,
        pc_var=pc_var,
        **plot_kwargs,
    )


def wrap_lon_180(da: xr.DataArray, lon_name: str = "lon") -> xr.DataArray:
    """Convert longitude from [0, 360) to [-180, 180) and sort."""
    lon = da[lon_name]
    if float(lon.max()) > 180:
        da = da.assign_coords({lon_name: ((lon + 180) % 360) - 180})
    return da.sortby(lon_name)


def as_latlon(da: xr.DataArray) -> xr.DataArray:
    """Return a lon-wrapped DataArray with lat/lon as trailing dimensions."""
    da = wrap_lon_180(da)
    if "lat" not in da.dims or "lon" not in da.dims:
        raise ValueError(f"Expected lat/lon dimensions, got {da.dims}")
    other_dims = [dim for dim in da.dims if dim not in ("lat", "lon")]
    return da.transpose(*other_dims, "lat", "lon")


def subset_latlon_domain(
    da: xr.DataArray,
    *,
    lat_bnds: Optional[Tuple[float, float]],
    lon_bnds: Optional[Tuple[float, float]],
) -> xr.DataArray:
    """Subset DataArray to a possibly wrapped lat/lon domain."""
    da = as_latlon(da)
    if lat_bnds is not None:
        da = da.sel(lat=slice(min(lat_bnds), max(lat_bnds)))
    if lon_bnds is not None:
        lon0, lon1 = lon_bnds
        if lon0 <= lon1:
            da = da.sel(lon=slice(lon0, lon1))
        else:
            da = xr.concat(
                [da.sel(lon=slice(lon0, 180)), da.sel(lon=slice(-180, lon1))],
                dim="lon",
            ).sortby("lon")
    return da


def calc_eof_svd(da: xr.DataArray, eof_num: int = 1) -> xr.DataArray:
    """Calculate an EOF pattern using latitude-weighted SVD."""
    da = as_latlon(da)
    if "time" not in da.dims:
        raise ValueError(f"EOF bootstrap requires a time dimension, got {da.dims}")

    da = da - da.mean("time")
    lat = da.lat.values
    weights = np.sqrt(np.cos(np.deg2rad(lat)))[:, None]
    weighted = da.values * weights[None, :, :]
    flat = weighted.reshape(weighted.shape[0], -1)
    valid = np.isfinite(flat).all(axis=0)
    if valid.sum() < eof_num:
        raise ValueError("Not enough valid grid cells for EOF bootstrap.")

    x = flat[:, valid]
    x = x - x.mean(axis=0, keepdims=True)
    _u, _s, vt = np.linalg.svd(x, full_matrices=False)
    eof_flat = np.full(flat.shape[1], np.nan, dtype=float)
    eof_flat[valid] = vt[eof_num - 1]
    eof_vals = eof_flat.reshape(da.sizes["lat"], da.sizes["lon"]) / weights
    return xr.DataArray(
        eof_vals,
        coords={"lat": da.lat, "lon": da.lon},
        dims=("lat", "lon"),
        name="eof",
    )


def align_and_rescale_eof(eof: xr.DataArray, target: xr.DataArray) -> xr.DataArray:
    """Align EOF sign/amplitude to a target EOF using least-squares projection."""
    eof = as_latlon(eof)
    target = as_latlon(target)
    if not (
        np.array_equal(eof.lon.values, target.lon.values)
        and np.array_equal(eof.lat.values, target.lat.values)
    ):
        target = target.interp(lon=eof.lon, lat=eof.lat)

    x = eof.values
    y = target.values
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() == 0:
        return eof
    denom = np.sum(x[mask] ** 2)
    if denom == 0 or not np.isfinite(denom):
        return eof
    scale = np.sum(x[mask] * y[mask]) / denom
    if not np.isfinite(scale):
        scale = 1.0
    return eof * scale


def bootstrap_diff_pvalue(
    *,
    case_anom: xr.DataArray,
    ref_anom: xr.DataArray,
    case_target: xr.DataArray,
    ref_target: xr.DataArray,
    eof_num: int = 1,
    n_boot: int = 300,
    seed: int = 42,
) -> xr.DataArray:
    """
    Bootstrap p-value for EOF pattern difference.

    At each grid cell, this tests whether the bootstrapped difference
    ``case_eof - ref_eof`` has a consistently positive or negative sign.
    """
    rng = np.random.default_rng(seed)
    case_anom = as_latlon(case_anom)
    ref_anom = as_latlon(ref_anom)
    case_target = as_latlon(case_target)
    ref_target = as_latlon(ref_target)
    valid = np.isfinite(case_target.values) & np.isfinite(ref_target.values)
    gt_count = np.zeros(case_target.shape, dtype=np.int32)
    lt_count = np.zeros(case_target.shape, dtype=np.int32)
    n_case = case_anom.sizes["time"]
    n_ref = ref_anom.sizes["time"]

    for iboot in range(n_boot):
        case_boot = case_anom.isel(time=rng.integers(0, n_case, size=n_case))
        ref_boot = ref_anom.isel(time=rng.integers(0, n_ref, size=n_ref))
        case_eof = align_and_rescale_eof(calc_eof_svd(case_boot, eof_num=eof_num), case_target)
        ref_eof = align_and_rescale_eof(calc_eof_svd(ref_boot, eof_num=eof_num), ref_target)
        diff = (case_eof - ref_eof).values
        gt_count += valid & (diff >= 0)
        lt_count += valid & (diff <= 0)
        if (iboot + 1) % 50 == 0 or (iboot + 1) == n_boot:
            print(f"  bootstrap {iboot + 1}/{n_boot}")

    pvals = 2.0 * np.minimum(
        (gt_count + 1) / (n_boot + 1),
        (lt_count + 1) / (n_boot + 1),
    )
    pvals = np.where(valid, np.clip(pvals, 0.0, 1.0), np.nan)
    return xr.DataArray(pvals, coords=case_target.coords, dims=case_target.dims, name="pvalue")


def calc_pattern_corr_rmse(
    ref: np.ndarray,
    test: np.ndarray,
    w2d: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Calculate pattern correlation and RMSE against a reference panel."""
    mask = np.isfinite(ref) & np.isfinite(test)
    if w2d is not None:
        mask = mask & np.isfinite(w2d)
    if mask.sum() == 0:
        return np.nan, np.nan

    ref1 = ref[mask]
    test1 = test[mask]
    wt = None if w2d is None else w2d[mask]
    if wt is None:
        ref_mean = ref1.mean()
        test_mean = test1.mean()
        ref_anom = ref1 - ref_mean
        test_anom = test1 - test_mean
        denom = np.sqrt(np.sum(ref_anom**2) * np.sum(test_anom**2))
        pcorr = np.nan if denom == 0 else np.sum(ref_anom * test_anom) / denom
        rmse = np.sqrt(np.mean((test1 - ref1) ** 2))
    else:
        ref_mean = np.average(ref1, weights=wt)
        test_mean = np.average(test1, weights=wt)
        ref_anom = ref1 - ref_mean
        test_anom = test1 - test_mean
        denom = np.sqrt(np.sum(wt * ref_anom**2) * np.sum(wt * test_anom**2))
        pcorr = np.nan if denom == 0 else np.sum(wt * ref_anom * test_anom) / denom
        rmse = np.sqrt(np.average((test1 - ref1) ** 2, weights=wt))
    return float(pcorr), float(rmse)


def plot_eof_difference_significance(
    *,
    obs_ds: xr.Dataset,
    model_results: Mapping[str, xr.Dataset],
    model_list: Sequence[str],
    mode: str,
    season: str,
    mode_info: Mapping[str, Any],
    field_var: str,
    fig_dir: str,
    unit_scale: float = 1.0,
    unit_label: str = "",
    obs_label: str = "Obs",
    config: Optional[Mapping[str, Any]] = None,
) -> tuple[Any, list[Any], list[tuple[str, xr.DataArray, Optional[xr.DataArray], Optional[float]]]]:
    """Plot EOF patterns with bootstrap stippling for case-reference differences."""
    cfg = {
        "ref_case": "OBS",
        "compare_cases": None,
        "sig_level": 0.05,
        "n_boot": 300,
        "random_seed": 42,
        "use_cache": True,
        "save_cache": True,
        "overwrite_cache": False,
        "compute_missing_cache": True,
        "cache_dir_name": "eof_diff_sig_cache",
        "dot_density": 2,
        "dot_size": 8,
        "pattern_levels": np.linspace(-5, 5, 11),
        "plot_field": "pattern",
        "cmap": "RdBu_r",
        "font_size": 18,
        "ncols": 2,
        "figsize_per_panel": (6, 4),
        "xtick_step": 20.0,
        "ytick_step": 10.0,
        "axis_label_size": None,
        "axis_label_pad": 4,
        "fig_prefix": None,
        "fig_format": "pdf",
        "dpi": 300,
        "ref_pattern_var": "eof",
        "model_pattern_var": "eof",
        "ref_frac_var": "frac",
        "model_frac_var": "frac",
        "pvalue_source": "bootstrap",
        "ref_pval_var": None,
        "model_pval_var": None,
        "pval_fallback_var": None,
    }
    if config is not None:
        cfg.update(dict(config))
    plot_field = str(cfg["plot_field"]).lower()
    if plot_field not in {"difference", "pattern"}:
        raise ValueError("config['plot_field'] must be 'difference' or 'pattern'.")
    pvalue_source = str(cfg["pvalue_source"]).lower()
    if pvalue_source not in {"auto", "dataset", "bootstrap"}:
        raise ValueError("config['pvalue_source'] must be 'auto', 'dataset', or 'bootstrap'.")

    ref_case = cfg["ref_case"]
    ref_is_obs = str(ref_case).lower() in {"obs", "observation", "observations", "reference"}
    ref_label = obs_label if ref_is_obs else str(ref_case)
    compare_cases = (
        [case for case in model_list if case in model_results and (ref_is_obs or case != ref_case)]
        if cfg["compare_cases"] is None
        else list(cfg["compare_cases"])
    )
    if not ref_is_obs and ref_case not in model_results:
        raise KeyError(f"ref_case={ref_case!r} is not OBS and is not in model_results.")
    if not compare_cases:
        raise ValueError("No model cases found to compare against reference case.")

    lat_bnds = mode_info.get("lat_bnds")
    lon_bnds = mode_info.get("lon_bnds")
    eof_num = int(mode_info.get("eof_num", 1))
    anom_key = f"{field_var}_anom"
    ref_ds = obs_ds if ref_is_obs else model_results[ref_case]
    if pvalue_source == "bootstrap" and anom_key not in ref_ds:
        raise KeyError(f"{anom_key!r} is required for EOF bootstrap but was not saved for {ref_label}.")

    ref_target = subset_latlon_domain(
        dataset_var(ref_ds, cfg["ref_pattern_var"]) * unit_scale,
        lat_bnds=lat_bnds,
        lon_bnds=lon_bnds,
    )
    ref_anom = None
    if anom_key in ref_ds:
        ref_anom = subset_latlon_domain(ref_ds[anom_key] * unit_scale, lat_bnds=lat_bnds, lon_bnds=lon_bnds)
    ref_frac = float(dataset_var(ref_ds, cfg["ref_frac_var"])) * 100 if cfg["ref_frac_var"] in ref_ds else None

    lat_plot = ref_target.lat.values
    lon_plot = ref_target.lon.values
    w = np.cos(np.deg2rad(lat_plot)).astype(float)
    w = w / np.nanmean(w)
    w2d = w[:, None] * np.ones((lat_plot.size, lon_plot.size), dtype=float)
    if plot_field == "difference":
        ref_plot = xr.zeros_like(ref_target)
    else:
        ref_plot = ref_target
    plot_records = [(ref_label, ref_plot, None, ref_frac)]
    stat_targets = [ref_target]
    cache_dir = os.path.join(fig_dir, cfg["cache_dir_name"])

    def cache_path(case_name: str, seed: int) -> str:
        safe = {
            "mode": str(mode).replace(" ", "_"),
            "season": str(season).replace(" ", "_"),
            "field": str(field_var).replace(" ", "_"),
            "case": str(case_name).replace(" ", "_"),
            "ref": str(ref_label).replace(" ", "_"),
            "ref_pattern": str(cfg["ref_pattern_var"]).replace(" ", "_"),
            "model_pattern": str(cfg["model_pattern_var"]).replace(" ", "_"),
        }
        fname = (
            f"{safe['mode']}_{safe['season']}_{safe['field']}_"
            f"{safe['case']}_vs_{safe['ref']}_eof{eof_num}_"
            f"{safe['ref_pattern']}_to_{safe['model_pattern']}_"
            f"nboot{int(cfg['n_boot'])}_seed{seed}_pval.nc"
        )
        return os.path.join(cache_dir, fname)

    def compatible(pval: xr.DataArray, case_name: str, seed: int, target: xr.DataArray) -> bool:
        expected_attrs = {
            "case": str(case_name),
            "reference_case": str(ref_label),
            "mode": str(mode),
            "season": str(season),
            "field_var": str(field_var),
            "ref_pattern_var": str(cfg["ref_pattern_var"]),
            "model_pattern_var": str(cfg["model_pattern_var"]),
            "eof_num": eof_num,
            "n_boot": int(cfg["n_boot"]),
            "seed": int(seed),
        }
        for key, expected in expected_attrs.items():
            if pval.attrs.get(key) != expected:
                return False
        return (
            pval.dims == target.dims
            and np.array_equal(pval.lon.values, target.lon.values)
            and np.array_equal(pval.lat.values, target.lat.values)
        )

    def get_dataset_pvalue(ds: xr.Dataset, pval_var: Optional[str], case_name: str) -> Optional[xr.DataArray]:
        candidates = [pval_var, cfg["pval_fallback_var"]]
        for candidate in candidates:
            if candidate and candidate in ds:
                pval = subset_latlon_domain(ds[candidate], lat_bnds=lat_bnds, lon_bnds=lon_bnds)
                if not (
                    np.array_equal(pval.lon.values, ref_target.lon.values)
                    and np.array_equal(pval.lat.values, ref_target.lat.values)
                ):
                    pval = pval.interp(lon=ref_target.lon, lat=ref_target.lat)
                print(f"Loading saved p-values for {case_name} from dataset variable {candidate!r}")
                return pval
        return None

    def get_pvalue(case_name: str, case_ds: xr.Dataset, case_anom: Optional[xr.DataArray], case_target: xr.DataArray, seed: int) -> Optional[xr.DataArray]:
        if pvalue_source in {"auto", "dataset"}:
            pval = get_dataset_pvalue(case_ds, cfg["model_pval_var"], case_name)
            if pval is not None:
                return pval
            if pvalue_source == "dataset":
                raise KeyError(
                    f"No saved p-value variable found for {case_name}. "
                    f"Tried {cfg['model_pval_var']!r} and {cfg['pval_fallback_var']!r}."
                )

        if ref_anom is None or case_anom is None:
            raise KeyError(f"{anom_key!r} is required for EOF bootstrap but was not saved for {case_name} or {ref_label}.")
        os.makedirs(cache_dir, exist_ok=True)
        path = cache_path(case_name, seed)
        if cfg["use_cache"] and os.path.exists(path) and not cfg["overwrite_cache"]:
            try:
                cached = xr.open_dataarray(path)
                pval_cached = cached.load()
                cached.close()
                if compatible(pval_cached, case_name, seed, case_target):
                    print(f"Loading cached EOF-difference p-values -> {path}")
                    return pval_cached
                print(f"Ignoring incompatible cached EOF-difference p-values -> {path}")
            except Exception as exc:
                print(f"Could not read cached EOF-difference p-values ({exc}); recomputing -> {path}")

        if cfg["use_cache"] and not cfg["overwrite_cache"] and not cfg["compute_missing_cache"]:
            print(f"No compatible cached EOF-difference p-values found for {case_name}; skipping recompute.")
            return None

        print(f"Computing EOF-difference p-values -> {path}")
        pval = bootstrap_diff_pvalue(
            case_anom=case_anom,
            ref_anom=ref_anom,
            case_target=case_target,
            ref_target=ref_target,
            eof_num=eof_num,
            n_boot=int(cfg["n_boot"]),
            seed=seed,
        )
        pval.attrs.update(
            {
                "description": f"Bootstrap p-value for EOF pattern difference: {case_name} - {ref_label}",
                "case": str(case_name),
                "reference_case": str(ref_label),
                "mode": str(mode),
                "season": str(season),
                "field_var": str(field_var),
                "ref_pattern_var": str(cfg["ref_pattern_var"]),
                "model_pattern_var": str(cfg["model_pattern_var"]),
                "eof_num": eof_num,
                "n_boot": int(cfg["n_boot"]),
                "seed": int(seed),
                "note": "Two-sided bootstrap sign test of EOF pattern difference at each grid cell.",
            }
        )
        if cfg["save_cache"]:
            pval.to_netcdf(path)
            print(f"Saved cached EOF-difference p-values -> {path}")
        return pval

    for case_name in compare_cases:
        case_ds = model_results[case_name]
        if pvalue_source == "bootstrap" and anom_key not in case_ds:
            raise KeyError(f"{anom_key!r} is required for EOF bootstrap but was not saved for {case_name}.")
        case_target = subset_latlon_domain(
            dataset_var(case_ds, cfg["model_pattern_var"]) * unit_scale,
            lat_bnds=lat_bnds,
            lon_bnds=lon_bnds,
        )
        case_anom = None
        if anom_key in case_ds:
            case_anom = subset_latlon_domain(case_ds[anom_key] * unit_scale, lat_bnds=lat_bnds, lon_bnds=lon_bnds)
        if not (
            np.array_equal(case_target.lon.values, ref_target.lon.values)
            and np.array_equal(case_target.lat.values, ref_target.lat.values)
        ):
            case_target = case_target.interp(lon=ref_target.lon, lat=ref_target.lat)
            if case_anom is not None:
                case_anom = case_anom.interp(lon=ref_target.lon, lat=ref_target.lat)
        seed = int(cfg["random_seed"]) + len(plot_records)
        pval = get_pvalue(case_name, case_ds, case_anom, case_target, seed)
        frac = float(dataset_var(case_ds, cfg["model_frac_var"])) * 100 if cfg["model_frac_var"] in case_ds else None
        plot_target = case_target - ref_target if plot_field == "difference" else case_target
        plot_records.append((case_name, plot_target, pval, frac))
        stat_targets.append(case_target)
        if pval is not None:
            sig_pct = float(np.nanmean(pval.values < float(cfg["sig_level"])) * 100.0)
            print(f"{case_name} - {ref_label}: significant grid cells p<{float(cfg['sig_level']):.2f}: {sig_pct:.1f}%")
        else:
            print(f"{case_name} - {ref_label}: no EOF-difference p-values available; no stippling.")

    levels = cfg["pattern_levels"]
    if levels is None:
        all_eof = np.concatenate([rec[1].values.ravel() for rec in plot_records])
        finite = np.isfinite(all_eof)
        vmax = float(np.nanpercentile(np.abs(all_eof[finite]), 98)) if np.any(finite) else 1.0
        if vmax == 0 or not np.isfinite(vmax):
            vmax = 1.0
        levels = np.linspace(-vmax, vmax, 11)

    ncols = int(cfg["ncols"])
    nrows = int(np.ceil(len(plot_records) / ncols))
    figsize = cfg["figsize_per_panel"]
    fig = plt.figure(figsize=(figsize[0] * ncols, figsize[1] * nrows + 0.9))
    proj = ccrs.PlateCarree(central_longitude=0.0)
    data_crs = ccrs.PlateCarree()
    axes = [fig.add_subplot(nrows, ncols, i + 1, projection=proj) for i in range(nrows * ncols)]
    fig.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.14, wspace=0.05, hspace=0.40)

    extent = [
        float(ref_target.lon.min()),
        float(ref_target.lon.max()),
        float(ref_target.lat.min()),
        float(ref_target.lat.max()),
    ]
    xtick_step = float(cfg["xtick_step"])
    ytick_step = float(cfg["ytick_step"])
    xticks = np.arange(
        np.ceil(extent[0] / xtick_step) * xtick_step,
        np.floor(extent[1] / xtick_step) * xtick_step + xtick_step * 0.5,
        xtick_step,
    )
    yticks = np.arange(
        np.ceil(extent[2] / ytick_step) * ytick_step,
        np.floor(extent[3] / ytick_step) * ytick_step + ytick_step * 0.5,
        ytick_step,
    )

    font_size = float(cfg["font_size"])
    axis_label_size = float(cfg["axis_label_size"]) if cfg["axis_label_size"] is not None else font_size * 0.85
    axis_label_pad = float(cfg["axis_label_pad"])
    panel_labels = list(string.ascii_lowercase)
    im = None
    for idx, ax in enumerate(axes):
        if idx >= len(plot_records):
            ax.set_visible(False)
            continue
        case_name, eof_field, pval, frac = plot_records[idx]
        im = ax.contourf(
            eof_field.lon.values,
            eof_field.lat.values,
            eof_field.values,
            levels=levels,
            cmap=cfg["cmap"],
            extend="both",
            transform=data_crs,
        )
        ax.contour(
            eof_field.lon.values,
            eof_field.lat.values,
            eof_field.values,
            levels=levels,
            colors="k",
            linewidths=0.3,
            alpha=0.65,
            transform=data_crs,
        )
        if pval is not None:
            add_sig_dots(
                ax,
                pval,
                sig_level=float(cfg["sig_level"]),
                dot_density=int(cfg["dot_density"]),
                dot_size=float(cfg["dot_size"]),
                transform=data_crs,
            )
        ax.coastlines(linewidth=0.6)
        ax.set_extent(extent, crs=data_crs)
        ax.set_xticks(xticks, crs=data_crs)
        ax.set_yticks(yticks, crs=data_crs)
        ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".0f"))
        ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".0f"))
        ax.tick_params(labelsize=font_size * 0.85)
        ax.gridlines(linewidth=0.4, color="gray", alpha=0.5, linestyle="--", draw_labels=False)
        if idx % ncols == 0:
            ax.set_ylabel("Latitude", fontsize=axis_label_size, labelpad=axis_label_pad)
        ax.set_xlabel("Longitude", fontsize=axis_label_size, labelpad=axis_label_pad)
        frac_str = f" ({frac:.1f}%)" if frac is not None and np.isfinite(frac) else ""
        if plot_field == "difference" and idx != 0:
            title_name = f"{case_name} - {ref_label}"
        else:
            title_name = case_name
        ax.set_title(f"({panel_labels[idx]}) {title_name}{frac_str}", fontsize=font_size, loc="left", pad=6)
        if idx != 0:
            pcorr, rmse = calc_pattern_corr_rmse(ref_target.values, stat_targets[idx].values, w2d=w2d)
            print(f"{case_name} vs {ref_label}: PCC={pcorr:.3f}, RMSE={rmse:.3f}")
            ax.text(
                0.05,
                0.05,
                f"PCC = {pcorr:.2f}\nRMSE = {rmse:.2f}",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=font_size * 0.85,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.6", alpha=0.5),
                zorder=5,
            )

    cbar = fig.colorbar(im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.10, aspect=45, ticks=levels)
    if plot_field == "difference":
        cbar_base_label = "EOF amplitude difference"
    else:
        cbar_base_label = "EOF amplitude"
    cbar.set_label(f"{cbar_base_label} ({unit_label})" if unit_label else cbar_base_label, fontsize=font_size)
    cbar.ax.tick_params(labelsize=font_size * 0.85)

    fig_prefix = cfg["fig_prefix"]
    if fig_prefix is None:
        if plot_field == "difference":
            fig_prefix = f"{mode}_{season}_eof_pattern_diff_sig_vs_{ref_label}_bootstrap_p{cfg['sig_level']}"
        else:
            fig_prefix = f"{mode}_{season}_eof_pattern_model_sig_vs_{ref_label}_bootstrap_p{cfg['sig_level']}"
    fig_path = os.path.join(fig_dir, f"{fig_prefix}.{cfg['fig_format']}")
    fig.savefig(fig_path, dpi=int(cfg["dpi"]), bbox_inches="tight", pad_inches=0.05)
    print(f"Saved -> {fig_path}")
    return fig, axes, plot_records
