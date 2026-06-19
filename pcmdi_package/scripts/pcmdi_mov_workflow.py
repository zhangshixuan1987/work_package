from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import xarray as xr

from pcmdi_mov_reader import EMOVDiagReader, ModeFileSpec
from movs_plotter import ExtrapropicalModeMapPlotter, MultimodelPCTimeSeriesPlotter


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
    model_maps = [dataset_var(ds, model_var, model_fallback) * unit_scale for _, ds in model_cases]
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
                model_pvals.append(dataset_var(ds, model_var, model_fallback))
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
