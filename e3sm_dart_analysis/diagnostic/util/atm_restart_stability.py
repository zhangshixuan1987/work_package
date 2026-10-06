"""Cached sanity checks for native EAM atmospheric restart files."""
from __future__ import annotations

import glob
import hashlib
import json
import os
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False


AVAILABLE_RESTART_VARIABLES = {
    "PS": {"level_dim": None, "suggested_range": (10000.0, 120000.0), "description": "surface pressure"},
    "T": {"level_dim": "lev", "suggested_range": (100.0, 400.0), "description": "air temperature"},
    "Q": {"level_dim": "lev", "suggested_range": (0.0, 0.1), "description": "specific humidity"},
    "U": {"level_dim": "lev", "suggested_range": (-250.0, 250.0), "description": "zonal wind"},
    "V": {"level_dim": "lev", "suggested_range": (-250.0, 250.0), "description": "meridional wind"},
}

AVAILABLE_METRICS = (
    "minimum",
    "maximum",
    "mean",
    "standard_deviation",
    "nan_count",
    "inf_count",
    "below_range_count",
    "above_range_count",
    "flagged_fraction",
    "out_of_range_fraction_map",
    "ensemble_mean",
    "ensemble_minimum",
    "ensemble_maximum",
    "rms_change",
    "maximum_absolute_change",
)


@dataclass(frozen=True)
class RestartVariableSpec:
    name: str
    level_index: Optional[int] = None
    level_dim: Optional[str] = None
    valid_min: Optional[float] = None
    valid_max: Optional[float] = None

    @property
    def label(self):
        return self.name if self.level_index is None else f"{self.name}_lev{self.level_index:02d}"


def discover_restart_files(case_dir, time_tags, ensemble_ids=None):
    """Return common EN## members and a time/member restart-file mapping."""
    case_dir = Path(case_dir)
    requested = set(ensemble_ids) if ensemble_ids is not None else None
    by_time = {}
    for tag in time_tags:
        patterns = (
            case_dir / "EN*" / "archive" / "rest" / tag / f"*.eam.i.{tag}.nc",
            case_dir / tag / f"*.eam.i.{tag}.nc",
        )
        files = []
        for pattern in patterns:
            files = sorted(glob.glob(str(pattern)))
            if files:
                break
        mapping = {}
        for path in files:
            match = re.search(r"(?:^|\.)(EN\d+)(?:\.|$)", Path(path).name)
            if match is None:
                continue
            member = match.group(1)
            if requested is None or member in requested:
                mapping[member] = path
        if not mapping:
            raise FileNotFoundError(f"No EAM restart files found for {tag} under {case_dir}")
        by_time[tag] = mapping
    common = sorted(set.intersection(*(set(mapping) for mapping in by_time.values())))
    if not common:
        raise ValueError("No ensemble members are common to every selected restart time")
    return common, by_time


def _extract_field(dataset, spec):
    if spec.name not in dataset:
        raise KeyError(f"{spec.name!r} not found; choose one of {list(dataset.data_vars)}")
    field = dataset[spec.name]
    if "time" in field.dims:
        field = field.isel(time=0, drop=True)
    if spec.level_index is not None:
        level_dim = spec.level_dim or next((name for name in ("lev", "plev", "level") if name in field.dims), None)
        if level_dim is None:
            raise ValueError(f"No vertical dimension found for {spec.label}")
        field = field.isel({level_dim: spec.level_index}, drop=True)
    if len(field.dims) != 1 or field.dims[0] not in {"ncol", "ncol_d"}:
        raise ValueError(f"{spec.label} did not reduce to ncol or ncol_d; remaining dims={field.dims}")
    return np.asarray(field.values, dtype=float), str(field.attrs.get("units", "")), field.dims[0]


def _request_hash(case_name, time_tags, requests, ensemble_ids):
    payload = {
        "case": case_name,
        "time_tags": list(time_tags),
        "requests": [asdict(item) for item in requests],
        "ensemble_ids": None if ensemble_ids is None else list(ensemble_ids),
    }
    return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:12]


def restart_stability_cache_path(cache_dir, case_name, time_tags, requests, ensemble_ids=None):
    key = _request_hash(case_name, time_tags, requests, ensemble_ids)
    return Path(cache_dir) / f"{case_name}_{time_tags[0]}_{time_tags[-1]}_{key}_restart_stability.nc"


def compute_restart_stability(case_dir, time_tags, requests, ensemble_ids=None):
    """Compute member statistics, spatial flag fractions, and inter-restart jumps."""
    requests = tuple(requests)
    time_tags = tuple(time_tags)
    if not requests:
        raise ValueError("At least one restart variable request is required")
    members, files = discover_restart_files(case_dir, time_tags, ensemble_ids=ensemble_ids)
    nr, nt, nm = len(requests), len(time_tags), len(members)
    shape = (nr, nt, nm)
    summary = {name: np.full(shape, np.nan) for name in ("minimum", "maximum", "mean", "standard_deviation", "flagged_fraction")}
    counts = {name: np.zeros(shape, dtype=np.int64) for name in ("nan_count", "inf_count", "below_range_count", "above_range_count", "total_count")}
    transitions = tuple(f"{a}_to_{b}" for a, b in zip(time_tags[:-1], time_tags[1:]))
    rms_change = np.full((nr, len(transitions), nm), np.nan)
    max_change = np.full_like(rms_change, np.nan)
    spatial_products = None
    lat = lon = None
    units = []

    for ir, spec in enumerate(requests):
        previous = {}
        request_units = ""
        for it, tag in enumerate(time_tags):
            flag_sum = valid_sum = None
            ens_sum = ens_count = ens_min = ens_max = None
            for im, member in enumerate(members):
                with xr.open_dataset(files[tag][member], decode_times=False) as dataset:
                    values, request_units, horizontal_dim = _extract_field(dataset, spec)
                    if lat is None:
                        suffix = "_d" if horizontal_dim == "ncol_d" else ""
                        lat = np.asarray(dataset[f"lat{suffix}"].values, dtype=float)
                        lon = np.asarray(dataset[f"lon{suffix}"].values, dtype=float)
                if spatial_products is None:
                    npoint = values.size
                    spatial_products = {
                        "out_of_range_fraction_map": np.full((nr, nt, npoint), np.nan),
                        "ensemble_mean": np.full((nr, nt, npoint), np.nan),
                        "ensemble_minimum": np.full((nr, nt, npoint), np.nan),
                        "ensemble_maximum": np.full((nr, nt, npoint), np.nan),
                    }
                if values.size != spatial_products["ensemble_mean"].shape[-1]:
                    raise ValueError("Restart grids differ across selected files")

                finite = np.isfinite(values)
                nan_mask = np.isnan(values)
                inf_mask = np.isinf(values)
                below = finite & (values < spec.valid_min) if spec.valid_min is not None else np.zeros(values.shape, bool)
                above = finite & (values > spec.valid_max) if spec.valid_max is not None else np.zeros(values.shape, bool)
                flagged = nan_mask | inf_mask | below | above
                finite_values = values[finite]
                if finite_values.size:
                    summary["minimum"][ir, it, im] = finite_values.min()
                    summary["maximum"][ir, it, im] = finite_values.max()
                    summary["mean"][ir, it, im] = finite_values.mean()
                    summary["standard_deviation"][ir, it, im] = finite_values.std()
                counts["nan_count"][ir, it, im] = nan_mask.sum()
                counts["inf_count"][ir, it, im] = inf_mask.sum()
                counts["below_range_count"][ir, it, im] = below.sum()
                counts["above_range_count"][ir, it, im] = above.sum()
                counts["total_count"][ir, it, im] = values.size
                summary["flagged_fraction"][ir, it, im] = flagged.mean()

                if flag_sum is None:
                    flag_sum = np.zeros(values.shape, dtype=np.int64)
                    valid_sum = np.zeros(values.shape, dtype=np.int64)
                    ens_sum = np.zeros(values.shape, dtype=float)
                    ens_count = np.zeros(values.shape, dtype=np.int64)
                    ens_min = np.full(values.shape, np.inf)
                    ens_max = np.full(values.shape, -np.inf)
                flag_sum += flagged
                valid_sum += 1
                ens_sum[finite] += values[finite]
                ens_count[finite] += 1
                ens_min[finite] = np.minimum(ens_min[finite], values[finite])
                ens_max[finite] = np.maximum(ens_max[finite], values[finite])

                if it > 0 and member in previous:
                    pair_valid = finite & np.isfinite(previous[member])
                    if pair_valid.any():
                        difference = values[pair_valid] - previous[member][pair_valid]
                        rms_change[ir, it - 1, im] = np.sqrt(np.mean(difference ** 2))
                        max_change[ir, it - 1, im] = np.max(np.abs(difference))
                previous[member] = values

            spatial_products["out_of_range_fraction_map"][ir, it] = flag_sum / valid_sum
            spatial_products["ensemble_mean"][ir, it] = np.divide(ens_sum, ens_count, out=np.full_like(ens_sum, np.nan), where=ens_count > 0)
            spatial_products["ensemble_minimum"][ir, it] = np.where(ens_count > 0, ens_min, np.nan)
            spatial_products["ensemble_maximum"][ir, it] = np.where(ens_count > 0, ens_max, np.nan)
        units.append(request_units)

    coords = {
        "request": [item.label for item in requests],
        "snapshot": list(time_tags),
        "member": members,
        "point": np.arange(lat.size),
        "lat": ("point", lat),
        "lon": ("point", lon),
        "transition": list(transitions),
        "variable": ("request", [item.name for item in requests]),
        "level_index": ("request", [-1 if item.level_index is None else item.level_index for item in requests]),
        "valid_min": ("request", [np.nan if item.valid_min is None else item.valid_min for item in requests]),
        "valid_max": ("request", [np.nan if item.valid_max is None else item.valid_max for item in requests]),
        "units": ("request", units),
    }
    dataset = xr.Dataset(coords=coords)
    for name, values in summary.items():
        dataset[name] = (("request", "snapshot", "member"), values)
    for name, values in counts.items():
        dataset[name] = (("request", "snapshot", "member"), values)
    for name, values in spatial_products.items():
        dataset[name] = (("request", "snapshot", "point"), values)
    dataset["rms_change"] = (("request", "transition", "member"), rms_change)
    dataset["maximum_absolute_change"] = (("request", "transition", "member"), max_change)
    dataset.attrs.update(case_name=Path(case_dir).name, source="native EAM restart files")
    return dataset


def cache_restart_stability(path, *, case_dir=None, time_tags=None, requests=None, ensemble_ids=None, force_compute=False):
    path = Path(path)
    if path.is_file() and not force_compute:
        return xr.load_dataset(path)
    if case_dir is None or time_tags is None or requests is None:
        raise ValueError("case_dir, time_tags, and requests are required for a missing cache")
    dataset = compute_restart_stability(case_dir, time_tags, requests, ensemble_ids=ensemble_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    dataset.to_netcdf(tmp)
    os.replace(tmp, path)
    return dataset


def plot_restart_summary(dataset, *, show=True):
    """Plot flagged fractions for every request, snapshot, and member."""
    requests = list(dataset.request.values.astype(str))
    snapshots = list(dataset.snapshot.values.astype(str))
    members = list(dataset.member.values.astype(str))
    fig, axes = plt.subplots(len(requests), 1, figsize=(max(7, len(snapshots) * 1.4), 2.7 * len(requests)), squeeze=False)
    for ir, request in enumerate(requests):
        axis = axes[ir, 0]
        values = dataset.flagged_fraction.sel(request=request).transpose("member", "snapshot")
        image = axis.imshow(values, vmin=0, vmax=1, cmap="magma", aspect="auto")
        axis.set_yticks(np.arange(len(members)), members)
        axis.set_xticks(np.arange(len(snapshots)), snapshots, rotation=30, ha="right")
        axis.set_title(f"{request}: fraction outside configured range or non-finite")
        fig.colorbar(image, ax=axis, label="flagged fraction")
    fig.tight_layout()
    if show:
        plt.show()
    return fig


def plot_flag_fraction_maps(dataset, requests=None, snapshots=None, *, show=True):
    """Plot the fraction of members flagged at each native-grid column."""
    request_names = list(dataset.request.values.astype(str)) if requests is None else list(requests)
    snapshot_names = list(dataset.snapshot.values.astype(str)) if snapshots is None else list(snapshots)
    projection = ccrs.PlateCarree() if HAS_CARTOPY else None
    fig, axes = plt.subplots(
        len(request_names), len(snapshot_names),
        figsize=(6 * len(snapshot_names), 3.4 * len(request_names)),
        squeeze=False,
        subplot_kw={"projection": projection} if HAS_CARTOPY else {},
    )
    scatter = None
    for i, request in enumerate(request_names):
        for j, snapshot in enumerate(snapshot_names):
            axis = axes[i, j]
            values = dataset.out_of_range_fraction_map.sel(request=request, snapshot=snapshot)
            kwargs = {"transform": ccrs.PlateCarree()} if HAS_CARTOPY else {}
            scatter = axis.scatter(dataset.lon, dataset.lat, c=values, s=2, vmin=0, vmax=1, cmap="magma", **kwargs)
            if HAS_CARTOPY:
                axis.set_global()
                axis.coastlines(linewidth=0.6)
            else:
                axis.set_xlim(-180, 180)
                axis.set_ylim(-90, 90)
            axis.set_title(f"{request} | {snapshot}")
    fig.colorbar(scatter, ax=axes.ravel().tolist(), orientation="horizontal", fraction=0.04, pad=0.07, label="fraction of members flagged")
    fig.suptitle("Restart stability flags")
    fig.subplots_adjust(bottom=0.14, top=0.90, hspace=0.22)
    if show:
        plt.show()
    return fig
