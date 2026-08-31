"""Configuration and dataset I/O for the online ML analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import xarray as xr
import yaml


def load_config(path: str | Path) -> dict[str, Any]:
    """Load YAML configuration and resolve paths relative to the config file."""
    path = Path(path).expanduser().resolve()
    with path.open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    config["config_path"] = str(path)
    config["config_dir"] = str(path.parent)
    return config


def _expand_files(spec: str | Sequence[str]) -> list[str]:
    patterns = [spec] if isinstance(spec, str) else list(spec)
    files: list[str] = []
    for pattern in patterns:
        expanded = Path(pattern).expanduser()
        if any(character in str(expanded) for character in "*?["):
            files.extend(str(item) for item in sorted(expanded.parent.glob(expanded.name)))
        elif expanded.exists():
            files.append(str(expanded))
    if not files:
        raise FileNotFoundError(f"No files matched: {patterns}")
    return files


def open_dataset(
    spec: str | Sequence[str],
    *,
    chunks: Mapping[str, int] | None = None,
    rename: Mapping[str, str] | None = None,
) -> xr.Dataset:
    """Open one or many files and apply optional dimension/variable renaming."""
    files = _expand_files(spec)
    if len(files) == 1:
        dataset = xr.open_dataset(files[0], chunks=chunks)
    else:
        dataset = xr.open_mfdataset(
            files,
            combine="by_coords",
            parallel=True,
            chunks=chunks,
        )
    valid_renames = {
        old: new for old, new in (rename or {}).items() if old in dataset and old != new
    }
    return dataset.rename(valid_renames)


def open_analysis_datasets(config: Mapping[str, Any]) -> dict[str, xr.Dataset]:
    """Open CTRL, UNET, UNETXTR, and reference datasets from configuration."""
    datasets: dict[str, xr.Dataset] = {}
    defaults = config.get("dataset_options", {})
    for name in ("CTRL", "UNET", "UNETXTR", "reference"):
        entry = config["datasets"][name]
        datasets[name] = open_dataset(
            entry["files"],
            chunks=entry.get("chunks", defaults.get("chunks")),
            rename={**defaults.get("rename", {}), **entry.get("rename", {})},
        )
    return datasets


def select_and_rename_variables(
    dataset: xr.Dataset,
    variable_map: Mapping[str, str],
) -> xr.Dataset:
    """Select configured source variables and rename them to canonical names."""
    missing = [source for source in variable_map.values() if source not in dataset]
    if missing:
        raise KeyError(f"Dataset is missing configured variables: {missing}")
    selected = dataset[list(variable_map.values())]
    return selected.rename({source: canonical for canonical, source in variable_map.items()})


def rename_available_variables(
    dataset: xr.Dataset,
    variable_map: Mapping[str, str],
) -> xr.Dataset:
    """Rename configured variables that exist, allowing optional tendencies."""
    renames = {
        source: canonical
        for canonical, source in variable_map.items()
        if source in dataset and source != canonical
    }
    return dataset.rename(renames)


def output_paths(config: Mapping[str, Any]) -> dict[str, Path]:
    """Create and return processed, figure, and table output directories."""
    root = Path(config["output_root"]).expanduser()
    paths = {
        "root": root,
        "processed": root / "processed",
        "figures": root / "figures",
        "tables": root / "tables",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def save_dataset(dataset: xr.Dataset, path: str | Path) -> Path:
    """Write a diagnostic dataset with compression for floating-point fields."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        name: {"zlib": True, "complevel": 2}
        for name, variable in dataset.data_vars.items()
        if variable.dtype.kind == "f"
    }
    dataset.to_netcdf(path, encoding=encoding)
    return path
