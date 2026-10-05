from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


@dataclass(frozen=True)
class CMIPMetricsFixConfig:
    """Configuration for patching selected CMIP mean-climate JSON metrics."""

    metrics_root: Path
    base_dataset: str
    update_dataset: str
    output_version: str
    output_root: Optional[Path] = None

    def dataset_parts(self, dataset: str) -> tuple[str, str, str]:
        parts = tuple(part.strip() for part in dataset.split("."))
        if len(parts) != 3 or any(not part for part in parts):
            raise ValueError(
                f"Invalid dataset name {dataset!r}; expected '<mip>.<exp>.<version>'."
            )
        return parts

    @property
    def base_dir(self) -> Path:
        mip, exp, version = self.dataset_parts(self.base_dataset)
        return Path(self.metrics_root) / mip / exp / version

    @property
    def update_dir(self) -> Path:
        mip, exp, version = self.dataset_parts(self.update_dataset)
        return Path(self.metrics_root) / mip / exp / version

    @property
    def output_dir(self) -> Path:
        mip, exp, _ = self.dataset_parts(self.base_dataset)
        root = Path(self.output_root) if self.output_root is not None else Path(self.metrics_root)
        return root / mip / exp / self.output_version


def _first_match(directory: Path, patterns: Iterable[str]) -> Optional[Path]:
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def find_metric_json(directory: Path, variable: str, version: str) -> Optional[Path]:
    """Find a variable JSON, supporting both dot and underscore version styles."""
    return _first_match(
        Path(directory),
        (
            f"{variable}*.{version}.json",
            f"{variable}*_{version}.json",
        ),
    )


def compare_json_files(old_file: str | Path, new_file: str | Path) -> None:
    """Print structural/value differences between two JSON files."""
    with Path(old_file).open() as handle:
        old_data = json.load(handle)
    with Path(new_file).open() as handle:
        new_data = json.load(handle)

    print(f"[INFO] Comparing:\n  OLD: {old_file}\n  NEW: {new_file}")

    def recurse_compare(left, right, path: str = "") -> None:
        if isinstance(left, dict) and isinstance(right, dict):
            for key in sorted(set(left) | set(right)):
                new_path = f"{path}/{key}" if path else str(key)
                if key not in left:
                    print(f"[ADDED]   {new_path}: {right[key]}")
                elif key not in right:
                    print(f"[REMOVED] {new_path}: {left[key]}")
                else:
                    recurse_compare(left[key], right[key], new_path)
            return
        if left != right:
            print(f"[CHANGED] {path}:\n   OLD: {left}\n   NEW: {right}")

    recurse_compare(old_data, new_data)


def update_cmip_metrics(
    config: CMIPMetricsFixConfig,
    variables: Iterable[str],
    *,
    dry_run: bool = True,
    compare_existing_output: bool = False,
) -> list[Path]:
    """Patch selected variables in a CMIP metrics set from a newer metrics set."""
    written: list[Path] = []
    _, _, base_version = config.dataset_parts(config.base_dataset)
    _, _, update_version = config.dataset_parts(config.update_dataset)

    for variable in variables:
        base_file = find_metric_json(config.base_dir, variable, base_version)
        update_file = find_metric_json(config.update_dir, variable, update_version)

        if base_file is None or update_file is None:
            print(f"[WARN] Missing file for {variable}: base={base_file}, update={update_file}")
            continue

        out_file = config.output_dir / base_file.name.replace(base_version, config.output_version)
        print(f"[INFO] Base file:   {base_file}")
        print(f"[INFO] Update from: {update_file}")
        print(f"[INFO] Output to:   {out_file}")

        if compare_existing_output and out_file.exists():
            compare_json_files(base_file, out_file)

        with base_file.open() as handle:
            base_data = json.load(handle)
        with update_file.open() as handle:
            update_data = json.load(handle)

        update_results = update_data.get("RESULTS", {})
        for model, model_data in base_data.get("RESULTS", {}).items():
            if model not in update_results:
                continue
            for obs_key, obs_data in model_data.items():
                if obs_key not in update_results[model] or not isinstance(obs_data, dict):
                    continue
                for obs_source, source_data in obs_data.items():
                    if obs_source == "source":
                        continue
                    update_source = update_results[model][obs_key].get(obs_source)
                    if not isinstance(source_data, dict) or not isinstance(update_source, dict):
                        continue
                    for region in source_data:
                        if region == "InputClimatologyFileName" or region not in update_source:
                            continue
                        source_data[region] = update_source[region]
                        print(f"[INFO] Updated: {model}/{obs_key}/{obs_source}/{region}")

        if dry_run:
            print("[dry-run] Not writing output.")
            continue

        out_file.parent.mkdir(parents=True, exist_ok=True)
        with out_file.open("w") as handle:
            json.dump(base_data, handle, indent=2, sort_keys=True)
        written.append(out_file)
        print(f"[SUCCESS] Wrote {out_file}")

    return written
