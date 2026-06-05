from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Mapping, Optional

try:
    from .metrics_group_merger import PCMDIRun
except ImportError:
    from metrics_group_merger import PCMDIRun


def _first_input_climatology_filename(data: Mapping) -> Optional[str]:
    if isinstance(data, dict):
        for key, value in data.items():
            if key == "InputClimatologyFileName":
                return str(value)
            found = _first_input_climatology_filename(value)
            if found:
                return found
    elif isinstance(data, list):
        for item in data:
            found = _first_input_climatology_filename(item)
            if found:
                return found
    return None


def _resolve_existing_path(
    path_or_name: str,
    search_roots: Iterable[str | Path],
    variable: Optional[str] = None,
) -> Optional[Path]:
    path = Path(path_or_name)
    if path.is_absolute() and path.exists():
        return path

    for root in search_roots:
        root = Path(root)
        candidates = [root / path_or_name, root / path.name]
        if variable:
            candidates.extend(
                [
                    root / variable / path.name,
                    root / "climo" / variable / path.name,
                    root / "model" / "clim" / variable / path.name,
                    root / "observations" / "clim" / variable / path.name,
                ]
            )
        candidates.append(root / "climo_ref" / path.name)
        for candidate in candidates:
            if candidate.exists():
                return candidate
    return None


def _find_local_climo_file(directory: Path, variable: str) -> Optional[Path]:
    patterns = (
        f"*.Amon.{variable}.*.nc",
        f"*{variable}*.nc",
    )
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def _find_metric_json(
    mean_climate_dir: Path,
    variable: str,
    model_name: str,
    case_id: str,
) -> Optional[Path]:
    patterns = (
        f"{variable}*.{model_name}.{case_id}.json",
        f"{variable}*{model_name}*{case_id}.json",
        f"{variable}*.{case_id}.json",
    )
    for pattern in patterns:
        matches = sorted(mean_climate_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


def build_experiments_from_pcmdi_runs(
    runs: Iterable[PCMDIRun],
    variables: Iterable[str],
    *,
    run_type: str = "model_vs_obs",
    search_roots: Iterable[str | Path] = (),
    search_roots_by_model: Optional[Mapping[str, Iterable[str | Path]]] = None,
    include_reference: bool = True,
    reference_name: str = "ERA5",
    strict: bool = False,
) -> dict[str, dict[str, Path]]:
    """Build a map-plot EXPERIMENTS dictionary from raw PCMDI metric JSONs."""
    experiments: dict[str, dict[str, Path]] = {}
    search_roots_by_model = search_roots_by_model or {}

    for run in runs:
        label = run.output_name or run.model_name
        mean_dir = run.metrics_data_dir(run_type=run_type) / "mean_climate"
        extra_roots = list(search_roots_by_model.get(label, ()))
        run_search_roots = [
            mean_dir,
            mean_dir / "climo",
            mean_dir / "climo_ref",
            *extra_roots,
            *search_roots,
        ]
        case_id = run.metrics_case_id or ""
        experiments.setdefault(label, {})

        for variable in variables:
            if variable == "tau_mag":
                continue
            metric_json = _find_metric_json(mean_dir, variable, run.model_name, case_id)
            if metric_json is None:
                message = f"[WARN] Missing metric JSON for {label}/{variable} in {mean_dir}"
                if strict:
                    raise FileNotFoundError(message)
                print(message)
                continue

            with metric_json.open() as handle:
                metric_data = json.load(handle)

            model_file_name = _first_input_climatology_filename(metric_data)
            if model_file_name:
                model_path = _resolve_existing_path(model_file_name, run_search_roots, variable)
                if model_path is None:
                    model_path = _find_local_climo_file(mean_dir / "climo", variable)
                if model_path is not None:
                    experiments[label][variable] = model_path
                else:
                    print(f"[WARN] Could not resolve {label}/{variable}: {model_file_name}")

            if include_reference:
                if variable in experiments.get(reference_name, {}):
                    continue
                ref_info = metric_data.get("References", {}).get("default", {})
                ref_path_text = ref_info.get("file_path") or ref_info.get("template")
                if ref_path_text:
                    ref_path = _resolve_existing_path(ref_path_text, run_search_roots, variable)
                    if ref_path is None:
                        ref_path = _find_local_climo_file(mean_dir / "climo_ref", variable)
                    if ref_path is not None:
                        experiments.setdefault(reference_name, {})[variable] = ref_path
                    else:
                        print(f"[WARN] Could not resolve {reference_name}/{variable}: {ref_path_text}")

    if include_reference and reference_name in experiments:
        return {
            reference_name: experiments[reference_name],
            **{name: paths for name, paths in experiments.items() if name != reference_name},
        }
    return experiments
