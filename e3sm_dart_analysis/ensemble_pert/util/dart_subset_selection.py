"""Helpers for selecting a representative E3SM-DART ensemble subset."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import xarray as xr


MEMBER_RE = re.compile(r"\.EN(\d{2})\.")


def configure_dask(
    *,
    enabled: bool = True,
    use_distributed: bool = False,
    n_workers: int = 4,
    threads_per_worker: int = 1,
    memory_limit: str = "4GB",
):
    """Configure Dask for this notebook and optionally start a local client."""
    if not enabled:
        return None

    try:
        import dask
    except ImportError:
        print("[WARN] dask is not available; falling back to serial execution.")
        return None

    dask.config.set(
        {
            "array.slicing.split_large_chunks": True,
            "distributed.worker.memory.target": 0.75,
            "distributed.worker.memory.spill": 0.85,
            "distributed.worker.memory.pause": 0.95,
        }
    )

    if not use_distributed:
        return None

    try:
        from dask.distributed import Client
    except ImportError:
        print("[WARN] dask.distributed is not available; using local scheduler.")
        return None

    client = Client(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        memory_limit=memory_limit,
    )
    print(client)
    return client


def files_by_member(paths: Iterable[Path]) -> dict[int, Path]:
    """Return a member-number keyed file map for DART restart files."""
    result = {}
    for path in paths:
        match = MEMBER_RE.search(path.name)
        if match:
            result[int(match.group(1))] = path
    return result


def discover_member_files(rest_dir: Path, date_tag: str) -> pd.DataFrame:
    """Discover paired EAM and ELM files in shared or per-member archives."""
    eam_paths = sorted(rest_dir.glob(f"*.eam.i.{date_tag}.nc"))
    elm_paths = sorted(rest_dir.glob(f"*.elm.r.{date_tag}.nc"))

    # Newer Perlmutter archives place each member under
    # ``ENxx/archive/rest/<date_tag>`` instead of using one shared directory.
    if not eam_paths and not elm_paths:
        member_rest_pattern = f"EN*/archive/rest/{date_tag}"
        eam_paths = sorted(
            path
            for member_rest_dir in rest_dir.glob(member_rest_pattern)
            for path in member_rest_dir.glob(f"*.eam.i.{date_tag}.nc")
        )
        elm_paths = sorted(
            path
            for member_rest_dir in rest_dir.glob(member_rest_pattern)
            for path in member_rest_dir.glob(f"*.elm.r.{date_tag}.nc")
        )

    eam_map = files_by_member(eam_paths)
    elm_map = files_by_member(elm_paths)
    common_members = sorted(set(eam_map) & set(elm_map))
    return pd.DataFrame(
        {
            "member": common_members,
            "eam_file": [str(eam_map[m]) for m in common_members],
            "elm_file": [str(elm_map[m]) for m in common_members],
        }
    )


def inventory(path: str | Path, candidates: Sequence[str]) -> pd.DataFrame:
    """Summarize candidate variables available in one NetCDF file."""
    rows = []
    with xr.open_dataset(path, decode_times=False) as ds:
        for name in candidates:
            if name in ds.variables:
                da = ds[name]
                rows.append(
                    {
                        "variable": name,
                        "dims": str(da.dims),
                        "shape": str(da.shape),
                        "dtype": str(da.dtype),
                        "units": da.attrs.get("units", ""),
                        "long_name": da.attrs.get("long_name", ""),
                    }
                )
    return pd.DataFrame(rows)


def stride_indexers(da: xr.DataArray, max_points: int) -> dict[str, slice]:
    """Build deterministic multidimensional strides capped near max_points."""
    shape = np.asarray(da.shape, dtype=int)
    if shape.size == 0 or np.prod(shape) <= max_points:
        return {}

    strides = np.ones(shape.size, dtype=int)
    while np.prod(np.ceil(shape / strides).astype(int)) > max_points:
        k = int(np.argmax(shape / strides))
        strides[k] += 1

    return {
        dim: slice(None, None, int(step))
        for dim, step in zip(da.dims, strides)
        if step > 1
    }


def _open_dataset_for_features(path: str | Path, *, use_dask: bool):
    kwargs = {
        "decode_times": False,
        "mask_and_scale": True,
        "cache": False,
    }
    if use_dask:
        kwargs["chunks"] = {}
    return xr.open_dataset(path, **kwargs)


def summarize_field(
    da: xr.DataArray,
    max_points: int,
    quantiles: Sequence[float],
) -> dict[str, float]:
    """Summarize a sampled field with robust finite-value filtering."""
    sampled = da.isel(stride_indexers(da, max_points))
    values = np.asarray(sampled.compute().values, dtype=np.float64).ravel()
    values = values[np.isfinite(values)]
    values = values[np.abs(values) < 1.0e30]

    if values.size == 0:
        output = {"mean": np.nan, "std": np.nan}
        output.update({f"q{int(q * 100):02d}": np.nan for q in quantiles})
        return output

    output = {
        "mean": float(values.mean()),
        "std": float(values.std(ddof=1)) if values.size > 1 else 0.0,
    }
    for q, value in zip(quantiles, np.quantile(values, quantiles)):
        output[f"q{int(q * 100):02d}"] = float(value)
    return output


def component_features(
    path: str | Path,
    component: str,
    variables: Sequence[str],
    *,
    max_points: int,
    quantiles: Sequence[float],
    use_dask: bool,
) -> dict[str, float]:
    """Extract compact feature summaries for one component restart file."""
    output = {}
    with _open_dataset_for_features(path, use_dask=use_dask) as ds:
        for variable in variables:
            if variable not in ds.variables:
                continue
            da = ds[variable]
            if not np.issubdtype(da.dtype, np.number):
                continue
            stats = summarize_field(da, max_points, quantiles)
            for statistic, value in stats.items():
                output[f"{component}:{variable}:{statistic}"] = value
    return output


def member_feature_row(
    row: Mapping[str, object],
    *,
    eam_variables: Sequence[str],
    elm_variables: Sequence[str],
    max_points: int,
    quantiles: Sequence[float],
    use_dask: bool,
) -> dict[str, float]:
    """Extract all EAM and ELM features for one paired member."""
    member = int(row["member"])
    feature_row = {"member": member}
    feature_row.update(
        component_features(
            row["eam_file"],
            "EAM",
            eam_variables,
            max_points=max_points,
            quantiles=quantiles,
            use_dask=use_dask,
        )
    )
    feature_row.update(
        component_features(
            row["elm_file"],
            "ELM",
            elm_variables,
            max_points=max_points,
            quantiles=quantiles,
            use_dask=use_dask,
        )
    )
    return feature_row


def compute_member_features(
    member_files: pd.DataFrame,
    *,
    eam_variables: Sequence[str],
    elm_variables: Sequence[str],
    max_points: int,
    quantiles: Sequence[float],
    use_dask: bool = True,
    task_retries: int = 1,
    min_unique: int = 10,
) -> pd.DataFrame:
    """Compute member-level features, using Dask delayed when available."""
    records = list(member_files.to_dict("records"))

    if use_dask:
        try:
            import dask
            from dask.diagnostics import ProgressBar
        except ImportError:
            print("[WARN] dask is unavailable; falling back to serial extraction.")
        else:
            tasks = [
                dask.delayed(member_feature_row)(
                    row,
                    eam_variables=eam_variables,
                    elm_variables=elm_variables,
                    max_points=max_points,
                    quantiles=quantiles,
                    use_dask=use_dask,
                )
                for row in records
            ]
            with ProgressBar():
                rows = dask.compute(*tasks, retries=task_retries)
            return _clean_feature_frame(rows, min_unique=min_unique)

    rows = []
    for j, row in enumerate(records):
        member = int(row["member"])
        print(f"EN{member:02d}: {j + 1}/{len(records)}")
        rows.append(
            member_feature_row(
                row,
                eam_variables=eam_variables,
                elm_variables=elm_variables,
                max_points=max_points,
                quantiles=quantiles,
                use_dask=False,
            )
        )
    return _clean_feature_frame(rows, min_unique=min_unique)


def _clean_feature_frame(
    rows: Sequence[Mapping[str, float]],
    *,
    min_unique: int = 10,
) -> pd.DataFrame:
    features_raw = pd.DataFrame(rows).set_index("member").sort_index()
    features_raw = features_raw.replace([np.inf, -np.inf], np.nan)
    valid = (
        features_raw.notna().all(axis=0)
        & (features_raw.std(axis=0, ddof=1) > 0)
        & (features_raw.nunique(axis=0) >= min_unique)
    )
    return features_raw.loc[:, valid]


def make_subset_scorer(
    z: np.ndarray,
    *,
    weight_mean: float,
    weight_std: float,
    weight_cov: float,
):
    """Create a subset objective scorer bound to a standardized feature matrix."""
    full_mean = z.mean(axis=0)
    full_std = z.std(axis=0, ddof=1)
    full_cov = np.cov(z, rowvar=False)
    full_cov_norm = np.linalg.norm(full_cov, ord="fro")
    if not np.isfinite(full_cov_norm) or full_cov_norm == 0:
        full_cov_norm = 1.0

    def subset_score(indices: Sequence[int]):
        x = z[indices]
        mean_error = float(np.mean((x.mean(axis=0) - full_mean) ** 2))
        std_error = float(np.mean((x.std(axis=0, ddof=1) - full_std) ** 2))
        cov_error = float(
            (
                np.linalg.norm(np.cov(x, rowvar=False) - full_cov, ord="fro")
                / full_cov_norm
            )
            ** 2
        )
        total = (
            weight_mean * mean_error
            + weight_std * std_error
            + weight_cov * cov_error
        )
        return total, {
            "total": total,
            "mean_error": mean_error,
            "std_error": std_error,
            "cov_error": cov_error,
        }

    diagnostics = {
        "full_mean": full_mean,
        "full_std": full_std,
        "full_cov": full_cov,
        "full_cov_norm": full_cov_norm,
    }
    return subset_score, diagnostics


def _score_random_batch(
    batch_seed: int,
    n_trials: int,
    n_members: int,
    subset_size: int,
    subset_score,
):
    rng = np.random.default_rng(batch_seed)
    scores = np.empty(n_trials, dtype=np.float64)
    best_indices = None
    best_score = np.inf
    best_parts = None

    for trial in range(n_trials):
        candidate = np.sort(
            rng.choice(n_members, size=subset_size, replace=False)
        )
        score, parts = subset_score(candidate)
        scores[trial] = score
        if score < best_score:
            best_indices = candidate.copy()
            best_score = score
            best_parts = parts

    return scores, best_indices, best_score, best_parts


def random_search_subsets(
    *,
    n_trials: int,
    n_members: int,
    subset_size: int,
    subset_score,
    random_seed: int,
    batch_size: int = 2_000,
    use_dask: bool = True,
):
    """Search random subset candidates in deterministic, Dask-friendly batches."""
    seed_sequence = np.random.SeedSequence(random_seed)
    batch_sizes = [
        min(batch_size, n_trials - start)
        for start in range(0, n_trials, batch_size)
    ]
    batch_seeds = [
        int(s.generate_state(1)[0]) for s in seed_sequence.spawn(len(batch_sizes))
    ]

    if use_dask:
        try:
            import dask
            from dask.diagnostics import ProgressBar
        except ImportError:
            print("[WARN] dask is unavailable; falling back to serial search.")
        else:
            tasks = [
                dask.delayed(_score_random_batch)(
                    seed,
                    size,
                    n_members,
                    subset_size,
                    subset_score,
                )
                for seed, size in zip(batch_seeds, batch_sizes)
            ]
            with ProgressBar():
                results = dask.compute(*tasks)
            return _combine_random_search_results(results)

    results = [
        _score_random_batch(seed, size, n_members, subset_size, subset_score)
        for seed, size in zip(batch_seeds, batch_sizes)
    ]
    return _combine_random_search_results(results)


def _combine_random_search_results(results):
    trial_scores = np.concatenate([result[0] for result in results])
    best_scores = np.asarray([result[2] for result in results])
    best_pos = int(np.argmin(best_scores))
    _, best_indices, best_score, best_parts = results[best_pos]
    return trial_scores, best_indices, best_score, best_parts


def improve_by_swaps(initial_indices, *, n_members: int, subset_score):
    """Greedily improve a subset by one-member swaps."""
    current = np.sort(initial_indices.copy())
    current_score, current_parts = subset_score(current)
    iterations = 0

    while True:
        selected = set(current.tolist())
        unselected = sorted(set(range(n_members)) - selected)
        replacement = None
        replacement_score = current_score
        replacement_parts = current_parts

        for outgoing in current:
            for incoming in unselected:
                candidate = np.array(
                    sorted((selected - {int(outgoing)}) | {int(incoming)}),
                    dtype=int,
                )
                score, parts = subset_score(candidate)
                if score < replacement_score - 1.0e-14:
                    replacement = candidate
                    replacement_score = score
                    replacement_parts = parts

        if replacement is None:
            break

        current = replacement
        current_score = replacement_score
        current_parts = replacement_parts
        iterations += 1

    return current, current_score, current_parts, iterations
