"""Shared Dask and lazy xarray helpers for diagnostic workflows."""

from pathlib import Path
import sys
from typing import Any, Mapping, Optional

CWD = Path.cwd().resolve()
if CWD.name == "diagnostic":
    DIAGNOSTIC_DIR = CWD
elif (CWD / "diagnostic").is_dir():
    DIAGNOSTIC_DIR = CWD / "diagnostic"
else:
    DIAGNOSTIC_DIR = next((p for p in [CWD, *CWD.parents] if p.name == "diagnostic"), CWD.parent)
if str(DIAGNOSTIC_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTIC_DIR))


DEFAULT_DASK_CHUNKS = {
    "time": 30,
    "lat": 90,
    "lon": 180,
}


def configure_dask(
    *,
    use_distributed: bool = False,
    n_workers: int = 4,
    threads_per_worker: int = 1,
    memory_limit: str = "4GB",
    n_threads: Optional[int] = None,
) -> Optional[Any]:
    """Configure Dask for notebook workflows and optionally start a local client.

    ``n_threads`` caps the threaded scheduler (``use_distributed=False``). Its
    default is one thread per core: on a many-core login node that keeps
    hundreds of read chunks in memory at once.
    """
    try:
        import dask
    except ImportError:
        print("[WARN] dask is not available; xarray will run with its default backend.")
        return None

    dask.config.set(
        {
            "array.slicing.split_large_chunks": True,
            "distributed.worker.memory.target": 0.75,
            "distributed.worker.memory.spill": 0.85,
            "distributed.worker.memory.pause": 0.95,
        }
    )
    if n_threads:
        dask.config.set(num_workers=n_threads)

    if not use_distributed:
        return None

    try:
        from dask.distributed import Client
    except ImportError:
        print("[WARN] dask.distributed is not available; using threaded scheduler only.")
        return None

    client = Client(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        memory_limit=memory_limit,
    )
    print(client)
    return client


def merge_chunks(
    chunks: Optional[Mapping[str, int]] = None,
    **overrides: int,
) -> dict:
    """Return default chunk sizes with optional per-workflow overrides."""
    merged = dict(DEFAULT_DASK_CHUNKS)
    if chunks:
        merged.update(chunks)
    merged.update({k: v for k, v in overrides.items() if v is not None})
    return merged


def _chunk_existing_dims(obj, chunks: Optional[Mapping[str, int]] = None):
    selected = {
        dim: size
        for dim, size in merge_chunks(chunks).items()
        if dim in obj.dims
    }
    if not selected:
        return obj
    if not hasattr(obj, "data_vars"):
        return obj.chunk({d: s for d, s in selected.items() if d in obj.dims})
    # Numeric data variables only: rechunking lazily decoded time bounds (cftime or
    # datetime64 ``time_bnds``) fails in recent xarray, and bounds are small anyway.
    numeric = {
        name: var.chunk({d: s for d, s in selected.items() if d in var.dims})
        for name, var in obj.data_vars.items()
        if var.dtype.kind not in "OmM" and set(var.dims) & set(selected)
    }
    return obj.assign(numeric)


def compute_all(objects: Mapping[Any, Any]) -> dict:
    """Compute a mapping of lazy xarray objects in one Dask call.

    Every source's reduction runs in a single graph, so file reads and
    reductions of different sources proceed in parallel instead of one
    ``.compute()`` after another. Values that are not Dask-backed pass through.
    """
    import dask

    keys = list(objects)
    return dict(zip(keys, dask.compute(*(objects[k] for k in keys))))


def open_dataset_lazy(path, *, chunks: Optional[Mapping[str, int]] = None, **kwargs):
    """Open one NetCDF file lazily and chunk only dimensions it actually has."""
    import xarray as xr

    ds = xr.open_dataset(path, chunks={}, **kwargs)
    return _chunk_existing_dims(ds, chunks)


def open_mfdataset_lazy(paths, *, chunks: Optional[Mapping[str, int]] = None, **kwargs):
    """Open many NetCDF files lazily and chunk only dimensions they actually have."""
    import xarray as xr

    defaults = {
        "combine": "by_coords",
        "parallel": True,
        "chunks": {},
    }
    defaults.update(kwargs)
    ds = xr.open_mfdataset(paths, **defaults)
    return _chunk_existing_dims(ds, chunks)


def persist_if_dask(obj, enabled: bool = True):
    """Persist dask-backed xarray objects when the caller will reuse them."""
    if not enabled:
        return obj
    try:
        return obj.persist()
    except Exception:
        return obj


def netcdf_encoding_for(ds, complevel: int = 1) -> dict:
    """Small compression defaults for robust NetCDF output without huge files."""
    return {
        name: {"zlib": True, "complevel": complevel}
        for name in ds.data_vars
    }
