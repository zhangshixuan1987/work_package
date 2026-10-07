"""Vectorized NaN-aware statistics.

numpy's nanquantile/nanpercentile fall back to a Python loop over every grid
cell for N-D input, which dominates runtime on lat/lon maps.
"""

from __future__ import annotations

import numpy as np
import xarray as xr


def nanquantile_numpy(values, q, axis=-1):
    """Linear-interpolation nanquantile matching np.nanquantile (q leading if 1-D)."""
    values = np.moveaxis(np.asarray(values, dtype=float), axis, -1)
    values = np.sort(values, axis=-1)  # NaNs sort to the end
    count = np.sum(~np.isnan(values), axis=-1)
    q_array = np.asarray(q, dtype=float)
    results = []
    for quantile in np.atleast_1d(q_array):
        position = np.clip(quantile * (count - 1), 0, None)
        low = np.floor(position).astype(int)
        high = np.minimum(low + 1, np.maximum(count - 1, 0))
        lower = np.take_along_axis(values, low[..., None], axis=-1)[..., 0]
        upper = np.take_along_axis(values, high[..., None], axis=-1)[..., 0]
        results.append(np.where(count > 0, lower + (position - low) * (upper - lower), np.nan))
    return results[0] if q_array.ndim == 0 else np.stack(results)


def nanquantile(data: xr.DataArray, q, dim: str) -> xr.DataArray:
    """Drop-in for ``data.quantile(q, dim, skipna=True)`` with the same output layout."""
    scalar = np.ndim(q) == 0
    q_values = float(q) if scalar else [float(value) for value in q]

    def reduce(values):
        result = nanquantile_numpy(values, q_values, axis=-1)
        return result if scalar else np.moveaxis(result, 0, -1)

    result = xr.apply_ufunc(
        reduce, data,
        input_core_dims=[[dim]],
        output_core_dims=[[] if scalar else ["quantile"]],
        dask="parallelized",
        output_dtypes=[float],
        dask_gufunc_kwargs=None if scalar else {"output_sizes": {"quantile": len(q_values)}},
        keep_attrs=False,
    )
    if scalar:
        return result.assign_coords(quantile=q_values)
    return result.assign_coords(quantile=q_values).transpose("quantile", ...)
