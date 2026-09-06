"""Per-case, per-variable Figure 3 checkpoints and diagnostic calculation."""

import json
from pathlib import Path
from uuid import uuid4

import numpy as np
import xarray as xr

from . import io, metrics

FIELDS = {
    'rmse_reduction', 'improvement_percent', 'regional_mean_reduction',
    'global_rmse', 'global_pcc', 'improved_grid_fraction', 'degraded_grid_fraction',
}


def signature(case, variable, case_dirs, post_subdir, period, analysis):
    """Exclude the selected variable/case lists so additions preserve old caches."""
    return json.dumps({
        'version': 1, 'case': case, 'variable': variable, 'period': period,
        'sources': {key: str(case_dirs[key] / post_subdir / f'{variable}_{period}.nc')
                    for key in ('CTRL', 'REF', case)},
        'start': analysis['start_date'], 'end': analysis['end_date'],
        'minimum_control_rmse': analysis['minimum_control_rmse'],
        'time_policy': 'exact common input timestamps before daily averaging',
    }, sort_keys=True)


def cache_valid(path, case, variable, expected_signature=None):
    """Read the complete small product to reject incomplete/corrupt checkpoints.

    Undefined correlations (constant series) and masked percentages are allowed.
    """
    if not Path(path).exists():
        return False
    try:
        with xr.open_dataset(path) as cached:
            if not FIELDS.issubset(cached.data_vars):
                return False
            if expected_signature is not None and cached.attrs.get('cache_signature') != expected_signature:
                return False
            selected = cached.sel(experiment=[case], variable=[variable]).load()
            return (selected.rmse_reduction.dims == ('experiment', 'variable', 'lat', 'lon')
                    and bool(np.isfinite(selected.rmse_reduction).any()))
    except (KeyError, OSError, ValueError, RuntimeError):
        return False


def legacy_valid(path, case, variable, case_dirs, period, analysis):
    """Reuse old annual caches only for the original full-period date selection.

    Old files lack a complete signature. This is a best-effort migration using
    their recorded case paths/period; callers can disable it or force recompute.
    """
    start, end = period.split('_')
    import calendar
    first = f'{start[:4]}-{start[4:]}-01'
    last = f'{end[:4]}-{end[4:]}-{calendar.monthrange(int(end[:4]), int(end[4:]))[1]:02d}'
    if analysis['start_date'] != first or analysis['end_date'] != last:
        return False
    if analysis['minimum_control_rmse'] != 1e-12:
        return False
    if not cache_valid(path, case, variable):
        return False
    with xr.open_dataset(path) as cached:
        return ('cache_signature' not in cached.attrs
                and cached.attrs.get('period') == period
                and cached.attrs.get('experiment') == case
                and cached.attrs.get('control_case') == str(case_dirs['CTRL'])
                and cached.attrs.get('reference_case') == str(case_dirs['REF']))


def save_atomic(dataset, path):
    """Publish only complete NetCDF products; interruption preserves prior cache."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f'.{path.stem}.{uuid4().hex}.tmp.nc')
    try:
        io.save_dataset(dataset, temporary)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def calculate(model, reference, control_rmse, area, minimum_control_rmse):
    """Compute the original seven Figure 3 statistics for one daily field."""
    model, reference = xr.align(model, reference, join='exact')
    ml = metrics.spatial_rmse(model, reference)
    reduction = control_rmse - ml
    valid = reduction.notnull()
    horizontal = ['lat', 'lon']
    daily_rmse = metrics.weighted_rmse(model, reference, area, horizontal)
    grid_pcc = metrics.temporal_pearson_correlation(model, reference)
    return xr.Dataset({
        'rmse_reduction': reduction,
        'improvement_percent': 100 * reduction / control_rmse.where(control_rmse > minimum_control_rmse),
        'regional_mean_reduction': metrics.regional_weighted_statistics(reduction, area, reference.lat),
        'global_rmse': daily_rmse.mean('time', skipna=True),
        'global_pcc': metrics.weighted_mean(grid_pcc, area, horizontal),
        'improved_grid_fraction': 100 * (reduction > 0).where(valid).mean(horizontal, skipna=True),
        'degraded_grid_fraction': 100 * (reduction < 0).where(valid).mean(horizontal, skipna=True),
    })


def load_products(pair_paths, cases, variables):
    """Assemble plotting data exclusively from selected small checkpoints."""
    products = []
    for case in cases:
        parts = []
        for variable in variables:
            with xr.open_dataset(pair_paths[case, variable]) as ds:
                parts.append(ds.load())
        products.append(xr.concat(parts, dim='variable', join='exact', combine_attrs='drop_conflicts'))
    return xr.concat(products, dim='experiment', join='exact', combine_attrs='drop_conflicts')
