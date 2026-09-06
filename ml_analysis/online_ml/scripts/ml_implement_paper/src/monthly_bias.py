"""Monthly common-pressure zonal bias diagnostics for Figure 6."""

import numpy as np
import xarray as xr

from . import preprocessing


def select_months(dataset, start, end):
    """Label monthly means by YYYYMM, using CF time-bound midpoints if present.

    E3SM may timestamp a monthly mean at the start of the following month.
    Integer month labels also permit matching different model calendars.
    """
    bounds_name = dataset.time.attrs.get('bounds') or dataset.time.attrs.get('climatology')
    if bounds_name in dataset:
        bounds = dataset[bounds_name]
        bound_dim = next(dim for dim in bounds.dims if dim != 'time')
        # Scalar arithmetic supports both numpy datetimes and cftime calendars.
        left = bounds.isel({bound_dim: 0}).values
        right = bounds.isel({bound_dim: 1}).values
        stamp = xr.DataArray([a + (b - a) / 2 for a, b in zip(left, right)], dims='time')
    else:
        stamp = dataset.time
    months = np.asarray(stamp.dt.year * 100 + stamp.dt.month, dtype=int)
    if len(np.unique(months)) != len(months):
        raise ValueError('Expected exactly one sample per month')
    days = np.asarray(stamp.dt.days_in_month, dtype=float)
    result = dataset.assign_coords(time=months, days_in_month=('time', days)).rename(time='month')
    first, last = int(start[:7].replace('-', '')), int(end[:7].replace('-', ''))
    result = result.sel(month=slice(first, last))
    if not result.sizes['month']:
        raise ValueError('No monthly samples in the requested period')
    return result


def pressure_fields(dataset, variables, levels, *, level_dim='lev', pressure_name=None,
                    hybrid_names=None):
    """Interpolate before differencing; reject unrecognized vertical coordinates."""
    # Postprocessed E3SM fields commonly already use plev in Pa.
    if 'plev' in dataset[variables[0]].dims:
        dataset = dataset.rename({'plev': '_source_level'})
        level_dim = '_source_level'
        if not dataset[level_dim].attrs.get('units'):
            # NCO monthly pressure products can omit plev units. Use the same
            # Pa/hPa magnitude convention as the existing preprocessing helper.
            dataset[level_dim].attrs['units'] = (
                'Pa' if float(dataset[level_dim].max()) > 2000 else 'hPa'
            )
        pressure_name = None
        hybrid_names = dict(hybrid_names or {'ps': 'PS', 'hyam': 'hyam', 'hybm': 'hybm', 'p0': 'P0'})
        # Any retained native coefficients do not describe pressure-level fields.
        dataset = dataset.drop_vars([hybrid_names[k] for k in ('hyam', 'hybm') if hybrid_names[k] in dataset])
    sample = dataset[variables[0]]
    if not {'month', 'lat', 'lon', level_dim}.issubset(sample.dims):
        raise ValueError('Use monthly fields on a regular lat/lon grid with a vertical dimension')
    hybrid_names = hybrid_names or {'ps': 'PS', 'hyam': 'hyam', 'hybm': 'hybm', 'p0': 'P0'}
    has_hybrid = all(hybrid_names[k] in dataset for k in ('ps', 'hyam', 'hybm'))
    if not pressure_name and not has_hybrid:
        units = dataset[level_dim].attrs.get('units', '').lower()
        if units not in {'pa', 'hpa', 'mb', 'mbar', 'millibar', 'millibars', 'pascal', 'pascals'}:
            raise ValueError('Pressure interpolation needs pressure units or PS/hyam/hybm metadata')
        if 'hybrid' in dataset[level_dim].attrs.get('standard_name', '').lower():
            raise ValueError('Hybrid levels require PS, hyam, and hybm')
    dataset = dataset.chunk({level_dim: -1}) if sample.chunks is not None else dataset
    result = preprocessing.interpolate_dataset_fields(
        dataset, variables, levels, level_dim=level_dim,
        pressure_name=pressure_name, hybrid_names=hybrid_names,
    )[list(variables)].assign_coords(plev=np.asarray(levels, dtype=float))
    # Mask below-ground targets when surface pressure is available.
    ps_name = hybrid_names['ps']
    if ps_name in dataset:
        result = result.where(result.plev <= preprocessing._pressure_to_hpa(dataset[ps_name]))
    result.plev.attrs.update(units='hPa', positive='down')
    return result


def zonal_bias(experiment, reference, variables=('U', 'V', 'T', 'Q')):
    """Mean paired finite experiment-minus-reference differences over longitude.

    Input is the regular 180x360 grid: equal longitude widths give equal weights
    within each latitude band. Q remains in its input units in the cache.
    """
    experiment, reference = xr.align(experiment, reference, join='exact')
    result = xr.Dataset()
    for name in variables:
        model, ref = experiment[name], reference[name]
        if model.attrs.get('units', '') != ref.attrs.get('units', ''):
            raise ValueError(f'{name}: experiment and reference units differ')
        difference = (model - ref).where(np.isfinite(model) & np.isfinite(ref))
        result[name] = difference.mean('lon', skipna=True)
        result[name].attrs.update(units=model.attrs.get('units', ''), long_name=f'{name} zonal-mean bias (experiment minus REF)')
    result.attrs['bias_definition'] = 'experiment minus REF; paired finite longitude mean'
    return result


def open_monthly_fields(case_dir, post_subdir, climo_subdir, period, variables,
                        start, end, *, hybrid_names=None, pressure_name=None, chunks=None):
    """Open monthly series, falling back to exact single-month climo products.

    Seasonal, annual, and multi-year climatologies are never substitutes for
    individual monthly samples. Use all fields from one mode to retain a
    consistent vertical grid and surface pressure within the case.
    """
    from pathlib import Path

    case_dir = Path(case_dir)
    directory = case_dir / post_subdir
    hybrid_names = hybrid_names or {'ps': 'PS', 'hyam': 'hyam', 'hybm': 'hybm', 'p0': 'P0'}
    files = [directory / f'{variable}_{period}.nc' for variable in variables]
    missing = [path for path in files if not path.exists()]
    mode = 'monthly time series'
    expected = None
    if missing:
        mode = 'single-month climatology files'
        climo_dir = case_dir / climo_subdir
        first, last = period.split('_')
        months = np.arange(np.datetime64(f'{first[:4]}-{first[4:]}', 'M'),
                           np.datetime64(f'{last[:4]}-{last[4:]}', 'M') + 1)
        expected = [int(str(month).replace('-', '')) for month in months]
        files = []
        for month in expected:
            pattern = f'*_{month % 100:02d}_{month}_{month}_climo.nc'
            matches = sorted(climo_dir.glob(pattern))
            if len(matches) != 1:
                raise FileNotFoundError(
                    'Missing monthly time-series fields: ' + ', '.join(map(str, missing))
                    + f'\nFallback needs exactly one single-month file: {climo_dir / pattern}'
                    + f' (found {len(matches)}). Supply these monthly products first.'
                )
            files.append(matches[0])
    else:
        ps_path = directory / f"{hybrid_names['ps']}_{period}.nc"
        if ps_path.exists() and ps_path not in files:
            files.append(ps_path)
    opened = []
    try:
        parts = []
        for index, path in enumerate(files):
            raw = xr.open_dataset(path, chunks=chunks, cache=False)
            opened.append(raw)
            if expected is not None:
                absent = set(variables) - set(raw.data_vars)
                if absent:
                    raise KeyError(f'{path} lacks required fields: {sorted(absent)}')
                years = raw.attrs.get('yrs_averaged')
                if years is not None and str(years).strip() != str(expected[index] // 100):
                    raise ValueError(f'{path}: expected a single-year monthly product, got yrs_averaged={years}')
                bounds_name = raw.time.attrs.get('climatology') or raw.time.attrs.get('bounds')
                if bounds_name in raw:
                    bounds = raw[bounds_name].load()
                    bound_dim = next(dim for dim in bounds.dims if dim != 'time')
                    left = bounds.isel({bound_dim: 0})
                    right = bounds.isel({bound_dim: 1})
                    current = expected[index]
                    following = current + 1 if current % 100 < 12 else (current // 100 + 1) * 100 + 1
                    if (left.size != 1 or right.size != 1
                            or int((left.dt.year * 100 + left.dt.month).item()) != current
                            or int((right.dt.year * 100 + right.dt.month).item()) != following
                            or int(left.dt.day.item()) != 1 or int(right.dt.day.item()) != 1):
                        raise ValueError(f'{path}: bounds do not span exactly the requested month')
                # Verify actual time/bounds, not just the filename.
                month_label = str(expected[index])
                check = select_months(raw, f'{month_label[:4]}-{month_label[4:]}-01',
                                      f'{month_label[:4]}-{month_label[4:]}-31')
                if raw.sizes.get('time') != 1 or check.sizes['month'] != 1:
                    raise ValueError(f'{path}: expected one monthly sample')
                part = check
            else:
                part = select_months(raw, start, end)
            needed = set(variables) | set(hybrid_names.values())
            if pressure_name:
                needed.add(pressure_name)
            parts.append(part[[name for name in part.data_vars if name in needed]])
        if expected is None:
            result = xr.merge(parts, join='exact', compat='no_conflicts')
        else:
            result = xr.concat(parts, dim='month', data_vars='minimal', coords='minimal',
                               compat='equals', join='exact')
            result = result.sel(month=slice(int(start[:7].replace('-', '')),
                                            int(end[:7].replace('-', ''))))
        result.attrs.update(monthly_input_mode=mode, monthly_input_files='\n'.join(map(str, files)))
        result.set_close(lambda: [dataset.close() for dataset in opened])
        return result
    except BaseException:
        for dataset in opened:
            dataset.close()
        raise


def open_era5_monthly(path, variables, start, end, *, chunks=None):
    """Load only requested months/fields from the processed ERA5 monthly file.

    This product uses the model's lat/lon grid and pressure levels in millibars.
    Normalize equivalent unit spellings; no wind or humidity values are rescaled.
    """
    raw = xr.open_dataset(path, cache=False)
    try:
        required = list(variables) + ['PS']
        absent = set(required) - set(raw.data_vars)
        if absent:
            raise KeyError(f'ERA5 lacks required fields: {sorted(absent)}')
        result = select_months(raw[required], start, end)
        # Subset the 41-year archive before constructing the Dask graph.
        if chunks is not None:
            result = result.chunk({('month' if dim == 'time' else dim): size
                                   for dim, size in chunks.items()})
        units = {
            'U': ({'m/s', 'ms**-1', 'ms^-1', 'ms-1'}, 'm/s'),
            'V': ({'m/s', 'ms**-1', 'ms^-1', 'ms-1'}, 'm/s'),
            'T': ({'k', 'kelvin'}, 'K'),
            'Q': ({'kg/kg', 'kgkg**-1', 'kgkg^-1', 'kgkg-1'}, 'kg/kg'),
            'PS': ({'pa', 'pascal', 'pascals'}, 'Pa'),
        }
        for name in required:
            accepted, canonical = units[name]
            unit = result[name].attrs.get('units', '').lower().replace(' ', '')
            if unit not in accepted:
                raise ValueError(f'ERA5 {name}: unsupported units {result[name].attrs.get("units")!r}')
            result[name].attrs['units'] = canonical
        result.attrs.update(monthly_input_mode='ERA5 monthly pressure-level reanalysis',
                            monthly_input_files=str(path), reference_label='ERA5')
        result.set_close(raw.close)
        return result
    except BaseException:
        raw.close()
        raise


def select_era5_pressure(dataset, variables, levels):
    """Select existing ERA5 standard pressure levels exactly, without interpolation."""
    pressure = preprocessing._pressure_to_hpa(dataset.lev.copy())
    missing = sorted(set(levels) - set(pressure.values))
    if missing:
        raise ValueError(f'ERA5 does not contain requested pressure levels [hPa]: {missing}')
    result = (dataset[list(variables)].assign_coords(lev=pressure)
              .sel(lev=list(levels)).rename(lev='plev'))
    result = result.where(result.plev <= preprocessing._pressure_to_hpa(dataset.PS))
    result.plev.attrs.update(units='hPa', positive='down')
    return result


def zonal_comparison(experiment, reference, variables=('U', 'V', 'T', 'Q')):
    """Retain paired zonal means as well as bias for cross-section statistics."""
    result = zonal_bias(experiment, reference, variables)
    for name in variables:
        model, ref = xr.align(experiment[name], reference[name], join='exact')
        valid = np.isfinite(model) & np.isfinite(ref)
        for suffix, field in [('model', model), ('reference', ref)]:
            mean = field.where(valid).mean('lon', skipna=True)
            mean.attrs.update(units=field.attrs.get('units', ''),
                              long_name=f'{name} paired zonal mean ({suffix})')
            result[f'{name}_{suffix}'] = mean
    return result


def cross_section_statistics(model, reference):
    """Cos(latitude)*dp weighted RMSE/PCC of one displayed mean cross section.

    Pressure-cell edges are adjacent-level midpoints, bounded by the selected
    minimum/maximum pressures. Statistics use common finite cells. PCC is NaN
    for a constant field or fewer than two valid cells.
    """
    model, reference = xr.align(model, reference, join='exact')
    model = model.transpose('plev', 'lat').sortby('plev')
    reference = reference.transpose('plev', 'lat').sortby('plev')
    pressure = np.asarray(model.plev.values, dtype=float)
    if len(pressure) < 2 or np.any(np.diff(pressure) <= 0):
        raise ValueError('Cross-section statistics require at least two distinct pressure levels')
    edges = np.r_[pressure[0], (pressure[:-1] + pressure[1:]) / 2, pressure[-1]]
    weights = np.diff(edges)[:, None] * np.cos(np.deg2rad(model.lat.values))[None, :]
    x, y = np.asarray(model, dtype=float), np.asarray(reference, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0)
    if not valid.any():
        return float('nan'), float('nan')
    x, y, w = x[valid], y[valid], weights[valid]
    w = w / w.sum()
    rmse = float(np.sqrt(np.sum(w * (x - y)**2)))
    dx, dy = x - np.sum(w*x), y - np.sum(w*y)
    denominator = np.sqrt(np.sum(w*dx**2) * np.sum(w*dy**2))
    pcc = float(np.clip(np.sum(w*dx*dy) / denominator, -1, 1)) if x.size > 1 and np.ptp(x) > 0 and np.ptp(y) > 0 and denominator > 0 else float('nan')
    return rmse, pcc
