from __future__ import annotations
"""Global-mean ocean time-series comparisons of alternate spin-ups vs Full-CPL
(MPAS-Analysis time series + surface fluxes): extraction to data/ caches and figures.

  FigureExtractSpec / FigureDataCollector / OCNDiagnosticsPlotter / PlotRequest
        <- timeseries_alternate/ocn_sfc_ts and tke_amoc (identical code)
  OHUFigureExtractSpec / OHUFigureDataCollector / OHUDiagnosticsPlotter
        <- timeseries_alternate/ohu_ts   (uses the v2 flux extractor: *V2 names below)
  OHCFigureDataCollectorConfig / OHCFigureDataCollector / OHCDiagnosticsPlotter
        <- timeseries_alternate/ohc_ts
"""
from dataclasses import dataclass
from typing import Optional, Tuple, List, Dict, Optional, Any, Sequence, Literal
import os, re, glob
import numpy as np
import xarray as xr
import cftime
import json
import pandas as pd
import math
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats, signal
from matplotlib.ticker import (
    FixedLocator, FixedFormatter, AutoMinorLocator,
    MultipleLocator, NullFormatter
)
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from util.experiments import build_exp_subdirs_from_run_dict
from util.ocn_timeseries import (
    OCNTimeSeriesExtractor,
    ExtractOptions,
    SeriesRequest,
    OCNRegistry,
    extract_all,
)
from util.flux_timeseries import (
    FluxTimeSeriesExtractor,
    FluxExtractOptions,
    FluxSeriesRequest,
    FluxRegistry,
    extract_all_flux,
)
from dataclasses import dataclass, replace
from typing import Optional, Tuple, List, Dict, Optional, Any, Literal,  Sequence
import string
from util.ohc_timeseries import (
    OHCExtractConfig,
    OHCSeriesRequest,
    extract_ohc,
    OHCRegistry,
    extract_all_ohc,
    OHCHovmollerRequest,
    extract_all_ohc_hovmoller,
)
from util.flux_timeseries import (
    FluxTimeSeriesExtractorV2,
    FluxExtractOptionsV2,
    FluxSeriesRequestV2,
    FluxRegistryV2,
    extract_all_flux_v2,
)
from dataclasses import dataclass, asdict
from typing import Optional, Sequence, Tuple, Union, Literal, Iterable, List, Dict, Optional, Any


@dataclass
class FigureExtractSpec:
    # OCN (mpas-analysis) side
    ocn_ts_base: str
    ocn_climo_len: int
    ocn_need_base: Sequence[str]

    # FLUX side
    flux_ts_base: str
    flux_climo_len: int
    flux_need_vars: Sequence[str]
    flux_need_base: Sequence[str] = ()
    flux_annual_mean: bool = False
    flux_annual_weighted: bool = False
    flux_coef_flx: float = 86400.0 * 365.0

    # per-area controls (passed through to extract_flux_data)
    flux_per_area_vars: Optional[Sequence[str]] = None
    flux_per_area_map: Optional[Dict[str, bool]] = None
    global_area_m2: Optional[float] = None
    mesh_path: Optional[str] = None

    mpas_analysis: bool = False


@dataclass
class OHUFigureExtractSpec:
    # ---------------- REQUIRED (no defaults) ----------------
    # FLUX side
    flux_ts_base: str
    flux_climo_len: int
    flux_need_vars: Sequence[str]

    # OHU side
    ohc_climo_len: int

    # ---------------- OPTIONAL (defaults) -------------------
    # FLUX options
    flux_need_base: Sequence[str] = ()
    flux_annual_mean: bool = False
    flux_annual_weighted: bool = False
    flux_mpas_analysis: bool = False

    flux_per_area_vars: Optional[Sequence[str]] = None
    flux_per_area_map: Optional[Dict[str, bool]] = None
    global_area_m2: Optional[float] = None
    mesh_path: Optional[str] = None

    # OHU options (OHU via OHCTimeSeriesUtil)
    ohc_need_vars: Sequence[str] = ()
    ohc_file_key: str = "mpasTimeSeriesOcean"
    ohc_ts_base: Optional[str] = None
    ohc_fdepth: Optional[str] = None

    ohc_anomaly: bool = True
    ohc_baseline_months: int = 12
    ohc_to_ZJ: bool = True
    ohc_use_mask_in_volume: bool = True

    ohc_annual_mean: bool = False
    ohc_annual_weighted: bool = True
    ohc_annual_drop_incomplete: bool = True
    ohc_annual_min_days: int = 360
    ohc_reassign_monthly_time: bool = True
    ohc_length_policy: str = "auto"
    ohc_mpas_analysis: bool = False


class FigureDataCollector:
    """Extracts (or loads cached) plot-ready series for each experiment (timeseries_alternate/ocn_sfc_ts, tke_amoc)."""

    def __init__(
        self,
        *,
        run_dict,
        exp_tags: Sequence[str],
        exp_type: str,
        ts_window: Tuple[int, int],
        data_dir: str,
        year_token: str,
        time_unit: str,
        calendar: str,
        region_key: str = "Global",
        reg_dim: str = "nOceanRegions",
        reg_var=None,
        engine: str = "netcdf4",
        chunks=None,
        use_mfdataset: bool = True,
        verbose: bool = False,
    ):
        self.run_dict = run_dict
        self.exp_tags = list(exp_tags)
        self.exp_type = str(exp_type)
        self.ts_window = tuple(ts_window)
        self.data_dir = str(data_dir)
        self.year_token = str(year_token)
        self.time_unit = str(time_unit)
        self.calendar = str(calendar)
        self.region_key = str(region_key)
        self.reg_dim = str(reg_dim)
        self.reg_var = reg_var
        self.engine = str(engine)
        self.chunks = chunks
        self.use_mfdataset = bool(use_mfdataset)
        self.verbose = bool(verbose)

    @staticmethod
    def _nc_sanitize_attrs(attrs: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        NetCDF attribute sanitizer:
          - bool / np.bool_ -> int (0/1)
          - numpy scalars -> python scalars
          - lists/tuples/dicts -> JSON string
          - None -> skip
        """
        if not attrs:
            return {}
        out: Dict[str, Any] = {}
        for k, v in attrs.items():
            if v is None:
                continue
            if isinstance(v, np.generic):
                v = v.item()
            if isinstance(v, (bool, np.bool_)):
                v = int(v)
            if isinstance(v, (list, tuple, dict)):
                v = json.dumps(v)
            if isinstance(v, (bytes, bytearray)):
                v = v.decode("utf-8", errors="replace")
            out[str(k)] = v
        return out

    @staticmethod
    def extract_ocn_data(
        *,
        run_dict,
        exp_tags,
        exp_type,
        ts_window,
        climo_len,
        ts_base,
        data_dir,
        year_token,
        time_unit,
        calendar,
        need_base=("Mass_net",),
        region_key="Global",
        reg_dim="nOceanRegions",
        reg_var=None,
        engine="netcdf4",
        chunks=None,
        use_mfdataset=True,
    ):
        exp_subdirs = build_exp_subdirs_from_run_dict(
            run_dict=run_dict,
            exp_tags=exp_tags,
            exp_type=exp_type,
            ts_window=ts_window,
            climo_len=climo_len,
            ts_base=ts_base,
        )
        exps = list(exp_subdirs.keys())

        registry = OCNRegistry()
        cfg = registry.build_extract_config(
            year_token=year_token,
            region_key=region_key,
            region_dim=reg_dim,
            region_names_var=reg_var,
        )

        var_dict = registry.get_var_dict()
        var_dict_small = {k: var_dict[k] for k in need_base if k in var_dict}

        for var in need_base:
            if var not in var_dict_small:
                raise KeyError(f"{var} not found in registry var_dict (check OCNRegistry keys).")

        data = extract_all(
            root_dir=data_dir,
            exps=exps,
            exp_subdirs=exp_subdirs,
            var_dict=var_dict_small,
            cfg=cfg,
            time_origin=time_unit,
            calendar=calendar,
            engine=engine,
            chunks=chunks,
            use_mfdataset=use_mfdataset,
        )

        if len(exps) != len(exp_tags):
            raise ValueError(
                f"Cannot remap experiment keys: len(exps)={len(exps)} != len(exp_tags)={len(exp_tags)}"
            )

        data = {tag: data[exp] for exp, tag in zip(exps, exp_tags)}

        print("\n[Extracted series summary]")
        for exp in exp_tags:
            print(f"  {exp}:")
            for k in var_dict_small.keys():
                arr = data[exp][k]
                print(
                    f"    {k:12s}  n={arr.size:6d}  "
                    f"min={np.nanmin(arr):.4g}  max={np.nanmax(arr):.4g}"
                )

        return data

    @staticmethod
    def extract_flux_data(
        *,
        run_dict,
        exp_tags,
        exp_type,
        ts_window,
        climo_len,
        ts_base,
        data_dir,
        year_token,
        time_unit,
        calendar,
        need_vars=("Fmass",),
        need_base=(),
        region_key="Global",
        reg_dim="nOceanRegions",
        reg_var=None,
        mpas_analysis=False,
        annual_mean=False,
        annual_weighted=False,
        coef_flx=86400.0 * 365.0,
        reassign_monthly_time=True,
        length_policy="auto",
        engine="netcdf4",
        chunks=None,
        verbose=False,
        per_area_vars=None,
        per_area_map=None,
        global_area_m2=None,
        mesh_path=None,
    ):
        exp_subdirs = build_exp_subdirs_from_run_dict(
            run_dict=run_dict,
            exp_tags=exp_tags,
            exp_type=exp_type,
            ts_window=ts_window,
            climo_len=climo_len,
            ts_base=ts_base,
            mpas_analysis=mpas_analysis,
        )
        exps = list(exp_subdirs.keys())

        registry = FluxRegistry()

        # Keep your current behavior here (derived registrations live here)
        # Net fluxes as sums of their components: MPAS netMassFlux / netEnergyFlux are 0 in the first month
        # after every restart (accumulator reset). Qnet reproduces netEnergyFlux exactly in all other months.
        registry.add_derived_terms("Fmass", terms=[(k, +coef_flx) for k in ("FW_rain", "FW_snow", "FW_evap", "FW_river", "FW_irun", "FW_si", "FW_frazil", "FW_li", "FW_berg")])
        registry.add_derived_terms("Qnet", terms=[(k, 1.0) for k in ("Qsw", "Qlw_dn", "Qlw_up", "Qlh", "Qsh", "Qsi", "Qfrazil", "Qsnow", "Qirunoff", "Qberg", "Qli")]
                                   + [(k, 1026.0 * 3996.0) for k in ("Qrain", "Qriver", "Qevap", "QsiT", "QbergT")])
        registry.add_derived_expr("Qatm", expr="Qsw + Qlw_dn - Qlw_up - Qlh - Qsh")
        registry.add_derived_terms("TKE", terms=[("KE", 1.0)])

        exemplar = (data_dir, exps[0], ts_base)

        cfg = registry.build_extract_config(
            year_token=year_token,
            region_key=region_key,
            region_dim=reg_dim,
            region_names_var=reg_var,
            ts_window=ts_window,
            annual_mean=annual_mean,
            annual_weighted=annual_weighted,
            calendar=calendar,
            reassign_monthly_time=reassign_monthly_time,
            length_policy=length_policy,
        )

        var_dict = registry.get_var_dict()
        derived = registry.get_derived()

        want_vars = list(need_vars) if need_vars is not None else list(need_base)
        want_base = [v for v in want_vars if v in var_dict]
        want_derived = [v for v in want_vars if v in derived]
        unknown = [v for v in want_vars if (v not in var_dict and v not in derived)]
        if unknown:
            raise KeyError(
                f"Unknown requested flux variables: {unknown}. Known derived: {list(derived.keys())}"
            )

        inferred_base = set()
        for dk in want_derived:
            rec = derived[dk]
            if "terms" in rec:
                for subkey, _coef in rec["terms"]:
                    inferred_base.add(subkey)

        base_to_extract = sorted(set(need_base) | set(want_base) | inferred_base)

        missing_base = [v for v in base_to_extract if v not in var_dict]
        if missing_base:
            raise KeyError(f"Missing in FluxRegistry var_dict: {missing_base}")

        var_dict_small = {k: var_dict[k] for k in base_to_extract}
        derived_small = {k: derived[k] for k in want_derived}

        # per-area selection
        per_area_set = set(per_area_vars or ())
        if per_area_map:
            for k, v in per_area_map.items():
                if v:
                    per_area_set.add(k)
                else:
                    per_area_set.discard(k)
        per_area_set = per_area_set.intersection(set(want_vars))
        per_area_vars_use = list(per_area_set)

        if per_area_vars_use:
            if global_area_m2 is None and mesh_path is None and exemplar is None:
                raise ValueError(
                    "per-area requested but no area source provided. "
                    "Provide one of: global_area_m2, mesh_path, or exemplar."
                )

        data = extract_all_flux(
            base_dir=data_dir,
            exps=exps,
            exp_subdirs=exp_subdirs,
            var_dict=var_dict_small,
            derived=derived_small,
            cfg=cfg,
            engine=engine,
            chunks=chunks,
            verbose=verbose,
            per_area_vars=per_area_vars_use,
            global_area_m2=global_area_m2,
            mesh_path=mesh_path,
            exemplar=exemplar,
        )

        if len(exps) != len(exp_tags):
            raise ValueError(
                f"Cannot remap experiment keys: len(exps)={len(exps)} != len(exp_tags)={len(exp_tags)}"
            )

        data = {tag: data[exp] for exp, tag in zip(exps, exp_tags)}

        if need_vars is not None:
            data = {exp: {k: data[exp][k] for k in want_vars} for exp in exp_tags}

        print("\n[Extracted FLUX series summary]")
        print(f"  per-area applied to: {per_area_vars_use if per_area_vars_use else 'None'}")
        for exp in exp_tags:
            print(f"  {exp}:")
            for k in want_vars:
                arr = data[exp][k]
                print(
                    f"    {k:12s}  n={arr.size:6d}  "
                    f"min={np.nanmin(arr):.4g}  max={np.nanmax(arr):.4g}"
                )

        return data

    @staticmethod
    def build_plot_data(ocn_data, flux_data, plot_vars, exp_tags):
        plot_data = {}
        for exp in exp_tags:
            plot_data[exp] = {}
            for v in plot_vars:
                if v in ocn_data.get(exp, {}):
                    plot_data[exp][v] = ocn_data[exp][v]
                elif v in flux_data.get(exp, {}):
                    plot_data[exp][v] = flux_data[exp][v]
                else:
                    raise KeyError(f"{v} not found for experiment '{exp}'")
        return plot_data

    def collect(
        self,
        *,
        spec: FigureExtractSpec,
        plot_vars: Sequence[str],
    ) -> Dict[str, Dict[str, np.ndarray]]:
        ocn_data = self.extract_ocn_data(
            run_dict=self.run_dict,
            exp_tags=self.exp_tags,
            exp_type=self.exp_type,
            ts_window=self.ts_window,
            climo_len=spec.ocn_climo_len,
            ts_base=spec.ocn_ts_base,
            data_dir=self.data_dir,
            year_token=self.year_token,
            time_unit=self.time_unit,
            calendar=self.calendar,
            need_base=tuple(spec.ocn_need_base),
            region_key=self.region_key,
            reg_dim=self.reg_dim,
            reg_var=self.reg_var,
            engine=self.engine,
            chunks=self.chunks,
            use_mfdataset=self.use_mfdataset,
        )

        flux_data = self.extract_flux_data(
            run_dict=self.run_dict,
            exp_tags=self.exp_tags,
            exp_type=self.exp_type,
            ts_window=self.ts_window,
            climo_len=spec.flux_climo_len,
            ts_base=spec.flux_ts_base,
            data_dir=self.data_dir,
            year_token=self.year_token,
            time_unit=self.time_unit,
            calendar=self.calendar,
            need_vars=tuple(spec.flux_need_vars),
            need_base=tuple(spec.flux_need_base),
            mpas_analysis=spec.mpas_analysis,
            annual_mean=spec.flux_annual_mean,
            annual_weighted=spec.flux_annual_weighted,
            coef_flx=spec.flux_coef_flx,
            engine=self.engine,
            chunks=self.chunks,
            verbose=self.verbose,
            per_area_vars=spec.flux_per_area_vars,
            per_area_map=spec.flux_per_area_map,
            global_area_m2=spec.global_area_m2,
            mesh_path=spec.mesh_path,
        )

        plot_data = self.build_plot_data(
            ocn_data=ocn_data,
            flux_data=flux_data,
            plot_vars=plot_vars,
            exp_tags=self.exp_tags,
        )

        if self.verbose:
            print("\n[FigureDataCollector] merged keys:")
            for exp in self.exp_tags:
                print(f"  {exp}: {list(plot_data[exp].keys())}")

        return plot_data

    def save_per_exp(
        self,
        *,
        out_nc_base: str,
        plot_data: Dict[str, Dict[str, np.ndarray]],
        plot_info: Dict[str, Dict[str, Any]],
        freq: Literal["monthly", "annual"],
        year_min_for_coords: int,
        year_max_for_meta: int,
        apply_scale: bool = True,
        attrs: Optional[Dict[str, Any]] = None,
        compress: bool = True,
    ) -> None:
        base = os.path.splitext(out_nc_base)[0]
        attrs_in = {} if attrs is None else dict(attrs)
        ts_var_keys = list(plot_info.keys())

        for exp in self.exp_tags:
            if exp not in plot_data:
                raise KeyError(f"[save_ts] missing exp in plot_data: {exp}")

            v0 = ts_var_keys[0]
            if v0 not in plot_data[exp]:
                raise KeyError(f"[save_ts] missing var '{v0}' in plot_data[{exp}]")
            arr0 = np.asarray(plot_data[exp][v0], float).reshape(-1)
            n = arr0.size

            if freq == "annual":
                dim = "year"
                year = np.arange(int(year_min_for_coords), int(year_min_for_coords) + n, dtype=int)
                coords = {"year": year}
            elif freq == "monthly":
                dim = "month"
                month = np.arange(n, dtype=int)
                time_year = float(year_min_for_coords) + month / 12.0
                coords = {"month": month, "time_year": (("month",), time_year)}
            else:
                raise ValueError(f"freq must be 'monthly' or 'annual', got {freq}")

            data_vars: Dict[str, xr.DataArray] = {}
            for var in ts_var_keys:
                if var not in plot_data[exp]:
                    raise KeyError(f"[save_ts] missing var '{var}' in plot_data[{exp}]")

                info = plot_info.get(var, {}) or {}
                units = info.get("unit", "")
                long_name = info.get("label", var)
                scale = float(info.get("scale", 1.0)) if apply_scale else 1.0

                arr = np.asarray(plot_data[exp][var], float).reshape(-1)
                if arr.size != n:
                    raise ValueError(f"[save_ts] length mismatch {exp}:{var} ({arr.size} vs {n})")

                data_vars[var] = xr.DataArray(
                    arr * scale,
                    dims=(dim,),
                    coords=coords,
                    attrs={"units": units, "long_name": long_name, "applied_scale": scale},
                )

            ds = xr.Dataset(data_vars)
            base_attrs = {
                "experiment": exp,
                "year_min": int(year_min_for_coords),
                "year_max_meta": int(year_max_for_meta),
                "freq": str(freq),
                "source_kind": "time_series",
            }
            ds.attrs.update(self._nc_sanitize_attrs(base_attrs))
            ds.attrs.update(self._nc_sanitize_attrs(attrs_in))

            fn = f"{base}_ts_{exp}.nc"
            encoding = None
            if compress:
                encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}

            ds.to_netcdf(fn, encoding=encoding)
            if self.verbose:
                print(f"[saved TS] {fn}")

    def load_per_exp(
        self,
        *,
        out_nc_base: str,
        ts_var_keys: Sequence[str],
        prefer_float64: bool = True,
        exp_order: Optional[Sequence[str]] = None,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        base = os.path.splitext(out_nc_base)[0]
        use_exps = list(exp_order) if exp_order is not None else self.exp_tags

        plot_data: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in use_exps:
            fn = f"{base}_ts_{exp}.nc"
            if not os.path.exists(fn):
                raise FileNotFoundError(f"[load_ts] missing file: {fn}")

            with xr.open_dataset(fn) as ds:
                plot_data[exp] = {}
                for var in ts_var_keys:
                    if var not in ds.data_vars:
                        raise KeyError(f"[load_ts] {fn} missing variable '{var}'")
                    vals = ds[var].values
                    arr = np.asarray(vals, dtype=float if prefer_float64 else vals.dtype).reshape(-1)
                    plot_data[exp][var] = arr

        if self.verbose:
            print(f"\n[FigureDataCollector] loaded from: {base}_ts_<exp>.nc")
            for exp in use_exps:
                print(f"  {exp}: {list(plot_data[exp].keys())}")

        return plot_data


class OHUFigureDataCollector(FigureDataCollector):
    """Collects OHC/OHU and surface-flux series (timeseries_alternate/ohu_ts); uses flux extractor v2."""

    @staticmethod
    def extract_ohc_ohu_data(
        *,
        run_dict,
        exp_tags,
        exp_type,
        ts_window,
        climo_len,
        ts_base,
        data_dir,
        year_token,
        calendar,
        region_key="Global",
        reg_dim="nOceanRegions",
        need_vars=(),
        file_key="mpasTimeSeriesOcean",
        mpas_analysis=False,
        fdepth="",
        anomaly=True,
        baseline_months=12,
        to_ZJ=True,
        use_mask_in_volume=True,
        annual_mean=False,
        annual_weighted=True,
        annual_drop_incomplete=True,
        annual_min_days=360,
        reassign_monthly_time=True,
        length_policy="auto",
        engine="netcdf4",
        chunks=None,
        verbose=False,
    ):
        if not need_vars:
            return {tag: {} for tag in exp_tags}

        # subdirs: same helper as you already use; ts_base is the OHCTimeSeriesOcean directory base
        exp_subdirs = build_exp_subdirs_from_run_dict(
            run_dict=run_dict,
            exp_tags=exp_tags,
            exp_type=exp_type,
            ts_window=ts_window,
            climo_len=climo_len,          # not used by this path; safe placeholder
            ts_base=ts_base,
            mpas_analysis=mpas_analysis,
        )
        exps = list(exp_subdirs.keys())

        registry = OHCRegistry()
        var_dict = registry.get_var_dict()

        # only pull what you asked for
        var_dict_small = {k: var_dict[k] for k in need_vars if k in var_dict}
        missing = [k for k in need_vars if k not in var_dict_small]
        if missing:
            raise KeyError(f"Unknown OHU/OHC keys in OHCRegistry: {missing}")

        region_index = registry.get_region_index(region_key)

        cfg = registry.build_extract_config(
            ts_window=ts_window,
            annual_mean=annual_mean,
            annual_weighted=annual_weighted,
            annual_drop_incomplete=annual_drop_incomplete,
            annual_min_days=annual_min_days,
            calendar=calendar,
            reassign_monthly_time=reassign_monthly_time,
            length_policy=length_policy,
        )
        
        if not fdepth:
            raise ValueError("ohc_fdepth is required for OHU/OHC extraction (path to depth.nc).")

        data = extract_all_ohc(
            root_dir=data_dir,
            exps=exps,
            exp_subdirs=exp_subdirs,
            var_dict=var_dict_small,
            cfg=cfg,
            region_index=region_index,
            region_dim=reg_dim,
            anomaly=anomaly,
            baseline_months=baseline_months,
            to_ZJ=to_ZJ,
            use_mask_in_volume=use_mask_in_volume,
            file_key=file_key,
            fdepth=fdepth,
            engine=engine,
            chunks=chunks,
            verbose=verbose,
        )

        # remap internal exp keys to exp_tags (same pattern as your other extractors)
        if len(exps) != len(exp_tags):
            raise ValueError(f"Cannot remap experiment keys: len(exps)={len(exps)} != len(exp_tags)={len(exp_tags)}")
        data = {tag: data[exp] for exp, tag in zip(exps, exp_tags)}

        if verbose:
            print("\n[Extracted OHC/OHU series summary]")
            for exp in exp_tags:
                print(f"  {exp}:")
                for k in var_dict_small.keys():
                    arr = data[exp][k]
                    print(f"    {k:12s}  n={arr.size:6d}  min={np.nanmin(arr):.4g}  max={np.nanmax(arr):.4g}")

        return data

    @staticmethod
    def extract_flux_data(
        *,
        run_dict,
        exp_tags,
        exp_type,
        ts_window,
        climo_len,
        ts_base,
        data_dir,
        year_token,
        time_unit,
        calendar,
        need_vars=("Fmass",),
        need_base=(),
        region_key="Global",
        reg_dim="nOceanRegions",
        reg_var=None,
        mpas_analysis=False,
        annual_mean=False,
        annual_weighted=False,
        reassign_monthly_time=True,
        length_policy="auto",
        engine="netcdf4",
        chunks=None,
        verbose=False,
        per_area_vars=None,
        per_area_map=None,
        global_area_m2=None,
        mesh_path=None,
    ):
        exp_subdirs = build_exp_subdirs_from_run_dict(
            run_dict=run_dict,
            exp_tags=exp_tags,
            exp_type=exp_type,
            ts_window=ts_window,
            climo_len=climo_len,
            ts_base=ts_base,
            mpas_analysis=mpas_analysis,
        )
        exps = list(exp_subdirs.keys())

        registry = FluxRegistryV2()

        # Keep your current behavior here (derived registrations live here)
        rho_sw = 1026.0  # density of salt water (kg/m^3)
        cp_sw = 3.996e3  # specific heat salt water
        mass_coef_flx = 86400.0 * 365.0 
        heat_coef_flx = rho_sw * cp_sw 
        # True derived only
        # Net fluxes as sums of their components: MPAS netMassFlux / netEnergyFlux are 0 in the first month
        # after every restart (accumulator reset). Qnet reproduces netEnergyFlux exactly in all other months.
        registry.add_derived_terms("Fmass", terms=[(k, +mass_coef_flx) for k in ("FW_rain", "FW_snow", "FW_evap", "FW_river", "FW_irun", "FW_si", "FW_frazil", "FW_li", "FW_berg")])
        registry.add_derived_terms("Qnet", terms=[(k, 1.0) for k in ("Qsw", "Qlw_dn", "Qlw_up", "Qlh", "Qsh", "Qsi", "Qfrazil", "Qsnow", "Qirunoff", "Qberg", "Qli")]
                                   + [(k, heat_coef_flx) for k in ("QrainT", "QriverT", "QevapT", "QsiT", "QbergT")])
        registry.add_derived_expr("Qatm", expr="Qsw + Qlw_dn - Qlw_up - Qlh - Qsh")
        
        #registry.add_derived_terms(
        #    "Qsource",
        #    terms=[
        #        ("Qsw",      +1.0),
        #        ("Qlw_dn",   +1.0),
        #        ("Qfrazil",  +1.0),
        #        ("QriverT",  +heat_coef_flx),
        #        ("QrainT",   +heat_coef_flx),
        #    ],
        #)
        
        #registry.add_derived_terms(
        #    "Qsink",
        #    terms=[
        #        ("Qlw_up",   +1.0),
        #        ("Qsh",      +1.0),
        #        ("Qlh",      +1.0),
        #        ("Qsi",      +1.0),
        #        ("Qsnow",    +1.0),
        #        ("Qirunoff", +1.0),
        #        ("Qberg",    +1.0),
        #        ("Qli",      +1.0),
        #        ("QbergT",   +heat_coef_flx),
        #        ("QevapT",   +heat_coef_flx),
        #        ("QsiT",     +heat_coef_flx),
        #    ],
        #)
        
        registry.add_derived_terms(
            "Qsource",
            terms=[
                ("Qsw",      +1.0),
                ("Qlw_dn",   +1.0)
            ],
        )
        
        registry.add_derived_terms(
            "Qsink",
            terms=[
                ("Qlw_up",   +1.0),
                ("Qsh",      +1.0),
                ("Qlh",      +1.0)
            ],
        )
        
        registry.add_derived_terms(
            "Qrad",
            terms=[("Qsw", +1.0), ("Qlw_dn", +1.0), ("Qlw_up", +1.0)],
        )
        
        registry.add_derived_terms(
            "Qsrf",
            terms=[("Qsw", +1.0), ("Qlw_dn", +1.0), ("Qlw_up", +1.0), ("Qsh", +1.0), ("Qlh", +1.0)],
        )

        exemplar = (data_dir, exps[0], ts_base)

        cfg = registry.build_extract_config(
            year_token=year_token,
            region_key=region_key,
            region_dim=reg_dim,
            region_names_var=reg_var,
            ts_window=ts_window,
            annual_mean=annual_mean,
            annual_weighted=annual_weighted,
            calendar=calendar,
            reassign_monthly_time=reassign_monthly_time,
            length_policy=length_policy,
        )

        var_dict = registry.get_var_dict()
        derived = registry.get_derived()

        want_vars = list(need_vars) if need_vars is not None else list(need_base)
        want_base = [v for v in want_vars if v in var_dict]
        want_derived = [v for v in want_vars if v in derived]
        unknown = [v for v in want_vars if (v not in var_dict and v not in derived)]
        if unknown:
            raise KeyError(
                f"Unknown requested flux variables: {unknown}. Known derived: {list(derived.keys())}"
            )

        inferred_base = set()
        for dk in want_derived:
            rec = derived[dk]
            if "terms" in rec:
                for subkey, _coef in rec["terms"]:
                    inferred_base.add(subkey)

        base_to_extract = sorted(set(need_base) | set(want_base) | inferred_base)

        missing_base = [v for v in base_to_extract if v not in var_dict]
        if missing_base:
            raise KeyError(f"Missing in FluxRegistry var_dict: {missing_base}")

        var_dict_small = {k: var_dict[k] for k in base_to_extract}
        derived_small = {k: derived[k] for k in want_derived}

        # per-area selection
        per_area_set = set(per_area_vars or ())
        if per_area_map:
            for k, v in per_area_map.items():
                if v:
                    per_area_set.add(k)
                else:
                    per_area_set.discard(k)
        per_area_set = per_area_set.intersection(set(want_vars))
        per_area_vars_use = list(per_area_set)

        if per_area_vars_use:
            if global_area_m2 is None and mesh_path is None and exemplar is None:
                raise ValueError(
                    "per-area requested but no area source provided. "
                    "Provide one of: global_area_m2, mesh_path, or exemplar."
                )

        data = extract_all_flux_v2(
            base_dir=data_dir,
            exps=exps,
            exp_subdirs=exp_subdirs,
            var_dict=var_dict_small,
            derived=derived_small,
            cfg=cfg,
            engine=engine,
            chunks=chunks,
            verbose=verbose,
            per_area_vars=per_area_vars_use,
            global_area_m2=global_area_m2,
            mesh_path=mesh_path,
            exemplar=exemplar,
        )

        if len(exps) != len(exp_tags):
            raise ValueError(
                f"Cannot remap experiment keys: len(exps)={len(exps)} != len(exp_tags)={len(exp_tags)}"
            )

        data = {tag: data[exp] for exp, tag in zip(exps, exp_tags)}

        if need_vars is not None:
            data = {exp: {k: data[exp][k] for k in want_vars} for exp in exp_tags}

        print("\n[Extracted FLUX series summary]")
        print(f"  per-area applied to: {per_area_vars_use if per_area_vars_use else 'None'}")
        for exp in exp_tags:
            print(f"  {exp}:")
            for k in want_vars:
                arr = data[exp][k]
                print(
                    f"    {k:12s}  n={arr.size:6d}  "
                    f"min={np.nanmin(arr):.4g}  max={np.nanmax(arr):.4g}"
                )

        return data

    @staticmethod
    def build_plot_data(ohc_data, flux_data, plot_vars, exp_tags):
        plot_data = {}
        for exp in exp_tags:
            plot_data[exp] = {}
            for v in plot_vars:
                if v in ohc_data.get(exp, {}):
                    plot_data[exp][v] = ohc_data[exp][v]
                elif v in flux_data.get(exp, {}):
                    plot_data[exp][v] = flux_data[exp][v]
                else:
                    raise KeyError(f"{v} not found for experiment '{exp}'")
        return plot_data

    def collect(
        self,
        *,
        spec: OHUFigureExtractSpec,
        plot_vars: Sequence[str],
    ) -> Dict[str, Dict[str, np.ndarray]]:
        ohc_data = self.extract_ohc_ohu_data(
            run_dict=self.run_dict,
            exp_tags=self.exp_tags,
            exp_type=self.exp_type,
            ts_window=self.ts_window,
            climo_len=spec.flux_climo_len,
            ts_base=spec.ohc_ts_base,
            data_dir=self.data_dir,
            year_token=self.year_token,
            calendar=self.calendar,
            region_key=self.region_key,
            reg_dim=self.reg_dim,
            need_vars=tuple(spec.ohc_need_vars),
            mpas_analysis=spec.ohc_mpas_analysis,
            file_key=spec.ohc_file_key,
            fdepth=spec.ohc_fdepth,
            anomaly=spec.ohc_anomaly,
            baseline_months=spec.ohc_baseline_months,
            to_ZJ=spec.ohc_to_ZJ,
            use_mask_in_volume=spec.ohc_use_mask_in_volume,
            annual_mean=spec.ohc_annual_mean,
            annual_weighted=spec.ohc_annual_weighted,
            annual_drop_incomplete=spec.ohc_annual_drop_incomplete,
            annual_min_days=spec.ohc_annual_min_days,
            reassign_monthly_time=spec.ohc_reassign_monthly_time,
            length_policy=spec.ohc_length_policy,
            engine=self.engine,
            chunks=self.chunks,
            verbose=self.verbose,
        )

        flux_data = self.extract_flux_data(
            run_dict=self.run_dict,
            exp_tags=self.exp_tags,
            exp_type=self.exp_type,
            ts_window=self.ts_window,
            climo_len=spec.flux_climo_len,
            ts_base=spec.flux_ts_base,
            data_dir=self.data_dir,
            year_token=self.year_token,
            time_unit=self.time_unit,
            calendar=self.calendar,
            need_vars=tuple(spec.flux_need_vars),
            need_base=tuple(spec.flux_need_base),
            mpas_analysis=spec.flux_mpas_analysis,
            annual_mean=spec.flux_annual_mean,
            annual_weighted=spec.flux_annual_weighted,
            engine=self.engine,
            chunks=self.chunks,
            verbose=self.verbose,
            per_area_vars=spec.flux_per_area_vars,
            per_area_map=spec.flux_per_area_map,
            global_area_m2=spec.global_area_m2,
            mesh_path=spec.mesh_path,
        )

        plot_data = self.build_plot_data(
            ohc_data=ohc_data,
            flux_data=flux_data,
            plot_vars=plot_vars,
            exp_tags=self.exp_tags,
        )

        if self.verbose:
            print("\n[FigureDataCollector] merged keys:")
            for exp in self.exp_tags:
                print(f"  {exp}: {list(plot_data[exp].keys())}")

        return plot_data

    def load_per_exp(
        self,
        *,
        out_nc_base: str,
        ts_var_keys: Sequence[str],
        prefer_float64: bool = True,
        exp_order: Optional[Sequence[str]] = None,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        base = os.path.splitext(out_nc_base)[0]
        use_exps = list(exp_order) if exp_order is not None else self.exp_tags

        plot_data: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in use_exps:
            fn = f"{base}_ts_{exp}.nc"
            if not os.path.exists(fn):
                raise FileNotFoundError(f"[load_ts] missing file: {fn}")

            with xr.open_dataset(fn) as ds:
                plot_data[exp] = {}
                for var in ts_var_keys:
                    if var not in ds.data_vars:
                        raise KeyError(f"[load_ts] {fn} missing variable '{var}'")
                    vals = ds[var].values
                    arr = np.asarray(vals, dtype=float if prefer_float64 else vals.dtype).reshape(-1)
                    plot_data[exp][var] = arr

        if self.verbose:
            print(f"\n[FigureDataCollector] loaded from: {base}_ts_<exp>.nc")
            for exp in use_exps:
                print(f"  {exp}: {list(plot_data[exp].keys())}")

        return plot_data


class OCNDiagnosticsPlotter:
    """Multi-experiment time-series panels with trends, std bands and bootstrap tests (timeseries_alternate/ocn_sfc_ts, tke_amoc)."""

    def __init__(
        self,
        *,
        year_min: int,
        year_max: int,
        annual_mean: bool = True,
        compute_anomaly: bool = False,
        baseline_months: int = 12,
        baseline_years: int = 1,

        running_mean: bool = False,
        running_mean_months: int = 0,
        running_mean_center: str = "center",   # center|trailing|leading
        running_mean_min_frac: float = 0.5,

        use_filter: bool = False,
        filt_cutoff_months: float = 24.0,
        filt_order: int = 4,
        filt_method: str = "pad",          # pad|gust
        filt_padlen: Optional[int] = None,
        filt_padtype: str = "odd",
        edge_fix_points: int = 0,

        trend_years: int = 30,
        trend_unit: str = "per_year",      # per_year|per_decade
        add_trend_ci: bool = True,
        ci_method: str = "hac",            # hac|glsar|neff
        hac_lag_years: int = 1,
        glsar_max_iter: int = 10,
    ):
        self.year_min = int(year_min)
        self.year_max = int(year_max)

        self.annual_mean = bool(annual_mean)
        self.compute_anomaly = bool(compute_anomaly)
        self.baseline_months = int(baseline_months)
        self.baseline_years = int(baseline_years)

        self.running_mean = bool(running_mean)
        self.running_mean_months = int(running_mean_months or 0)
        self.running_mean_center = str(running_mean_center).lower()
        self.running_mean_min_frac = float(running_mean_min_frac)

        self.use_filter = bool(use_filter)
        self.filt_cutoff_months = float(filt_cutoff_months)
        self.filt_order = int(filt_order)
        self.filt_method = str(filt_method)
        self.filt_padlen = None if filt_padlen is None else int(filt_padlen)
        self.filt_padtype = str(filt_padtype)
        self.edge_fix_points = int(edge_fix_points)

        self.trend_years = int(trend_years)
        self.trend_unit = str(trend_unit)
        self.add_trend_ci = bool(add_trend_ci)
        self.ci_method = str(ci_method)
        self.hac_lag_years = int(hac_lag_years)
        self.glsar_max_iter = int(glsar_max_iter)

        # pre-design butter if needed
        self._butter_b = self._butter_a = None
        if self.use_filter:
            fc = 1.0 / self.filt_cutoff_months
            wn = max(min(fc / 0.5, 0.999999), 1e-6)
            self._butter_b, self._butter_a = signal.butter(self.filt_order, wn, btype="low")

    @staticmethod
    def to_annual_mean(x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, float).reshape(-1)
        n_years = len(x) // 12
        if n_years <= 0:
            return np.asarray([], float)
        return x[: n_years * 12].reshape(n_years, 12).mean(axis=1)

    @staticmethod
    def _nan_running_mean(x, win, center="center", min_valid_frac=0.5):
        x = np.asarray(x, float).reshape(-1)
        n = len(x)
        if win <= 1 or n == 0:
            return x.copy()

        valid = np.isfinite(x).astype(float)
        x0 = np.where(np.isfinite(x), x, 0.0)
        k = np.ones(int(win), float)

        if center == "trailing":
            s = np.convolve(x0, k, mode="full")[win - 1 : win - 1 + n]
            c = np.convolve(valid, k, mode="full")[win - 1 : win - 1 + n]
        elif center == "leading":
            s = np.convolve(x0, k, mode="full")[:n]
            c = np.convolve(valid, k, mode="full")[:n]
        else:
            s = np.convolve(x0, k, mode="same")
            c = np.convolve(valid, k, mode="same")

        with np.errstate(invalid="ignore", divide="ignore"):
            y = s / c
        min_count = max(1, int(np.ceil(min_valid_frac * win)))
        y[c < min_count] = np.nan
        return y

    @staticmethod
    def _rolling_nanstd(x, win):
        x = np.asarray(x, float).reshape(-1)
        n = len(x)
        if n == 0:
            return x.copy()
        win = max(1, int(win))
        if win % 2 == 0:
            win += 1
        half = win // 2
        out = np.empty(n, dtype=float)
        for i in range(n):
            a = max(0, i - half)
            b = min(n, i + half + 1)
            out[i] = np.nanstd(x[a:b], ddof=1) if (b - a) > 1 else np.nanstd(x[a:b])
        return out

    def _to_anomaly(self, x: np.ndarray, *, samples_per_year: int) -> np.ndarray:
        x = np.asarray(x, float).reshape(-1)
        n0 = self.baseline_years if samples_per_year == 1 else self.baseline_months
        if not (0 < n0 <= len(x)):
            return x - np.nanmean(x)
        clim = float(np.nanmean(x[:n0]))
        return x - clim

    def _lowpass(self, x: np.ndarray) -> np.ndarray:
        if not self.use_filter:
            return np.asarray(x, float)

        x = np.asarray(x, float).reshape(-1)
        n = len(x)
        if n < 3:
            return x.copy()

        if self._butter_b is None or self._butter_a is None:
            fc = 1.0 / self.filt_cutoff_months
            wn = max(min(fc / 0.5, 0.999999), 1e-6)
            self._butter_b, self._butter_a = signal.butter(self.filt_order, wn, btype="low")

        b, a = self._butter_b, self._butter_a

        if self.filt_method == "gust":
            y = signal.filtfilt(b, a, x, method="gust")
        else:
            default_pad = 3 * (max(len(a), len(b)) - 1)
            cutoff_pad = int(max(1, round(2 * self.filt_cutoff_months)))
            padlen = self.filt_padlen if self.filt_padlen is not None else max(default_pad, cutoff_pad)
            padlen = min(padlen, max(1, n - 1))
            if padlen >= n:
                y = signal.filtfilt(b, a, x, method="gust")
            else:
                y = signal.filtfilt(b, a, x, method="pad", padtype=self.filt_padtype, padlen=padlen)

        k = self.edge_fix_points
        if k > 0 and n >= 2 * k + 1:
            y[:k] = y[k]
            y[-k:] = y[-k - 1]
        return y

    def _last_n_year_stats(self, raw: np.ndarray, filt: np.ndarray) -> Dict[str, Any]:
        """
        Compute last-N-year summary + trend diagnostics.
        (This is your implementation; left intact except formatting.)
        """
        raw = np.asarray(raw, float).reshape(-1)
        filt = np.asarray(filt, float).reshape(-1)
        n = int(min(raw.size, filt.size))

        out: Dict[str, Any] = {
            "n_total": n,
            "window_years": int(self.trend_years),
            "trend_unit": str(self.trend_unit),
            "ci_method": str(self.ci_method),
            "add_trend_ci": bool(self.add_trend_ci),
        }

        if n == 0:
            out.update(
                {
                    "n_used": 0,
                    "start_idx": None,
                    "end_idx": None,
                    "raw_mean": np.nan,
                    "raw_std": np.nan,
                    "filt_mean": np.nan,
                    "filt_std": np.nan,
                    "trend_slope_raw": np.nan,
                    "trend_slope_filt": np.nan,
                    "trend_ci95_raw": (np.nan, np.nan),
                    "trend_ci95_filt": (np.nan, np.nan),
                }
            )
            return out

        spyr = 1 if self.annual_mean else 12

        win = max(2, int(round(self.trend_years * spyr)))
        i1 = n
        i0 = max(0, n - win)
        out["start_idx"] = int(i0)
        out["end_idx"] = int(i1 - 1)

        raw_w = raw[i0:i1]
        filt_w = filt[i0:i1]

        def _basic_stats(x: np.ndarray) -> Tuple[float, float, int]:
            m = np.isfinite(x)
            if m.sum() == 0:
                return (np.nan, np.nan, 0)
            xx = x[m]
            std = float(np.nanstd(xx, ddof=1)) if xx.size > 1 else float(np.nanstd(xx))
            return (float(np.nanmean(xx)), std, int(xx.size))

        raw_mean, raw_std, n_raw = _basic_stats(raw_w)
        filt_mean, filt_std, n_filt = _basic_stats(filt_w)

        out.update(
            {
                "raw_mean": raw_mean,
                "raw_std": raw_std,
                "filt_mean": filt_mean,
                "filt_std": filt_std,
                "n_used_raw": n_raw,
                "n_used_filt": n_filt,
            }
        )

        def _trend_and_ci(y: np.ndarray) -> Tuple[float, Tuple[float, float], Dict[str, Any]]:
            diag: Dict[str, Any] = {}
            m = np.isfinite(y)
            if m.sum() < 3:
                return (np.nan, (np.nan, np.nan), {"reason": "too_few_points"})

            yy = y[m]
            t = np.arange(yy.size, dtype=float)

            tt = t - t.mean()
            slope = float(np.dot(tt, yy - yy.mean()) / np.dot(tt, tt))
            diag["n"] = int(yy.size)
            diag["slope_per_sample"] = slope

            if not self.add_trend_ci:
                return (slope, (np.nan, np.nan), diag)

            method = str(self.ci_method).lower()

            try:
                import statsmodels.api as sm
                X = sm.add_constant(t)

                if method == "hac":
                    maxlags = max(1, int(round(self.hac_lag_years * spyr)))
                    res = sm.OLS(yy, X).fit(cov_type="HAC", cov_kwds={"maxlags": maxlags})
                    se = float(res.bse[1])
                    diag.update({"method_used": "hac", "hac_maxlags": maxlags})

                elif method == "glsar":
                    model = sm.GLSAR(yy, X, rho=1)
                    res = model.iterative_fit(maxiter=int(self.glsar_max_iter))
                    se = float(res.bse[1])
                    diag.update({"method_used": "glsar", "glsar_maxiter": int(self.glsar_max_iter)})

                else:
                    raise ImportError("fallback_to_neff")

                lo = slope - 1.96 * se
                hi = slope + 1.96 * se
                diag["se_slope"] = se
                return (slope, (float(lo), float(hi)), diag)

            except Exception:
                diag["method_used"] = "neff"

                yhat = (yy.mean() + slope * (t - t.mean()))
                r = yy - yhat

                if r.size < 3:
                    return (slope, (np.nan, np.nan), {**diag, "reason": "too_few_for_neff"})

                r0 = r[:-1]
                r1 = r[1:]
                denom = np.sqrt(np.dot(r0, r0) * np.dot(r1, r1))
                rho1 = float(np.dot(r0, r1) / denom) if denom > 0 else 0.0
                rho1 = float(np.clip(rho1, -0.99, 0.99))
                diag["rho1"] = rho1

                neff = yy.size * (1.0 - rho1) / (1.0 + rho1)
                neff = float(np.clip(neff, 3.0, yy.size))
                diag["neff"] = neff

                sigma2 = float(np.nanvar(r, ddof=1)) if r.size > 1 else float(np.nanvar(r))
                sxx = float(np.dot(tt[: yy.size], tt[: yy.size]))
                if sxx <= 0:
                    return (slope, (np.nan, np.nan), {**diag, "reason": "degenerate_time"})

                infl = (yy.size / neff)
                se = np.sqrt(sigma2 * infl / sxx)
                diag["se_slope"] = float(se)

                lo = slope - 1.96 * se
                hi = slope + 1.96 * se
                return (slope, (float(lo), float(hi)), diag)

        slope_raw, ci_raw, diag_raw = _trend_and_ci(raw_w)
        slope_filt, ci_filt, diag_filt = _trend_and_ci(filt_w)

        per_year_factor = float(spyr)
        unit_factor = per_year_factor * 10.0 if str(self.trend_unit).lower() == "per_decade" else per_year_factor

        def _scale_slope_and_ci(slope, ci):
            if not np.isfinite(slope):
                return (np.nan, (np.nan, np.nan))
            slo = slope * unit_factor
            clo, chi = ci
            if np.isfinite(clo) and np.isfinite(chi):
                return (float(slo), (float(clo * unit_factor), float(chi * unit_factor)))
            return (float(slo), (np.nan, np.nan))

        slope_raw_u, ci_raw_u = _scale_slope_and_ci(slope_raw, ci_raw)
        slope_filt_u, ci_filt_u = _scale_slope_and_ci(slope_filt, ci_filt)

        out.update(
            {
                "trend_slope_raw": slope_raw_u,
                "trend_slope_filt": slope_filt_u,
                "trend_ci95_raw": ci_raw_u,
                "trend_ci95_filt": ci_filt_u,
                "trend_diag_raw": diag_raw,
                "trend_diag_filt": diag_filt,
            }
        )

        out["n_used"] = int(max(out.get("n_used_raw", 0), out.get("n_used_filt", 0)))
        return out

    @staticmethod
    def _ols_slope(y: np.ndarray) -> float:
        """OLS slope of y vs t=0..n-1 (per sample)."""
        y = np.asarray(y, float).reshape(-1)
        m = np.isfinite(y)
        y = y[m]
        n = y.size
        if n < 3:
            return np.nan
        t = np.arange(n, dtype=float)
        tt = t - t.mean()
        denom = float(np.dot(tt, tt))
        if denom <= 0:
            return np.nan
        return float(np.dot(tt, y - y.mean()) / denom)

    @staticmethod
    def _block_bootstrap_samples(
        x: np.ndarray,
        *,
        block_len: int,
        n_boot: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Moving block bootstrap for a 1D series x (assumed regularly spaced).
        Returns array of shape (n_boot, n) containing bootstrap resamples.
        """
        x = np.asarray(x, float).reshape(-1)
        n = x.size
        if n == 0:
            return np.empty((0, 0), dtype=float)
        block_len = int(max(1, min(block_len, n)))
        n_blocks = int(math.ceil(n / block_len))

        # start indices for blocks (circular MBB to avoid edge bias)
        starts = rng.integers(0, n, size=(n_boot, n_blocks), endpoint=False)

        out = np.empty((n_boot, n_blocks * block_len), dtype=float)
        for b in range(n_blocks):
            s = starts[:, b]
            idx = (s[:, None] + np.arange(block_len)[None, :]) % n
            out[:, b * block_len : (b + 1) * block_len] = x[idx]

        return out[:, :n]

    @staticmethod
    def _two_sided_p_from_boot(boot: np.ndarray, obs: float) -> float:
        """
        Two-sided p-value from bootstrap distribution under H0: statistic == 0,
        using sign test around 0 (standard in paired-diff bootstrap).
        """
        boot = np.asarray(boot, float).reshape(-1)
        boot = boot[np.isfinite(boot)]
        if boot.size == 0 or not np.isfinite(obs):
            return np.nan
        # p = 2*min(P(T<=0), P(T>=0)) for bootstrap distribution of the statistic on differences
        p_lo = float(np.mean(boot <= 0.0))
        p_hi = float(np.mean(boot >= 0.0))
        p = 2.0 * min(p_lo, p_hi)
        return float(min(max(p, 0.0), 1.0))

    @staticmethod
    def _boot_ci(boot: np.ndarray, lo=2.5, hi=97.5) -> Tuple[float, float]:
        boot = np.asarray(boot, float).reshape(-1)
        boot = boot[np.isfinite(boot)]
        if boot.size == 0:
            return (np.nan, np.nan)
        return (float(np.percentile(boot, lo)), float(np.percentile(boot, hi)))

    def _paired_block_bootstrap_tests(
        self,
        y_ref: np.ndarray,
        y_alt: np.ndarray,
        *,
        block_len: int,
        n_boot: int,
        seed: Optional[int],
        stat: str,                   # "mean" | "trend"
        samples_per_year: int,        # 1 (annual) or 12 (monthly)
        trend_unit: str,              # per_year|per_decade
    ) -> Dict[str, Any]:
        """
        Paired block bootstrap on the *difference series* d = alt - ref.

        Returns:
          - obs_diff
          - p
          - ci95 (percentile)
          - boot_mean/boot_trend not returned (kept small)
        """
        y_ref = np.asarray(y_ref, float).reshape(-1)
        y_alt = np.asarray(y_alt, float).reshape(-1)

        n = int(min(y_ref.size, y_alt.size))
        y_ref = y_ref[:n]
        y_alt = y_alt[:n]

        m = np.isfinite(y_ref) & np.isfinite(y_alt)
        d = (y_alt[m] - y_ref[m]).astype(float)
        n_eff = d.size

        out: Dict[str, Any] = {
            "n_pair": int(n_eff),
            "block_len": int(block_len),
            "n_boot": int(n_boot),
            "stat": str(stat),
        }

        if n_eff < 8:
            out.update({"obs_diff": np.nan, "p": np.nan, "ci95": (np.nan, np.nan), "reason": "too_few_points"})
            return out

        rng = np.random.default_rng(None if seed is None else int(seed))

        if stat == "mean":
            obs = float(np.mean(d))
            bb = self._block_bootstrap_samples(d, block_len=block_len, n_boot=n_boot, rng=rng)
            boot_stat = np.nanmean(bb, axis=1)

        elif stat == "trend":
            # slope of d vs time (per sample), then convert to per_year/per_decade
            obs_slope = self._ols_slope(d)
            obs = float(obs_slope)

            bb = self._block_bootstrap_samples(d, block_len=block_len, n_boot=n_boot, rng=rng)
            boot_stat = np.array([self._ols_slope(bb[i]) for i in range(bb.shape[0])], dtype=float)

        else:
            raise ValueError(f"Unknown stat='{stat}' (expected 'mean' or 'trend').")

        # Convert trend units if requested
        if stat == "trend":
            per_year_factor = float(samples_per_year)
            unit_factor = per_year_factor * 10.0 if str(trend_unit).lower() == "per_decade" else per_year_factor
            obs = float(obs * unit_factor)
            boot_stat = boot_stat * unit_factor

        p = self._two_sided_p_from_boot(boot_stat, obs)
        ci = self._boot_ci(boot_stat, 2.5, 97.5)

        out.update({"obs_diff": obs, "p": p, "ci95": ci})
        return out

    @staticmethod
    def _json_sanitize(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, tuple):
            return [OCNDiagnosticsPlotter._json_sanitize(x) for x in obj]
        if isinstance(obj, dict):
            return {str(k): OCNDiagnosticsPlotter._json_sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [OCNDiagnosticsPlotter._json_sanitize(x) for x in obj]
        return obj

    def save_stats(
        self,
        *,
        out_base: str,
        stats_by_exp: Dict[str, Dict[str, Any]],
        fmt: str = "json",   # "json" | "csv" | "both"
    ) -> None:
        base = os.path.splitext(out_base)[0]
        fmt = str(fmt).lower()

        if fmt in ("json", "both"):
            fn = f"{base}_stats.json"
            payload = self._json_sanitize(stats_by_exp)
            with open(fn, "w") as f:
                json.dump(payload, f, indent=2)

        if fmt in ("csv", "both"):
            rows = []
            for exp, st in stats_by_exp.items():
                r = dict(st)

                for k in ("trend_ci95_raw", "trend_ci95_filt"):
                    if k in r and isinstance(r[k], (tuple, list)) and len(r[k]) == 2:
                        r[f"{k}_lo"] = r[k][0]
                        r[f"{k}_hi"] = r[k][1]

                # flatten sig CIs if present
                for k in ("bb_mean_raw_vs_ref", "bb_mean_filt_vs_ref", "bb_trend_raw_vs_ref", "bb_trend_filt_vs_ref"):
                    if k in r and isinstance(r[k], dict) and "ci95" in r[k]:
                        lo, hi = r[k]["ci95"]
                        r[f"{k}_ci95_lo"] = lo
                        r[f"{k}_ci95_hi"] = hi

                r.pop("trend_diag_raw", None)
                r.pop("trend_diag_filt", None)
                rows.append(self._json_sanitize(r))

            fn = f"{base}_stats.csv"
            pd.DataFrame(rows).to_csv(fn, index=False)

    def load_stats(self, *, out_base: str) -> Dict[str, Dict[str, Any]]:
        base = os.path.splitext(out_base)[0]
        fn = f"{base}_stats.json"
        if not os.path.exists(fn):
            raise FileNotFoundError(f"Missing stats file: {fn}")
        with open(fn, "r") as f:
            return json.load(f)

    def _set_time_axis(
        self,
        ax,
        *,
        step_yrs: int,
        xlim_years: Optional[Tuple[int, int]],
        ref_yrs: Optional[List[int]],
        fontsize: int,
        labels_zero_based: bool,
    ):
        if self.annual_mean:
            ny = self.year_max - self.year_min + 1
            ax.set_xlim(0, ny - 1)
            step = max(1, int(step_yrs))
            xt = np.arange(0, ny, step)
            labels = [str(int(t if labels_zero_based else (self.year_min + t))) for t in xt]
            ax.xaxis.set_major_locator(FixedLocator(xt))
            ax.xaxis.set_major_formatter(FixedFormatter(labels))
            ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))

            if ref_yrs:
                for y in ref_yrs:
                    xpos = int(y) - self.year_min
                    if 0 <= xpos < ny:
                        ax.axvline(xpos, color="k", ls="--", lw=1.2, alpha=0.6)

            if xlim_years:
                a, b = xlim_years
                ax.set_xlim(a - self.year_min, b - self.year_min)

        else:
            nm = (self.year_max - self.year_min + 1) * 12
            ax.set_xlim(0, nm)
            step = max(1, int(step_yrs))
            xt = np.arange(0, nm + 1, 12 * step)
            labels = [str(int((t / 12) if labels_zero_based else (self.year_min + t / 12.0))) for t in xt]
            ax.xaxis.set_major_locator(FixedLocator(xt))
            ax.xaxis.set_major_formatter(FixedFormatter(labels))
            ax.xaxis.set_minor_locator(AutoMinorLocator(n=4))

            if ref_yrs:
                for y in ref_yrs:
                    xpos = (int(y) - self.year_min) * 12
                    if 0 <= xpos <= nm:
                        ax.axvline(xpos, color="k", ls="--", lw=1.2, alpha=0.6)

            if xlim_years:
                a, b = xlim_years
                ax.set_xlim((a - self.year_min) * 12, (b - self.year_min) * 12)

    @staticmethod
    def _legend(ax, labels, colors, fontsize, *, loc=None, ncol=None, y=0.995):
        handles = [Line2D([0], [0], color=c, lw=4.0) for c in colors[:len(labels)]]
        if ncol is None:
            ncol = len(labels)

        return ax.legend(
            handles, labels,
            ncol=ncol,
            loc=loc if loc else "upper center",
            bbox_to_anchor=(0.05, y),
            fontsize=int(fontsize * 0.85),
            frameon=False,
            handlelength=2.2,
            handletextpad=0.6,
            columnspacing=1.4,
            borderaxespad=0.0,
            labelspacing=0.2,
        )

    def _shade_segments(
        self,
        ax,
        *,
        exp: str,
        segments,
        color: str,
        alpha: float = 0.12,
        hatch: str = "///",
        zorder: int = 0,
    ):
        if not segments:
            return

        for seg in segments:
            if not (isinstance(seg, (list, tuple)) and len(seg) == 2):
                raise ValueError(f"Bad segment for {exp}: {seg} (expected (y0, y1))")
            y0, y1 = seg
            y0, y1 = int(y0), int(y1)
            if y1 < y0:
                y0, y1 = y1, y0

            if self.annual_mean:
                x0 = y0 - self.year_min
                x1 = y1 - self.year_min + 1
            else:
                x0 = (y0 - self.year_min) * 12
                x1 = (y1 - self.year_min + 1) * 12

            # rasterized: PDF hatch patterns have a fixed tile size that does not shrink when the
            # figure is scaled in LaTeX (coarse, thick hatching); a bitmap scales with the figure
            ax.axvspan(
                x0, x1,
                facecolor="none",
                edgecolor=color,
                hatch=hatch,
                lw=0.0,
                alpha=alpha,
                zorder=zorder,
                rasterized=True,
            )

    def plot_on_axis(
        self,
        ax,
        data: Dict[str, Dict[str, np.ndarray]],
        *,
        var_key: str,
        exp_order: List[str],
        exp_labels: Optional[List[str]] = None,
        colors: Optional[List[str]] = None,
        req: Optional[PlotRequest] = None,
        fontsize: int = 18,
        step_yrs: int = 25,
        ref_yrs: Optional[List[int]] = None,
        xlim_years: Optional[Tuple[int, int]] = None,
        labels_zero_based: bool = False,
        show_raw: bool = False,
        show_raw_faint: bool = False,
        std_band: bool = True,
        std_center: str = "filtered",
        std_source: str = "raw",
        std_mode: str = "global",
        std_window_months: int = 12,
        std_alpha: float = 0.15,
        legend_loc: Optional[str] = None,
        exp_shas: Optional[Dict[str, List[Tuple[int, int]]]] = None,
        sha_alpha: float = 0.10,

        # stats saving
        save_stats_base: Optional[str] = None,
        stats_format: str = "json",

        # significance test controls (BLOCK BOOTSTRAP)
        ref_exp: str = "Full-CPL",
        sig_alpha: float = 0.05,
        bb_n: int = 2000,
        bb_block_years: int = 5,
        bb_seed: Optional[int] = 12345,
    ):
        if req is None:
            req = PlotRequest(key=var_key, vname=var_key, vunit="", scale=1.0)

        exp_labels = exp_labels or exp_order
        colors = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

        ymins, ymaxs = [], []
        used_colors = []

        stats_by_exp: Dict[str, Dict[str, Any]] = {}

        # Keep last-window slices for bootstrap tests (raw/filt)
        raw_win_by_exp: Dict[str, np.ndarray] = {}
        filt_win_by_exp: Dict[str, np.ndarray] = {}

        for i, exp in enumerate(exp_order):
            if exp not in data or var_key not in data[exp]:
                raise KeyError(f"Missing data[{exp}][{var_key}]")

            col = colors[i % len(colors)]
            used_colors.append(col)

            raw0 = np.asarray(data[exp][var_key], float).reshape(-1) * float(req.scale)

            if self.annual_mean:
                raw = self.to_annual_mean(raw0)
                samples_per_year = 1
            else:
                raw = raw0
                samples_per_year = 12

            if (not self.annual_mean) and self.running_mean and self.running_mean_months > 1:
                raw = self._nan_running_mean(
                    raw,
                    self.running_mean_months,
                    center=self.running_mean_center,
                    min_valid_frac=self.running_mean_min_frac,
                )

            if self.compute_anomaly:
                raw = self._to_anomaly(raw, samples_per_year=samples_per_year)

            filt = self._lowpass(raw) if self.use_filter else raw.copy()

            st = self._last_n_year_stats(raw, filt)
            st = dict(st)
            st.update({
                "exp": exp,
                "var_key": var_key,
                "vname": req.vname,
                "vunit": req.vunit,
                "annual_mean": bool(self.annual_mean),
                "compute_anomaly": bool(self.compute_anomaly),
                "use_filter": bool(self.use_filter),
            })
            stats_by_exp[exp] = st

            # last-window slices
            i0 = int(st["start_idx"]) if st.get("start_idx") is not None else 0
            i1 = int(st["end_idx"]) + 1 if st.get("end_idx") is not None else raw.size
            raw_win_by_exp[exp] = np.asarray(raw[i0:i1], float)
            filt_win_by_exp[exp] = np.asarray(filt[i0:i1], float)

            # shading behind
            sha = None if exp_shas is None else exp_shas.get(exp, None)
            if sha:
                if isinstance(sha, dict):
                    segments = sha.get("segments", None)
                    hatch = sha.get("hatch", "///")
                    alpha = sha.get("alpha", sha_alpha)
                    color_use = sha.get("color", col)
                else:
                    segments = sha
                    hatch = "///"
                    alpha = sha_alpha
                    color_use = col

                if segments:
                    self._shade_segments(ax, exp=exp, segments=segments,
                                        color=color_use, hatch=hatch, alpha=alpha, zorder=-50)

            # std band
            if std_band and raw.size:
                center = filt if (std_center.lower() == "filtered" and self.use_filter) else raw
                if std_source.lower() == "residual" and self.use_filter:
                    base = raw - filt
                elif std_source.lower() == "filtered" and self.use_filter:
                    base = filt
                else:
                    base = raw

                if std_mode.lower() == "rolling":
                    win = std_window_months if (not self.annual_mean) else max(1, std_window_months // 12)
                    sigma = self._rolling_nanstd(base, win)
                else:
                    s = float(np.nanstd(base, ddof=1)) if base.size > 1 else float(np.nanstd(base))
                    sigma = np.full_like(center, s, dtype=float)

                ax.fill_between(np.arange(center.size), center - sigma, center + sigma,
                                color=col, alpha=std_alpha, linewidth=0, zorder=-5)

                ymins.append(np.nanmin(center - sigma))
                ymaxs.append(np.nanmax(center + sigma))

            # plot
            if self.use_filter:
                if show_raw:
                    ax.plot(raw, color=col, lw=1.0, alpha=0.7, ls="--", zorder=5)
                elif show_raw_faint:
                    ax.plot(raw, color=col, lw=0.5, alpha=0.2, zorder=5)
                ax.plot(filt, color=col, lw=2.0, alpha=1.0, zorder=5)
                ysrc = filt
            else:
                ax.plot(raw, color=col, lw=2.0, alpha=1.0, zorder=5)
                ysrc = raw

            if ysrc.size:
                ymins.append(np.nanmin(ysrc))
                ymaxs.append(np.nanmax(ysrc))

        # axis formatting
        self._set_time_axis(ax, step_yrs=step_yrs, xlim_years=xlim_years, ref_yrs=ref_yrs,
                            fontsize=fontsize, labels_zero_based=labels_zero_based)

        title_default = req.title or (req.vname if not self.compute_anomaly else f"{req.vname} Anomaly (Δ{req.vname})")
        ylab_default = req.y_label or (f"{req.vname} ({req.vunit})" if req.vunit else req.vname)

        if req.show_title:
            ax.set_title(title_default, fontsize=fontsize, loc="left", pad=8)

        ax.set_xlabel("Model Time (year)", fontsize=int(fontsize * 0.95) if req.show_xlabel else 0)
        if not req.show_xlabel:
            ax.set_xlabel("")
        ax.tick_params(axis="x", labelbottom=bool(req.show_xticks))
        if not req.show_xticks:
            ax.tick_params(axis="x", which="both", length=0)

        if req.show_ylabel:
            ax.set_ylabel(ylab_default, fontsize=int(fontsize * 0.95))
        else:
            ax.set_ylabel("")

        if req.ylim:
            ax.set_ylim(*req.ylim)
        elif ymins and ymaxs:
            lo, hi = float(np.nanmin(ymins)), float(np.nanmax(ymaxs))
            pad = 0.05 * max(1e-12, hi - lo)
            ax.set_ylim(lo - pad, hi + pad)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ylow, yhigh = ax.get_ylim()
        if ylow < 0 < yhigh:
            ax.axhline(0.0, ls="--", lw=1.0, color="k", alpha=0.8)

        ax.grid(True, alpha=0.1, color="grey", linewidth=0.01)
        ax.tick_params(labelsize=int(fontsize * 0.85), length=4, width=1)

        if req.panel_label:
            x0, y0 = req.label_loc
            ax.text(x0, y0, req.panel_label, transform=ax.transAxes,
                    va="top", ha="left", fontsize=int(fontsize * 0.9))

        loc_use = legend_loc if legend_loc is not None else req.legend_loc
        if loc_use:
            self._legend(ax, exp_labels, used_colors, fontsize=fontsize, loc=loc_use, ncol=len(exp_labels))

        # ---------------- NEW: BLOCK-BOOTSTRAP significance vs reference ----------------
        ref_exp = str(ref_exp)
        alpha = float(sig_alpha)
        spyr = 1 if self.annual_mean else 12

        # choose block length in *samples* (years -> samples)
        bb_block_len = int(max(2, round(float(bb_block_years) * spyr)))

        if ref_exp in stats_by_exp:
            for exp in exp_order:
                if exp == ref_exp:
                    continue

                # mean difference on last-window RAW
                bb_mean_raw = self._paired_block_bootstrap_tests(
                    raw_win_by_exp.get(ref_exp, np.array([])),
                    raw_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=bb_seed,
                    stat="mean",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                # mean difference on last-window FILTERED
                bb_mean_flt = self._paired_block_bootstrap_tests(
                    filt_win_by_exp.get(ref_exp, np.array([])),
                    filt_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=None if bb_seed is None else int(bb_seed) + 1,
                    stat="mean",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                # trend difference on last-window RAW (trend of difference series)
                bb_trend_raw = self._paired_block_bootstrap_tests(
                    raw_win_by_exp.get(ref_exp, np.array([])),
                    raw_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=None if bb_seed is None else int(bb_seed) + 2,
                    stat="trend",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                # trend difference on last-window FILTERED
                bb_trend_flt = self._paired_block_bootstrap_tests(
                    filt_win_by_exp.get(ref_exp, np.array([])),
                    filt_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=None if bb_seed is None else int(bb_seed) + 3,
                    stat="trend",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                stats_by_exp[exp].update({
                    "ref_exp": ref_exp,
                    "sig_alpha": alpha,
                    "sig_method": "paired_block_bootstrap",
                    "bb_n": int(bb_n),
                    "bb_block_years": int(bb_block_years),
                    "bb_block_len_samples": int(bb_block_len),

                    # store full test dicts (sanitized on save)
                    "bb_mean_raw_vs_ref": bb_mean_raw,
                    "bb_mean_filt_vs_ref": bb_mean_flt,
                    "bb_trend_raw_vs_ref": bb_trend_raw,
                    "bb_trend_filt_vs_ref": bb_trend_flt,

                    # convenient booleans
                    "signif_mean_raw_vs_ref": (bool(np.isfinite(bb_mean_raw.get("p", np.nan)) and bb_mean_raw["p"] < alpha)),
                    "signif_mean_filt_vs_ref": (bool(np.isfinite(bb_mean_flt.get("p", np.nan)) and bb_mean_flt["p"] < alpha)),
                    "signif_trend_raw_vs_ref": (bool(np.isfinite(bb_trend_raw.get("p", np.nan)) and bb_trend_raw["p"] < alpha)),
                    "signif_trend_filt_vs_ref": (bool(np.isfinite(bb_trend_flt.get("p", np.nan)) and bb_trend_flt["p"] < alpha)),
                })

        # save stats if requested
        if save_stats_base:
            self.save_stats(out_base=save_stats_base, stats_by_exp=stats_by_exp, fmt=stats_format)

        return stats_by_exp

    def plot_from_extracted(
        self,
        data: Dict[str, Dict[str, np.ndarray]],
        *,
        var_key: str,
        exp_order: List[str],
        exp_labels: Optional[List[str]] = None,
        colors: Optional[List[str]] = None,
        req: Optional[PlotRequest] = None,
        outname: Optional[str] = None,
        figsize: Tuple[int, int] = (24, 6),
        fontsize: int = 28,
        step_yrs: int = 25,
        ref_yrs: Optional[List[int]] = None,
        xlim_years: Optional[Tuple[int, int]] = None,
        labels_zero_based: bool = False,
        show_raw: bool = False,
        show_raw_faint: bool = False,
        std_band: bool = True,
        std_center: str = "filtered",
        std_source: str = "raw",
        std_mode: str = "global",
        std_window_months: int = 12,
        std_alpha: float = 0.15,
        legend_loc: str = "upper left",
        exp_shas: Optional[Dict[str, List[Tuple[int, int]]]] = None,
        sha_alpha: float = 0.10,
        # (optional) pass-through for sig config if desired
        ref_exp: str = "Full-CPL",
        sig_alpha: float = 0.05,
        bb_n: int = 5000,
        bb_block_years: int = 5,
        bb_seed: Optional[int] = 20240207,
        save_stats_base: Optional[str] = None,
        stats_format: str = "json",
    ):
        if req is None:
            req = PlotRequest(key=var_key, vname=var_key, vunit="", scale=1.0)

        fig, ax = plt.subplots(figsize=figsize)

        stats_list = self.plot_on_axis(
            ax=ax,
            data=data,
            var_key=var_key,
            exp_order=exp_order,
            exp_labels=exp_labels,
            colors=colors,
            req=req,
            fontsize=fontsize,
            step_yrs=step_yrs,
            ref_yrs=ref_yrs,
            xlim_years=xlim_years,
            labels_zero_based=labels_zero_based,
            show_raw=show_raw,
            show_raw_faint=show_raw_faint,
            std_band=std_band,
            std_center=std_center,
            std_source=std_source,
            std_mode=std_mode,
            std_window_months=std_window_months,
            std_alpha=std_alpha,
            legend_loc=legend_loc,
            exp_shas=exp_shas,
            sha_alpha=sha_alpha,
            ref_exp=ref_exp,
            sig_alpha=sig_alpha,
            bb_n=bb_n,
            bb_block_years=bb_block_years,
            bb_seed=bb_seed,
            save_stats_base=save_stats_base,
            stats_format=stats_format,
        )

        fig.set_constrained_layout(False)
        fig.subplots_adjust(wspace=0.18, hspace=0.25)

        if outname:
            fig.savefig(outname, dpi=600, bbox_inches="tight")

        return fig, ax, stats_list


class OHUDiagnosticsPlotter(OCNDiagnosticsPlotter):
    """Raw-vs-filtered heat uptake panels with Butterworth low-pass (timeseries_alternate/ohu_ts)."""

    def __init__(
        self,
        *,
        year_min: int,
        year_max: int,
        annual_mean: bool = True,
        compute_anomaly: bool = False,
        baseline_months: int = 12,
        baseline_years: int = 1,

        running_mean: bool = False,
        running_mean_months: int = 0,
        running_mean_center: str = "center",   # center|trailing|leading
        running_mean_min_frac: float = 0.5,

        use_filter: bool = False,
        filt_cutoff_months: float = 24.0,      # used when annual_mean=False OR when filt_cutoff_years is None
        filt_cutoff_years: Optional[float] = None,  # NEW: used when annual_mean=True
        filt_order: int = 4,
        filt_method: str = "pad",          # pad|gust
        filt_padlen: Optional[int] = None,
        filt_padtype: str = "odd",
        edge_fix_points: int = 0,

        trend_years: int = 30,
        trend_unit: str = "per_year",      # per_year|per_decade
        add_trend_ci: bool = True,
        ci_method: str = "hac",            # hac|glsar|neff
        hac_lag_years: int = 1,
        glsar_max_iter: int = 10,
    ):
        self.year_min = int(year_min)
        self.year_max = int(year_max)

        self.annual_mean = bool(annual_mean)
        self.compute_anomaly = bool(compute_anomaly)
        self.baseline_months = int(baseline_months)
        self.baseline_years = int(baseline_years)

        self.running_mean = bool(running_mean)
        self.running_mean_months = int(running_mean_months or 0)
        self.running_mean_center = str(running_mean_center).lower()
        self.running_mean_min_frac = float(running_mean_min_frac)

        self.use_filter = bool(use_filter)
        self.filt_cutoff_months = float(filt_cutoff_months)
        self.filt_cutoff_years = None if filt_cutoff_years is None else float(filt_cutoff_years)
        self.filt_order = int(filt_order)
        self.filt_method = str(filt_method)
        self.filt_padlen = None if filt_padlen is None else int(filt_padlen)
        self.filt_padtype = str(filt_padtype)
        self.edge_fix_points = int(edge_fix_points)

        self.trend_years = int(trend_years)
        self.trend_unit = str(trend_unit)
        self.add_trend_ci = bool(add_trend_ci)
        self.ci_method = str(ci_method)
        self.hac_lag_years = int(hac_lag_years)
        self.glsar_max_iter = int(glsar_max_iter)

        # butter cache (depends on sample rate + cutoff)
        self._butter_b = None
        self._butter_a = None
        self._butter_spyr = None
        self._butter_cutoff_years = None

    def _effective_cutoff_years(self, samples_per_year: int) -> float:
        """
        If annual_mean=True and filt_cutoff_years is provided, use that directly.
        Otherwise convert filt_cutoff_months -> years.
        """
        if self.annual_mean and (self.filt_cutoff_years is not None):
            return float(self.filt_cutoff_years)
        return float(self.filt_cutoff_months) / 12.0

    def _ensure_butter(self, samples_per_year: int) -> None:
        if not self.use_filter:
            return

        spyr = int(samples_per_year)
        cutoff_years = self._effective_cutoff_years(spyr)

        if (
            self._butter_b is not None
            and self._butter_a is not None
            and self._butter_spyr == spyr
            and self._butter_cutoff_years == cutoff_years
        ):
            return

        # Wn = fc/fN = (1/T)/(0.5*fs) = 2/(T*fs)
        wn = 2.0 / (max(cutoff_years, 1e-12) * max(spyr, 1))
        wn = float(np.clip(wn, 1e-6, 0.999999))

        self._butter_b, self._butter_a = signal.butter(self.filt_order, wn, btype="low")
        self._butter_spyr = spyr
        self._butter_cutoff_years = cutoff_years

    def _lowpass(self, x: np.ndarray, *, samples_per_year: int) -> np.ndarray:
        if not self.use_filter:
            return np.asarray(x, float)

        x = np.asarray(x, float).reshape(-1)
        n = len(x)
        if n < 3:
            return x.copy()

        self._ensure_butter(samples_per_year)

        b, a = self._butter_b, self._butter_a
        if b is None or a is None:
            return x.copy()

        if self.filt_method == "gust":
            y = signal.filtfilt(b, a, x, method="gust")
        else:
            default_pad = 3 * (max(len(a), len(b)) - 1)

            cutoff_years = self._effective_cutoff_years(samples_per_year)
            cutoff_samples = int(max(1, round(2.0 * cutoff_years * samples_per_year)))

            padlen = self.filt_padlen if self.filt_padlen is not None else max(default_pad, cutoff_samples)
            padlen = min(padlen, max(1, n - 1))

            if padlen >= n:
                y = signal.filtfilt(b, a, x, method="gust")
            else:
                y = signal.filtfilt(b, a, x, method="pad", padtype=self.filt_padtype, padlen=padlen)

        k = self.edge_fix_points
        if k > 0 and n >= 2 * k + 1:
            y[:k] = y[k]
            y[-k:] = y[-k - 1]
        return y

    @staticmethod
    def _json_sanitize(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, tuple):
            return [OHUDiagnosticsPlotter._json_sanitize(x) for x in obj]
        if isinstance(obj, dict):
            return {str(k): OHUDiagnosticsPlotter._json_sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [OHUDiagnosticsPlotter._json_sanitize(x) for x in obj]
        return obj

    @staticmethod
    def _legend(ax, labels, colors, fontsize, *, loc=None, ncol=None):
        handles = [Line2D([0], [0], color=c, lw=4.0) for c in colors[:len(labels)]]
        if ncol is None:
            ncol = len(labels)
        if "upper" in loc:
            y = 0.995
        elif "lower" in loc:
            y = 0.05
        else:
            y = 0.5
        return ax.legend(
            handles, labels,
            ncol=ncol,
            loc=loc if loc else "upper center",
            bbox_to_anchor=(0.05, y),
            fontsize=int(fontsize * 0.85),
            frameon=False,
            handlelength=2.2,
            handletextpad=0.6,
            columnspacing=1.4,
            borderaxespad=0.0,
            labelspacing=0.2,
        )

    def plot_on_axis(
        self,
        ax,
        data: Dict[str, Dict[str, np.ndarray]],
        *,
        var_key: str,
        exp_order: List[str],
        exp_labels: Optional[List[str]] = None,
        colors: Optional[List[str]] = None,
        req: Optional[PlotRequest] = None,
        fontsize: int = 18,
        step_yrs: int = 25,
        ref_yrs: Optional[List[int]] = None,
        xlim_years: Optional[Tuple[int, int]] = None,
        labels_zero_based: bool = False,
        show_raw: bool = False,
        show_raw_faint: bool = False,
        std_band: bool = True,
        std_center: str = "filtered",
        std_source: str = "raw",
        std_mode: str = "global",
        std_window_months: int = 12,
        std_alpha: float = 0.15,
        legend_loc: Optional[str] = None,
        exp_shas: Optional[Dict[str, List[Tuple[int, int]]]] = None,
        sha_alpha: float = 0.10,

        save_stats_base: Optional[str] = None,
        stats_format: str = "json",

        ref_exp: str = "Full-CPL",
        sig_alpha: float = 0.05,
        bb_n: int = 2000,
        bb_block_years: int = 5,
        bb_seed: Optional[int] = 12345,
    ):
        if req is None:
            req = PlotRequest(key=var_key, vname=var_key, vunit="", scale=1.0)

        exp_labels = exp_labels or exp_order
        colors = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]        
        ltstr,rtstr = req.title

        ymins, ymaxs = [], []
        used_colors = []
        stats_by_exp: Dict[str, Dict[str, Any]] = {}

        raw_win_by_exp: Dict[str, np.ndarray] = {}
        filt_win_by_exp: Dict[str, np.ndarray] = {}

        for i, exp in enumerate(exp_order):
            if exp not in data or var_key not in data[exp]:
                raise KeyError(f"Missing data[{exp}][{var_key}]")

            col = colors[i % len(colors)]
            used_colors.append(col)

            raw0 = np.asarray(data[exp][var_key], float).reshape(-1) * float(req.scale)

            if self.annual_mean:
                raw = self.to_annual_mean(raw0)
                samples_per_year = 1
            else:
                raw = raw0
                samples_per_year = 12

            if (not self.annual_mean) and self.running_mean and self.running_mean_months > 1:
                raw = self._nan_running_mean(
                    raw,
                    self.running_mean_months,
                    center=self.running_mean_center,
                    min_valid_frac=self.running_mean_min_frac,
                )

            if self.compute_anomaly:
                raw = self._to_anomaly(raw, samples_per_year=samples_per_year)

            filt = self._lowpass(raw, samples_per_year=samples_per_year) if self.use_filter else raw.copy()

            st = self._last_n_year_stats(raw, filt)
            st = dict(st)
            st.update({
                "exp": exp,
                "var_key": var_key,
                "vname": req.vname,
                "vunit": req.vunit,
                "annual_mean": bool(self.annual_mean),
                "compute_anomaly": bool(self.compute_anomaly),
                "use_filter": bool(self.use_filter),
            })
            stats_by_exp[exp] = st

            i0 = int(st["start_idx"]) if st.get("start_idx") is not None else 0
            i1 = int(st["end_idx"]) + 1 if st.get("end_idx") is not None else raw.size
            raw_win_by_exp[exp] = np.asarray(raw[i0:i1], float)
            filt_win_by_exp[exp] = np.asarray(filt[i0:i1], float)

            sha = None if exp_shas is None else exp_shas.get(exp, None)
            if sha:
                if isinstance(sha, dict):
                    segments = sha.get("segments", None)
                    hatch = sha.get("hatch", "///")
                    alpha = sha.get("alpha", sha_alpha)
                    color_use = sha.get("color", col)
                else:
                    segments = sha
                    hatch = "///"
                    alpha = sha_alpha
                    color_use = col

                if segments:
                    self._shade_segments(ax, exp=exp, segments=segments,
                                         color=color_use, hatch=hatch, alpha=alpha, zorder=-50)

            if std_band and raw.size:
                center = filt if (std_center.lower() == "filtered" and self.use_filter) else raw
                if std_source.lower() == "residual" and self.use_filter:
                    base = raw - filt
                elif std_source.lower() == "filtered" and self.use_filter:
                    base = filt
                else:
                    base = raw

                if std_mode.lower() == "rolling":
                    win = std_window_months if (not self.annual_mean) else max(1, std_window_months // 12)
                    sigma = self._rolling_nanstd(base, win)
                else:
                    s = float(np.nanstd(base, ddof=1)) if base.size > 1 else float(np.nanstd(base))
                    sigma = np.full_like(center, s, dtype=float)

                ax.fill_between(np.arange(center.size), center - sigma, center + sigma,
                                color=col, alpha=std_alpha, linewidth=0, zorder=-5)

                ymins.append(np.nanmin(center - sigma))
                ymaxs.append(np.nanmax(center + sigma))

            if self.use_filter:
                if show_raw:
                    ax.plot(raw, color=col, lw=1.0, alpha=0.7, ls="--", zorder=5)
                elif show_raw_faint:
                    ax.plot(raw, color=col, lw=0.5, alpha=0.2, zorder=5)
                ax.plot(filt, color=col, lw=2.0, alpha=1.0, zorder=5)
                ysrc = filt
            else:
                ax.plot(raw, color=col, lw=2.0, alpha=1.0, zorder=5)
                ysrc = raw

            if ysrc.size:
                ymins.append(np.nanmin(ysrc))
                ymaxs.append(np.nanmax(ysrc))

        self._set_time_axis(ax, step_yrs=step_yrs, xlim_years=xlim_years, ref_yrs=ref_yrs,
                            fontsize=fontsize, labels_zero_based=labels_zero_based)

        title_default = req.title or (req.vname if not self.compute_anomaly else f"{req.vname} Anomaly (Δ{req.vname})")
        ylab_default = req.y_label or (f"{req.vname} ({req.vunit})" if req.vunit else req.vname)

        if req.show_title:
            if ltstr is not None and rtstr is not None:
                ax.set_title(f"{req.panel_label} {ltstr}", fontsize=fontsize, loc="left", pad=8)
                ax.set_title(rtstr, fontsize=fontsize, loc="right", pad=8)
            else:
                ax.set_title(f"{req.panel_label} {title_default}", fontsize=fontsize, loc="left", pad=8)

        ax.set_xlabel("Model Time (year)", fontsize=int(fontsize * 0.95) if req.show_xlabel else 0)
        if not req.show_xlabel:
            ax.set_xlabel("")
        ax.tick_params(axis="x", labelbottom=bool(req.show_xticks))
        if not req.show_xticks:
            ax.tick_params(axis="x", which="both", length=0)

        if req.show_ylabel:
            ax.set_ylabel(ylab_default, fontsize=int(fontsize * 0.95))
        else:
            ax.set_ylabel("")

        if req.ylim:
            ax.set_ylim(*req.ylim)
        elif ymins and ymaxs:
            lo, hi = float(np.nanmin(ymins)), float(np.nanmax(ymaxs))
            pad = 0.05 * max(1e-12, hi - lo)
            ax.set_ylim(lo - pad, hi + pad)

        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ylow, yhigh = ax.get_ylim()
        if ylow < 0 < yhigh:
            ax.axhline(0.0, ls="--", lw=1.0, color="k", alpha=0.8)

        ax.grid(True, alpha=0.1, color="grey", linewidth=0.01)
        ax.tick_params(labelsize=int(fontsize * 0.85), length=4, width=1)

        loc_use = legend_loc if legend_loc is not None else req.legend_loc
        if loc_use:
            self._legend(ax, exp_labels, used_colors, fontsize=fontsize, loc=loc_use, ncol=len(exp_labels))

        # ---- block-bootstrap significance vs reference (unchanged behavior) ----
        ref_exp = str(ref_exp)
        alpha = float(sig_alpha)
        spyr = 1 if self.annual_mean else 12
        bb_block_len = int(max(2, round(float(bb_block_years) * spyr)))

        if ref_exp in stats_by_exp:
            for exp in exp_order:
                if exp == ref_exp:
                    continue

                bb_mean_raw = self._paired_block_bootstrap_tests(
                    raw_win_by_exp.get(ref_exp, np.array([])),
                    raw_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=bb_seed,
                    stat="mean",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                bb_mean_flt = self._paired_block_bootstrap_tests(
                    filt_win_by_exp.get(ref_exp, np.array([])),
                    filt_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=None if bb_seed is None else int(bb_seed) + 1,
                    stat="mean",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                bb_trend_raw = self._paired_block_bootstrap_tests(
                    raw_win_by_exp.get(ref_exp, np.array([])),
                    raw_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=None if bb_seed is None else int(bb_seed) + 2,
                    stat="trend",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                bb_trend_flt = self._paired_block_bootstrap_tests(
                    filt_win_by_exp.get(ref_exp, np.array([])),
                    filt_win_by_exp.get(exp, np.array([])),
                    block_len=bb_block_len,
                    n_boot=int(bb_n),
                    seed=None if bb_seed is None else int(bb_seed) + 3,
                    stat="trend",
                    samples_per_year=spyr,
                    trend_unit=self.trend_unit,
                )

                stats_by_exp[exp].update({
                    "ref_exp": ref_exp,
                    "sig_alpha": alpha,
                    "sig_method": "paired_block_bootstrap",
                    "bb_n": int(bb_n),
                    "bb_block_years": int(bb_block_years),
                    "bb_block_len_samples": int(bb_block_len),

                    "bb_mean_raw_vs_ref": bb_mean_raw,
                    "bb_mean_filt_vs_ref": bb_mean_flt,
                    "bb_trend_raw_vs_ref": bb_trend_raw,
                    "bb_trend_filt_vs_ref": bb_trend_flt,

                    "signif_mean_raw_vs_ref": bool(np.isfinite(bb_mean_raw.get("p", np.nan)) and bb_mean_raw["p"] < alpha),
                    "signif_mean_filt_vs_ref": bool(np.isfinite(bb_mean_flt.get("p", np.nan)) and bb_mean_flt["p"] < alpha),
                    "signif_trend_raw_vs_ref": bool(np.isfinite(bb_trend_raw.get("p", np.nan)) and bb_trend_raw["p"] < alpha),
                    "signif_trend_filt_vs_ref": bool(np.isfinite(bb_trend_flt.get("p", np.nan)) and bb_trend_flt["p"] < alpha),
                })

        if save_stats_base:
            self.save_stats(out_base=save_stats_base, stats_by_exp=stats_by_exp, fmt=stats_format)

        return stats_by_exp

    def plot_raw_vs_filtered_from_extracted(
        self,
        data: Dict[str, Dict[str, np.ndarray]],
        *,
        var_key: str,
        exp_order: List[str],
        exp_labels: Optional[List[str]] = None,
        colors: Optional[List[str]] = None,
        req: Optional[PlotRequest] = None,
        outname: Optional[str] = None,
        figsize: Tuple[int, int] = (22, 6),
        fontsize: int = 28,
        step_yrs: int = 25,
        ref_yrs: Optional[List[int]] = None,
        xlim_years: Optional[Tuple[int, int]] = None,
        labels_zero_based: bool = False,
        show_raw_faint: bool = True,
        std_band: bool = True,
        std_source: str = "raw",
        std_mode: str = "global",
        std_window_months: int = 12,
        std_alpha: float = 0.15,
        legend_loc: str = "upper left",
        legend_panels: str = "both",   # "both" | "left" | "right" | "none"
        wspace: float = 0.18,

        exp_shas: Optional[Dict[str, List[Tuple[int, int]]]] = None,
        sha_alpha: float = 0.10,

        ref_exp: str = "Full-CPL",
        sig_alpha: float = 0.05,
        bb_n: int = 5000,
        bb_block_years: int = 5,
        bb_seed: Optional[int] = 20240207,

        save_stats_base: Optional[str] = None,   # writes *_raw_stats.* and *_filt_stats.*
        stats_format: str = "json",
    ):
        if req is None:
            req = PlotRequest(key=var_key, vname=var_key, vunit="", scale=1.0)

        fig, (axL, axR) = plt.subplots(1, 2, figsize=figsize, sharey=False)
        fig.subplots_adjust(wspace=wspace)

        def _legend_for(panel: str) -> Optional[str]:
            if legend_panels == "none":
                return None
            if legend_panels == "both":
                return legend_loc
            if legend_panels == "left":
                return legend_loc if panel == "left" else None
            if legend_panels == "right":
                return legend_loc if panel == "right" else None
            return legend_loc

        # preserve & restore
        _use_filter0 = bool(self.use_filter)
        _butter_b0, _butter_a0 = self._butter_b, self._butter_a
        _butter_spyr0, _butter_cut0 = self._butter_spyr, self._butter_cutoff_years

        try:
            # ---- Left: RAW ----
            self.use_filter = False
            reqL = replace(req)
            reqL.title = (f"{req.vname}","(Annual Mean)")
            stats_raw = self.plot_on_axis(
                ax=axL,
                data=data,
                var_key=var_key,
                exp_order=exp_order,
                exp_labels=exp_labels,
                colors=colors,
                req=reqL,
                fontsize=fontsize,
                step_yrs=step_yrs,
                ref_yrs=ref_yrs,
                xlim_years=xlim_years,
                labels_zero_based=labels_zero_based,
                show_raw=False,
                show_raw_faint=False,
                std_band=std_band,
                std_center="raw",
                std_source=std_source,
                std_mode=std_mode,
                std_window_months=std_window_months,
                std_alpha=std_alpha,
                legend_loc=_legend_for("left"),
                exp_shas=exp_shas,
                sha_alpha=sha_alpha,
                save_stats_base=(None if not save_stats_base else f"{save_stats_base}_raw"),
                stats_format=stats_format,
                ref_exp=ref_exp,
                sig_alpha=sig_alpha,
                bb_n=bb_n,
                bb_block_years=bb_block_years,
                bb_seed=bb_seed,
            )

            # ---- Right: FILTERED ----
            self.use_filter = True
            reqR = replace(req)

            # Make the title explicitly "10-year low-pass" when annual_mean=True and filt_cutoff_years is set
            if self.annual_mean:
                cutoff_yrs = self.filt_cutoff_years if self.filt_cutoff_years is not None else (self.filt_cutoff_months / 12.0)
                reqR.title = (f"{req.vname}",f"({cutoff_yrs:g}-yr Filter)")
            else:
                reqR.title = (f"{req.vname}",f"({self.filt_cutoff_months:g}-month Filter)")

            stats_filt = self.plot_on_axis(
                ax=axR,
                data=data,
                var_key=var_key,
                exp_order=exp_order,
                exp_labels=exp_labels,
                colors=colors,
                req=reqR,
                fontsize=fontsize,
                step_yrs=step_yrs,
                ref_yrs=ref_yrs,
                xlim_years=xlim_years,
                labels_zero_based=labels_zero_based,
                show_raw=False,
                show_raw_faint=bool(show_raw_faint),
                std_band=std_band,
                std_center="filtered",
                std_source=std_source,
                std_mode=std_mode,
                std_window_months=std_window_months,
                std_alpha=std_alpha,
                legend_loc=_legend_for("right"),
                exp_shas=exp_shas,
                sha_alpha=sha_alpha,
                save_stats_base=(None if not save_stats_base else f"{save_stats_base}_filt"),
                stats_format=stats_format,
                ref_exp=ref_exp,
                sig_alpha=sig_alpha,
                bb_n=bb_n,
                bb_block_years=bb_block_years,
                bb_seed=bb_seed,
            )

            fig.set_constrained_layout(False)

            if outname:
                fig.savefig(outname, dpi=600, bbox_inches="tight")

            return fig, (axL, axR), (stats_raw, stats_filt)

        finally:
            self.use_filter = _use_filter0
            self._butter_b, self._butter_a = _butter_b0, _butter_a0
            self._butter_spyr, self._butter_cutoff_years = _butter_spyr0, _butter_cut0


class OHCDiagnosticsPlotter(OCNDiagnosticsPlotter):
    """OHC time series + Hovmoller figure and stats JSON (timeseries_alternate/ohc_ts)."""

    def _to_anomaly(self, x: np.ndarray, *, samples_per_year: int) -> np.ndarray:
        x = np.asarray(x, float).reshape(-1)
        if samples_per_year == 1:
            n0 = self.baseline_years
        else:
            n0 = self.baseline_months
        if not (0 < n0 <= len(x)):
            return x - np.nanmean(x)
        clim = float(np.nanmean(x[:n0]))
        return x - clim

    def _last_n_year_stats(self, raw: np.ndarray, filt: np.ndarray) -> Dict[str, Any]:
        """
        Compute last-N-year summary + trend diagnostics.
    
        Parameters
        ----------
        raw, filt : np.ndarray
            1D series (annual if self.annual_mean else monthly), same length.
    
        Returns
        -------
        stats : dict
            Keys include:
              - n_total, n_used, window_years, start_idx, end_idx
              - raw_mean, raw_std, filt_mean, filt_std
              - trend_slope_{raw|filt}  (per year or per decade depending on self.trend_unit)
              - trend_ci95_{raw|filt}   ((lo, hi) in same units), if enabled and possible
              - method used, and some supporting diagnostics
        """
        raw = np.asarray(raw, float).reshape(-1)
        filt = np.asarray(filt, float).reshape(-1)
        n = int(min(raw.size, filt.size))
    
        out: Dict[str, Any] = {
            "n_total": n,
            "window_years": int(self.trend_years),
            "trend_unit": str(self.trend_unit),
            "ci_method": str(self.ci_method),
            "add_trend_ci": bool(self.add_trend_ci),
        }
    
        if n == 0:
            out.update(
                {
                    "n_used": 0,
                    "start_idx": None,
                    "end_idx": None,
                    "raw_mean": np.nan,
                    "raw_std": np.nan,
                    "filt_mean": np.nan,
                    "filt_std": np.nan,
                    "trend_slope_raw": np.nan,
                    "trend_slope_filt": np.nan,
                    "trend_ci95_raw": (np.nan, np.nan),
                    "trend_ci95_filt": (np.nan, np.nan),
                }
            )
            return out
    
        # samples-per-year depends on whether we're plotting annual or monthly
        spyr = 1 if self.annual_mean else 12
    
        # choose last window length in samples
        win = max(2, int(round(self.trend_years * spyr)))
        i1 = n  # exclusive
        i0 = max(0, n - win)
        out["start_idx"] = int(i0)
        out["end_idx"] = int(i1 - 1)
    
        raw_w = raw[i0:i1]
        filt_w = filt[i0:i1]
    
        # mask NaNs consistently per series
        def _basic_stats(x: np.ndarray) -> Tuple[float, float, int]:
            m = np.isfinite(x)
            if m.sum() == 0:
                return (np.nan, np.nan, 0)
            xx = x[m]
            std = float(np.nanstd(xx, ddof=1)) if xx.size > 1 else float(np.nanstd(xx))
            return (float(np.nanmean(xx)), std, int(xx.size))
    
        raw_mean, raw_std, n_raw = _basic_stats(raw_w)
        filt_mean, filt_std, n_filt = _basic_stats(filt_w)
    
        out.update(
            {
                "raw_mean": raw_mean,
                "raw_std": raw_std,
                "filt_mean": filt_mean,
                "filt_std": filt_std,
                "n_used_raw": n_raw,
                "n_used_filt": n_filt,
            }
        )
    
        # trend helper
        def _trend_and_ci(y: np.ndarray) -> Tuple[float, Tuple[float, float], Dict[str, Any]]:
            """
            Returns slope (in units per sample), CI (lo,hi) for slope, and diagnostics dict.
            Slope later converted to per-year/per-decade as requested.
            """
            diag: Dict[str, Any] = {}
            m = np.isfinite(y)
            if m.sum() < 3:
                return (np.nan, (np.nan, np.nan), {"reason": "too_few_points"})
    
            yy = y[m]
            t = np.arange(yy.size, dtype=float)
    
            # OLS slope (per sample index)
            # slope = cov(t,yy)/var(t)
            tt = t - t.mean()
            slope = float(np.dot(tt, yy - yy.mean()) / np.dot(tt, tt))
            diag["n"] = int(yy.size)
            diag["slope_per_sample"] = slope
    
            if not self.add_trend_ci:
                return (slope, (np.nan, np.nan), diag)
    
            method = str(self.ci_method).lower()
    
            # --- Try statsmodels-based methods if available ---
            try:
                import statsmodels.api as sm
                X = sm.add_constant(t)
    
                if method == "hac":
                    # Newey-West HAC
                    # lag in samples: hac_lag_years * spyr
                    maxlags = max(1, int(round(self.hac_lag_years * spyr)))
                    res = sm.OLS(yy, X).fit(cov_type="HAC", cov_kwds={"maxlags": maxlags})
                    se = float(res.bse[1])
                    diag.update({"method_used": "hac", "hac_maxlags": maxlags})
    
                elif method == "glsar":
                    # GLSAR AR(1)
                    # You can generalize order, but AR(1) is usually enough here
                    model = sm.GLSAR(yy, X, rho=1)
                    res = model.iterative_fit(maxiter=int(self.glsar_max_iter))
                    se = float(res.bse[1])
                    diag.update({"method_used": "glsar", "glsar_maxiter": int(self.glsar_max_iter)})
    
                else:
                    raise ImportError("fallback_to_neff")
    
                # 95% normal approx CI
                lo = slope - 1.96 * se
                hi = slope + 1.96 * se
                diag["se_slope"] = se
                return (slope, (float(lo), float(hi)), diag)
    
            except Exception:
                # --- NEFF fallback (lag-1 autocorr approximation) ---
                diag["method_used"] = "neff"
    
                # residuals from OLS fit
                yhat = (yy.mean() + slope * (t - t.mean()))
                r = yy - yhat
    
                # lag-1 autocorr (guarded)
                if r.size < 3:
                    return (slope, (np.nan, np.nan), {**diag, "reason": "too_few_for_neff"})
    
                r0 = r[:-1]
                r1 = r[1:]
                denom = np.sqrt(np.dot(r0, r0) * np.dot(r1, r1))
                rho1 = float(np.dot(r0, r1) / denom) if denom > 0 else 0.0
                rho1 = float(np.clip(rho1, -0.99, 0.99))
                diag["rho1"] = rho1
    
                neff = yy.size * (1.0 - rho1) / (1.0 + rho1)
                neff = float(np.clip(neff, 3.0, yy.size))
                diag["neff"] = neff
    
                # approximate SE(slope) under iid with neff
                # Var(slope) = sigma^2 / sum((t - mean)^2)
                sigma2 = float(np.nanvar(r, ddof=1)) if r.size > 1 else float(np.nanvar(r))
                sxx = float(np.dot(tt[: yy.size], tt[: yy.size]))  # same as sum((t-mean)^2)
                if sxx <= 0:
                    return (slope, (np.nan, np.nan), {**diag, "reason": "degenerate_time"})
    
                # scale sigma^2 by (n/neff) to inflate uncertainty
                infl = (yy.size / neff)
                se = np.sqrt(sigma2 * infl / sxx)
                diag["se_slope"] = float(se)
    
                lo = slope - 1.96 * se
                hi = slope + 1.96 * se
                return (slope, (float(lo), float(hi)), diag)
    
        # compute slopes per-sample on raw and filt window
        slope_raw, ci_raw, diag_raw = _trend_and_ci(raw_w)
        slope_filt, ci_filt, diag_filt = _trend_and_ci(filt_w)
    
        # convert slope units: per-sample -> per-year (or per-decade)
        # per-sample corresponds to 1 time step (month or year)
        per_year_factor = float(spyr)  # slope per month * 12 = per year; per year * 1 = per year
        if str(self.trend_unit).lower() == "per_decade":
            unit_factor = per_year_factor * 10.0
        else:
            unit_factor = per_year_factor
    
        def _scale_slope_and_ci(slope, ci):
            if not np.isfinite(slope):
                return (np.nan, (np.nan, np.nan))
            slo = slope * unit_factor
            clo, chi = ci
            if np.isfinite(clo) and np.isfinite(chi):
                return (float(slo), (float(clo * unit_factor), float(chi * unit_factor)))
            return (float(slo), (np.nan, np.nan))
    
        slope_raw_u, ci_raw_u = _scale_slope_and_ci(slope_raw, ci_raw)
        slope_filt_u, ci_filt_u = _scale_slope_and_ci(slope_filt, ci_filt)
    
        out.update(
            {
                "trend_slope_raw": slope_raw_u,
                "trend_slope_filt": slope_filt_u,
                "trend_ci95_raw": ci_raw_u,
                "trend_ci95_filt": ci_filt_u,
                "trend_diag_raw": diag_raw,
                "trend_diag_filt": diag_filt,
            }
        )
    
        # a convenience "n_used" for the window (max of the two)
        out["n_used"] = int(max(out.get("n_used_raw", 0), out.get("n_used_filt", 0)))
        return out

    @staticmethod
    def _json_sanitize(obj):
        if isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, tuple):
            return [OHCDiagnosticsPlotter._json_sanitize(x) for x in obj]
        if isinstance(obj, dict):
            return {str(k): OHCDiagnosticsPlotter._json_sanitize(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [OHCDiagnosticsPlotter._json_sanitize(x) for x in obj]
        return obj

    @staticmethod
    def _json_sanitize_deep(obj):
        """JSON-safe conversion for nested dict/list + numpy scalars/arrays."""
        if obj is None:
            return None
        if isinstance(obj, (str, int, float, bool)):
            return obj
        if isinstance(obj, (np.generic,)):
            return obj.item()
        if isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        if isinstance(obj, (tuple, list)):
            return [OHCDiagnosticsPlotter._json_sanitize_deep(x) for x in obj]
        if isinstance(obj, dict):
            return {str(k): OHCDiagnosticsPlotter._json_sanitize_deep(v) for k, v in obj.items()}
        return str(obj)

    @classmethod
    def save_figure_stats_json(
        cls,
        *,
        out_stats_json: str,
        stats_by_var: Dict[str, Dict[str, Any]],
        meta: Optional[Dict[str, Any]] = None,
        indent: int = 2,
        verbose: bool = True,
    ) -> None:
        """
        Write one stats JSON for a multi-panel figure.

        Structure:
          {
            "meta": {...},
            "stats": {
              "<var_key>": {
                "<exp>": { ...stats fields... },
                ...
              },
              ...
            }
          }
        """
        if not out_stats_json:
            return

        payload = {
            "meta": cls._json_sanitize_deep(meta or {}),
            "stats": cls._json_sanitize_deep(stats_by_var),
        }
        os.makedirs(os.path.dirname(out_stats_json) or ".", exist_ok=True)
        with open(out_stats_json, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=indent, sort_keys=True)
        if verbose:
            print("[saved STATS]", out_stats_json)

    @staticmethod
    def _legend(ax, labels, colors, fontsize, *, loc=None, ncol=None, y=0.995):
        """
        Compact, inside-panel legend (1 row, NOT stretched).
        No transforms, no expand.
        """
        handles = [Line2D([0], [0], color=c, lw=4.0)
                   for c in colors[:len(labels)]]
        if ncol is None:
            ncol = len(labels)
    
        return ax.legend(
            handles, labels,
            ncol=ncol,
            loc=loc if loc else "upper center",
            bbox_to_anchor=(0.05, y),     # anchor center-top
            fontsize=int(fontsize * 0.90),
            frameon=False,
            handlelength=2.2,
            handletextpad=0.6,
            columnspacing=1.4,           # tighten spacing
            borderaxespad=0.0,
            labelspacing=0.2,
        )

    def _shade_segments(
        self,
        ax,
        *,
        exp: str,
        segments,
        color: str,
        alpha: float = 0.12,
        hatch: str, 
        zorder: int = 0,
    ):
        """
        Shade year-ranges for one experiment on the current axis.
    
        segments: list of (start_year, end_year) in *model years* (same convention as your plots).
                  Assumed inclusive.
        """
        if not segments:
            return
    
        for seg in segments:
            if not (isinstance(seg, (list, tuple)) and len(seg) == 2):
                raise ValueError(f"Bad segment for {exp}: {seg} (expected (y0, y1))")
            y0, y1 = seg
            y0, y1 = int(y0), int(y1)
            if y1 < y0:
                y0, y1 = y1, y0
    
            if self.annual_mean:
                # x is 0..ny-1 where x=0 corresponds to year_min
                x0 = y0 - self.year_min
                x1 = y1 - self.year_min + 1  # +1 to make end inclusive in axvspan
            else:
                # monthly: x is months since year_min
                x0 = (y0 - self.year_min) * 12
                x1 = (y1 - self.year_min + 1) * 12
    
            #ax.axvspan(x0, x1, facecolor=color, alpha=alpha, lw=0.0, zorder=zorder)
            # rasterized: PDF hatch patterns have a fixed tile size that does not shrink when the
            # figure is scaled in LaTeX (coarse, thick hatching); a bitmap scales with the figure
            ax.axvspan(
                x0, x1,
                facecolor="none",
                edgecolor=color,   # use gray/red/etc
                hatch=hatch,
                lw=0.0,
                alpha=alpha,
                zorder=zorder,
                rasterized=True,
            )

    def plot_on_axis(
        self,
        ax,
        data: Dict[str, Dict[str, np.ndarray]],
        *,
        var_key: str,
        exp_order: List[str],
        exp_labels: Optional[List[str]] = None,
        colors: Optional[List[str]] = None,
        req: Optional[PlotRequest] = None,
        fontsize: int = 18,
        step_yrs: int = 25,
        ref_yrs: Optional[List[int]] = None,
        xlim_years: Optional[Tuple[int, int]] = None,
        labels_zero_based: bool = False,
        show_raw: bool = False,
        show_raw_faint: bool = False,
        std_band: bool = True,
        std_center: str = "filtered",   # filtered|raw
        std_source: str = "raw",        # raw|filtered|residual
        std_mode: str = "global",       # global|rolling
        std_window_months: int = 12,
        std_alpha: float = 0.20,
        legend_loc: Optional[str] = None,  # None => no legend
        exp_shas: Optional[Dict[str, List[Tuple[int,int]]]] = None,
        sha_alpha: float = 0.10, 
    ):
        if req is None:
            req = PlotRequest(key=var_key, vname=var_key, vunit="", scale=1.0)

        exp_labels = exp_labels or exp_order
        colors = colors or ["black", "#cc0000", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]

        ymins, ymaxs = [], []
        used_colors = []
        stats_list = []

        for i, exp in enumerate(exp_order):
            if exp not in data or var_key not in data[exp]:
                raise KeyError(f"Missing data[{exp}][{var_key}]")

            col = colors[i % len(colors)]
            used_colors.append(col)

            raw0 = np.asarray(data[exp][var_key], float).reshape(-1) * float(req.scale)

            # -------- choose sampling & preprocess --------
            if self.annual_mean:
                raw = self.to_annual_mean(raw0)
                samples_per_year = 1
            else:
                raw = raw0
                samples_per_year = 12

            if (not self.annual_mean) and self.running_mean and self.running_mean_months > 1:
                raw = self._nan_running_mean(
                    raw,
                    self.running_mean_months,
                    center=self.running_mean_center,
                    min_valid_frac=self.running_mean_min_frac,
                )

            if self.compute_anomaly:
                raw = self._to_anomaly(raw, samples_per_year=samples_per_year)

            filt = self._lowpass(raw) if self.use_filter else raw.copy()

            st = self._last_n_year_stats(raw, filt)
            stats_list.append(st)

            # --- shading FIRST (behind everything) ---
            sha = None if exp_shas is None else exp_shas.get(exp, None)
            if sha:
                # --- parse flexible schema ---
                if isinstance(sha, dict):
                    segments = sha.get("segments", None)
                    hatch = sha.get("hatch", "///")
                    alpha = sha.get("alpha", sha_alpha)
                    color_use = sha.get("color", col)
                else:
                    # backward compatibility: list of segments
                    segments = sha
                    hatch = "///"
                    alpha = sha_alpha
                    color_use = col
            
                if segments:
                    self._shade_segments(
                        ax,
                        exp=exp,
                        segments=segments,
                        color=color_use,
                        hatch = hatch,
                        alpha=alpha,
                        zorder=-50,
                    )

            # -------- std band (optional) --------
            if std_band and raw.size:
                center = filt if (std_center.lower() == "filtered" and self.use_filter) else raw
                if std_source.lower() == "residual" and self.use_filter:
                    base = raw - filt
                elif std_source.lower() == "filtered" and self.use_filter:
                    base = filt
                else:
                    base = raw

                if std_mode.lower() == "rolling":
                    win = std_window_months if (not self.annual_mean) else max(1, std_window_months // 12)
                    sigma = self._rolling_nanstd(base, win)
                else:
                    s = float(np.nanstd(base, ddof=1)) if base.size > 1 else float(np.nanstd(base))
                    sigma = np.full_like(center, s, dtype=float)

                ax.fill_between(
                    np.arange(center.size),
                    center - sigma,
                    center + sigma,
                    color=col,
                    alpha=std_alpha,
                    linewidth=0,
                    zorder=-5,
                )
                ymins.append(np.nanmin(center - sigma))
                ymaxs.append(np.nanmax(center + sigma))
                
            # -------- plot lines --------
            if self.use_filter:
                if show_raw:
                    ax.plot(raw, color=col, lw=1.0, alpha=0.7, ls="--",zorder=5)
                elif show_raw_faint:
                    ax.plot(raw, color=col, lw=0.5, alpha=0.2, zorder=5)
                ax.plot(filt, color=col, lw=2.0, alpha=1.0, zorder=5)
                ysrc = filt
            else:
                ax.plot(raw, color=col, lw=2.0, alpha=1.0, zorder=5)
                ysrc = raw

            if ysrc.size:
                ymins.append(np.nanmin(ysrc))
                ymaxs.append(np.nanmax(ysrc))

        # x axis
        self._set_time_axis(
            ax,
            step_yrs=step_yrs,
            xlim_years=xlim_years,
            ref_yrs=ref_yrs,
            fontsize=fontsize,
            labels_zero_based=labels_zero_based,
        )

        # title/labels
        title_default = req.title or (req.vname if not self.compute_anomaly else f"{req.vname} Anomaly (Δ{req.vname})")
        ylab_default = req.y_label or (f"{req.vname} ({req.vunit})" if req.vunit else req.vname)

        if req.show_title:
            ax.set_title(title_default, fontsize=fontsize, loc="left", pad=8)

        # x label + tick labels (avoid double-setting xlabel)
        if req.show_xlabel:
            ax.set_xlabel("Model Time (year)", fontsize=int(fontsize * 0.95))
        else:
            ax.set_xlabel("")

        ax.tick_params(axis="x", labelbottom=bool(req.show_xticks))
        if not req.show_xticks:
            ax.tick_params(axis="x", which="both", length=0)
            
        if req.show_ylabel:
            ax.set_ylabel(ylab_default, fontsize=int(fontsize * 0.95))
        else:
            ax.set_ylabel("")

        # y limits
        if req.ylim:
            ax.set_ylim(*req.ylim)
        elif ymins and ymaxs:
            lo, hi = float(np.nanmin(ymins)), float(np.nanmax(ymaxs))
            pad = 0.05 * max(1e-12, hi - lo)
            ax.set_ylim(lo - pad, hi + pad)
            
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

        # zero line & grid
        ylow, yhigh = ax.get_ylim()
        if ylow < 0 < yhigh:
            ax.axhline(0.0, ls="--", lw=1.0, color="k", alpha=0.8)

        # panel label
        if req.panel_label:
            x0, y0 = req.label_loc
            ax.text(x0, y0, req.panel_label, transform=ax.transAxes,
                    va="top", ha="left", fontsize=int(fontsize * 0.9))

        # --- x ticks visibility (safe) ---
        ax.tick_params(axis="x", labelbottom=bool(req.show_xticks))
        if not req.show_xticks:
            ax.tick_params(axis="x", which="both", length=0)
            
        ax.tick_params(labelsize=int(fontsize * 0.9), length=4, width=1)
        
        ax.grid(True, alpha=0.1, color="grey", linewidth=0.01)

        # --- legend behavior (safe) ---
        # Priority:
        #   1) explicit legend_loc argument to plot_on_axis()
        #   2) req.legend_loc
        #   3) None (no legend)
        loc_use = legend_loc if legend_loc is not None else req.legend_loc
        if loc_use:
            self._legend(ax, exp_labels, used_colors, fontsize=fontsize, loc=loc_use, ncol=len(exp_labels))
            
        return stats_list

    def plot_from_extracted(
        self,
        data: Dict[str, Dict[str, np.ndarray]],
        *,
        var_key: str,
        exp_order: List[str],
        exp_labels: Optional[List[str]] = None,
        colors: Optional[List[str]] = None,
        req: Optional[PlotRequest] = None,
        outname: Optional[str] = None,
        figsize: Tuple[int, int] = (24, 6),
        fontsize: int = 28,
        step_yrs: int = 25,
        ref_yrs: Optional[List[int]] = None,
        xlim_years: Optional[Tuple[int, int]] = None,
        labels_zero_based: bool = False,
        show_raw: bool = False,
        show_raw_faint: bool = False,
        std_band: bool = True,
        std_center: str = "filtered",
        std_source: str = "raw",
        std_mode: str = "global",
        std_window_months: int = 12,
        std_alpha: float = 0.20,
        legend_loc: str = "upper left",
        exp_shas: Optional[Dict[str, List[Tuple[int,int]]]] = None,
        sha_alpha: float = 0.10, 
    ):
        if req is None:
            req = PlotRequest(key=var_key, vname=var_key, vunit="", scale=1.0)

        fig, ax = plt.subplots(figsize=figsize)

        stats_list = self.plot_on_axis(
            ax=ax,
            data=data,
            var_key=var_key,
            exp_order=exp_order,
            exp_labels=exp_labels,
            colors=colors,
            req=req,
            fontsize=fontsize,
            step_yrs=step_yrs,
            ref_yrs=ref_yrs,
            xlim_years=xlim_years,
            labels_zero_based=labels_zero_based,
            show_raw=show_raw,
            show_raw_faint=show_raw_faint,
            std_band=std_band,
            std_center=std_center,
            std_source=std_source,
            std_mode=std_mode,
            std_window_months=std_window_months,
            std_alpha=std_alpha,
            legend_loc=legend_loc,
            exp_shas=exp_shas,
            sha_alpha=sha_alpha
        )
        
        fig.set_constrained_layout(False)   # disable layout engine
        
        fig.subplots_adjust(
            wspace=0.18,   # horizontal space between panels
            hspace=0.25    # vertical space between panels
        )
        if outname:
            fig.savefig(outname, dpi=600, bbox_inches="tight")

        return fig, ax, stats_list

    def _infer_time_lev_dims(self, obj, time_dim=None, lev_dim=None):
        """
        Infer names of time and vertical dimensions from an xarray DataArray/Dataset variable.
        Works for:
          - time dims: time, Time, xtime, etc.
          - vertical dims: depth, z, nVertLevels, etc.
        """
        import numpy as np
    
        dims = tuple(obj.dims)
    
        # --- time dim ---
        if time_dim is None:
            # Look for datetime-like coordinate
            for d in dims:
                coord = obj.coords.get(d, None)
                if coord is not None:
                    try:
                        if np.issubdtype(coord.dtype, np.datetime64) or hasattr(coord, "dt"):
                            time_dim = d
                            break
                    except Exception:
                        pass
    
            # Common fallbacks
            if time_dim is None:
                for cand in ("time", "Time", "xtime", "t"):
                    if cand in dims:
                        time_dim = cand
                        break
    
            # Last resort: name contains 'time'
            if time_dim is None:
                for d in dims:
                    if "time" in d.lower():
                        time_dim = d
                        break
    
        if time_dim is None:
            raise KeyError(f"Could not infer time dimension from dims={dims}")
    
        # --- vertical dim ---
        if lev_dim is None:
            for cand in ("depth", "z", "nVertLevels", "lev", "level", "layer", "nz"):
                if cand in dims:
                    lev_dim = cand
                    break
    
            # Last resort: any other non-time dim with size > 1
            if lev_dim is None:
                for d in dims:
                    if d != time_dim and obj.sizes.get(d, 0) > 1:
                        lev_dim = d
                        break
    
        if lev_dim is None:
            raise KeyError(f"Could not infer vertical dimension from dims={dims}")
    
        return time_dim, lev_dim

    def _to_annual_mean_hov(self, da: xr.DataArray, time_dim: str) -> xr.DataArray:
        """
        Convert monthly (time, depth) to annual mean.
        Drops remainder months if nt not divisible by 12.
        Output has dim 'year' and coordinate year = [year_min, year_min+1, ...].
        """
        nt = int(da.sizes[time_dim])
        ny = nt // 12
        if ny <= 0:
            raise ValueError(f"Not enough samples for annual mean: nt={nt}")
    
        da0 = da.isel({time_dim: slice(0, ny * 12)})
        ann = da0.coarsen({time_dim: 12}, boundary="trim").mean()
        ann = ann.rename({time_dim: "year"})
        ann = ann.assign_coords(year=np.arange(self.year_min, self.year_min + ny))
        return ann

    def plot_ohc_hov_ts_figure(
        self,
        *,
        plotter,                       # OHCDiagnosticsPlotter (TS) instance (often `self`)
        hov_data,                      # dict: hov_data[exp] -> xr.DataArray or xr.Dataset
        ts_plot_data,                  # dict: ts_plot_data[exp][var_key] -> 1D numpy
        exp_order,                     # e.g. ["Full-CPL", "AltBG1", "AltBG2"]
        exp_labels,                    # e.g. ["Full-CPL", "FC-FOSI-S1", "FC-FOSI-S2"]
        colors,                        # e.g. ["black", "#cc0000", "#377eb8"]
        ts_var_keys,                   # ["OHC_0_bottom", "OHC_0_700", ...]
        ts_plot_info,                  # your plot_info dict
        exp_shas=None,                 # optional shading segments (TS + HOV, if desired)
        sha_alpha=0.05,

        # hov settings
        hov_key="OHC_layer",           # if hov_data[exp] is Dataset, use this key
        cmap="RdBu_r",
        vmin=None,
        vmax=None,                     # if None -> symmetric based on data
        add_colorbar=True,

        # layout
        ref_yrs=None,                  # e.g. [221] for TS panels
        step_yrs=25,
        fontsize=18,
        figsize=(18, 10),
        outpath=None,

        # NEW: stats output
        out_stats_json: Optional[str] = None,      # default: derived from outpath if provided
        stats_meta: Optional[Dict[str, Any]] = None,

        # optional: vertical lines on hov panels (years)
        hov_vlines=None,               # dict exp -> list of years, or None
        hov_spin_year=None,            # e.g. 221 to match your TS

        ts_hspace=0.35,
        ts_wspace=0.20,
        hov_hspace=0.20,
        hov_wspace=0.25,

        hov_y_label=None,
        hov_y_label_unit=None,
        use_panel_labels=False,

        std_band=True,
        std_center="filtered",
        std_source="raw",
        std_mode="global",
    ):
        # -------- default stats path --------
        if out_stats_json is None and outpath:
            base = os.path.splitext(str(outpath))[0]
            out_stats_json = f"{base}_stats.json"

        # -----------------------
        # 1) Collect hov arrays
        # -----------------------
        hov_list = []
        for exp in exp_order:
            obj = hov_data[exp]
            da = obj[hov_key] if hasattr(obj, "data_vars") else obj
            hov_list.append(da)

        # -----------------------
        # 2) Symmetric color scale
        # -----------------------
        if vmin is None or vmax is None:
            mx = 0.0
            for da in hov_list:
                a = np.asarray(da.values, float)
                mx = max(mx, float(np.nanmax(np.abs(a))) if a.size else 0.0)
            mx = max(mx, 1e-12)
            if vmin is None:
                vmin = -mx
            if vmax is None:
                vmax = mx

        # -----------------------
        # 3) Layout: master grid (Option A)
        # -----------------------
        n_exp = len(exp_order)
        n_ts  = len(ts_var_keys)

        ts_ncols = 2
        ts_nrows = int(np.ceil(n_ts / ts_ncols))

        fig = plt.figure(figsize=figsize, constrained_layout=False)

        # master grid: n_exp hov panels + 1 colorbar column
        gs = fig.add_gridspec(
            nrows=1 + 1 + ts_nrows + 1,   # hov + top pad + TS + bottom pad
            ncols=n_exp + 1,
            height_ratios=[
                1.05,        # Hov row
                0.45,        # top pad
                *([1.0] * ts_nrows),
                0.35,        # bottom pad
            ],
            width_ratios=[1.0] * n_exp + [0.05],
            hspace=0.0,
            wspace=hov_wspace,
        )

        # TS block spans FULL width (INCLUDING cbar column)
        # uses symmetric small spacer columns so TS panels stay equal width
        ts_gs = gs[2 : 2 + ts_nrows, :].subgridspec(
            nrows=ts_nrows,
            ncols=ts_ncols + 2,
            height_ratios=[1.05] + [1.0] * (ts_nrows - 1),
            width_ratios=[0.025] + [1.0] * ts_ncols + [0.025],
            hspace=ts_hspace,
            wspace=ts_wspace,
        )

        # -----------------------
        # 4) Plot Hov panels
        # -----------------------
        def _edges_from_centers(c: np.ndarray) -> np.ndarray:
            c = np.asarray(c, float)
            if c.size == 1:
                return np.array([c[0] - 0.5, c[0] + 0.5], dtype=float)
            e = np.empty(c.size + 1, dtype=float)
            e[1:-1] = 0.5 * (c[1:] + c[:-1])
            e[0] = c[0] - (e[1] - c[0])
            e[-1] = c[-1] + (c[-1] - e[-2])
            return e

        im = None
        for j, (exp, lab, da) in enumerate(zip(exp_order, exp_labels, hov_list)):
            ax = fig.add_subplot(gs[0, j])

            time_dim, z_dim = self._infer_time_lev_dims(da)

            # --- choose annual vs monthly for Hov ---
            if self.annual_mean:
                da2 = self._to_annual_mean_hov(da, time_dim)
                x_dim = "year"
            else:
                da2 = da
                x_dim = time_dim

            da2 = da2.transpose(x_dim, z_dim)

            z = np.asarray(da2[z_dim].values, float)

            depth_units = str((getattr(da2[z_dim], "attrs", {}) or {}).get("units", "")).strip().lower()
            is_meters = depth_units.startswith("m")  # "m", "meter", "meters", ...
            z_plot = z / 1000.0 if is_meters else z.copy()
            z_plot_unit = "km" if is_meters else (depth_units if depth_units else "depth")

            if z_plot.size >= 2 and np.nanmedian(np.diff(z_plot)) < 0:
                z_plot = z_plot[::-1]
                da2 = da2.isel({z_dim: slice(None, None, -1)})

            if z_plot.size and np.nanmean(z_plot) < 0:
                z_plot = np.abs(z_plot)

            nx = int(da2.sizes[x_dim])
            x_edges = np.arange(nx + 1, dtype=float)
            y_edges = _edges_from_centers(z_plot)

            Z = np.asarray(da2.values, float).T  # (z, x)

            x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
            y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

            if Z.shape != (y_centers.size, x_centers.size):
                raise ValueError(
                    f"[hov] shape mismatch exp={exp}: "
                    f"Z{Z.shape} vs centers {(y_centers.size, x_centers.size)} "
                    f"(dims x={x_dim}, z={z_dim})"
                )

            levels = np.linspace(vmin, vmax, 41)
            im = ax.contourf(
                x_centers, y_centers, Z,
                levels=levels,
                cmap=cmap,
                extend="both",
            )

            ax.set_title(f"({chr(ord('a') + j)}) {lab}", loc="left", fontsize=int(fontsize), pad=6)
            ax.set_title(r"$\Delta$OHC (layer)", loc="right", fontsize=int(fontsize), pad=6)
            ax.set_xlabel("Model Time (year)", fontsize=int(fontsize * 0.95))

            if j == 0:
                ax.set_ylabel(f"Depth ({z_plot_unit})", fontsize=int(fontsize * 0.95))
            else:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)

            # --- hatched shading for exp_shas windows on HOV (optional) ---
            sha = None if exp_shas is None else exp_shas.get(exp, None)
            if sha:
                if isinstance(sha, dict):
                    segments  = sha.get("segments", None)
                    # "hov_alpha"/"hov_color" override "alpha"/"color" on the Hovmöller only
                    # (e.g. stronger hatching so it stays visible over the red/blue colormap)
                    alpha     = sha.get("hov_alpha", sha.get("alpha", sha_alpha))
                    hatch     = sha.get("hatch", "///")
                    color_use = sha.get("hov_color", sha.get("color", "red"))
                else:
                    segments  = sha
                    alpha     = sha_alpha
                    hatch     = "///"
                    color_use = "red"

                if segments:
                    for (y0, y1) in segments:
                        y0, y1 = sorted((int(y0), int(y1)))

                        if self.annual_mean:
                            x0 = y0 - self.year_min
                            x1 = (y1 - self.year_min) + 1
                        else:
                            x0 = (y0 - self.year_min) * 12
                            x1 = ((y1 - self.year_min) + 1) * 12

                        x0c = max(0, x0)
                        x1c = min(nx, x1)

                        if x1c > x0c:
                            # rasterized: PDF hatch patterns have a fixed tile size that does not shrink when the
                            # figure is scaled in LaTeX (coarse, thick hatching); a bitmap scales with the figure
                            ax.axvspan(
                                x0c, x1c,
                                facecolor="none",
                                edgecolor=color_use,
                                hatch=hatch,
                                linewidth=0.0,
                                alpha=alpha,
                                zorder=5,
                                rasterized=True,
                            )

            ax.invert_yaxis()
            ax.tick_params(labelsize=int(fontsize * 0.90))

            # -----------------------
            # x ticks + vertical lines
            # -----------------------
            if self.annual_mean:
                step = max(1, int(step_yrs * 2))
                xt = np.arange(0, nx + 1, step)
                ax.set_xticks(xt)
                ax.set_xticklabels([str(self.year_min + int(k)) for k in xt])

                if hov_vlines is not None:
                    for y in hov_vlines.get(exp, []):
                        idx = int(y) - self.year_min
                        if 0 <= idx < nx:
                            ax.axvline(idx, color="r", lw=1.0, ls="--", alpha=0.8)

                if hov_spin_year is not None:
                    idx = int(hov_spin_year) - self.year_min
                    if 0 <= idx < nx:
                        ax.axvline(idx, color="k", lw=1.0, ls="--", alpha=0.9)

            else:
                step = max(1, int(step_yrs * 2))
                xt = np.arange(0, nx, 12 * step)
                ax.set_xticks(xt)
                ax.set_xticklabels([str(self.year_min + int(k // 12)) for k in xt])

                if hov_vlines is not None:
                    for y in hov_vlines.get(exp, []):
                        idx = (int(y) - self.year_min) * 12
                        if 0 <= idx < nx:
                            ax.axvline(idx, color="r", lw=1.0, ls="--", alpha=0.8)

                if hov_spin_year is not None:
                    idx = (int(hov_spin_year) - self.year_min) * 12
                    if 0 <= idx < nx:
                        ax.axvline(idx, color="k", lw=1.0, ls="--", alpha=0.9)

        # -----------------------
        # Colorbar
        # -----------------------
        if add_colorbar and im is not None:
            cax = fig.add_subplot(gs[0, -1])
            cb = fig.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=int(fontsize * 0.90))
            units = hov_list[0].attrs.get("units", "")
            cb.set_label(r"$\Delta$OHC per layer" + (f" [{units}]" if units else ""),
                         fontsize=int(fontsize * 0.95))

        # -----------------------
        # 5) Plot TS panels using ts_gs
        # -----------------------
        def _ts_spec(rr, cc):
            return ts_gs[rr, 0:2] if cc == 0 else ts_gs[rr, 2:4]

        # NEW: collect stats across vars/experiments
        stats_by_var: Dict[str, Dict[str, Any]] = {}

        for i, var_key in enumerate(ts_var_keys):
            rr = i // ts_ncols
            cc = i % ts_ncols

            ax = fig.add_subplot(_ts_spec(rr, cc))

            vinfo = ts_plot_info[var_key]
            is_last_row = (rr == ts_nrows - 1)
            show_xticks = is_last_row

            panel_label = f"({chr(ord('a') + n_exp + i)})" if use_panel_labels else None

            req = PlotRequest(
                key=var_key,
                vname=vinfo["short"],
                vunit=vinfo["unit"],
                scale=vinfo.get("scale", 1.0),
                ylim=vinfo.get("ylim", None),
                title=vinfo.get("label", None),
                y_label=fr"{vinfo['short']} ({vinfo['unit']})",
                panel_label=panel_label,
                show_xlabel=show_xticks,
                show_xticks=show_xticks,
                show_ylabel=(cc == 0),
                show_title=True,
                legend_loc="upper left",
            )

            stats_list = plotter.plot_on_axis(
                ax=ax,
                data=ts_plot_data,
                var_key=var_key,
                exp_order=exp_order,
                exp_labels=exp_labels,
                colors=colors,
                req=req,
                fontsize=fontsize,
                step_yrs=step_yrs,
                ref_yrs=ref_yrs,
                std_band=std_band,
                std_center=std_center,
                std_source=std_source,
                std_mode=std_mode,
                exp_shas=exp_shas,
                sha_alpha=sha_alpha,
            )

            # stats_list aligns with exp_order
            stats_by_var[var_key] = {exp: stats_list[ii] for ii, exp in enumerate(exp_order)}

        # hide unused TS panels
        nslots = ts_nrows * ts_ncols
        for k in range(n_ts, nslots):
            rr = k // ts_ncols
            cc = k % ts_ncols
            fig.add_subplot(_ts_spec(rr, cc)).axis("off")

        # save figure
        if outpath:
            os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
            fig.savefig(outpath, dpi=600, bbox_inches="tight")

        # write stats JSON
        if out_stats_json:
            meta = {} if stats_meta is None else dict(stats_meta)
            meta.setdefault("exp_order", list(exp_order))
            meta.setdefault("exp_labels", list(exp_labels))
            meta.setdefault("ts_var_keys", list(ts_var_keys))
            meta.setdefault("annual_mean", bool(plotter.annual_mean))
            meta.setdefault("compute_anomaly", bool(plotter.compute_anomaly))
            meta.setdefault("trend_years", int(plotter.trend_years))
            meta.setdefault("trend_unit", str(plotter.trend_unit))
            meta.setdefault("ci_method", str(plotter.ci_method))
            meta.setdefault("add_trend_ci", bool(plotter.add_trend_ci))
            meta.setdefault("year_min", int(plotter.year_min))
            meta.setdefault("year_max", int(plotter.year_max))

            self.save_figure_stats_json(
                out_stats_json=out_stats_json,
                stats_by_var=stats_by_var,
                meta=meta,
                verbose=True,
            )

        return fig, stats_by_var


@dataclass
class PlotRequest:
    """
    Plot settings for one variable (e.g., SST).
    """
    key: str
    vname: str
    vunit: str
    scale: float = 1.0
    ylim: Optional[Tuple[float, float]] = None
    title: Optional[str] = None
    y_label: Optional[str] = None

    # multipanel-friendly controls
    panel_label: Optional[str] = None
    show_xlabel: bool = True
    show_ylabel: bool = True
    show_xticks: bool = True
    show_title: bool = True
    legend_loc: Optional[str] = None
    label_loc: Tuple[float, float] = (0.01, 0.98)


@dataclass
class OHCFigureDataCollectorConfig:
    """
    Defaults that are common across OHC series + Hov extractions.
    You can override per-call, but this mirrors the FigureDataCollector idea
    where config lives on the collector.
    """
    # IO
    file_key: str = "mpasTimeSeriesOcean"
    engine: str = "netcdf4"
    chunks: Optional[dict] = None
    verbose: bool = False

    # time handling (passed into registry.build_extract_config)
    calendar: str = "noleap"
    reassign_monthly_time: bool = True
    length_policy: str = "auto"

    # annualization defaults (series + hov)
    annual_mean: bool = False
    annual_weighted: bool = True
    annual_drop_incomplete: bool = True
    annual_min_days: int = 360


class OHCFigureDataCollector:
    """
    FigureDataCollector-style wrapper for OHC:
      - Collect 1D OHC/OHU series
      - Collect OHC Hovmöller matrices
      - Provide a build_plot_data merge helper
      - Save/load TS + HOV per experiment
      - Save/load a full "bundle" (TS + HOV + meta.json)

    It expects your existing build_exp_subdirs_from_run_dict(...) to define the
    experiment file subdirectories.
    """

    def __init__(
        self,
        *,
        run_dict: Dict[str, Any],
        exp_tags: Sequence[str],
        exp_type: str,
        ts_window: Tuple[int, int],
        climo_len: int,
        ts_base: str,
        data_dir: str,
        mpas_depth: str,  # fdepth path (required for OHC util)
        cfg: Optional[OHCFigureDataCollectorConfig] = None,
    ):
        self.run_dict = run_dict
        self.exp_tags = list(exp_tags)
        self.exp_type = str(exp_type)
        self.ts_window = tuple(ts_window)
        self.climo_len = int(climo_len)
        self.ts_base = str(ts_base)
        self.data_dir = str(data_dir)
        self.mpas_depth = str(mpas_depth)

        self.cfg = cfg or OHCFigureDataCollectorConfig()

        # registry is cheap; keep one
        self.registry = OHCRegistry()

        if not self.mpas_depth or (not os.path.exists(self.mpas_depth)):
            raise FileNotFoundError(f"[OHC] mpas_depth (fdepth) not found: {self.mpas_depth}")

        # build subdirs once
        self.exp_subdirs = build_exp_subdirs_from_run_dict(
            run_dict=self.run_dict,
            exp_tags=self.exp_tags,
            exp_type=self.exp_type,
            ts_window=self.ts_window,
            climo_len=self.climo_len,
            ts_base=self.ts_base,
            mpas_analysis=True,
        )
        self.exps = list(self.exp_subdirs.keys())

        if len(self.exps) != len(self.exp_tags):
            raise ValueError(
                f"Cannot remap experiment keys: len(exps)={len(self.exps)} != len(exp_tags)={len(self.exp_tags)}"
            )

    # -------------------------
    # small utilities
    # -------------------------
    def _remap_exps_to_tags(self, obj_by_exp: Dict[str, Any]) -> Dict[str, Any]:
        return {tag: obj_by_exp[exp] for exp, tag in zip(self.exps, self.exp_tags)}

    def _resolve_region_index(self, *, region_key: str, region_index: Optional[int]) -> int:
        if region_index is None:
            region_index = self.registry.get_region_index(region_key)
        return int(region_index)

    # -------------------------
    # 1D series collection
    # -------------------------
    def collect_ohc_series(
        self,
        *,
        need_vars: Sequence[str] = ("OHC_0_700",),
        region_key: str = "Global",
        reg_dim: str = "nOceanRegions",
        region_index: Optional[int] = None,
        # anomaly/unit controls
        l_anom: bool = True,
        to_ZJ: bool = True,
        rho: float = 1026.0,
        cp: float = 3996.0,
        baseline_months: int = 12,
        # geometry conventions
        use_mask_in_volume: bool = True,
        # annualization overrides (defaults from cfg if None)
        annual_mean: Optional[bool] = None,
        annual_weighted: Optional[bool] = None,
        annual_drop_incomplete: Optional[bool] = None,
        annual_min_days: Optional[int] = None,
        # symmetry args (unused but accepted)
        year_token: Optional[str] = None,
        time_unit: Optional[str] = None,
        use_mfdataset: bool = True,
        verbose: Optional[bool] = None,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Returns:
          data[tag][var] -> 1D numpy array (monthly or annual)
        """
        v = self.cfg.verbose if verbose is None else bool(verbose)

        ridx = self._resolve_region_index(region_key=region_key, region_index=region_index)

        a_mean = self.cfg.annual_mean if annual_mean is None else bool(annual_mean)
        a_wgt  = self.cfg.annual_weighted if annual_weighted is None else bool(annual_weighted)
        a_drop = self.cfg.annual_drop_incomplete if annual_drop_incomplete is None else bool(annual_drop_incomplete)
        a_min  = self.cfg.annual_min_days if annual_min_days is None else int(annual_min_days)

        cfg = self.registry.build_extract_config(
            ts_window=self.ts_window,
            annual_mean=a_mean,
            annual_weighted=a_wgt,
            annual_drop_incomplete=a_drop,
            annual_min_days=a_min,
            calendar=self.cfg.calendar,
            reassign_monthly_time=self.cfg.reassign_monthly_time,
            length_policy=self.cfg.length_policy,
        )

        var_dict = self.registry.get_var_dict()
        var_dict_small = {k: var_dict[k] for k in need_vars if k in var_dict}
        missing = [k for k in need_vars if k not in var_dict_small]
        if missing:
            raise KeyError(
                f"Requested OHC vars not found in registry: {missing}\n"
                f"Known: {list(var_dict.keys())}"
            )

        if v:
            print(f"[OHC] region_key='{region_key}' -> region_index={ridx}")
            print(f"[OHC] fdepth='{self.mpas_depth}' file_key='{self.cfg.file_key}' annual_mean={a_mean}")
            print(f"[OHC] use_mask_in_volume={use_mask_in_volume}")

        data_by_exp = extract_all_ohc(
            root_dir=self.data_dir,
            exps=self.exps,
            exp_subdirs=self.exp_subdirs,
            var_dict=var_dict_small,
            cfg=cfg,
            region_index=ridx,
            region_dim=str(reg_dim),
            anomaly=bool(l_anom),
            baseline_months=int(baseline_months),
            to_ZJ=bool(to_ZJ),
            rho=float(rho),
            cp=float(cp),
            use_mask_in_volume=bool(use_mask_in_volume),
            file_key=self.cfg.file_key,
            fdepth=self.mpas_depth,
            engine=self.cfg.engine,
            chunks=self.cfg.chunks,
            verbose=v,
        )

        data = self._remap_exps_to_tags(data_by_exp)

        if v:
            print("\n[Extracted OHC series summary]")
            print(f"  region_key='{region_key}' -> region_index={ridx}")
            print(f"  fdepth='{self.mpas_depth}'")
            print(f"  use_mask_in_volume={use_mask_in_volume}")
            for tag in self.exp_tags:
                print(f"  {tag}:")
                for k in var_dict_small.keys():
                    arr = data[tag][k]
                    print(f"    {k:16s}  n={arr.size:6d}  min={np.nanmin(arr):.4g}  max={np.nanmax(arr):.4g}")

        return data

    # -------------------------
    # Hovmöller collection
    # -------------------------
    def collect_ohc_hov(
        self,
        *,
        region_key: str = "Global",
        reg_dim: str = "nOceanRegionsTmp",
        region_index: Optional[int] = None,
        # annualization overrides (defaults from cfg if None)
        annual_mean: Optional[bool] = True,
        annual_weighted: Optional[bool] = None,
        annual_drop_incomplete: Optional[bool] = None,
        annual_min_days: Optional[int] = None,
        # hov request knobs
        mode: str = "total",
        to_units: str = "1e22 J",
        anomaly: bool = True,
        baseline_mode: str = "calendar_year",
        first_N: int = 12,
        depth_coord: str = "bottom",
        return_cumulative: bool = False,
        from_top: bool = True,
        # conventions / closure controls
        use_mask_in_area: bool = False,
        normalize_by_surface_area: bool = False,
        window: Optional[Tuple[float, float]] = None,
        apply_window_overlap: bool = False,
        # physics
        rho: float = 1026.0,
        cp: float = 3996.0,
        verbose: Optional[bool] = None,
    ) -> Dict[str, Union[xr.DataArray, xr.Dataset]]:
        """
        Returns:
          hov[tag] -> DataArray or Dataset

        If return_cumulative=False:
          Dataset with OHC_layer(time, depth)

        If return_cumulative=True:
          Dataset with OHC_layer and OHC_cumulative
        """
        v = self.cfg.verbose if verbose is None else bool(verbose)

        ridx = self._resolve_region_index(region_key=region_key, region_index=region_index)

        a_mean = (True if annual_mean is None else bool(annual_mean))
        a_wgt  = self.cfg.annual_weighted if annual_weighted is None else bool(annual_weighted)
        a_drop = self.cfg.annual_drop_incomplete if annual_drop_incomplete is None else bool(annual_drop_incomplete)
        a_min  = self.cfg.annual_min_days if annual_min_days is None else int(annual_min_days)

        cfg = self.registry.build_extract_config(
            ts_window=self.ts_window,
            annual_mean=a_mean,
            annual_weighted=a_wgt,
            annual_drop_incomplete=a_drop,
            annual_min_days=a_min,
            calendar=self.cfg.calendar,
            reassign_monthly_time=self.cfg.reassign_monthly_time,
            length_policy=self.cfg.length_policy,
        )

        req = OHCHovmollerRequest(
            mode=str(mode),
            anomaly=bool(anomaly),
            baseline_mode=str(baseline_mode),
            first_N=int(first_N),
            to_units=str(to_units),
            rho=float(rho),
            cp=float(cp),
            region_index=int(ridx),
            region_dim=str(reg_dim),
            depth_coord=str(depth_coord),
            return_cumulative=bool(return_cumulative),
            from_top=bool(from_top),
            normalize_by_surface_area=bool(normalize_by_surface_area),
            window=window,
            apply_window_overlap=bool(apply_window_overlap),
            use_mask_in_area=bool(use_mask_in_area),
        )

        hov_by_exp = extract_all_ohc_hovmoller(
            root_dir=self.data_dir,
            exps=self.exps,
            exp_subdirs=self.exp_subdirs,
            cfg=cfg,
            req=req,
            file_key=self.cfg.file_key,
            fdepth=self.mpas_depth,
            engine=self.cfg.engine,
            chunks=self.cfg.chunks,
            verbose=v,
        )

        hov = self._remap_exps_to_tags(hov_by_exp)

        if v:
            print("\n[HOV extracted]")
            for tag in self.exp_tags:
                obj = hov[tag]
                if isinstance(obj, xr.Dataset):
                    print(f"  {tag}: Dataset vars={list(obj.data_vars)}")
                    da0 = obj["OHC_layer"]
                else:
                    print(f"  {tag}: DataArray name={obj.name}")
                    da0 = obj
                print(f"    dims={da0.dims} shape={tuple(da0.shape)} units={da0.attrs.get('units','')}")
                if "depth" in da0.dims:
                    dmin = float(np.nanmin(da0["depth"].values))
                    dmax = float(np.nanmax(da0["depth"].values))
                    print(f"    depth: [{dmin:.3g}, {dmax:.3g}] {da0['depth'].attrs.get('units','')}")

        return hov

    # -------------------------
    # plot data merge helper
    # -------------------------
    @staticmethod
    def build_plot_data(
        *sources: Dict[str, Dict[str, np.ndarray]],
        plot_vars: Sequence[str],
        exp_tags: Sequence[str],
    ) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Merge multiple extracted-data dictionaries into one plot_data dict.

        Each source must look like:
          source[exp][var] -> 1D array

        Later sources are only used if var not found in earlier sources.
        """
        plot_data: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in exp_tags:
            plot_data[exp] = {}
            for v in plot_vars:
                for src in sources:
                    if not src:
                        continue
                    if v in src.get(exp, {}):
                        plot_data[exp][v] = src[exp][v]
                        break
                else:
                    raise KeyError(f"{v} not found for experiment '{exp}' in any provided source")
        return plot_data

    # ==========================================================
    # NetCDF attribute sanitizer (self-contained)
    # ==========================================================
    @staticmethod
    def _nc_sanitize_attrs(attrs: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        NetCDF-safe attrs:
          - None: skip
          - bool/np.bool_: int (0/1)
          - numpy scalar: python scalar
          - list/tuple/dict: JSON string
          - bytes: decode
          - otherwise: str(...)
        """
        if not attrs:
            return {}

        out: Dict[str, Any] = {}
        for k, v in attrs.items():
            if v is None:
                continue
            if isinstance(v, np.generic):
                v = v.item()
            if isinstance(v, (bool, np.bool_)):
                v = int(v)
            if isinstance(v, (list, tuple, dict)):
                v = json.dumps(v)
            if isinstance(v, (bytes, bytearray)):
                v = v.decode("utf-8", errors="replace")
            if isinstance(v, (str, int, float)):
                out[str(k)] = v
            else:
                out[str(k)] = str(v)
        return out

    # ==========================================================
    # TS save/load (one file per experiment)
    # ==========================================================
    def save_per_exp(
        self,
        *,
        out_nc_base: str,
        plot_data: Dict[str, Dict[str, np.ndarray]],
        plot_info: Dict[str, Dict[str, Any]],
        freq: Literal["monthly", "annual"],
        year_min_for_coords: int,
        year_max_for_meta: int,
        apply_scale: bool = True,
        attrs: Optional[Dict[str, Any]] = None,
        compress: bool = True,
    ) -> None:
        base = os.path.splitext(out_nc_base)[0]
        attrs_in = {} if attrs is None else dict(attrs)
        ts_var_keys = list(plot_info.keys())

        for exp in self.exp_tags:
            if exp not in plot_data:
                raise KeyError(f"[save_ts] missing exp in plot_data: {exp}")

            v0 = ts_var_keys[0]
            if v0 not in plot_data[exp]:
                raise KeyError(f"[save_ts] missing var '{v0}' in plot_data[{exp}]")
            arr0 = np.asarray(plot_data[exp][v0], float).reshape(-1)
            n = arr0.size

            if freq == "annual":
                dim = "year"
                year = np.arange(int(year_min_for_coords), int(year_min_for_coords) + n, dtype=int)
                coords = {"year": year, "time_year": (("year",), year.astype(float))}
            elif freq == "monthly":
                dim = "month"
                month = np.arange(n, dtype=int)
                time_year = float(year_min_for_coords) + month / 12.0
                coords = {"month": month, "time_year": (("month",), time_year)}
            else:
                raise ValueError(f"freq must be 'monthly' or 'annual', got {freq}")

            data_vars: Dict[str, xr.DataArray] = {}
            for var in ts_var_keys:
                if var not in plot_data[exp]:
                    raise KeyError(f"[save_ts] missing var '{var}' in plot_data[{exp}]")

                info = plot_info.get(var, {}) or {}
                units = info.get("unit", "")
                long_name = info.get("label", var)
                scale = float(info.get("scale", 1.0)) if apply_scale else 1.0

                arr = np.asarray(plot_data[exp][var], float).reshape(-1)
                if arr.size != n:
                    raise ValueError(f"[save_ts] length mismatch {exp}:{var} ({arr.size} vs {n})")

                data_vars[var] = xr.DataArray(
                    arr * scale,
                    dims=(dim,),
                    coords=coords,
                    attrs={"units": units, "long_name": long_name, "applied_scale": scale},
                )

            ds = xr.Dataset(data_vars)
            base_attrs = {
                "experiment": exp,
                "year_min": int(year_min_for_coords),
                "year_max_meta": int(year_max_for_meta),
                "freq": str(freq),
                "source_kind": "time_series",
            }
            ds.attrs.update(self._nc_sanitize_attrs(base_attrs))
            ds.attrs.update(self._nc_sanitize_attrs(attrs_in))

            fn = f"{base}_ts_{exp}.nc"
            encoding = None
            if compress:
                encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}

            ds.to_netcdf(fn, encoding=encoding)
            if self.cfg.verbose:
                print(f"[saved TS] {fn}")

    def load_per_exp(
        self,
        *,
        out_nc_base: str,
        ts_var_keys: Sequence[str],
        prefer_float64: bool = True,
        exp_order: Optional[Sequence[str]] = None,
    ) -> Dict[str, Dict[str, np.ndarray]]:
        base = os.path.splitext(out_nc_base)[0]
        use_exps = list(exp_order) if exp_order is not None else list(self.exp_tags)

        plot_data: Dict[str, Dict[str, np.ndarray]] = {}
        for exp in use_exps:
            fn = f"{base}_ts_{exp}.nc"
            if not os.path.exists(fn):
                raise FileNotFoundError(f"[load_ts] missing file: {fn}")

            with xr.open_dataset(fn) as ds:
                plot_data[exp] = {}
                for var in ts_var_keys:
                    if var not in ds.data_vars:
                        raise KeyError(f"[load_ts] {fn} missing variable '{var}'")
                    vals = ds[var].values
                    arr = np.asarray(vals, dtype=float if prefer_float64 else vals.dtype).reshape(-1)
                    plot_data[exp][var] = arr

        if self.cfg.verbose:
            print(f"\n[OHCDataCollector] loaded TS from: {base}_ts_<exp>.nc")
            for exp in use_exps:
                print(f"  {exp}: {list(plot_data[exp].keys())}")

        return plot_data

    # ==========================================================
    # HOV save/load (one file per experiment)
    # ==========================================================
    def save_hov_per_exp(
        self,
        *,
        out_nc_base: str,
        hov: Dict[str, Union[xr.DataArray, xr.Dataset]],
        hov_key: str = "OHC_layer",
        attrs: Optional[Dict[str, Any]] = None,
        compress: bool = True,
    ) -> None:
        """
        Writes:
          {base}_hov_<exp>.nc

        Accepts hov[tag] as xr.Dataset or xr.DataArray.
        Always saves as Dataset (DataArray -> .to_dataset()).
        """
        base = os.path.splitext(out_nc_base)[0]
        attrs_in = {} if attrs is None else dict(attrs)

        for exp in self.exp_tags:
            if exp not in hov:
                raise KeyError(f"[save_hov] missing exp in hov: {exp}")

            obj = hov[exp]
            if isinstance(obj, xr.DataArray):
                ds = obj.to_dataset(name=obj.name or hov_key)
            elif isinstance(obj, xr.Dataset):
                ds = obj
            else:
                raise TypeError(f"[save_hov] hov[{exp}] must be DataArray or Dataset, got {type(obj)}")

            # ensure expected key exists if possible
            if hov_key not in ds.data_vars and len(ds.data_vars) == 1:
                only = list(ds.data_vars)[0]
                ds = ds.rename({only: hov_key})

            ds = ds.copy()
            base_attrs = {"experiment": exp, "source_kind": "hovmoller"}
            ds.attrs.update(self._nc_sanitize_attrs(base_attrs))
            ds.attrs.update(self._nc_sanitize_attrs(attrs_in))

            fn = f"{base}_hov_{exp}.nc"
            encoding = None
            if compress:
                encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
            ds.to_netcdf(fn, encoding=encoding)

            if self.cfg.verbose:
                print(f"[saved HOV] {fn}")

    def load_hov_per_exp(
        self,
        *,
        out_nc_base: str,
        exp_order: Optional[Sequence[str]] = None,
        hov_key: str = "OHC_layer",
        as_dataset: bool = True,
    ) -> Dict[str, Union[xr.Dataset, xr.DataArray]]:
        """
        Reads:
          {base}_hov_<exp>.nc

        If as_dataset=True (default): returns Dataset per exp.
        If as_dataset=False: returns DataArray if dataset contains only hov_key, else Dataset.
        """
        base = os.path.splitext(out_nc_base)[0]
        use_exps = list(exp_order) if exp_order is not None else list(self.exp_tags)

        hov: Dict[str, Union[xr.Dataset, xr.DataArray]] = {}
        for exp in use_exps:
            fn = f"{base}_hov_{exp}.nc"
            if not os.path.exists(fn):
                raise FileNotFoundError(f"[load_hov] missing file: {fn}")

            ds = xr.open_dataset(fn)
            if hov_key not in ds.data_vars and len(ds.data_vars) == 1:
                only = list(ds.data_vars)[0]
                ds = ds.rename({only: hov_key})

            if as_dataset:
                hov[exp] = ds
            else:
                if hov_key in ds.data_vars and len(ds.data_vars) == 1:
                    hov[exp] = ds[hov_key]
                    ds.close()
                else:
                    hov[exp] = ds

        if self.cfg.verbose:
            print(f"\n[OHCDataCollector] loaded HOV from: {base}_hov_<exp>.nc")

        return hov

    # ==========================================================
    # Bundle save/load (TS + HOV + meta.json)
    # ==========================================================
    def save_bundle(
        self,
        *,
        out_base: str,
        ts_plot_data: Optional[Dict[str, Dict[str, np.ndarray]]] = None,
        ts_plot_info: Optional[Dict[str, Dict[str, Any]]] = None,
        ts_freq: Optional[Literal["monthly", "annual"]] = None,
        ts_year_min_for_coords: Optional[int] = None,
        ts_year_max_for_meta: Optional[int] = None,
        hov: Optional[Dict[str, Union[xr.DataArray, xr.Dataset]]] = None,
        hov_key: str = "OHC_layer",
        meta: Optional[Dict[str, Any]] = None,
        compress: bool = True,
    ) -> None:
        """
        Writes (as available):
          {out_base}_ts_<exp>.nc
          {out_base}_hov_<exp>.nc
          {out_base}_meta.json
        """
        base = os.path.splitext(out_base)[0]
        os.makedirs(os.path.dirname(base) or ".", exist_ok=True)

        if ts_plot_data is not None:
            if ts_plot_info is None or ts_freq is None or ts_year_min_for_coords is None or ts_year_max_for_meta is None:
                raise ValueError(
                    "[save_bundle] TS requested but missing ts_plot_info/ts_freq/ts_year_min_for_coords/ts_year_max_for_meta"
                )

            self.save_per_exp(
                out_nc_base=base,
                plot_data=ts_plot_data,
                plot_info=ts_plot_info,
                freq=ts_freq,
                year_min_for_coords=ts_year_min_for_coords,
                year_max_for_meta=ts_year_max_for_meta,
                apply_scale=True,
                attrs={"bundle": "OHCFigureDataCollector"},
                compress=compress,
            )

        if hov is not None:
            self.save_hov_per_exp(
                out_nc_base=base,
                hov=hov,
                hov_key=hov_key,
                attrs={"bundle": "OHCFigureDataCollector"},
                compress=compress,
            )

        meta_out: Dict[str, Any] = {} if meta is None else dict(meta)
        meta_out.setdefault(
            "collector",
            {
                "class": "OHCFigureDataCollector",
                "cfg": asdict(self.cfg),
                "ts_window": list(self.ts_window),
                "exp_tags": list(self.exp_tags),
                "exp_type": self.exp_type,
                "data_dir": self.data_dir,
                "ts_base": self.ts_base,
                "climo_len": self.climo_len,
                "mpas_depth": self.mpas_depth,
            },
        )
        meta_out.setdefault("hov_key", hov_key)

        # include plot_info so a load can infer ts_var_keys
        if ts_plot_info is not None:
            meta_out.setdefault("ts_plot_info", ts_plot_info)
        if ts_freq is not None:
            meta_out.setdefault("ts_freq", str(ts_freq))
        if ts_year_min_for_coords is not None:
            meta_out.setdefault("ts_year_min_for_coords", int(ts_year_min_for_coords))
        if ts_year_max_for_meta is not None:
            meta_out.setdefault("ts_year_max_for_meta", int(ts_year_max_for_meta))

        with open(f"{base}_meta.json", "w", encoding="utf-8") as f:
            json.dump(meta_out, f, indent=2, sort_keys=True)

        if self.cfg.verbose:
            print(f"[saved META] {base}_meta.json")

    def load_bundle(
        self,
        *,
        out_base: str,
        ts_var_keys: Optional[Sequence[str]] = None,
        load_ts: bool = True,
        load_hov: bool = True,
        hov_as_dataset: bool = True,
        hov_key: str = "OHC_layer",
    ) -> Dict[str, Any]:
        """
        Returns dict with keys:
          - "ts_plot_data" (optional)
          - "hov" (optional)
          - "meta" (optional)
        """
        base = os.path.splitext(out_base)[0]
        out: Dict[str, Any] = {"meta": {}}

        meta_fn = f"{base}_meta.json"
        if os.path.exists(meta_fn):
            with open(meta_fn, "r", encoding="utf-8") as f:
                out["meta"] = json.load(f)

        if load_ts:
            if ts_var_keys is None:
                if "ts_plot_info" in out["meta"]:
                    ts_var_keys = list(out["meta"]["ts_plot_info"].keys())
                else:
                    raise ValueError("[load_bundle] ts_var_keys must be provided (or present in meta['ts_plot_info']).")

            out["ts_plot_data"] = self.load_per_exp(
                out_nc_base=base,
                ts_var_keys=ts_var_keys,
                prefer_float64=True,
            )

        if load_hov:
            hk = out["meta"].get("hov_key", hov_key)
            out["hov"] = self.load_hov_per_exp(
                out_nc_base=base,
                as_dataset=hov_as_dataset,
                hov_key=hk,
            )

        return out
