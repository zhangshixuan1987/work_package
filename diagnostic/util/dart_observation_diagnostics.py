import os
import glob
from typing import Dict, Tuple, List, Optional

import numpy as np
import xarray as xr
import xcdat as xc

class DartObsDiagReader:
    """
    Lightweight reader for DART diagnostics (obs_diag & obs_seq).

    Required config
    ---------------
    {
        "data_path": "/compyfs/.../v3_dart_cda_scratch",
        "path_template": "{data_path}/{run_id}/archive/{group_key}/dart_diagnostics/{diag_set}",
        "file_templates": {
            "obs_diag": "{run_id}.dart.e.eam_{diag_set}_output.{period}.nc",
            "obs_seq":  "{name}.dart.e.eam_{diag_set}_final.{period}-*.nc",
        },
        # optional:
        "region_dict": {...}  # supports [(lat_lo,lat_hi), (lon_lo,lon_hi)] or + (z_lo,z_hi)
    }

    Notes
    -----
    - Resolver tries: (primary) -> (swap {name}<->{run_id}) -> (wildcard {period}).
    - For obs_diag files, use extract_metrics_data(...).
    - For obs_seq files, use make_obs_seq_datadict(...).
    """

    # ------------------------------------------------------------------ #
    # Init / config
    # ------------------------------------------------------------------ #
    def __init__(self, config: dict):
        if "file_templates" in config:
            required = ("data_path", "path_template", "file_templates")
            missing = [k for k in required if k not in config]
            if missing:
                raise ValueError(f"Missing required config keys: {missing}")
            if not isinstance(config["file_templates"], dict):
                raise ValueError("'file_templates' must be a dict keyed by diag_set names.")
        else:
            required = ("path_template", "file_template")
            missing = [k for k in required if k not in config]
            if missing:
                raise ValueError(f"Missing required config keys: {missing}")

        self.config = dict(config)
        self.dask_chunks = dict(self.config.get("dask_chunks", {}))
        self.return_dask = bool(self.config.get("return_dask", False))
        self.persist_arrays = bool(self.config.get("persist_arrays", False))

        # Region tuples: [(lat_lo,lat_hi), (lon_lo,lon_hi)] or [(lat_lo,lat_hi),(lon_lo,lon_hi),(z_lo,z_hi)]
        default_regions = {
            'GB': {"name": "Global",              "region": [(-90,  90), (0,   360), (0, 120000)]},
            'NH': {"name": "Northern Hemisphere", "region": [( 20,  90), (0,   360), (0, 120000)]},
            'SH': {"name": "Southern Hemisphere", "region": [(-90,  20), (0,   360), (0, 120000)]},
            'TP': {"name": "Tropics",             "region": [(-20,  20), (0,   360), (0, 120000)]},
            'NA': {"name": "North America",       "region": [( 25,  55), (235, 295), (0, 120000)]},
        }
        self.region_dict = self.config.get("region_dict", default_regions)

    def _chunk_existing_dims(self, ds):
        chunks = {dim: size for dim, size in self.dask_chunks.items() if dim in ds.dims}
        return ds.chunk(chunks) if chunks else ds

    def _open_mfdataset(self, rpath, *, decode_times=False):
        ds = xc.open_mfdataset(rpath, decode_times=decode_times, chunks={})
        return self._chunk_existing_dims(ds)

    def _maybe_return_array(self, arr):
        if self.return_dask and hasattr(arr, "dims"):
            return arr.persist() if self.persist_arrays else arr
        if hasattr(arr, "compute"):
            return arr.compute().values
        return np.asarray(arr)

    # ------------------------------------------------------------------ #
    # Unified resolver
    # ------------------------------------------------------------------ #
    def resolve_file_path(self, exp_subrun: dict, diag_set: str) -> Tuple[str, str]:
        """
        Resolve (path, file_or_glob) for a given diag_set using config["file_templates"][diag_set].

        exp_subrun must provide: run_id, name, period, and group or group_key.
        Returns (path, file_or_glob); the second item may contain '*' if wildcarding is used.
        """
        if diag_set not in self.config["file_templates"]:
            raise KeyError(
                f"diag_set '{diag_set}' not in file_templates. "
                f"Available: {list(self.config['file_templates'].keys())}"
            )

        data_path  = self.config["data_path"]
        path_tmpl  = self.config["path_template"]
        file_tmpl  = self.config["file_templates"][diag_set]

        run_id    = exp_subrun.get("run_id", "")
        name      = exp_subrun.get("name", "")
        period    = exp_subrun.get("period", "")
        group_key = exp_subrun.get("group", exp_subrun.get("group_key", ""))

        if not (run_id and group_key and period):
            missing = [k for k in ("run_id", "group/group_key", "period")
                       if (k == "run_id" and not run_id)
                       or (k == "group/group_key" and not group_key)
                       or (k == "period" and not period)]
            raise ValueError(f"resolve_file_path() missing required fields: {missing}; got subrun={exp_subrun}")

        tokens = {
            "data_path": data_path,
            "run_id":    run_id,
            "name":      name,
            "group_key": group_key,
            "diag_set":  diag_set,
            "period":    period,
        }

        path = path_tmpl.format(**tokens)

        tried = []

        # 1) Primary (allow glob straight away)
        primary = file_tmpl.format(**tokens)
        p_full = os.path.join(path, primary)
        tried.append(p_full)
        if glob.glob(p_full):
            print(f"[RESOLVE] {diag_set}: {p_full}")
            return path, primary

        # 2) Swap {run_id} <-> {name} and try again
        alt = None
        if "{run_id}" in file_tmpl and name and name != run_id:
            alt = file_tmpl.replace("{run_id}", "{name}").format(**tokens)
        elif "{name}" in file_tmpl and run_id and name and name != run_id:
            alt = file_tmpl.replace("{name}", "{run_id}").format(**tokens)

        if alt:
            a_full = os.path.join(path, alt)
            tried.append(a_full)
            if glob.glob(a_full):
                print(f"[RESOLVE] {diag_set} (swap): {a_full}")
                return path, alt

        # 3) Wildcard period
        if "{period}" in file_tmpl:
            wc = file_tmpl.replace("{period}", "*").format(**tokens)
            w_full = os.path.join(path, wc)
            tried.append(w_full)
            if glob.glob(w_full):
                print(f"[RESOLVE] {diag_set} (wildcard): {w_full}")
                return path, wc

            if alt:
                swapped = (file_tmpl.replace("{run_id}", "{name}")
                           if "{run_id}" in file_tmpl else
                           file_tmpl.replace("{name}", "{run_id}"))
                wc2 = swapped.replace("{period}", "*").format(**tokens)
                w2_full = os.path.join(path, wc2)
                tried.append(w2_full)
                if glob.glob(w2_full):
                    print(f"[RESOLVE] {diag_set} (swap+wildcard): {w2_full}")
                    return path, wc2

        print("[WARN] No matching files for", diag_set, "Tried:")
        for t in tried:
            print("  ", t)
        return path, primary

    def resolve_dart_file_path(
        self, exp_info: dict, exp: str,
        diag_set: str, date: str
    ) -> Tuple[str, str]:
        """Resolve legacy obs_diag config entries used by the analysis_atm_da notebooks."""
        diag_name = exp_info.get(diag_set, diag_set)
        run_name = exp_info.get("run", exp_info.get("run_id", exp))
        case_name = exp_info.get("key", exp_info.get("group", exp_info.get("group_key", "")))

        path_template = exp_info.get("path_template", self.config["path_template"])
        path = (
            path_template
            .replace("%(RUNNAME)", run_name)
            .replace("%(CASENAME)", case_name)
            .replace("%(DIAG)", diag_name)
        )
        file = (
            self.config["file_template"]
            .replace("%(RUNNAME)", run_name)
            .replace("%(KEY)", self.config.get("diag_key", "obs_diag_output"))
            .replace("%(TIME)", date)
        )
        return path, file

    @staticmethod
    def define_region(regnam: str = "global") -> Tuple[Tuple[float, float], Tuple[float, float]]:
        reg_dict = {
            "global": [(-90, 90), (-180, 180)],
            "Northern Hemisphere": [(20, 90), (-180, 180)],
            "Southern Hemisphere": [(-90, -20), (-180, 180)],
            "Tropics": [(-20, 20), (-180, 180)],
            "North America": [(25, 55), (-125, -65)],
        }
        if regnam not in reg_dict:
            raise KeyError(f"Region {regnam!r} is not defined. Available regions: {list(reg_dict)}")
        return reg_dict[regnam]

    # ------------------------------------------------------------------ #
    # obs_diag helpers
    # ------------------------------------------------------------------ #
    def build_ts_var_dict(
        self, var_key: str = None, name: str = None,
        y1axis: Optional[List[float]] = None, y2axis: Optional[List[float]] = None
    ) -> dict:
        """
        Build variable-entry dict for time/level diagnostics.
        Provides corrected keys and legacy typos for back-compat.
        """
        var_key = var_key or "RADIOSONDE_U"
        name = name or "RADIOSONDE_U_WIND_COMPONENT"

        y1_default = [0, 10] if y1axis is None else y1axis
        y2_default = [0, 100] if y2axis is None else y2axis

        entry = {
            "name": name,
            "lev_type": "pressure",
            "CopySpread": "totalspread",
            "CopyRMSE": "rmse",
            "CopyBias": "bias",
            "CopyNposs": "Nposs",
            "CopyNused": "Nused",
            "type1": "guess",
            "type2": "VPguess",
            "type3": "guess_RankHist",

            # Correct keys
            "y1axis": y1_default,
            "y2axis": [0, 100],
            "y1axis0": y2_default,
            "y2axis0": [0, 100],

            # Legacy aliases
            "y1aix": y1_default,
            "y2aix": [0, 100],
            "y1aix0": y2_default,
            "y2aix0": [0, 100],
        }
        return {var_key: entry}

    def create_lev_str(self, lev, levp, lev_type) -> List[str]:
        """Create human-readable layer labels from mid-levels and edges (defensive)."""
        levstr = []
        n = len(lev)
        for i in range(n):
            if i + 1 >= len(levp):
                break
            if lev_type == "pressure":
                levstr.append(f"{int(levp[i+1])}-{int(levp[i])} hPa")
            elif lev_type == "height":
                levstr.append(f"{levp[i]}-{levp[i+1]} m")
            elif lev_type == "model":
                levstr.append(f"{levp[i]}-{levp[i+1]} layer")
            else:
                raise ValueError(f"Unknown lev_type: {lev_type}")
        return levstr

    # ------------------------------------------------------------------ #
    # obs_diag reader
    # ------------------------------------------------------------------ #
    def read_dart_obs_diag(
        self, regnam: str, var: str, dtype: str, var_dict: Dict[str, str],
        date: str, path: str, file: str
    ) -> Tuple:
        """
        Read a DART obs_diag file (or glob) and extract arrays.
        """
        rpath = os.path.join(path, file)
        print(f"Reading file(s): {rpath}")

        def _opt(ds, name):
            return ds[name].values if name in ds.variables else np.array([])

        def _first_idx(tup):
            arr = tup[0] if isinstance(tup, tuple) else np.asarray(tup)
            return int(arr[0]) if arr.size else None

        ds = self._open_mfdataset(rpath, decode_times=False)
        close_ds = not self.return_dask
        try:
            time_var = ds["time"]
            time = time_var.values
            self._last_time_unit = time_var.attrs.get("units")

            mlevel       = _opt(ds, "mlevel")
            mlevel_edges = _opt(ds, "mlevel_edges")
            plevel       = _opt(ds, "plevel")
            plevel_edges = _opt(ds, "plevel_edges")
            hlevel       = _opt(ds, "hlevel")
            hlevel_edges = _opt(ds, "hlevel_edges")
            rank_bins    = _opt(ds, "rank_bins")  # optional

            # region name -> index (fallback to 0)
            if "region_names" in ds.variables:
                region_names = np.array([
                    ch.decode("utf-8").strip() if isinstance(ch, (bytes, bytearray)) else str(ch).strip()
                    for ch in ds["region_names"].values
                ])
                matches = np.where(region_names == regnam)[0]
                if matches.size == 0:
                    raise ValueError(
                        f"Region '{regnam}' not found. Available: {', '.join(map(str, region_names))}"
                    )
                ind_reg = int(matches[0])
            else:
                ind_reg = 0

            CopyMetaData = np.array([
                ch.decode("utf-8").strip() if isinstance(ch, (bytes, bytearray)) else str(ch).strip()
                for ch in ds["CopyMetaData"].values
            ])
            ind_vars = _first_idx(np.where(CopyMetaData == var_dict["CopySpread"]))
            ind_rmse = _first_idx(np.where(CopyMetaData == var_dict["CopyRMSE"]))
            ind_bias = _first_idx(np.where(CopyMetaData == var_dict.get("CopyBias", "bias")))
            ind_npos = _first_idx(np.where(CopyMetaData == var_dict["CopyNposs"]))
            ind_nuse = _first_idx(np.where(CopyMetaData == var_dict["CopyNused"]))

            sprd = rmse = bias = npos = nuse = hrank = np.array([])

            if dtype == "guess":
                varname = f"{var_dict['name']}_guess"
                if varname in ds.variables:
                    base = ds[varname]  # [time, copy, level, region]
                    sprd = self._maybe_return_array(base[:, ind_vars, :, ind_reg])
                    rmse = self._maybe_return_array(base[:, ind_rmse, :, ind_reg])
                    if ind_bias is not None:
                        bias = self._maybe_return_array(base[:, ind_bias, :, ind_reg])
                    npos = self._maybe_return_array(base[:, ind_npos, :, ind_reg])
                    nuse = self._maybe_return_array(base[:, ind_nuse, :, ind_reg])

            elif dtype == "VPguess":
                varname = f"{var_dict['name']}_{var_dict['type2']}"
                if varname in ds.variables:
                    base = ds[varname]  # [copy, level, region]
                    sprd = self._maybe_return_array(base[ind_vars, :, ind_reg])
                    rmse = self._maybe_return_array(base[ind_rmse, :, ind_reg])
                    if ind_bias is not None:
                        bias = self._maybe_return_array(base[ind_bias, :, ind_reg])
                    npos = self._maybe_return_array(base[ind_npos, :, ind_reg])
                    nuse = self._maybe_return_array(base[ind_nuse, :, ind_reg])

            elif dtype == "guess_RankHist":
                varname = f"{var_dict['name']}_{var_dict['type3']}"
                if varname in ds.variables:
                    hrank = self._maybe_return_array(ds[varname][0, :, :, ind_reg])
        finally:
            if close_ds:
                ds.close()

        return (time, plevel, plevel_edges, mlevel, mlevel_edges,
                hlevel, hlevel_edges, sprd, rmse, bias, npos, nuse, hrank)

    # ------------------------------------------------------------------ #
    # obs_diag -> metrics dict
    # ------------------------------------------------------------------ #
    def extract_metrics_data(
        self, var: str, var_dict: Dict[str, str],
        dtype: str, regnam: str, diag_set: str,
        exp_dict: Dict[str, dict],
    ) -> Tuple[Dict, List[str]]:
        """
        Assemble obs_diag metrics for experiments in exp_dict.

        exp_dict: { "<exp_name>": <subrun_dict> } with run_id, name, period, group/group_key
        regnam: short region key ("GB","NH",...) mapped to the long name in file coords.
        """
        data_dict: Dict[str, dict] = {}
        levstr: List[str] = []

        region_long = self.region_dict.get(regnam, {}).get("name", regnam)

        for exp, subrun in exp_dict.items():
            date = subrun["period"]
            if "file_templates" in self.config:
                path, file = self.resolve_file_path(subrun, diag_set=diag_set)
                left = date.split("-")[0]
                time_unit = f"days since {left[:4]}-{left[4:6]}-{left[6:8]} 00:00:00"
            else:
                path, file = self.resolve_dart_file_path(subrun, exp, diag_set, date)
                time_unit = f"days since {date[:4]}-{date[4:6]}-{date[6:8]}"

            (time, plev, plev_edges, mlev, mlev_edges,
             hlev, hlev_edges, sprd, rmse, bias, npos, nuse, hrank) = self.read_dart_obs_diag(
                region_long, var, dtype, var_dict, date, path, file
            )
            time_unit = getattr(self, "_last_time_unit", None) or time_unit

            # rejection (%)
            if hasattr(npos, "size") and hasattr(nuse, "size") and npos.size > 0 and nuse.size > 0:
                if isinstance(npos, xr.DataArray) or isinstance(nuse, xr.DataArray):
                    rejection = xr.where(npos > 0, 100.0 - (nuse * 100.0 / npos), np.nan)
                else:
                    rejection = np.where(npos > 0, 100.0 - (nuse * 100.0 / npos), np.nan)
            else:
                rejection = np.array([])

            # choose levels / labels
            lev_type = var_dict.get("lev_type", "pressure")
            lev_map = {"pressure": (plev, plev_edges),
                       "height":   (hlev, hlev_edges),
                       "model":    (mlev, mlev_edges)}
            if lev_type not in lev_map:
                raise ValueError(f"Invalid level type: {lev_type}")
            lev, levp = lev_map[lev_type]
            levstr = self.create_lev_str(lev, levp, lev_type)

            rel_time = np.asarray(time) - time[0] if len(time) else np.array([])

            data_dict[exp] = {
                "time": rel_time,
                "time_unit": time_unit,
                "rmse": rmse,
                "bias": bias,
                "spread": sprd,
                "nused": nuse,
                "npos": npos,
                "rejection": rejection,
                "histrank": hrank,
                "rmse_str": "RMSE",
                "bias_str": "Bias",
                "spread_str": "Total Spread",
                "nused_str": "Number Assimilated",
                "npos_str": "Number Processed",
                "rejection_str": "Data Rejection(%)",
                "lev": lev,
                "levp": levp,
                "period": date,
            }

        return data_dict, levstr

    # ------------------------------------------------------------------ #
    # shared helpers
    # ------------------------------------------------------------------ #
    def _norm_str(self, arr):
        arr = np.asarray(arr)
        flat = arr.ravel()
        out = [(x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x)).strip()
               for x in flat]
        return np.array(out).reshape(arr.shape)

    def _first_index(self, bool_or_inds):
        """Return a scalar index from np.where result or boolean array."""
        if isinstance(bool_or_inds, tuple):
            inds = bool_or_inds[0]
        else:
            inds = np.asarray(bool_or_inds).nonzero()[0]
        if inds.size == 0:
            return None
        return int(inds[0])

    def _unpack_region(self, region_val):
        """
        Accepts region as [(lat_lo,lat_hi), (lon_lo,lon_hi)] or
        [(lat_lo,lat_hi), (lon_lo,lon_hi), (z_lo,z_hi)].
        Returns (lat_lo, lat_hi, lon_lo, lon_hi, z_lo, z_hi or None).
        """
        (lat_lo, lat_hi), (lon_lo, lon_hi) = region_val[:2]
        if len(region_val) >= 3 and region_val[2] is not None:
            z_lo, z_hi = region_val[2]
        else:
            z_lo = z_hi = None
        return lat_lo, lat_hi, lon_lo, lon_hi, z_lo, z_hi

    def _resolve_obs_type_selection(self, obs_filter, ObsTypesMetaData, ObsTypes) -> np.ndarray:
        """
        Return a sorted unique array of selected type CODES (1-based, as in ObsTypes).
        obs_filter supports:
          - None / "all"
          - "group:Conventional" or "group:Satellite"
          - "regex:<pattern>"
          - single str type name or list[str]
          - single int code or list[int]
        """
        names = np.asarray(ObsTypesMetaData).astype(str).ravel()
        codes = np.asarray(ObsTypes).astype(int).ravel()

        if obs_filter is None or (isinstance(obs_filter, str) and obs_filter.lower() == "all"):
            return np.unique(codes)

        # numeric code(s)
        if isinstance(obs_filter, int):
            return np.unique(np.array([obs_filter], dtype=int))
        if isinstance(obs_filter, (list, tuple)) and all(isinstance(x, int) for x in obs_filter):
            return np.unique(np.array(obs_filter, dtype=int))

        # group:*  (by friendly name sets)
        if isinstance(obs_filter, str) and obs_filter.lower().startswith("group:"):
            gname = obs_filter.split(":", 1)[1].strip()
            group = self.extract_obs_group().get(gname, {})
            name_pool = set(group.get("plevel", [])) | set(group.get("surface", [])) | set(group.get("hlevel", []))
            if not name_pool:
                return np.array([], dtype=int)
            mask = np.isin(names, list(name_pool))
            return np.unique(codes[mask])

        # regex:*
        if isinstance(obs_filter, str) and obs_filter.startswith("regex:"):
            import re
            patt = re.compile(obs_filter.split(":", 1)[1])
            mask = np.array([bool(patt.search(n)) for n in names])
            return np.unique(codes[mask])

        # string name(s)
        if isinstance(obs_filter, str):
            mask = (names == obs_filter)
            return np.unique(codes[mask])

        if isinstance(obs_filter, (list, tuple)) and all(isinstance(x, str) for x in obs_filter):
            mask = np.isin(names, np.array(obs_filter, dtype=str))
            return np.unique(codes[mask])

        # fallback: nothing selected
        return np.array([], dtype=int)

    # ------------------------------------------------------------------ #
    # obs_seq reader
    # ------------------------------------------------------------------ #
    def read_dart_obs_seq(
        self,
        var: str,
        var_dict: dict,
        date: str,
        path: str,
        file: str,
    ):
        """
        Read DART obs_seq netCDFs (possibly multi-file with globs) and return:
        (time, lat, lon, lev, which_vert,
         ObsTypesMetaData, ObsTypes, ObsIndex, obs_type, obs_keys, qc, observations)
        """
        rpath = os.path.join(path, file)
        print(f"[READ obs_seq] {rpath}")

        with self._open_mfdataset(rpath, decode_times=False) as dr:
            # packed location or separate vars
            if "location" in dr.variables:
                location = dr["location"].compute().values
                lon = location[:, 0]
                lat = location[:, 1]
                lev = location[:, 2]
            else:
                lon = dr["lon"].compute().values if "lon" in dr.variables else dr["longitude"].compute().values
                lat = dr["lat"].compute().values if "lat" in dr.variables else dr["latitude"].compute().values
                lev = dr["lev"].compute().values if "lev" in dr.variables else dr.get(
                    "level", xr.DataArray(np.full_like(lon, np.nan))
                ).compute().values

            time       = dr["time"].compute().values
            which_vert = dr["which_vert"].compute().values if "which_vert" in dr.variables else dr.get(
                "WhichVert", xr.DataArray(np.full_like(time, -2))
            ).compute().values

            # metadata / channels
            CopyMetaData = self._norm_str(dr["CopyMetaData"].values)
            observations = dr["observations"].compute().values
            QCMetaData   = self._norm_str(dr["QCMetaData"].values)
            qc           = dr["qc"].compute().values

            copy_name = var_dict.get("CopyString", "observation")
            qc_name   = var_dict.get("QCString", "DART quality control")
            ind_copy  = self._first_index(np.where(CopyMetaData == copy_name))
            ind_qc    = self._first_index(np.where(QCMetaData   == qc_name))
            if ind_copy is None:
                raise KeyError(f"Copy channel '{copy_name}' not found in CopyMetaData: {CopyMetaData.tolist()}")
            if ind_qc is None:
                raise KeyError(f"QC channel '{qc_name}' not found in QCMetaData: {QCMetaData.tolist()}")

            # observation type bookkeeping
            ObsTypesMetaData = self._norm_str(dr["ObsTypesMetaData"].values)
            ObsTypes         = dr["ObsTypes"].values
            ObsIndex         = dr["ObsIndex"].compute().values
            obs_type         = dr["obs_type"].compute().values
            obs_keys         = dr["obs_keys"].compute().values

            # region filter (lon/lat always; z if finite)
            reg_xyz = np.array(str(var_dict.get("region", "0 360 -90 90 -1e36 1e36")).split(), dtype=float)
            lon_min, lon_max, lat_min, lat_max, z_min, z_max = reg_xyz

            keep = ((lon_min <= lon) & (lon <= lon_max) &
                    (lat_min <= lat) & (lat <= lat_max))

            if np.isfinite(z_min) and np.isfinite(z_max):
                keep &= ((z_min <= lev) & (lev <= z_max))

            ind_loc = np.where(keep)[0]
            if ind_loc.size == 0:
                raise ValueError("No observations within requested lon/lat/z window: "
                                 f"[{lon_min}, {lon_max}] x [{lat_min}, {lat_max}] x [{z_min}, {z_max}]")

            # slice consistently
            time       = time[ind_loc]
            lon        = lon[ind_loc]
            lat        = lat[ind_loc]
            lev        = lev[ind_loc]
            which_vert = which_vert[ind_loc]
            ObsIndex   = ObsIndex[ind_loc]
            obs_type   = obs_type[ind_loc]
            obs_keys   = obs_keys[ind_loc]
            observations = observations[ind_loc, ind_copy]
            qc           = qc[ind_loc, ind_qc]

        return (time, lat, lon, lev, which_vert,
                ObsTypesMetaData, ObsTypes, ObsIndex, obs_type, obs_keys, qc, observations)

    # ------------------------------------------------------------------ #
    # obs_seq one-call pipeline (returns legacy plotting data_dict)
    # ------------------------------------------------------------------ #
    def make_obs_seq_datadict(
        self,
        exp_name: str,
        da_run_subrun: dict,
        group_key: str,
        date_stamp: str,
        reg_key: str = "GB",
        var: str = "all",
        var_dict_entry: Optional[dict] = None,
        obs_filter: Optional[object] = None,   # flexible selection (see helper)
    ) -> Dict[str, Dict[str, Dict[str, list]]]:
        """
        Build the legacy obs_seq data_dict using a single DA sub-run.

        obs_filter: None/"all", str name, list[str], int code, list[int],
                    "group:Conventional"/"group:Satellite", or "regex:<pattern>"
        """
        # --- defaults (region comes from region_dict) ---
        if var_dict_entry is None:
            var_dict_entry = {
                "ObsTypeString": "all",
                "CopyString": "observation",
                "QCString": "DART quality control",
                "verbose": False,
            }
        else:
            var_dict_entry = dict(var_dict_entry)

        # --- region override (lat/lon from region_dict; z from var_dict_entry or region_dict or defaults) ---
        if hasattr(self, "region_dict") and reg_key in self.region_dict:
            lat_lo, lat_hi, lon_lo, lon_hi, z_lo, z_hi = self._unpack_region(self.region_dict[reg_key]["region"])
            try:
                # Use z from incoming var_dict_entry if present
                zmin, zmax = map(float, str(var_dict_entry.get("region", "")).split()[4:6])
            except Exception:
                if z_lo is not None and z_hi is not None:
                    zmin, zmax = float(z_lo), float(z_hi)
                else:
                    zmin, zmax = (0.0, 120000.0)
            var_dict_entry["region"] = f"{lon_lo} {lon_hi} {lat_lo} {lat_hi} {zmin} {zmax}"

        # --- resolver tokens ---
        subrun = dict(da_run_subrun)
        subrun["group"]  = group_key
        subrun["period"] = date_stamp
        if "name" not in subrun or not subrun["name"]:
            subrun["name"] = exp_name

        # --- resolve files ---
        path, file_glob = self.resolve_file_path(subrun, diag_set="obs_seq")
        print(f"[OBS-SEQ] files: {os.path.join(path, file_glob)}")

        # --- read ---
        (time, lat, lon, lev, which_vert,
         ObsTypesMetaData, ObsTypes, ObsIndex, obs_type, obs_keys, qc, observations) = self.read_dart_obs_seq(
            var, var_dict_entry, date_stamp, path, file_glob
        )

        # --- select obs types if requested (by codes) ---
        sel_codes = self._resolve_obs_type_selection(obs_filter, ObsTypesMetaData, ObsTypes)
        if sel_codes.size > 0:
            keep = np.isin(np.asarray(obs_type).astype(int).ravel(), sel_codes)
            if not np.any(keep):
                return {}
            idx = np.where(keep)[0]
            time       = np.asarray(time)[idx]
            lat        = np.asarray(lat)[idx]
            lon        = np.asarray(lon)[idx]
            lev        = np.asarray(lev)[idx]
            which_vert = np.asarray(which_vert)[idx]
            ObsIndex   = np.asarray(ObsIndex)[idx]
            obs_type   = np.asarray(obs_type)[idx]
            obs_keys   = np.asarray(obs_keys)[idx]
            observations = np.asarray(observations)[idx]
            qc           = np.asarray(qc)[idx]

        # --- proceed with grouping ---
        obs_type = np.asarray(obs_type).astype(int).ravel()
        unique_types = np.unique(obs_type)
        valid_mask = (unique_types >= 1) & (unique_types <= len(ObsTypesMetaData))
        unique_types = unique_types[valid_mask]
        if unique_types.size == 0:
            return {}

        ObsTypesFiltered         = np.asarray(ObsTypes)[unique_types - 1]
        ObsTypesMetaDataFiltered = np.asarray(ObsTypesMetaData)[unique_types - 1]

        # --- diagnostic printout: show extracted obs types -------------------
        print(f"[INFO] Extracted {len(ObsTypesMetaDataFiltered)} obs types from file:")
        for i, name in enumerate(ObsTypesMetaDataFiltered):
            code = int(ObsTypesFiltered[i]) if i < len(ObsTypesFiltered) else None
            print(f"    {i+1:3d}: {name:40s} (code={code})")

        obs_group = self.extract_obs_group()
        data_dict: Dict[str, Dict[str, Dict[str, list]]] = {}

        for grp in obs_group:
            group_spec = obs_group[grp]
            for i, otyp in enumerate(ObsTypesMetaDataFiltered):
                if (otyp in group_spec.get("plevel", [])
                        or otyp in group_spec.get("surface", [])
                        or otyp in group_spec.get("hlevel", [])):
                    subtyp = otyp.split("_")[0].lower().capitalize()
                    dd_grp = data_dict.setdefault(grp, {}).setdefault(subtyp, {
                        "lat": [], "lon": [], "time": [], "lev": [], "obs": []
                    })
                    tcode = ObsTypesFiltered[i]
                    inds = np.where(obs_type == int(tcode))
                    dd_grp["lat" ].append(np.asarray(lat)[inds])
                    dd_grp["lon" ].append(np.asarray(lon)[inds])
                    dd_grp["time"].append(np.asarray(time)[inds])
                    dd_grp["lev" ].append(np.asarray(lev)[inds])
                    dd_grp["obs" ].append(np.asarray(observations)[inds])

        return data_dict

    # ------------------------------------------------------------------ #
    # Obs group mapping (as in your original code)
    # ------------------------------------------------------------------ #
    @staticmethod
    def extract_obs_group() -> dict:
        return {
            "Conventional": {
                "plevel": [
                    "TEMPERATURE", "SPECIFIC_HUMIDITY", "PRESSURE",
                    "RADIOSONDE_U_WIND_COMPONENT", "RADIOSONDE_V_WIND_COMPONENT",
                    "RADIOSONDE_GEOPOTENTIAL_HGT", "RADIOSONDE_TEMPERATURE",
                    "RADIOSONDE_SPECIFIC_HUMIDITY", "DROPSONDE_TEMPERATURE",
                    "DROPSONDE_U_WIND_COMPONENT", "DROPSONDE_V_WIND_COMPONENT",
                    "DROPSONDE_SPECIFIC_HUMIDITY", "AIRCRAFT_U_WIND_COMPONENT",
                    "AIRCRAFT_V_WIND_COMPONENT", "AIRCRAFT_TEMPERATURE",
                    "AIRCRAFT_SPECIFIC_HUMIDITY", "ACARS_U_WIND_COMPONENT",
                    "ACARS_V_WIND_COMPONENT", "ACARS_TEMPERATURE",
                    "ACARS_SPECIFIC_HUMIDITY",
                ],
                "surface": [
                    "RADIOSONDE_SURFACE_PRESSURE", "DROPSONDE_SURFACE_PRESSURE",
                    "RADIOSONDE_SURFACE_ALTIMETER", "DROPSONDE_SURFACE_ALTIMETER",
                    "METAR_ALTIMETER", "MESONET_SURFACE_ALTIMETER",
                    "MARINE_SFC_U_WIND_COMPONENT", "MARINE_SFC_V_WIND_COMPONENT",
                    "MARINE_SFC_TEMPERATURE", "MARINE_SFC_SPECIFIC_HUMIDITY",
                    "MARINE_SFC_PRESSURE", "LAND_SFC_U_WIND_COMPONENT",
                    "LAND_SFC_V_WIND_COMPONENT", "LAND_SFC_TEMPERATURE",
                    "LAND_SFC_SPECIFIC_HUMIDITY", "LAND_SFC_PRESSURE",
                    "MARINE_SFC_ALTIMETER", "LAND_SFC_ALTIMETER",
                ],
            },
            "Satellite": {
                "hlevel": ["GPSRO_REFRACTIVITY"],
                "plevel": [
                    "SAT_TEMPERATURE", "SAT_TEMPERATURE_ELECTRON",
                    "SAT_TEMPERATURE_ION", "SAT_DENSITY_NEUTRAL_O3P", "SAT_DENSITY_NEUTRAL_O2",
                    "SAT_DENSITY_NEUTRAL_N2", "SAT_DENSITY_NEUTRAL_N4S", "SAT_DENSITY_NEUTRAL_NO",
                    "SAT_DENSITY_NEUTRAL_N2D", "SAT_DENSITY_NEUTRAL_N2P", "SAT_DENSITY_NEUTRAL_H",
                    "SAT_DENSITY_NEUTRAL_HE", "SAT_DENSITY_NEUTRAL_CO2", "SAT_DENSITY_NEUTRAL_O1D",
                    "SAT_DENSITY_ION_O4SP", "SAT_DENSITY_ION_O2P", "SAT_DENSITY_ION_N2P",
                    "SAT_DENSITY_ION_NP", "SAT_DENSITY_ION_O2DP", "SAT_DENSITY_ION_O2PP",
                    "SAT_DENSITY_ION_HP", "SAT_DENSITY_ION_HEP", "SAT_DENSITY_ION_E",
                    "SAT_VELOCITY_U", "SAT_DENSITY_ION_NOP", "SAT_VELOCITY_V", "SAT_VELOCITY_W",
                    "SAT_VELOCITY_U_ION", "SAT_VELOCITY_V_ION", "SAT_VELOCITY_W_ION",
                    "SAT_VELOCITY_VERTICAL_O3P", "SAT_VELOCITY_VERTICAL_O2",
                    "SAT_VELOCITY_VERTICAL_N2", "SAT_VELOCITY_VERTICAL_N4S",
                    "SAT_VELOCITY_VERTICAL_NO", "SAT_F107", "SAT_RHO", "GPS_PROFILE",
                    "COSMIC_ELECTRON_DENSITY", "GND_GPS_VTEC", "CHAMP_DENSITY",
                    "MIDAS_TEC", "SSUSI_O_N2_RATIO", "GPS_VTEC_EXTRAP", "SABER_TEMPERATURE",
                    "AURAMLS_TEMPERATURE", "SAT_U_WIND_COMPONENT", "SAT_V_WIND_COMPONENT",
                    "ATOV_TEMPERATURE", "AIRS_TEMPERATURE", "AIRS_SPECIFIC_HUMIDITY",
                    "GPS_PRECIPITABLE_WATER", "CIMMS_AMV_U_WIND_COMPONENT",
                    "CIMMS_AMV_V_WIND_COMPONENT",
                ],
                "surface": [
                    "VADWND_U_WIND_COMPONENT", "VADWND_V_WIND_COMPONENT",
                ],
            },
        }
