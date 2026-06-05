import os
import glob
from typing import Dict, List, Tuple, Optional

import numpy as np
import xarray as xr

class ENSODiagReader:
    """
    Helper to read ENSO diagnostics (perf/proc/tel) across
    historical vs future and ensemble members.

    Usage:
        reader = ENSODiagReader(...)
        da_hist = reader.load("ENSO_perf", "enso_amplitude", period="hist")
        da_fut  = reader.load("ENSO_perf", "enso_amplitude", period="future")
    """

    def __init__(
        self,
        data_dir: str,
        model: str,
        groups: List[str],
        period_list: List[Tuple[int, int]],
        nens: List[int],
        enso_groups: Optional[Dict[str, List[str]]] = None,
        file_suffix_map: Optional[Dict[str, str]] = None,
        members: Optional[List[int]] = None,
        verbose: bool = False,
    ):
        self.data_dir = data_dir
        self.model = model
        self.groups = groups
        self.period_list = period_list
        self.nens = nens

        # Use provided maps if given; otherwise fall back to internal defaults
        self.enso_groups = enso_groups if enso_groups is not None else self.get_enso_var()
        self.file_suffix_map = (
            file_suffix_map if file_suffix_map is not None else self.get_file_suffix()
        )

        # Ensemble member IDs; if not given, infer from directory listing
        self.members = members  # e.g. [51, 91, 101, ...]
        self._cached_member_dirs: Dict[str, List[str]] = {}  # (group_key) -> [paths]

        # Verbosity
        self.verbose = verbose

    # ---------------------------------------------------------
    # Dictionary getters (FLEXIBLE)
    # ---------------------------------------------------------
    def get_enso_var(self) -> Dict[str, List[str]]:
        """Return mapping: ENSO group -> list of variable names."""
        return {
            "ENSO_perf": [
                "pr_lat_rmse", "pr_lon_rmse", "sst_lon_rmse", "taux_lon_rmse",
                "enso_amplitude", "enso_duration", "enso_seasonality",
                "enso_sst_diversity_mode1", "enso_sst_diversity_mode2",
                "enso_sst_lon_rmse", "enso_sst_skewness", "enso_sst_ts_rmse",
                "seasonal_pr_lat_rmse", "seasonal_pr_lon_rmse",
                "seasonal_sst_lon_rmse", "seasonal_taux_lon_rmse",
            ],

            "ENSO_proc": [
                "sst_lon_rmse", "taux_lon_rmse", "enso_amplitude",
                "enso_dsst_oce_mode1", "enso_dsst_oce_mode2",
                "enso_fb_ssh_sst", "enso_fb_sst_taux", "enso_fb_sst_thf",
                "enso_fb_taux_ssh", "enso_seasonality",
                "enso_sst_lon_rmse", "enso_sst_skewness",
            ],

            "ENSO_tel": [
                "enso_amplitude", "enso_pr_map_djf", "enso_pr_map_jja",
                "enso_seasonality", "enso_sst_lon_rmse",
                "enso_sst_map_djf", "enso_sst_map_jja",
            ],
        }

    def get_file_suffix(self) -> Dict[str, str]:
        """Return mapping: logical variable name → file suffix used in filenames."""
        return {
            # ENSO_perf
            "pr_lat_rmse":              "BiasPrLatRmse",
            "pr_lon_rmse":              "BiasPrLonRmse",
            "sst_lon_rmse":             "BiasSstLonRmse",
            "taux_lon_rmse":            "BiasTauxLonRmse",
            "enso_amplitude":           "EnsoAmpl",
            "enso_duration":            "EnsoDuration",
            "enso_seasonality":         "EnsoSeasonality",
            "enso_sst_diversity_mode1": "EnsoSstDiversity_1",
            "enso_sst_diversity_mode2": "EnsoSstDiversity_2",
            "enso_sst_lon_rmse":        "EnsoSstLonRmse",
            "enso_sst_skewness":        "EnsoSstSkew",
            "enso_sst_ts_rmse":         "EnsoSstTsRmse",
            "seasonal_pr_lat_rmse":     "SeasonalPrLatRmse",
            "seasonal_pr_lon_rmse":     "SeasonalPrLonRmse",
            "seasonal_sst_lon_rmse":    "SeasonalSstLonRmse",
            "seasonal_taux_lon_rmse":   "SeasonalTauxLonRmse",

            # ENSO_proc
            "enso_dsst_oce_mode1":      "EnsoDeltaSstOceMode1",
            "enso_dsst_oce_mode2":      "EnsoDeltaSstOceMode2",
            "enso_fb_ssh_sst":          "EnsoFbSshSst",
            "enso_fb_sst_taux":         "EnsoFbSstTaux",
            "enso_fb_sst_thf":          "EnsoFbSstThf",
            "enso_fb_taux_ssh":         "EnsoFbTauxSsh",

            # ENSO_tel
            "enso_pr_map_djf":          "EnsoPrMapDJF",
            "enso_pr_map_jja":          "EnsoPrMapJJA",
            "enso_sst_map_djf":         "EnsoSstMapDJF",
            "enso_sst_map_jja":         "EnsoSstMapJJA",
        }

    # ------------------------------ small helpers ------------------------------

    def available_groups(self) -> List[str]:
        """Return list of available ENSO groups."""
        return list(self.enso_groups.keys())

    def available_vars(self, enso_group: str) -> List[str]:
        """Return list of variable names for a given ENSO group."""
        if enso_group not in self.enso_groups:
            raise ValueError(
                f"Unknown ENSO group '{enso_group}'. "
                f"Must be one of {list(self.enso_groups.keys())}"
            )
        return self.enso_groups[enso_group]

    # ------------------------------ core helpers ------------------------------
    def _unify_longitude_name(
        self,
        ds: xr.Dataset,
        coord_candidates=("longitude", "lon", "LONGITUDE", "LON"),
        unified_name="longitude",
    ) -> xr.Dataset:
        lon_name = None
        for cand in coord_candidates:
            if cand in ds.coords:
                lon_name = cand
                break

        if lon_name is None:
            return ds

        if lon_name != unified_name:
            ds = ds.rename({lon_name: unified_name})

        return ds

    def _get_period_index(self, period: str) -> int:
        if period not in self.groups:
            raise ValueError(f"period must be one of {self.groups}, got {period}")
        return self.groups.index(period)

    def _get_enso_case_name(self, period: str) -> str:
        idx = self._get_period_index(period)
        start, end = self.period_list[idx]
        return f"ENSO_{start}-{end}"

    def _list_member_dirs(self, period: str) -> List[str]:
        """
        Return the list of member directories for a given period, e.g.
        [..., '<DATA_DIR>/hist/v3.LR.historical_0051', ...].
        """
        key = period
        if key in self._cached_member_dirs:
            return self._cached_member_dirs[key]

        base_dir = os.path.join(self.data_dir, period)

        # match v3.LR.historical_****
        pattern = os.path.join(base_dir, f"{self.model}_*")
        dirs = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))

        # If user specified explicit members, filter accordingly
        if self.members is not None:
            keep = []
            for d in dirs:
                mstr = os.path.basename(d).split("_")[-1]
                try:
                    mid = int(mstr)
                except ValueError:
                    continue
                if mid in self.members:
                    keep.append(d)
            dirs = keep

        # Optionally trim to NENS[period_index]
        idx = self._get_period_index(period)
        if len(dirs) > self.nens[idx]:
            dirs = dirs[: self.nens[idx]]

        self._cached_member_dirs[key] = dirs
        return dirs

    def _find_nc_file(
        self,
        member_dir: str,
        enso_group: str,
        suffix: str,
    ) -> str:
        """
        Find the NetCDF file in the enso_group directory whose name ends with
        '_<suffix>.nc'. We keep it flexible w.r.t. the date stamp, etc.
        """
        enso_root = os.path.join(
            member_dir,
            "pcmdi_diags",
            "model_vs_obs",
            "metrics_data",
            "enso_metric",
            enso_group,
        )

        if not os.path.isdir(enso_root):
            raise FileNotFoundError(f"ENSO directory not found: {enso_root}")

        pattern = os.path.join(enso_root, f"*_{suffix}.nc")
        matches = sorted(glob.glob(pattern))
        print(pattern)
        if not matches:
            raise FileNotFoundError(f"No file matching *_{suffix}.nc in {enso_root}")
        if len(matches) > 1:
            # In practice you likely have only one; if multiple, take the last.
            return matches[-1]

        if self.verbose:
            print(f"enso metrics file found: {matches[0]}")

        return matches[0]

    def _choose_default_var(
        self,
        ds: xr.Dataset,
        candidates: List[str],
        member_str: str,
    ) -> str:
        """
        Heuristic to choose a model diagnostic variable when nc_var is not provided.
        Prefer:
          1) vars containing member_str,
          2) 1D longitude vars over 2D maps,
          3) first candidate as a final fallback.
        """
        # prefer vars that contain this member id in the name
        cand_member = [v for v in candidates if member_str in v]

        # prefer 1D longitude vars (no latitude) for amplitude-like metrics
        cand_lon = [
            v for v in cand_member
            if ("longitude" in ds[v].dims and "latitude" not in ds[v].dims)
        ]
        if len(cand_lon) == 1:
            return cand_lon[0]
        if len(cand_lon) > 1:
            return cand_lon[0]

        if len(cand_member) == 1:
            return cand_member[0]
        if len(cand_member) > 1:
            return cand_member[0]

        # if all else fails, just take the first candidate
        return candidates[0]

    def _extract_obs(
        self,
        ds: xr.Dataset,
        base_var: str = "sstStd_lon",
        ref_tag: str = "ERA-Interim",
    ) -> xr.DataArray:
        """
        Load the observational ENSO amplitude zonal std curve
        (e.g., sstStd_lon__ERA-Interim(longitude)) from one file.

        Parameters
        ----------
        base_var : str
            Base pattern for the variable name, e.g., "sstStd_lon".
        ref_tag : str
            Substring that identifies the obs variable name,
            e.g., "ERA-Interim".

        Returns
        -------
        obs_da : DataArray(longitude)
        """
        if base_var is None or ref_tag is None:
            raise ValueError(
                "base_var and ref_tag must be non-None when extracting observations."
            )

        # require BOTH base pattern and ref_tag
        candidates = [
            name for name, var in ds.data_vars.items()
            if (base_var in name and ref_tag in name)
        ]
        if not candidates:
            raise RuntimeError(
                f"No obs vars containing '{base_var}' and '{ref_tag}' "
                f"found in dataset variables."
            )
            
        chosen = candidates[0]          # e.g. sstStd_lon__ERA-Interim
        obs_da = ds[chosen].squeeze()
        # rename to drop the observation name 
        obs_da = obs_da.rename(base_var)
        
        ds.close()
        return obs_da

    # ------------------------------ public API ------------------------------
    def load(
        self,
        enso_group: str,
        var_name: str,
        period: str = "hist",
        nc_var: Optional[str] = None,
        ref_tag: Optional[str] = None,
    ) -> xr.DataArray:
        """
        Load a given ENSO diagnostic for all members of one period.

        Parameters
        ----------
        enso_group : {"ENSO_perf", "ENSO_proc", "ENSO_tel"}
        var_name   : variable key from self.enso_groups[enso_group], e.g. "enso_amplitude"
        period     : "hist" or "future"
        nc_var     : optional variable name inside the NetCDF file.
                     If None, the first data_var is used.

        Returns
        -------
        da_model, da_obs : xarray.DataArray
            With dims: member + (whatever dims the metric has).
            Coordinates include 'member' (int) and 'member_str'.
        """
        # Sanity checks
        if enso_group not in self.enso_groups:
            raise ValueError(
                f"Unknown ENSO group {enso_group}. "
                f"Must be one of {list(self.enso_groups.keys())}"
            )
        if var_name not in self.enso_groups[enso_group]:
            raise ValueError(
                f"var_name '{var_name}' not in enso_groups['{enso_group}']"
            )
        if var_name not in self.file_suffix_map:
            raise KeyError(
                f"No file suffix mapping for '{var_name}'. "
                f"Add it to the file_suffix_map."
            )

        suffix = self.file_suffix_map[var_name]
        member_dirs = self._list_member_dirs(period)

        da_model = {}
        da_obs = {}
        for i,mdir in enumerate(member_dirs):
            base = os.path.basename(mdir)  # e.g. v3.LR.historical_0051
            mstr = "{:02d}".format(i)
            
            nc_path = self._find_nc_file(mdir, enso_group, suffix)
            ds = xr.open_dataset(nc_path,decode_times=False)

            data_vars = list(ds.data_vars)
            # drop bounds_* helpers
            candidates = [v for v in data_vars if not v.startswith("bounds_")]

            if nc_var is not None:
                # nc_var is a BASE PATTERN like "sstStd_lon"
                # -> for model, require pattern + member id
                cand_pattern = [v for v in candidates if nc_var in v]
                cand_member  = [v for v in cand_pattern if mstr in v]

                if len(cand_member) >= 1:
                    data_var = cand_member[0]   # e.g. sstStd_lon__v3-LR_0051
                elif len(cand_pattern) >= 1:
                    # fallback: first match on pattern
                    data_var = cand_pattern[0]
                else:
                    # fallback to generic heuristic
                    data_var = self._choose_default_var(ds, candidates, mstr)
            else:
                data_var = self._choose_default_var(ds, candidates, mstr)

            # introduce member coordinate
            mod_da = ds[data_var].squeeze()

            # rename to drop the observation name 
            mod_da = mod_da.rename(nc_var)
        
            da_model[mstr] = mod_da
            
            # extract observation or reference vars 
            do = self._extract_obs(
                ds,
                base_var=nc_var,
                ref_tag=ref_tag
            )
            da_obs[mstr] = do
            
            ds.close()

        if not da_obs or not da_model:
            raise RuntimeError(
                f"No data loaded for {enso_group}/{var_name}/{period}"
            )
            
        return da_model, da_obs

    def load_metric_data(
        self,
        enso_group: str,
        var_name: str,
        nc_var: Optional[str] = None,
        ref_dict: Optional[dict] = None,
        period_list: Optional[list] = None,
    ) -> xr.Dataset:
        """
        Convenience wrapper: load one or more periods and return a Dataset
        with an extra 'period' dimension.

        Parameters
        ----------
        period_list : sequence of str, optional
            Period tags to pass to `self.load`, e.g. ["hist"] or
            ["hist", "future"]. If None, defaults to ["hist", "future"].

        Returns
        -------
        ds_model, ds_obs : xr.Dataset
            Both with variable:
            - 'metric' : dims (period, member, ...)
        """
        # Decide which period tags to load
        if period_list is None:
            periods = list(self.groups)
        else:
            periods = period_list

        dm_list = {}
        do_list = {}

        for per in periods:
            print(f"processing period: {per}")
            dm, do = self.load(
                enso_group, var_name, period=per, nc_var=nc_var, ref_tag=ref_dict[per]
            )
            dm_list[per] = dm 
            do_list[per] = do 
            
        return dm_list, do_list
        
    def combine_members_to_array(self, member_dict, sample_dim=None):
        """
        Combine member DataArrays into a single DataArray with a new 'member' dim.
    
        Parameters
        ----------
        member_dict : dict
            Mapping like {'00': DataArray, '01': DataArray', ...}
    
        Returns
        -------
        da_model : xarray.DataArray
            dims: ('member', *original_dims)
            coords:
              - 'member'     : numeric or string member index
              - 'member_str' : string label for each member
              - all original coords preserved for non-'member' dims
        """
        if not member_dict:
            raise ValueError("member_dict is empty.")
    
        # Use first DataArray as template
        template = next(iter(member_dict.values()))
        template = template.squeeze()
    
        das = []
    
        for mem_key, da in member_dict.items():
            da = da.squeeze()
    
            # Optional: sanity check dims/sizes match template
            if da.dims != template.dims:
                raise ValueError(
                    f"Member '{mem_key}' has dims {da.dims}, expected {template.dims}"
                )
            for d in template.dims:
                if da.sizes[d] != template.sizes[d]:
                    raise ValueError(
                        f"Member '{mem_key}' has size {da.sizes[d]} on dim '{d}', "
                        f"expected {template.sizes[d]}"
                    )
    
            # Try to use an integer member index if possible, else keep as string
            try:
                mid = int(mem_key)
            except (TypeError, ValueError):
                mid = mem_key
    
            # Add member dimension and coords
            da = da.expand_dims({"member": [mid]})
            da = da.assign_coords(member=("member", [mid]))
            da = da.assign_coords(member_str=("member", [str(mem_key)]))
    
            das.append(da)
    
        # Concatenate along member
        da_model = xr.concat(das, dim="member")
    
        # If member is numeric, sort by member index
        if np.issubdtype(da_model["member"].dtype, np.number):
            order = np.argsort(da_model["member"].values)
            da_model = da_model.isel(member=order)
    
        return da_model

    def pool_members_to_samples(self, member_dict, sample_dim=None):
        """
        Pool member DataArrays along `sample_dim` into a single DataArray with
        a new 'sample' dimension, while keeping all other dims unchanged.
    
        Parameters
        ----------
        member_dict : dict
            Mapping like {'00': DataArray, '01': DataArray, ...}
        sample_dim : str or None, optional
            Name of the dimension to be treated as samples (e.g. 'years').
            If None, and the DataArray is 1D, that single dim is used.
    
        Returns
        -------
        pooled : xarray.DataArray
            dims: ('sample', *other_dims)
            coords:
              - 'sample': integer index
              - 'member': ('sample',) coord indicating which member each sample came from
              - all coords for other_dims copied from the first DataArray
        """
        if not member_dict:
            raise ValueError("member_dict is empty.")
    
        # Use first DataArray as template
        template = next(iter(member_dict.values()))
    
        # Infer sample_dim if not provided
        if sample_dim is None:
            if len(template.dims) == 1:
                sample_dim = template.dims[0]
            else:
                raise ValueError(
                    f"sample_dim is None, but DataArray has multiple dims {template.dims}. "
                    "Please specify sample_dim explicitly."
                )
    
        if sample_dim not in template.dims:
            raise ValueError(
                f"sample_dim='{sample_dim}' not found in dims {template.dims}"
            )
    
        # All dims, with sample_dim excluded
        other_dims = [d for d in template.dims if d != sample_dim]
    
        # Check that non-sample dims are consistent across members
        template_other_sizes = {d: template.sizes[d] for d in other_dims}
        for mem, da in member_dict.items():
            if set(da.dims) != set(template.dims):
                raise ValueError(f"Member {mem} has different dims: {da.dims}")
            for d in other_dims:
                if da.sizes[d] != template_other_sizes[d]:
                    raise ValueError(
                        f"Member {mem} has different size along dim '{d}': "
                        f"{da.sizes[d]} vs {template_other_sizes[d]}"
                    )
    
        vals = []
        members = []
    
        for mem, da in member_dict.items():
            # Put sample_dim first to make concatenation easy
            da_t = da.transpose(sample_dim, *other_dims)
            v = da_t.values  # shape: (n_sample, *other_shapes)
            n = v.shape[0]
            vals.append(v)
            members.extend([mem] * n)
    
        # Concatenate along sample axis (axis=0)
        all_vals = np.concatenate(vals, axis=0)
    
        # Build coords for other dims from template
        other_coords = {
            d: template.coords[d] for d in other_dims if d in template.coords
        }
    
        pooled = xr.DataArray(
            all_vals,
            dims=("sample", *other_dims),
            coords={
                "sample": np.arange(all_vals.shape[0]),
                "member": ("sample", members),
                **other_coords,
            },
            name=template.name,
            attrs=template.attrs,
        )
        return pooled

    def validate_constant_observation(
        self,
        obs_dict,
        ref_group="hist",
        ref_member="00",
        sample_dim=None,
        use_allclose=False,
        rtol=1e-8,
        atol=0.0,
        pool_ensemble=True,
        
    ):
        """
        Validate that observation DataArrays are identical across all periods/members,
        then return a DataArray shaped consistently with the model by pooling along
        the specified sample dimension.
    
        Parameters
        ----------
        self : object
            Object that provides `pool_members_to_samples` as a method.
            If you use this as a standalone function, remove `self` and call the
            function version of `pool_members_to_samples` instead.
        obs_dict : dict
            Nested dict: obs_dict[period][member] -> xarray.DataArray.
        ref_group : str, optional
            Period key used for the reference, e.g. 'hist'.
        ref_member : str, optional
            Member key used for the reference, e.g. '00'.
        sample_dim : str or None, optional
            Dimension to treat as samples when pooling.
            - If None and the DataArray is 1D, the single dim is inferred.
            - If None and the DataArray has multiple dims, `pool_members_to_samples`
              will raise a ValueError asking for an explicit `sample_dim`.
        use_allclose : bool, optional
            If True, use np.allclose instead of np.array_equal.
        rtol, atol : float, optional
            Tolerances passed to np.allclose if use_allclose=True.
    
        Returns
        -------
        out : xarray.DataArray
            Observation DataArray with dims ('sample', *other_dims), created via
            `pool_members_to_samples` so that its shape/dim ordering is consistent
            with the pooled model arrays.
    
        Raises
        ------
        KeyError
            If the reference group/member is missing.
        ValueError
            If any entry differs from the reference.
        """
        # --- get reference ---
        if ref_group not in obs_dict:
            raise KeyError(
                f"Reference group '{ref_group}' not found in obs_dict keys: "
                f"{list(obs_dict.keys())}"
            )
    
        if ref_member not in obs_dict[ref_group]:
            raise KeyError(
                f"Reference member '{ref_member}' not found in obs_dict['{ref_group}'] "
                f"keys: {list(obs_dict[ref_group].keys())}"
            )
    
        ref = obs_dict[ref_group][ref_member]
    
        # --- validate all obs equal to ref (by values) ---
        for period, members in obs_dict.items():
            for member, dx in members.items():
                if use_allclose:
                    same = np.allclose(
                        dx.values, ref.values,
                        rtol=rtol, atol=atol, equal_nan=True
                    )
                else:
                    same = np.array_equal(dx.values, ref.values)
    
                if not same:
                    raise ValueError(
                        f"Observation `obs_dict` differs at period={period}, member={member} "
                        f"(values are not identical to reference {ref_group}/{ref_member})."
                    )
    
        # --- make sure returned obs has pooled/sample shape consistent with model ---
        # Use a fake 1-member dict so we can reuse pool_members_to_samples logic
        if pool_ensemble:
            out = self.pool_members_to_samples({"00": ref}, sample_dim=sample_dim)
        else:
            out = ref 
    
        return out