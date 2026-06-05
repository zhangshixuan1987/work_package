from collections.abc import MutableMapping
from copy import deepcopy
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from logger import _setup_child_logger
from utils import (
    get_highlight_models,
    shift_row_to_bottom,
)

logger = _setup_child_logger(__name__)


class ClimMetricsMerger:
    def __init__(
        self,
        parameter,
        reference_alias: Optional[Dict[str, str]] = None,
        units_all: Optional[Dict[str, str]] = None,
    ):
        self.model_names = parameter.get("test_highlight_models") or parameter["test_name"]
        self.ref_group = parameter["ref_group"]
        self.test_group = parameter["test_group"]
        self.mean_group1_name = parameter.get("mean_group1_name", self.ref_group)
        self.mean_group2_name = parameter.get("mean_group2_name", self.test_group)
        self.unit_check = parameter.get("unit_check", True)
        self.error_norm = parameter.get("error_norm", "default")
        self.test_model_only = parameter.get("test_model_only", False)
        self.show_mean_columns = parameter.get(
            "show_mean_columns",
            parameter.get("plot_mean_groups", True),
        )
        
        # set aliases/units, honoring user-provided dicts if given
        if reference_alias is None:
            self.reference_alias = {
                "ceres_ebaf_toa_v4.1": "ceres_ebaf_v4_1",
                "ceres_ebaf_toa_v4.0": "ceres_ebaf_v4_0",
                "ceres_ebaf_toa_v2.8": "ceres_ebaf_v2_8",
                "ceres_ebaf_surface_v4.1": "ceres_ebaf_v4_1",
                "ceres_ebaf_surface_v4.0": "ceres_ebaf_v4_0",
                "ceres_ebaf_surface_v2.8": "ceres_ebaf_v2_8",
                "CERES-EBAF-4-1": "ceres_ebaf_v4_1",
                "CERES-EBAF-4-0": "ceres_ebaf_v4_0",
                "CERES-EBAF-2-8": "ceres_ebaf_v2_8",
                "GPCP_v2.3": "GPCP_v2_3",
                "GPCP_v2.2": "GPCP_v2_2",
                "GPCP_v3.2": "GPCP_v3_2",
                "GPCP-2-3": "GPCP_v2_3",
                "GPCP-2-2": "GPCP_v2_2",
                "GPCP-3-2": "GPCP_v3_2",
                "NOAA_20C": "NOAA-20C",
                "ERA-INT": "ERA-Interim",
                "ERA-5": "ERA5",
            }
        else:
            self.reference_alias = reference_alias

        if units_all is None:
            self.units_all = {
                "prw": "[kg m$^{-2}$]",
                "pr": "[mm d$^{-1}$]",
                "prsn": "[mm d$^{-1}$]",
                "prc": "[mm d$^{-1}$]",
                "hfls": "[W m$^{-2}$]",
                "hfss": "[W m$^{-2}$]",
                "clivi": "[kg $m^{-2}$]",
                "clwvi": "[kg $m^{-2}$]",
                "psl": "[Pa]",
                "rlds": "[W m$^{-2}$]",
                "rldscs": "[W $m^{-2}$]",
                "evspsbl": "[kg m$^{-2} s^{-1}$]",
                "rtmt": "[W m$^{-2}$]",
                "rsdt": "[W m$^{-2}$]",
                "rlus": "[W m$^{-2}$]",
                "rluscs": "[W m$^{-2}$]",
                "rlut": "[W m$^{-2}$]",
                "rlutcs": "[W m$^{-2}$]",
                "rsds": "[W m$^{-2}$]",
                "rsdscs": "[W m$^{-2}$]",
                "rstcre": "[W m$^{-2}$]",
                "rltcre": "[W m$^{-2}$]",
                "rsus": "[W m$^{-2}$]",
                "rsuscs": "[W m$^{-2}$]",
                "rsut": "[W m$^{-2}$]",
                "rsutcs": "[W m$^{-2}$]",
                "ts": "[K]",
                "tas": "[K]",
                "tauu": "[Pa]",
                "tauv": "[Pa]",
                "zg-500": "[m]",
                "ta-200": "[K]",
                "sfcWind": "[m s$^{-1}$]",
                "ta-850": "[K]",
                "ua-200": "[m s$^{-1}$]",
                "ua-850": "[m s$^{-1}$]",
                "va-200": "[m s$^{-1}$]",
                "va-850": "[m s$^{-1}$]",
                "uas": "[m s$^{-1}$]",
                "vas": "[m s$^{-1}$]",
                "tasmin": "[K]",
                "tasmax": "[K]",
                "clt": "[%]",
            }
        else:
            self.units_all = units_all

    def _check_references(
        self,
        data_dict: MutableMapping[str, Optional[List[str]]],
    ) -> MutableMapping[str, Optional[List[str]]]:

        for key, values in data_dict.items():
            if isinstance(values, list):
                data_dict[key] = [self.reference_alias.get(val, val) for val in values]
            elif values is not None:
                data_dict[key] = self.reference_alias.get(values, values)
            else:
                logger.warning(f"Reference for key '{key}' is None — skipping.")
                continue

        logger.debug(
            f"Checked references for {len(data_dict)} keys "
            f"({sum(v is None for v in data_dict.values())} None values skipped)"
        )

        return data_dict

    def _check_regions(self, data_lib, refr_lib):
        shared_regions = [
            region for region in data_lib.regions if region in refr_lib.regions
        ]

        for lib in [refr_lib, data_lib]:
            for stat in lib.df_dict:
                for season in lib.df_dict[stat]:
                    lib.df_dict[stat][season] = {
                        region: lib.df_dict[stat][season][region]
                        for region in shared_regions
                        if region in lib.df_dict[stat][season]
                    }

        data_lib.regions = shared_regions
        refr_lib.regions = shared_regions

        return data_lib, refr_lib

    @staticmethod
    def _prune_empty_dfs(lib):
        for stat in lib.df_dict:
            for season in lib.df_dict[stat]:
                lib.df_dict[stat][season] = {
                    region: df
                    for region, df in lib.df_dict[stat][season].items()
                    if isinstance(df, pd.DataFrame)
                    and not df.empty
                    and not df.isna().all().all()
                }
        return lib

    @staticmethod
    def _safe_merge_libs(lib1, lib2):
        """
        Merge two data libraries with nested dicts of DataFrames,
        gracefully handling missing or inconsistent keys.
        """
        merged = deepcopy(lib1)  # Avoid modifying original

        for stat in lib2.df_dict:
            if stat not in merged.df_dict:
                merged.df_dict[stat] = {}

            for season in lib2.df_dict[stat]:
                if season not in merged.df_dict[stat]:
                    merged.df_dict[stat][season] = {}

                for region, df2 in lib2.df_dict[stat][season].items():
                    df1 = merged.df_dict[stat][season].get(region)

                    valid_dfs = []
                    for df in (df1, df2):
                        if (
                            isinstance(df, pd.DataFrame)
                            and not df.empty
                            and not df.isna().all().all()
                        ):
                            df_clean = df.dropna(axis=1, how="all")
                            if not df_clean.empty and df_clean.shape[1] > 0:
                                valid_dfs.append(df_clean)

                    merged_df = (
                        pd.concat(valid_dfs, ignore_index=True, sort=False)
                        if valid_dfs
                        else pd.DataFrame()
                    )
                    merged.df_dict[stat][season][region] = merged_df

        return merged

    def _check_units(self, data_lib, verbose: bool = False):
        # Identify common variables and handle aliases like 'rt' or 'rmt'
        common_vars = [var for var in getattr(data_lib, "var_list", []) if var in self.units_all]
        if "rtmt" not in common_vars and any(
            var in getattr(data_lib, "var_list", []) for var in ["rt", "rmt"]
        ):
            common_vars.append("rtmt")

        # Collect units for these variables
        common_unts = [self.units_all[var] for var in common_vars if var in self.units_all]

        # Filter and correct reference list
        new_var_ref_dict = {}
        for var, ref in getattr(data_lib, "var_ref_dict", {}).items():
            if var in common_vars:
                new_var_ref_dict[var] = ref
            elif var in ["rt", "rmt"]:
                new_var_ref_dict["rtmt"] = ref
                if verbose:
                    logger.info(f"Alias {var} mapped to 'rtmt' in references.")

        if new_var_ref_dict:
            data_lib.var_ref_dict = self._check_references(new_var_ref_dict)

        # Clean DataFrames
        for stat, seasons in data_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in regions.items():
                    df = pd.DataFrame(df).copy()
                    # Handle aliases
                    if "rt" in df.columns:
                        df["rtmt"] = df["rt"]
                    elif "rmt" in df.columns:
                        df["rtmt"] = df["rmt"]

                    # Drop irrelevant variables
                    drop_cols = [
                        var for var in df.columns[3:] if var not in common_vars
                    ]
                    if drop_cols and verbose:
                        logger.info(
                            f"Dropping variables in {stat}/{season}/{region}: {drop_cols}"
                        )
                    df = df.drop(columns=drop_cols, errors="ignore")
                    data_lib.df_dict[stat][season][region] = df

        logger.debug(f"Setting data_lib.var_list={common_vars}")
        data_lib.var_list = common_vars
        data_lib.var_unit_list = common_unts

        return data_lib

    def _filter_regions(self, ref_lib, model_lib):
        # keep order consistent: (ref_lib, model_lib) in and out
        model_lib, ref_lib = self._check_regions(model_lib, ref_lib)
        return ref_lib, model_lib
    
    def _normalize_references(self, data_lib):
        ref_attr = None
        if hasattr(data_lib, "references"):
            ref_attr = "references"
        elif hasattr(data_lib, "References"):
            ref_attr = "References"

        if ref_attr is not None:
            refs = getattr(data_lib, ref_attr)
            if isinstance(refs, dict):
                setattr(data_lib, ref_attr, self._check_references(refs))

        return data_lib

    def _merge_and_standardize_units(self, ref_lib, model_lib):
        # Prune/normalize model lib
        cleaned_model_lib = self._prune_empty_dfs(model_lib)
        cleaned_model_lib = self._normalize_references(cleaned_model_lib)
        if hasattr(cleaned_model_lib, "var_list"):
            logger.debug(f"cleaned_model_lib.var_list: {cleaned_model_lib.var_list}")

        if self.test_model_only:
            ref_lib = None 
            
        # Prepare reference lib
        cleaned_ref_lib = ref_lib
        if ref_lib is not None:
            cleaned_ref_lib = self._prune_empty_dfs(ref_lib)
            if hasattr(cleaned_ref_lib, "var_list"):
                logger.debug(f"cleaned_ref_lib.var_list: {cleaned_ref_lib.var_list}")
                var_set_model = set(getattr(cleaned_model_lib, "var_list", []))
                var_set_ref = set(getattr(cleaned_ref_lib, "var_list", []))
                logger.debug(
                    f"Var list sizes - model: {len(var_set_model)}, ref: {len(var_set_ref)}"
                )
                logger.debug(
                    f"Var list diff (model - ref): {var_set_model - var_set_ref}"
                )
                logger.debug(
                    f"Var list diff (ref - model): {var_set_ref - var_set_model}"
                )
            if self.unit_check:
                cleaned_ref_lib = self._check_units(cleaned_ref_lib)

            # normalize regions (only keep shared)
            cleaned_ref_lib, cleaned_model_lib = self._filter_regions(cleaned_ref_lib, cleaned_model_lib)

            merged_lib = self._safe_merge_libs(cleaned_ref_lib, cleaned_model_lib)
        else:
            merged_lib = cleaned_model_lib

        if hasattr(merged_lib, "var_list"):
            logger.debug(f"merged_lib.var_list: {getattr(merged_lib, 'var_list', [])}")

        # Standardize units after merging
        if self.unit_check:
            merged_lib = self._check_units(merged_lib)

        if hasattr(merged_lib, "var_list"):
            logger.debug(
                f"Post-unit-check merged_lib.var_list: {getattr(merged_lib, 'var_list', [])}"
            )
        return cleaned_ref_lib, cleaned_model_lib, merged_lib

    def _highlight_and_sort_models(self, merged_lib):
        if merged_lib is None:
            raise ValueError("merged_lib is None")

        all_highlights: List[str] = []
        for stat, seasons in merged_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in regions.items():
                    df = pd.DataFrame(df)
                    highlight_models = get_highlight_models(
                        df.get("model", []), self.model_names
                    )
                    for model in highlight_models:
                        if model not in all_highlights:
                            all_highlights.append(model)
                    for model in highlight_models:
                        for idx in df[df["model"] == model].index:
                            df = shift_row_to_bottom(df, idx)
                    merged_lib.df_dict[stat][season][region] = df.fillna(np.nan)

        return all_highlights, merged_lib

    def _add_group_means(self, data_lib, mean_lib, mean_name, overwrite=True):
        """
        For each (stat, season, region):
          - Compute a SINGLE global mean across all numeric columns from mean_lib (excluding its mean rows).
          - Append that one-row mean (model == mean_name) to data_lib.
          - If overwrite=True, remove any existing mean row in data_lib before appending.
          - If overwrite=False and a mean row already exists in data_lib, leave it as-is (skip recompute).

        Returns
        -------
        data_lib, mean_model_list
            mean_model_list includes all contributing models from mean_lib plus the mean_name itself.
        """
        contributing = set()

        for stat, seasons in data_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in list(regions.items()):
                    df = pd.DataFrame(df)

                    mf_raw = mean_lib.df_dict.get(stat, {}).get(season, {}).get(region, None)
                    mf = pd.DataFrame(mf_raw) if mf_raw is not None else pd.DataFrame()

                    # Must have 'model' on both sides to proceed
                    if "model" not in df.columns or "model" not in mf.columns:
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    df["model"] = df["model"].astype(str)
                    mf["model"] = mf["model"].astype(str)

                    # Respect overwrite flag
                    if not overwrite and (df["model"] == mean_name).any():
                        data_lib.df_dict[stat][season][region] = df
                        continue
                    if overwrite and (df["model"] == mean_name).any():
                        df = df[df["model"] != mean_name]

                    # Compute global mean from mean_lib (no grouping)
                    base = mf[mf["model"] != mean_name]
                    if base.empty:
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    # track contributors for this slice
                    contributing.update(base["model"].astype(str).unique())

                    num_cols = base.select_dtypes(include=[np.number]).columns
                    if len(num_cols) == 0:
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    mean_vals = base[num_cols].mean(skipna=True)
                    mean_df = pd.DataFrame([mean_vals])

                    # label and align to data_lib's schema
                    mean_df["model"] = mean_name
                    cols = ["model"] + [c for c in df.columns if c != "model"]
                    mean_df = mean_df.reindex(columns=cols, fill_value=np.nan)

                    # append the single mean row
                    df_out = pd.concat([df, mean_df], ignore_index=True, sort=False)
                    data_lib.df_dict[stat][season][region] = df_out

        # add mean_name itself to the list of contributors
        contributing.add(mean_name)

        mean_model_list = sorted(contributing)
        return data_lib, mean_model_list

    def merge(
        self,
        ref_lib,
        model_lib,
        clim_vars: Optional[List[str]] = None,
        clim_regions: Optional[List[str]] = None,
    ):
        # 1) merge + unit standardization
        ref_lib, model_lib, merged_lib = self._merge_and_standardize_units(ref_lib, model_lib)

        # 2) highlight and sort (use the same var name consistently)
        e3sm_model_list, merged_lib = self._highlight_and_sort_models(merged_lib)

        # 3) Variables (preserve original behavior unless filters provided)
        var_list = list(getattr(merged_lib, "var_list", []))
        var_unit_list = list(getattr(merged_lib, "var_unit_list", []))
        if clim_vars is not None:
            name_to_unit = dict(zip(getattr(merged_lib, "var_list", []),
                                    getattr(merged_lib, "var_unit_list", [])))
            missing = [v for v in clim_vars if v not in name_to_unit]
            if missing:
                logger.warning("[mean_climate] Requested variables not found and will be skipped: %s", missing)
            var_list = [v for v in clim_vars if v in name_to_unit]
            var_unit_list = [name_to_unit[v] for v in var_list]

        # 4) Regions (preserve order)
        regions = list(getattr(merged_lib, "regions", []))
        if clim_regions is not None:
            missing_r = [r for r in clim_regions if r not in getattr(merged_lib, "regions", [])]
            if missing_r:
                logger.warning("[mean_climate] Requested regions not found and will be skipped: %s", missing_r)
            regions = [r for r in clim_regions if r in getattr(merged_lib, "regions", [])]

        # 5) data used for normalization for metrics plot
        norm_lib = ref_lib if self.error_norm == "reference" else merged_lib

        # 6) Append group means computed from ref_lib / model_lib into merged_lib
        mean_model_list: List[str] = []   # names of the mean entries to highlight (e.g., ["CMIP (mean)", "E3SMv3-LE (mean)"])
        ref_model_list: List[str] = []    # contributors used for reference mean + mean_name (per your _add_group_means)
        test_model_list: List[str] = []   # contributors used for test mean + mean_name

        if self.show_mean_columns and ref_lib is not None and getattr(self, "mean_group1_name", None):
            merged_lib, ref_model_list = self._add_group_means(
                merged_lib, ref_lib, mean_name=self.mean_group1_name, overwrite=True
            )
            if self.mean_group1_name in ref_model_list:
                mean_model_list.append(self.mean_group1_name)

        if self.show_mean_columns and model_lib is not None and getattr(self, "mean_group2_name", None):
            merged_lib, test_model_list = self._add_group_means(
                merged_lib, model_lib, mean_name=self.mean_group2_name, overwrite=True
            )
            if self.mean_group2_name in test_model_list:
                mean_model_list.append(self.mean_group2_name)

        # (optional) de-dup mean list in case names repeat
        if mean_model_list:
            mean_model_list = list(dict.fromkeys(mean_model_list))

        return (
            merged_lib,
            norm_lib,
            var_list,
            var_unit_list,
            regions,
            ref_model_list,
            test_model_list,
            e3sm_model_list,
            mean_model_list,
        )
