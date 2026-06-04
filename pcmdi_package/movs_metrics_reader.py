import glob
import json
import os
import re

import numpy as np
import pandas as pd
from pcmdi_metrics.utils import sort_human

from typing import Optional, Union, List   

from logger import _setup_child_logger
from utils import (
    find_latest_file_list,
    get_highlight_models,
    shift_row_to_bottom,
)

logger = _setup_child_logger(__name__)


class MoVsMetricsReader:
    def __init__(self, parameter):
        self.ref_path = parameter["ref_path"]
        self.ref_name = parameter["ref_name"]
        self.test_path = parameter["test_path"]
        
        self.mips = parameter["test_mip"]
        self.tests = parameter["test_name"]
        self.caseids = parameter["test_id"]

        self.group = parameter["movs_group"]
        self.movs_mode = parameter["movs_mode"]
        self.movc_obs = parameter["movc_obs"]
        self.mova_obs = parameter["mova_obs"]

        self.diag_vars = parameter["diag_vars"]
        self.var_pattern = re.compile(r"var_mode_(\w+)\.EOF\d+\..*\.json$")
        self.time_pattern = re.compile(r"\.v(\d{8})\.json$")

        self.exclude_models = parameter.get("exclude_models", [])
        self.error_norm = parameter.get("error_norm", "default")
        self.test_model_only = parameter.get("test_model_only", False)
        
        self.ref_group = parameter["ref_group"]
        self.test_group = parameter["test_group"]
        self.mean_group1_name = parameter.get("mean_group1_name", self.ref_group)
        self.mean_group2_name = parameter.get("mean_group2_name", self.test_group)
        self.test_combined = parameter.get("test_combined", False)

        self.ref_key = self.ref_model = self.ref_version = None
        if self.ref_name is not None:
            if not isinstance(self.ref_name, str):
                raise TypeError("ref_name must be a string if provided")
            name = self.ref_name.strip()
            parts = name.split(".")
            if len(parts) != 3 or any(p.strip() == "" for p in parts):
                raise ValueError(
                    f"Invalid ref_name format '{self.ref_name}'. Expected '<group>.<model>.<version>'."
                )
            self.ref_key, self.ref_model, self.ref_version = (p.strip() for p in parts)
            
        if self.test_combined: 
            self.test_key = self.test_model = self.test_version = None
            self.test_name = parameter["test_name"]
            if self.test_name is not None:
                if not isinstance(self.test_name, str):
                    raise TypeError("test_name must be a string if provided")
                name = self.test_name.strip()
                parts = name.split(".")
                if len(parts) != 3 or any(p.strip() == "" for p in parts):
                    raise ValueError(
                        f"Invalid ref_name format '{self.test_name}'. Expected '<group>.<model>.<version>'."
                    )
                self.test_key, self.test_model, self.test_version = (p.strip() for p in parts)

    @staticmethod
    def _prune_empty_dfs(lib):
        """
        If lib is {key: DataFrame}, drop empty/NaN-only frames.
        If lib is nested dict (e.g., mode -> model -> runs...), drop keys whose values are falsy/empty dicts.
        """
        if lib is None:
            return None
        # Dict of DataFrames
        if all(isinstance(v, (pd.DataFrame, type(None), dict)) for v in lib.values()):
            pruned = {}
            for k, v in lib.items():
                if isinstance(v, pd.DataFrame):
                    if not v.empty and not v.isna().all().all():
                        pruned[k] = v
                elif isinstance(v, dict):
                    if v:  # keep non-empty
                        pruned[k] = v
                elif v is not None:
                    pruned[k] = v
            return pruned
        return lib

    def _highlight_and_sort_models(self, merged_lib):
        if merged_lib is None:
            raise ValueError("merged_lib is None")

        all_highlights = set()
        for stat, df in list(merged_lib.items()):
            df = pd.DataFrame(df)
            highlight_models = get_highlight_models(df.get("model", []), self.tests)
            all_highlights.update(highlight_models)
            for model in highlight_models:
                for idx in df[df["model"] == model].index:
                    df = shift_row_to_bottom(df, idx)
            merged_lib[stat] = df.fillna(np.nan)

        return sorted(all_highlights), merged_lib
    
    def _exclude_models(
            self,
            data_lib: dict,
            model_list: Optional[Union[List[str], str]] = None,
            *,
            verbose: bool = False
        ):
        if not model_list:
            return data_lib
        if isinstance(model_list, str):
            model_list = [model_list]
        model_set = {str(m) for m in model_list}

        for stat, table in list(data_lib.items()):
            if not isinstance(table, (pd.DataFrame, dict, list, tuple)):
                continue
            df = table if isinstance(table, pd.DataFrame) else pd.DataFrame(table)
            if "model" not in df.columns:
                if verbose:
                    logger.info(f"[exclude_models] No 'model' column in '{stat}'; skipping.")
                data_lib[stat] = df
                continue
            before = len(df)
            df = df[~df["model"].astype(str).isin(model_set)].reset_index(drop=True)
            if verbose:
                removed = before - len(df)
                if removed:
                    logger.info(f"[exclude_models] ({stat}) removed {removed} rows.")
            data_lib[stat] = df
        return data_lib

    def _get_ref_files(self, ref_path, ref_key, ref_model, ref_version):
        current_dir: str = os.path.abspath(os.getcwd())
        pattern: str = os.path.join(
            ref_path, ref_key, ref_model, ref_version, "*/*/var_mode_*.json"
        )
        logger.info(f"From {current_dir}, checking for reference files matching {pattern}")
        matching_files = glob.glob(pattern)
        logger.info(f"Found {len(matching_files)} matching files")
        for file_name in matching_files:
            logger.debug(f"  - {file_name}")
        if not matching_files:
            logger.warning(
                f"[MoVsMetricsReader]: No matching files found for pattern: {pattern}. "
                "Ensure the path and pattern are correct."
            )
        return matching_files

    def _load_movs_files(self, file_lists):
        json_lib = {}
        for mode in self.movs_mode:
            eof = {"PSA1": "EOF2", "NPO": "EOF2", "NPGO": "EOF2", "PSA2": "EOF3"}.get(mode, "EOF1")
            for json_file in file_lists:
                if mode in json_file and eof in json_file:
                    try:
                        with open(json_file, "r") as fj:
                            data = json.load(fj)
                            json_lib[mode] = data.get("RESULTS", {})
                    except (FileNotFoundError, json.JSONDecodeError) as e:
                        logger.warning(f"[MoVsMetricsReader]: Could not load {json_file}: {e}")
                    break
        return json_lib
    
    def _pick_case(self, d: dict, key: str):
        # try exact, UPPER, then lower
        if key in d: return d[key]
        ku = key.upper()
        if ku in d: return d[ku]
        kl = key.lower()
        if kl in d: return d[kl]
        return None
    
    def _movs_dict_to_df(self, movs_dict, group, stat):
        full_list = ["NAM", "NAO", "PNA", "NPO", "SAM", "PSA1", "PSA2", "PDO", "NPGO", "AMO"]
        MODE_TO_EOF = {
            "NAM": "eof1",
            "NAO": "eof1",
            "PNA": "eof2",
            "NPO": "eof3",
            "SAM": "eof1",
            "PSA1": "eof2",
            "PSA2": "eof3",
            "PDO": "eof1",
            "NPGO": "eof2",
            "AMO": "eof1",
        }
        MODE_TO_SEASONS = {"PDO": ["monthly"], "NPGO": ["monthly"], "AMO": ["yearly"]}
        DEFAULT_SEASONS = ["DJF", "MAM", "JJA", "SON"]
        season_label_map = {"yearly": "ANN", "monthly": "MON"}

        mode_list = [m for m in self.movs_mode if m in movs_dict]
        if not mode_list:
            if movs_dict:
                fallback = next(iter(movs_dict))
                mode_list = [fallback]
                logger.warning(
                    f"[MoVsMetricsReader]: No requested modes found; falling back to first available: '{fallback}'."
                )
            else:
                logger.warning("[MoVsMetricsReader]: movs_dict is empty — nothing to process.")
                return pd.DataFrame(columns=["model", "num_runs"]), []

        models = sorted((movs_dict.get(mode_list[0]) or {}).keys())
        df = pd.DataFrame({"model": models, "num_runs": np.nan})
        mode_season_list = []

        for mode in self.movs_mode:
            seasons = MODE_TO_SEASONS.get(mode, DEFAULT_SEASONS)
            group_for_mode = group if group == "cbf" else MODE_TO_EOF.get(mode, "eof1")

            for season in seasons:
                seastr = season_label_map.get(season.lower(), season.upper())
                col_name = f"{mode}_{seastr}"
                df[col_name] = np.nan
                mode_season_list.append(col_name)

                for idx, model in enumerate(models):
                    value = np.nan
                    num_runs = 0

                    if mode in movs_dict and model in movs_dict[mode]:
                        runs = sort_human(list(movs_dict[mode][model].keys()))
                        stat_values = []
                        for run in runs:
                            try:
                                scope = movs_dict[mode][model][run]["defaultReference"][mode][season]
                                grp = self._pick_case(scope, group_for_mode)
                                if grp is None:
                                    continue
                                run_stat = grp.get(stat, np.nan)
                                stat_values.append(run_stat)
                            except KeyError:
                                continue

                        if stat_values:
                            # keep only finite, non-null numbers
                            finite_vals = [x for x in stat_values if pd.notna(x) and np.isfinite(x)]
                            if finite_vals:
                                value = float(np.mean(finite_vals))
                                num_runs = len(finite_vals)
        
                    df.at[idx, col_name] = value
                    if np.isnan(df.at[idx, "num_runs"]):
                        df.at[idx, "num_runs"] = num_runs
                    elif num_runs > 0:
                        df.at[idx, "num_runs"] = max(df.at[idx, "num_runs"], num_runs)

        return df, mode_season_list
    
    def _add_group_means(self, merged_lib: dict, base_lib: dict, mean_name: str, overwrite: bool = True):
        """
        Append a single mean row (model == mean_name) to merged_lib[stat], computed
        over numeric columns from base_lib[stat]. Both args are dict[str, DataFrame].
        Returns (merged_lib, contributors_list).
        """
        contributors = set()

        for stat, mdf in merged_lib.items():
            if stat not in base_lib:
                continue
            mdf = pd.DataFrame(mdf)
            bdf = pd.DataFrame(base_lib[stat])

            if "model" not in mdf.columns or "model" not in bdf.columns:
                merged_lib[stat] = mdf
                continue

            base = bdf[bdf["model"].astype(str) != str(mean_name)]
            if base.empty:
                merged_lib[stat] = mdf
                continue

            contributors.update(base["model"].astype(str).unique())

            num_cols = base.select_dtypes(include=[np.number]).columns
            if len(num_cols) == 0:
                merged_lib[stat] = mdf
                continue

            mean_vals = base[num_cols].mean(skipna=True)
            mean_row = pd.DataFrame([mean_vals])
            mean_row["model"] = mean_name

            cols = ["model"] + [c for c in mdf.columns if c != "model"]
            mean_row = mean_row.reindex(columns=cols, fill_value=np.nan)

            if overwrite and "model" in mdf.columns:
                mdf = mdf[mdf["model"].astype(str) != str(mean_name)]

            merged_lib[stat] = pd.concat([mdf, mean_row], ignore_index=True, sort=False)

        if mean_name:
            contributors.add(mean_name)
        return merged_lib, sorted(contributors)
    
    def collect_metrics(self):
        """
        Build merged_lib using the same logic as the old working version:
          1. Start from reference (cmip_lib / ref_lib) per stat.
          2. Append each test model one by one.
          3. Return merged_lib and related info.
        """
        # --- Read reference and test data ---
        ref_lib = {}
        if self.ref_name is not None:
            ref_files = self._get_ref_files(self.ref_path, self.ref_key, self.ref_model, self.ref_version)
            if not ref_files:
                raise FileNotFoundError("No reference MoVs metric files found.")
            ref_lib = self._load_movs_files(ref_files)
            ref_lib = self._prune_empty_dfs(ref_lib)
            
        if self.test_combined and self.test_name is not None:
            test_files = self._get_ref_files(self.test_path, self.test_key, self.test_model, self.test_version)
            if not test_files:
                raise FileNotFoundError("No combined model MoVs metric files found.")
            test_lib = self._load_movs_files(test_files)
            test_lib = self._prune_empty_dfs(test_lib)
            
        merged_lib = {}
        ref_df_lib = {}
        test_df_lib = {}
        mode_season_list = []
            
        # --- Build per-stat merged tables ---
        for stat, _ in self.diag_vars.items():
            # Reference first
            if ref_lib:
                ref_df, mode_season_ref = self._movs_dict_to_df(ref_lib, self.group, stat)
                merge_df = ref_df.copy(deep=True)
                ref_df_lib[stat] = ref_df
                mode_season_list.extend(mode_season_ref)
            else:
                merge_df = pd.DataFrame(columns=["model", "num_runs"])
                ref_df_lib[stat] = pd.DataFrame(columns=["model", "num_runs"])

            if self.test_combined: 
                test_df, mode_season_tst = self._movs_dict_to_df(test_lib, self.group, stat)
                merge_df = pd.concat([merge_df, test_df], ignore_index=True)
                
                # accumulate mode/season names
                for m in mode_season_tst:
                    if m not in mode_season_list:
                        mode_season_list.append(m)
                        
                test_df_lib[stat] = test_df

            else:
                test_all = []
                # Loop through each test model exactly as old version
                for test_mip, model_name, case_id in zip(self.mips, self.tests, self.caseids):
                    model_path = self.test_path.replace("put_model_here", model_name)
                    all_model_files = []
                    for mode in self.movs_mode:
                        obs_str = self.movc_obs if mode in ["PDO", "NPGO", "AMO"] else self.mova_obs
                        model_files = find_latest_file_list(
                            path=f"{model_path}/{mode}/{obs_str}",
                            file_pattern=f"var_mode_*{test_mip}*{case_id}.json",
                            var_pattern=self.var_pattern,
                            time_pattern=self.time_pattern,
                        )
                        if model_files:
                            all_model_files.extend(model_files)

                    all_model_files = list(dict.fromkeys(all_model_files))
                    if not all_model_files or not os.path.exists(all_model_files[0]):
                        raise FileNotFoundError(f"No Synthetic MoVs Metrics Data For {model_name}, aborting.")

                    logger.info(f"Found Synthetic MoVs Metrics for {model_name}, reading…")
                    model_lib = self._load_movs_files(all_model_files)  # {mode: RESULTS}

                    # Normalize to per-model structure (same as old code)
                    model_lib = {
                        mode: {model_name: next(iter(model_data.values()))}
                        for mode, model_data in model_lib.items()
                        if model_data
                    }

                    model_df, mode_season_tst = self._movs_dict_to_df(model_lib, self.group, stat)
                    merge_df = pd.concat([merge_df, model_df], ignore_index=True)
                    test_all.append(model_df)

                    # accumulate mode/season names
                    for m in mode_season_tst:
                        if m not in mode_season_list:
                            mode_season_list.append(m)

                test_df_lib[stat] = pd.concat(test_all, ignore_index=True) if test_all else pd.DataFrame()
                
            merged_lib[stat] = merge_df.reindex(columns=["model"] + [c for c in merge_df.columns if c != "model"])

        # --- Exclusions (same as before) ---
        if self.exclude_models:
            merged_lib = self._exclude_models(merged_lib, self.exclude_models)
            if ref_df_lib:
                ref_df_lib = self._exclude_models(ref_df_lib, self.exclude_models)
            if test_df_lib:
                test_df_lib = self._exclude_models(test_df_lib, self.exclude_models)

        # --- Normalization source ---
        norm_lib = ref_df_lib if (self.error_norm == "reference" and ref_df_lib) else merged_lib

        # --- Highlight and reorder models ---
        e3sm_model_list, merged_lib = self._highlight_and_sort_models(merged_lib)

        # --- Group means ---
        mean_model_list, ref_model_list, test_model_list = [], [], []
        if ref_df_lib and getattr(self, "mean_group1_name", None):
            merged_lib, ref_model_list = self._add_group_means(merged_lib, ref_df_lib, self.mean_group1_name, True)
            if self.mean_group1_name in ref_model_list:
                mean_model_list.append(self.mean_group1_name)
        if test_df_lib and getattr(self, "mean_group2_name", None):
            merged_lib, test_model_list = self._add_group_means(merged_lib, test_df_lib, self.mean_group2_name, True)
            if self.mean_group2_name in test_model_list:
                mean_model_list.append(self.mean_group2_name)

        if mean_model_list:
            mean_model_list = list(dict.fromkeys(mean_model_list))
        if mode_season_list:
            mode_season_list = list(dict.fromkeys(mode_season_list))

        return (
            merged_lib,
            norm_lib,
            mode_season_list,
            ref_model_list,
            test_model_list,
            e3sm_model_list,
            mean_model_list,
        )
