import glob
import json
import os
import re

import numpy as np
import pandas as pd
from pcmdi_metrics.graphics import Metrics

from typing import Optional, Dict, List, Union

from logger import _setup_child_logger

from utils import find_latest_file_list

logger = _setup_child_logger(__name__)


class ClimMetricsReader:
    def __init__(self, parameter):
        """
        Initialize the climate metrics collector.

        Args:
            parameter (dict): Contains path, model info, and identifiers.
        """
        self.ref_path = parameter["ref_path"]
        self.ref_name = parameter["ref_name"]
        self.test_path = parameter["test_path"]
        self.mips = parameter["test_mip"]
        self.tests = parameter["test_name"]
        self.caseids = parameter["test_id"]
        self.var_pattern = re.compile(r"^([A-Za-z0-9\-]+)\.")
        self.time_pattern = re.compile(r"\.v(\d{8})\.json$")
        self.exclude_models = parameter.get("exclude_models", [])
        self.exclude_vars = parameter.get("exclude_vars", {})
        self.test_combined = parameter.get("test_combined", False)

        if self.ref_name is not None:
            if not isinstance(self.ref_name, str):
                raise TypeError("ref_name must be a string if provided")
            name = self.ref_name.strip()
            parts = name.split(".")
            if len(parts) != 3 or any(p.strip() == "" for p in parts):
                raise ValueError(
                    f"Invalid ref_name format '{self.ref_name}'. "
                    "Expected '<group>.<model>.<version>'."
                )
            self.ref_group, self.ref_model, self.ref_version = (p.strip() for p in parts)
            
        if self.test_combined: 
            self.test_name = parameter["test_name"]
            if self.test_name is not None:
                if not isinstance(self.test_name, str):
                    raise TypeError("test_name must be a string if provided")
                name = self.test_name.strip()
                parts = name.split(".")
                if len(parts) != 3 or any(p.strip() == "" for p in parts):
                    raise ValueError(
                        f"Invalid test_name format '{self.test_name}'. "
                        "Expected '<group>.<model>.<version>'."
                    )
                self.test_group, self.test_model, self.test_version = (p.strip() for p in parts)

    def _load_clim_metrics_from_files(self, file_paths):
        """
        Loads and processes synthetic climate metric data from JSON files.

        Parameters:
            file_paths (list): List of file paths to load.

        Returns:
            Metrics: Processed Metrics object.
        """
        logger.info(f"file_paths= ({len(file_paths)})")
        for i, fp in enumerate(file_paths):
            logger.info(f"{i}. {fp}")
        """
        FAILURE:

        list(results_dict_var["RESULTS"][model_list[0]]["default"][run_list[0]].keys())
        IndexError: list index out of range

        cat /lcrc/group/e3sm/public_html/diagnostic_output/ac.forsyth2/zppy_pr719_output/unique_id_21/v3.LR.amip_0101/pcmdi_diags/model_vs_obs/metrics_data/mean_climate/rlus.2.5x2.5.e3sm.amip.v3-LR_0101.v20250725.json

        "RESULTS": {
            "v3-LR": {
                "default": {
                    "source": "ceres_ebaf_v4_1"
                }
            }
        },

        SYNTHETIC PLOTS ERROR #1: "source" is supposed to be a dictionary itself, even though mean_climate job completed successfully!
        """
        lib = Metrics(file_paths)
        
        return lib

    def _load_combined_metrics(
            self,
            ref_path,
            ref_group, 
            ref_model, 
            ref_version 
    ):
        ref_dir = os.path.join(
            ref_path,
            ref_group, 
            ref_model, 
            ref_version 
        )

        ref_files = sorted(
            glob.glob(os.path.join(ref_dir, f"*.{ref_version}.json"))
        )
        if not ref_files:
            raise FileNotFoundError(f"No reference metrics found in: {ref_dir}")

        logger.info(f"Loading reference metrics from {len(ref_files)} files...")
        ref_lib = self._load_clim_metrics_from_files(ref_files)
        return ref_lib

    def _process_test_model(self, mip_name, test_name, case_id):
        test_key = mip_name.split(".")[1] if "." in mip_name else mip_name
        test_path = self.test_path.replace("put_model_here", test_name)

        model_files = find_latest_file_list(
            path=test_path,
            file_pattern=f"*.{case_id}.json",
            var_pattern=self.var_pattern,
            time_pattern=self.time_pattern,
        )

        if not model_files or not os.path.exists(model_files[0]):
            raise FileNotFoundError(
                f"No synthetic mean climate metrics found for model: {test_name}"
            )

        logger.info(
            f"Reading metrics for model: {test_name} from {len(model_files)} file(s)..."
        )

        valid_model_files = []

        for file_path in model_files:
            try:
                with open(file_path, "r") as f:
                    data = json.load(f)

                results = data.get("RESULTS", {})
                modified = False

                for model, model_data in results.items():
                    if test_key in model_data:
                        model_data["default"] = model_data.pop(test_key)
                        modified = True

                if modified:
                    with open(file_path, "w") as f:
                        json.dump(data, f, indent=2)
                    logger.info(f"Updated file: {file_path}")

                valid_model_files.append(file_path)

            except (FileNotFoundError, json.JSONDecodeError) as e:
                logger.info(f"Warning: Could not load {file_path}: {e}")
                
        # after the file loop
        if not valid_model_files:
            # fall back to original list if none were writable/parsable but still exist
            valid_model_files = [fp for fp in model_files if os.path.exists(fp)]
            if not valid_model_files:
                raise RuntimeError(f"No usable metric files after validation for model: {test_name}")

        # Load metrics from valid files
        model_lib = self._load_clim_metrics_from_files(valid_model_files)
        
        # Standardize model name in metric DataFrames
        for stat, seasons in model_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in regions.items():
                    df = pd.DataFrame(df)
                    if "model" in df.columns:
                        # ensure string type, then standardize
                        df["model"] = str(test_name)
                    model_lib.df_dict[stat][season][region] = df
                    
        return model_lib
    
    def _check_badvals(
            self, 
            data_lib, 
            var_model: Optional[Dict[str, List[str]]] = None, 
            *, 
            verbose: bool = False
        ):
        """
        For each DataFrame in data_lib.df_dict[stat][season][region], set selected
        variable columns to NaN for specified models.
        """
        if var_model is None:
            var_model = {
                "E3SM-1-0":     ["ta-850"],
                "E3SM-1-1-ECA": ["ta-850"],
                "CIESM":        ["pr"],
                "KIOST-ESM":    ["zg-500"],
            }

        target_models = set(var_model.keys())

        for stat, seasons in data_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, table in regions.items():
                    df = table if isinstance(table, pd.DataFrame) else pd.DataFrame(table)
                    df = df.copy()

                    if "model" not in df.columns:
                        if verbose:
                            logger.warning(f"[badvals] Missing 'model' column for ({stat}, {season}, {region}); skipping.")
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    if not set(df["model"].astype(str).unique()).intersection(target_models):
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    cols_to_masks: Dict[str, np.ndarray] = {}

                    for model_name, cols in var_model.items():
                        model_mask = (df["model"].astype(str) == model_name).to_numpy()
                        if not model_mask.any():
                            continue
                        valid_cols = [c for c in cols if c in df.columns]
                        for col in valid_cols:
                            if col not in cols_to_masks:
                                cols_to_masks[col] = np.zeros(len(df), dtype=bool)
                            cols_to_masks[col] |= model_mask

                    nulled_counts = {}
                    for col, mask in cols_to_masks.items():
                        if not mask.any():
                            continue
                        if pd.api.types.is_integer_dtype(df[col].dtype):
                            df[col] = df[col].astype("float64")
                        before = df.loc[mask, col].notna().sum()
                        df.loc[mask, col] = np.nan
                        after = df.loc[mask, col].notna().sum()
                        nulled_counts[col] = before - after

                    if verbose and nulled_counts:
                        total = sum(nulled_counts.values())
                        detail = ", ".join(f"{k}:{v}" for k, v in sorted(nulled_counts.items()))
                        logger.info(f"[badvals] ({stat}, {season}, {region}) nulled {total} -> {detail}")

                    data_lib.df_dict[stat][season][region] = df

        return data_lib
    
    def _exclude_models(
            self, 
            data_lib, 
            model_list: Optional[Union[List[str], str]] = None,
            *, 
            verbose: bool = False
        ):
        """
        Exclude rows where df['model'] is in model_list.
        """
        if not model_list:  # handles None, [], ""
            return data_lib
        if isinstance(model_list, str):
            model_list = [model_list]

        for stat, seasons in data_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, table in regions.items():
                    df = table if isinstance(table, pd.DataFrame) else pd.DataFrame(table)
                    if "model" not in df.columns:
                        if verbose:
                            logger.info(f"[exclude_models] No 'model' column in ({stat}, {season}, {region}); skipping.")
                        data_lib.df_dict[stat][season][region] = df
                        continue
                    before = len(df)
                    df["model"] = df["model"].astype(str)
                    df = df[~df["model"].isin(model_list)].reset_index(drop=True)
                    if verbose:
                        removed = before - len(df)
                        if removed:
                            logger.info(f"[exclude_models] ({stat}, {season}, {region}) removed {removed} rows.")
                    data_lib.df_dict[stat][season][region] = df
        return data_lib
    
    def collect(self):
        # --- reference (optional) ---
        ref_lib = None
        if self.ref_name is not None:
            ref_lib = self._load_combined_metrics(
                self.ref_path,
                self.ref_group,
                self.ref_model,
                self.ref_version,
            )
            if ref_lib is None:
                logger.warning("ref_name provided but _load_combined_metrics() returned None.")
            elif hasattr(ref_lib, "var_list"):
                logger.debug(f"[ClimMetricsReader] ref_lib.var_list: {ref_lib.var_list}")
                
            # exclude bad values in specific models (treat {} as None to use defaults)
            ref_lib = self._check_badvals(ref_lib, var_model=(self.exclude_vars or None))
            
            # exclude specific models if needed
            ref_lib = self._exclude_models(ref_lib, model_list=self.exclude_models)

        # --- All models ---
        all_lib = None
        if self.test_combined and self.test_name is not None:
            all_lib = self._load_combined_metrics(
                self.test_path,
                self.test_group,
                self.test_model,
                self.test_version,
            )
            if all_lib is None:
                logger.warning("test_name provided but _load_combined_metrics() returned None.")
            elif hasattr(all_lib, "var_list"):
                logger.debug(f"[ClimMetricsReader] all_lib.var_list: {all_lib.var_list}")
                
            # exclude bad values in specific models (treat {} as None to use defaults)
            all_lib = self._check_badvals(all_lib, var_model=(self.exclude_vars or None))
            
            # exclude specific models if needed
            all_lib = self._exclude_models(all_lib, model_list=self.exclude_models)
        else:
            all_names = []
            for i, (mip_name, test_name, case_id) in enumerate(
                zip(self.mips, self.tests, self.caseids), start=1
            ):
                logger.debug(f"Processing model {i}: mip_name={mip_name}, test_name={test_name}")
                model_lib = self._process_test_model(mip_name, test_name, case_id)

                if model_lib is None:
                    raise RuntimeError(f"_process_test_model returned None for model '{test_name}'")

                if hasattr(model_lib, "var_list"):
                    logger.debug(f"[ClimMetricsReader] model_lib.var_list ({test_name}): {model_lib.var_list}")

                all_lib = model_lib if all_lib is None else all_lib.merge(model_lib)
                all_names.append(test_name)

            logger.info(f"Merging all model metrics: {all_names}")

            if all_lib is not None and hasattr(all_lib, "var_list"):
                logger.debug(f"[ClimMetricsReader] all_lib.var_list: {all_lib.var_list}")

            # exclude bad values in specific models (treat {} as None to use defaults)
            all_lib = self._check_badvals(all_lib, var_model=(self.exclude_vars or None))

            # exclude specific models if needed
            all_lib = self._exclude_models(all_lib, model_list=self.exclude_models)
        
        return ref_lib, all_lib
