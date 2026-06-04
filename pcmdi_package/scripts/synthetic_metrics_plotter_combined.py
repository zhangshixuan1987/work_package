import os
import re
from collections import OrderedDict
from typing import Any, Dict, List, Iterable, Optional, Tuple, Union

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import PathPatch
from matplotlib.lines import Line2D

from pcmdi_metrics.enso.lib import enso_portrait_plot
from pcmdi_metrics.graphics import (
    normalize_by_median,
    parallel_coordinate_plot,
    portrait_plot,
)

from logger import _setup_child_logger

from clim_metrics_reader import (
    ClimMetricsReader,
)

from clim_metrics_merger import (
    ClimMetricsMerger,
)

from enso_metrics_reader import (
    EnsoMetricsReader,
)
from movs_metrics_reader import (
    MoVsMetricsReader,
)
from utils import (
    get_highlight_models,
    realign_cbar_and_legend,
    drop_vars,
    archive_data
)

logger = _setup_child_logger(__name__)


class SyntheticMetricsPlotter:
    def __init__(
        self,
        test_group: str = None,
        test_clim_dir: Optional[str] = None,
        test_clim_set: Optional[str] = None,
        test_movs_dir: Optional[str] = None,
        test_movs_set: Optional[str] = None,
        test_enso_dir: Optional[str] = None,
        test_enso_set: Optional[str] = None,
        figure_format: str,
        metric_dict: Dict[str, Any],
        save_data: bool,
        base_test_input_path: str,
        # --- everything after here has defaults ---
        results_dir: Optional[str] = None,
        # Reference data 
        ref_group: str = None,
        ref_clim_dir: Optional[str] = None,
        ref_clim_set: Optional[str] = None,
        ref_movs_dir: Optional[str] = None,
        ref_movs_set: Optional[str] = None,
        ref_enso_dir: Optional[str] = None,
        ref_enso_set: Optional[str] = None,
        # Mean-climate viewer
        clim_viewer: bool = True,
        clim_vars: Optional[Union[List[str], str]] = None,
        clim_regions: Optional[Union[List[str], str]] = None,
        # Atmosphere modes (MOVA)
        mova_viewer: bool = True,
        mova_modes: Optional[Union[List[str], str]] = None,
        mova_obs: str = None, 
        # Coupled modes (MOVC)
        movc_viewer: bool = True,
        movc_modes: Optional[Union[List[str], str]] = None,
        movc_obs: Optional[str] = None,
        # ENSO viewer
        enso_viewer: bool = True,
        # Setup for visualization
        test_model_only: bool = False,
        movs_group: Optional[str] = None,
        exclude_vars: Optional[Dict[str, Any]] = None, 
        mean_group1_name: Optional[str] = None,
        mean_group2_name: Optional[str] = None,
        exclude_models: Optional[Union[List[str], str]] = None,  
        error_norm: Optional[str] = None, 
        extra_groups_name: Optional[Union[List[str], str]] = None,
    ):
        # Core
        self.figure_format = figure_format
        self.metric_dict = metric_dict
        self.save_data = bool(save_data)
        self.base_test_input_path = base_test_input_path
        self.results_dir = results_dir or "."
        
        # Test data 
        self.test_group = test_group if test_group is not None else "E3SM"
        self.test_clim_dir = test_clim_dir
        self.test_clim_set = test_clim_set
        self.test_movs_dir = test_movs_dir
        self.test_movs_set = test_movs_set
        self.test_enso_dir = test_enso_dir
        self.test_enso_set = test_enso_set
        
        # Reference data 
        self.ref_group = ref_group if ref_group is not None else "CMIP"
        self.ref_clim_dir = ref_clim_dir
        self.ref_clim_set = ref_clim_set
        self.ref_movs_dir = ref_movs_dir
        self.ref_movs_set = ref_movs_set
        self.ref_enso_dir = ref_enso_dir
        self.ref_enso_set = ref_enso_set

        # Mean climate
        self.clim_viewer = bool(clim_viewer)
        self.clim_vars = self._to_list(clim_vars)  # [] => all available
        self.clim_regions = self._to_list(clim_regions)  # [] => all regions
        self.test_model_only = test_model_only
        self.exclude_vars = exclude_vars or {}
        self.exclude_models = self._to_list(exclude_models)
        self.error_norm = error_norm if error_norm is not None else "default"
        self.movs_group = movs_group if movs_group is not None else "cbf"
            
        # MOVA
        self.mova_viewer = bool(mova_viewer)
        self.mova_modes = self._to_list(mova_modes) if self.mova_viewer else []
        self.mova_obs = mova_obs or ("NOAA-20C" if self.mova_viewer else None)
        self.movc_viewer = bool(movc_viewer)
        self.movc_modes = self._to_list(movc_modes) if self.movc_viewer else []
        self.movc_obs = movc_obs or ("HadISST" if self.movc_viewer else None)
        self.movs_viewer = self.mova_viewer or self.movc_viewer

        # ENSO
        self.enso_viewer = bool(enso_viewer)
        
        # Group Mean 
        self.mean_group1_name = mean_group1_name if mean_group1_name is not None else ref_group
        self.mean_group2_name = mean_group2_name if mean_group2_name is not None else test_group
        self.extra_groups_name = self._to_list(extra_groups_name)

        # Final bundle for downstream readers/builders
        self.parameter: Dict[str, Any] = self._initialize_parameter()

    # ---------- helpers ----------
    def _initialize_parameter(self):
        out_dir = os.path.join(self.results_dir, "ERROR_metric")
        os.makedirs(out_dir, exist_ok=True)
        
        param = OrderedDict({
            "save_data": self.save_data,
            "out_dir": out_dir,
            "test_group": self.test_group,
            "ref_group" : self.ref_group,
            "mean_group1_name": self.mean_group1_name,
            "mean_group2_name": self.mean_group2_name,
            "extra_groups_name": self.extra_groups_name, 
            "movs_group": self.movs_group,
            "test_model_only": self.test_model_only,
        })
        return param

    def generate(self, figure_sets=None, debug=False) -> None:
        logger.info("Generating synthetic metrics plots ...")
        if figure_sets is None: 
            tasks = [
                (self.clim_viewer, "mean_climate", self._handle_mean_climate),
                (self.movs_viewer, "variability_modes", self._handle_variability_modes),
                (self.enso_viewer, "enso_metric", self._handle_enso_metric),
            ]
        else:
            tasks = [] 
            if "mean_climate" in figure_sets:
                tasks.append((self.clim_viewer, "mean_climate", self._handle_mean_climate))
            if "variability" in figure_sets:
                tasks.append((self.movs_viewer, "variability_modes", self._handle_variability_modes))
            if "enso" in figure_sets:
                tasks.append((self.enso_viewer, "enso_metric", self._handle_enso_metric))

        at_least_one_success = False
        failures = []

        for enabled, metric, handler in tasks:
            self.parameter["test_path"] = self.base_test_input_path.replace(
                "%(group_type)", metric
            )
            self.parameter["diag_vars"] = self.metric_dict[metric]
            if not enabled:
                continue
            logger.info("Processing metric: %s", metric)
            
            if debug: 
                handler(metric)
                at_least_one_success = True
            else:
                try:
                    handler(metric)
                    at_least_one_success = True
                except Exception as e:
                    logger.error("Failed to handle metric=%s: %s", metric, e, exc_info=True)
                    failures.append(metric)

        if not at_least_one_success:
            raise RuntimeError("No synthetic metrics plots could be generated.")

        if failures:
            logger.warning("Completed with partial failures: %s", ", ".join(failures))

    def _handle_mean_climate(self, metric: str) -> None:
        logger.info("Handling mean climate…")
        
        self.parameter.update(
            {"ref_group": self.ref_group, 
             "ref_path": self.ref_clim_dir, 
             "ref_name": self.ref_clim_set,
             "exclude_vars": self.exclude_vars,
             "exclude_models": self.exclude_models,
             "error_norm": self.error_norm,
             "unit_check": True
            }
        )
        
        # sanity check ---
        mips  = self.parameter.get("test_mip", [])
        tests = self.parameter.get("test_name", [])
        caseids = self.parameter.get("test_id", [])
        # Guard: strict 1:1 and non-empty
        if not (len(mips) == len(tests) == len(caseids) > 0):
            logger.error("[mean_climate] test_name/model_name/case_id not aligned or empty.")
            return
        
        # --- collect metrics data ---
        collector = ClimMetricsReader(self.parameter)
        ref_lib, test_lib = collector.collect()
        
        # --- Prepare data for plotting function ---
        merger = ClimMetricsMerger(self.parameter)
        merge_lib, norm_lib, var_list, var_unit_list, regions, ref_list, test_list, e3sm_list, mean_list = merger.merge(
            ref_lib=ref_lib,
            model_lib=test_lib,
            clim_vars=self.clim_vars,
            clim_regions=self.clim_regions,
        )
        
        # Use the same `metric` variable as before (assuming it's defined in scope)
        for stat, diag_ in self.metric_dict[metric].items():
            logger.debug(f"[mean_climate] Running plot driver: stat={stat}")
            # Keep the exact positional calling convention you had before
            mean_climate_plot_driver(
                metric,
                stat,
                diag_,
                regions,
                var_list,
                var_unit_list,
                merge_lib.df_dict[stat],
                norm_lib.df_dict[stat],
                ref_list,
                test_list,
                e3sm_list,
                mean_list,
                self.parameter["test_group"],
                self.parameter["ref_group"],
                self.parameter["extra_groups_name"],
                self.parameter["save_data"],
                self.parameter["out_dir"],
                self.figure_format,
            )

    def _handle_variability_modes(self, metric: str) -> None:
        logger.info("Handling modes variability …")

        # Combine atmospheric and coupled modes (already lists)
        modes_list = (self.mova_modes or []) + (self.movc_modes or [])

        if not modes_list:
            logger.warning(
                "[variability_modes] No modes specified; skipping variability mode plots."
            )
            return

        if (self.movc_viewer and not self.movc_obs) or (self.mova_viewer and not self.mova_obs):
            logger.error("[variability_modes] Missing or empty reference data for atm_mode or cpl_mode viewers.")
            return
        
        # Update parameters for reader
        self.parameter.update(
            {"ref_group": self.ref_group, 
             "ref_path": self.ref_movs_dir, 
             "ref_name": self.ref_movs_set,
             "movs_mode": modes_list,
             "mova_obs" : self.mova_obs,
             "movc_obs" : self.movc_obs,
             "exclude_models": self.exclude_models,
             "error_norm": self.error_norm
            }
        )

        # sanity check ---
        mips  = self.parameter.get("test_mip", [])
        tests = self.parameter.get("test_name", [])
        caseids = self.parameter.get("test_id", [])
        # Guard: strict 1:1 and non-empty
        if not (len(mips) == len(tests) == len(caseids) > 0):
            logger.error("[variability_modes] test_name/model_name/case_id not aligned or empty.")
            return
        
        # Collect metrics
        reader = MoVsMetricsReader(self.parameter)
        merge_lib, norm_lib, mode_season_list, ref_list, test_list, e3sm_list, mean_list = reader.collect_metrics() 
        
        # Ensure metric exists in dictionary
        if metric not in self.metric_dict:
            logger.error(
                f"[variability_modes] Metric '{metric}' not found in metric_dict keys={list(self.metric_dict.keys())}"
            )
            return

        # Loop through stats and plot
        for stat, diag_ in self.metric_dict[metric].items():
            if stat not in merge_lib:
                logger.warning(
                    f"[variability_modes] stat='{stat}' not found in merge_lib; available={list(merge_lib.keys())}"
                )
                continue

            logger.debug(f"[variability_modes] Running plot driver for stat={stat}")
            variability_modes_plot_driver(
                metric,
                stat,
                diag_,
                mode_season_list,
                merge_lib[stat],
                norm_lib[stat],
                ref_list,
                test_list,
                e3sm_list,
                mean_list,
                self.parameter["test_group"],
                self.parameter["ref_group"],
                self.parameter["extra_groups_name"],
                self.parameter["save_data"],
                self.parameter["out_dir"],
                self.figure_format,
            )
              
    def _handle_enso_metric(self, metric: str, debug: bool = False) -> None:
        logger.info("Handling ENSO metrics…")
        
        # Update parameters for reader
        self.parameter.update(
            {"ref_group": self.ref_group, 
             "ref_path": self.ref_enso_dir, 
             "ref_name": self.ref_enso_set,
             "exclude_models": self.exclude_models,
             "error_norm": self.error_norm
            }
        )
        
        # --- Validate metric entry ---
        if metric not in self.metric_dict or not isinstance(
            self.metric_dict[metric], dict
        ):
            logger.error(
                f"[enso] metric_dict['{metric}'] missing or not a dict. Available: {list(self.metric_dict.keys())}"
            )
            return
        
        # sanity check ---
        mips  = self.parameter.get("test_mip", [])
        tests = self.parameter.get("test_name", [])
        caseids = self.parameter.get("test_id", [])
        # Guard: strict 1:1 and non-empty
        if not (len(mips) == len(tests) == len(caseids) > 0):
            logger.error("[enso] test_name/model_name/case_id not aligned or empty.")
            return
        
        # Initialize metrics reader 
        reader = EnsoMetricsReader(self.parameter)
        
        # --- Main loop over stats ---
        for stat in self.metric_dict[metric].keys():
            try:
                enso_mips,diag_dict,dict_json_path = reader.run(stat)
            except Exception as e:
                logger.exception(f"[enso] Reader failed for stat='{stat}': {e}")
                continue

            if not dict_json_path:
                logger.warning(
                    f"[enso] Reader returned empty path for stat='{stat}'. Skipping plot."
                )
                continue

            try:
                enso_plot_driver(
                    metric, 
                    stat, 
                    dict_json_path,
                    diag_dict,
                    enso_mips,
                    self.parameter["out_dir"],
                    self.figure_format
                )
                logger.debug(f"[enso] Plotted stat='{stat}' successfully.")
            except Exception as e:
                logger.exception(f"[enso] Plot driver failed for stat='{stat}': {e}")
                
                
def mean_climate_plot_driver(
    metric,
    stat,
    diag_dict,
    regions,
    var_list,
    var_unit_list,
    df_dict,
    norm_dict,
    ref_model_list,
    test_model_list,
    e3sm_model_list,    
    mean_model_list,
    test_group,
    ref_group,
    meanx_model_list,
    save_data,
    out_path,
    fig_format,
):
    """Driver Function for the mean climate metrics plot"""
    mout_name = f"{test_group}_vs_{ref_group}"
    
    for mtype in diag_dict["type"]:
        for region in regions:
            do_plot = region in diag_dict["region"]

            if do_plot and mtype == "portrait":
                logger.info(
                    "Processing Portrait  Plots for {} {} {}....".format(
                        metric, region, stat
                    )
                )

                data_nor = dict()
                run_list = []
                for season in diag_dict["season"]:
                    # drop data if all is NaNs
                    data_dict, var_names, var_units = drop_vars(
                        df_dict[season][region].copy(),
                        var_list.copy(),
                        var_unit_list.copy(),
                    )

                    ref_dict, var_names, var_units = drop_vars(
                        norm_dict[season][region].copy(),
                        var_list.copy(),
                        var_unit_list.copy(),
                    )
                    
                    logger.debug(
                        f"var_names={var_names} derived from var_list={var_list}."
                    )
                    
                    logger.debug(f"Available columns: {data_dict.columns.tolist()}")
                    
                    run_list = []                        
                    for name in data_dict['model']:
                        if name not in run_list:
                            run_list.append(name)
                            
                    # --- guard: ensure all requested columns exist ---
                    missing_in_data = [v for v in var_names if v not in data_dict.columns]
                    missing_in_ref  = [v for v in var_names if v not in ref_dict.columns]
                    if missing_in_data or missing_in_ref:
                        parts = []
                        if missing_in_data:
                            parts.append(f"missing in data_dict: {missing_in_data}")
                        if missing_in_ref:
                            parts.append(f"missing in ref_dict: {missing_in_ref}")
                        msg = "; ".join(parts)
                        logger.error("Column mismatch for var_names -> " + msg)
                        raise KeyError(msg)

                    try:
                        if stat == "cor_xy":
                            data_nor[season] = data_dict.loc[:, var_names].to_numpy().T
                        else:
                            data_nor[season] = normalize_by_median(
                                data_dict.loc[:, var_names].to_numpy().T, 
                                data_median=ref_dict.loc[:, var_names].to_numpy().T, 
                                axis=1
                            )
                        # save data if requested (save normalized or raw-for-cor_xy)
                        if save_data:
                            outdir = os.path.join(out_path, metric, region)
                            os.makedirs(outdir, exist_ok=True)

                            outdic = data_dict.drop(columns=["model_run"], errors="ignore").copy()
                            outdic.loc[:, var_names] = data_nor[season].T

                            archive_data(
                                region,
                                stat,
                                season,
                                outdic,          # <-- save the normalized dataframe
                                mout_name,
                                var_names,
                                var_units,
                                outdir,
                            )
                    except KeyError as e:
                        logger.error(f"KeyError on var_names={var_names}: {e}")
                        raise

                outdir = os.path.join(out_path, metric)
                os.makedirs(outdir, exist_ok=True)

                portrait_metric_plot(
                    region,
                    stat,
                    metric,
                    data_nor,
                    diag_dict["name"],
                    var_names,
                    var_units,
                    run_list,
                    ref_model_list,
                    test_model_list,
                    e3sm_model_list,    
                    outdir,
                    fig_format,
                    mean_list=mean_model_list,
                    meanx_list=meanx_model_list,
                    base_fontsize=20,
                    base_figsize=(50, 20),
                    base_legend_lw=1.5,
                )

            elif do_plot and mtype == "parcoord":
                logger.info(
                    "Processing Parallel Coordinate Plots for {} {} {}....".format(
                        metric, region, stat
                    )
                )
                for season in diag_dict["season"]:
                    if season in df_dict.keys():
                        # drop data if all is NaNs
                        data_dict, var_names, var_units = drop_vars(
                            df_dict[season][region].copy(),
                            var_list.copy(),
                            var_unit_list.copy(),
                        )
                        
                        ref_dict, var_names, var_units = drop_vars(
                            norm_dict[season][region].copy(),
                            var_list.copy(),
                            var_unit_list.copy(),
                        )
                                
                        if save_data:
                            outdir = os.path.join(out_path, metric, region)
                            os.makedirs(outdir, exist_ok=True)
                            outdic = data_dict.drop(columns=["model_run"], errors="ignore").copy()
                            archive_data(
                                region, stat, season, outdic, mout_name,
                                var_names, var_units, outdir,
                            )
                            
                        run_list = []                        
                        for name in data_dict['model']:
                            if name not in run_list:
                                run_list.append(name)
                                
                        outdir = os.path.join(out_path, metric)
                        os.makedirs(outdir, exist_ok=True)
                        parcoord_metric_plot(
                            region, 
                            stat, 
                            metric, 
                            data_dict,
                            diag_dict["name"], 
                            var_names, 
                            var_units, 
                            run_list,
                            ref_model_list,
                            test_model_list,
                            e3sm_model_list, 
                            outdir, 
                            fig_format,
                            model_group=test_group,
                            ref_group=ref_group,
                            mean_list=mean_model_list,
                            meanx_list=meanx_model_list,
                            base_fontsize=24,
                            base_figsize=(50, 20),
                            base_legend_lw=1.5,
                        )
    return

def variability_modes_plot_driver(
    metric,
    stat,
    metric_dict,
    mode_season_list,
    df_dict,
    norm_dict,
    ref_model_list,
    test_model_list,
    e3sm_model_list,
    mean_model_list,
    test_group,
    ref_group,
    meanx_model_list,
    save_data,
    out_path,
    fig_format,
):
    """Driver Function for the modes variability metrics plot"""
    season = "all"
    mout_name = f"{test_group}_vs_{ref_group}"
    
    # ---- model list & exclusions
    run_list: List[str] = df_dict["model"].to_list()
        
    # Name of plot variable 
    stat_name = metric_dict["name"]

    # drop data if all is NaNs
    ref_dict, var_names, var_units = drop_vars(
        norm_dict.loc[:, mode_season_list].copy(), 
        mode_season_list.copy(), 
        None
    )
    
    data_dict, var_names, var_units = drop_vars(
        df_dict.loc[:, var_names].copy(), 
        var_names.copy(), 
        None
    )
    
    # Loop and plot 
    for mtype in metric_dict["type"]:
        if mtype == "portrait":
            
            logger.info("Processing Portrait  Plots for {} {}....".format(metric, stat))
            
            if stat not in ["stdv_pc_ratio_to_obs"]:
                data_nor = normalize_by_median(
                    data_dict.to_numpy().T, 
                    data_median=ref_dict.to_numpy().T, 
                    axis=1
                )
            else:
                data_nor = data_dict.to_numpy().T
                
            if save_data:
                outdir = os.path.join(out_path, metric)
                os.makedirs(outdir, exist_ok=True)
                data_dict[mode_season_list] = data_nor.T
                archive_data(
                    metric,stat,season,data_dict,
                    mout_name,mode_season_list,
                    None,outdir,
                )
                
            portrait_metric_plot(
                    metric,
                    stat,
                    season,
                    data_nor,
                    stat_name,
                    var_names,
                    var_units,
                    run_list,
                    ref_model_list,
                    test_model_list,
                    e3sm_model_list,    
                    out_path,
                    fig_format,
                    mean_list=mean_model_list,
                    meanx_list=meanx_model_list,
                    base_fontsize=20,
                    base_figsize=(50, 20),
                    base_legend_lw=1.5,
                )

        elif mtype == "parcoord":
            
            logger.info(
                "Processing Parallel Coordinate Plots for {} {}....".format(
                    metric, stat
                )
            )
            
            if save_data:
                outdir = os.path.join(out_path, metric)
                os.makedirs(outdir, exist_ok=True)
                archive_data(
                    metric,stat,season,data_dict,
                    mout_name,mode_season_list,
                    None,outdir,
                )
                
            parcoord_metric_plot(
                metric,
                stat,
                season,
                data_dict,
                stat_name,
                var_names,
                var_units,
                run_list,
                ref_model_list,
                test_model_list,
                e3sm_model_list, 
                out_path,
                fig_format,
                model_group=test_group,
                ref_group=ref_group,
                mean_list=mean_model_list,
                meanx_list=meanx_model_list,
                base_fontsize=20,
                base_figsize=(50, 20),
                base_legend_lw=1.5,
            )

    return

def enso_plot_driver(
        metric, 
        stat, 
        dict_json_path,
        metric_dict, 
        enso_mips,
        out_dir,
        fig_format,
        reduced_set=True,
        sort_y_names=True, 
        show_proj_means=False, 
        show_ref_row=True, 
        show_alt_obs_rows=False
):
    """
    Driver function to plot ENSO metrics based on specified type (e.g., portrait).
    """
    metrics_collections = metric_dict["collection"]
    tests, models, caseid = map(list, zip(*enso_mips))
    for mtype in metric_dict["type"]:
        if mtype == "portrait":
            logger.info(f"Processing Portrait Plots for {metric} {stat}...")

            list_project = models
            list_obs: List[object] = (
                []
            )  # fill in if observational references are needed
            outdir = os.path.join(out_dir, metric)
            os.makedirs(outdir, exist_ok=True)

            outfile = f"{metric}_{stat}_portrait.{fig_format}"
            figure_name = os.path.join(outdir, outfile)
            
            fig, ref_info_dict = enso_portrait_plot(
                metrics_collections, 
                list_project, 
                list_obs, 
                dict_json_path, 
                figure_name=figure_name, 
                reduced_set=reduced_set,
                sort_y_names=sort_y_names, 
                show_proj_means=show_proj_means, 
                show_ref_row=show_ref_row, 
                show_alt_obs_rows=show_alt_obs_rows
            )

    return

def portrait_metric_plot(
    region,
    stat,
    group,
    data_dict,
    stat_name,
    var_list,
    unit_list,
    run_list,
    ref_list,
    test_list,
    e3sm_list,
    out_path,
    fig_format,
    mean_list=None,
    meanx_list=None,
    show_unit=False,
    base_fontsize=20,
    base_figsize=(50, 20),
    base_legend_lw=1.5,
    box_as_square=True,
    missing_color="white",
    logo_rect=[0, 0, 0, 0],
    logo_off=True,
):
    
    # === Figure scaling setup ===
    fscale = len(var_list) / 35.0
    fscale = max(0.5, min(fscale, 1.5))
    fontsize = base_fontsize
    figsize = (base_figsize[0], base_figsize[1] * fscale)
    legend_box_xy = (1.025, 0.98)
    legend_box_size = 4 * fscale
    legend_lw = base_legend_lw * fscale
    shrink = 0.8 * fscale
    pad = 0.015
    legend_fontsize = fontsize * 0.8

    # --- SMALL GUARD A: basic inputs present ---
    if not var_list:
        logger.warning("[Portrait]: No variables to plot (var_list empty); returning.")
        return

    var_plot = var_list.copy()
    if show_unit and unit_list is not None:
        for i, unit in enumerate(unit_list):
            var_plot[i] = f"{var_list[i]}{unit}"  # keep your original format

    if run_list is None:
        logger.warning("[Portrait]: No models to plot (run_list empty); returning.")
        return

    # Normalize possibly-None lists to empty lists for safe iteration
    mean_list = mean_list or []
    e3sm_list = e3sm_list or [] 
    test_list = test_list or []
    ref_list = ref_list or []
    meanx_list = meanx_list or []
    
    if not test_list:
        logger.warning("[Portrait]: No test models to plot (test_list empty).")
    if not ref_list:
        logger.warning("[Portrait]: No reference models to plot (ref_list empty).")
    if not e3sm_list:
        logger.warning("[Portrait]: No e3sm models to plot (e3sm_list empty).")
    if not mean_list:
        logger.warning("[Portrait]: No groupd mean models to plot (mean_list empty).")
    if not meanx_list: 
        logger.warning("[Portrait]: No extra mean models to plot (meanx_list empty).")
    
    # --- Build data array ---
    if group == "mean_climate":
        required = ["djf", "mam", "jja", "son"]
        missing = [k for k in required if (k not in data_dict or data_dict[k] is None)]
        if missing:
            logger.warning("[Portrait]: Missing seasonal arrays for %s; returning. Missing=%s", group, missing)
            return
        try:
            arrs = [np.asarray(data_dict[k]) for k in required]
            if any(a.size == 0 for a in arrs):
                logger.warning("[Portrait]: One or more seasonal arrays are empty; returning.")
                return
            # Expect shape (season, nvar, nmodel) OR (season, something, ...)
            data_all_nor = np.stack(arrs, axis=0)
        except Exception as e:
            logger.warning("[Portrait]: Failed to stack seasonal arrays: %s; returning.", e)
            return
        legend_on = True
        legend_labels = ["DJF", "MAM", "JJA", "SON"]
        title = f"{region} — {group} ({stat_name})"
    else:
        data_all_nor = np.asarray(data_dict)
        if data_all_nor.size == 0:
            logger.warning("[Portrait]: Input data array is empty; returning.")
            return
        legend_on = False
        legend_labels = []
        title = f"{region} ({stat_name})"
        
    # --- Construct highlight model list ---
    highlight_models = [] 
    for m in e3sm_list:
        if m in run_list and m not in highlight_models:
            highlight_models.append(m)
    for m in mean_list:
        if m in run_list and m not in highlight_models:
            highlight_models.append(m)
    for m in meanx_list: 
        if m in run_list and m not in highlight_models:
            highlight_models.append(m)
    
    # --- REORDER COLUMNS (minimal, but critical for the separator lines) ---
    # Keep original order for each bucket: regular -> highlighted(non-mean) -> means
    run_index = {m: i for i, m in enumerate(run_list)}
    hi_set = set(highlight_models)
    mean_set = set(mean_list)

    regular = [m for m in run_list if m not in hi_set]
    hi_non_means = [m for m in run_list if (m in hi_set and m not in mean_set)]
    means = [m for m in run_list if m in mean_set and m in hi_set]

    new_run_list = regular + hi_non_means + means
    reindex = [run_index[m] for m in new_run_list]

    # Reindex along the last axis (models live on the x-axis)
    data_all_nor = data_all_nor[..., reindex]
    model_plot = new_run_list

    # Colors for x tick labels
    label_colors = []
    for m in model_plot:
        if m in mean_set:
            label_colors.append("#FC5A50")  # red for means 
        elif m in hi_set:
            label_colors.append("#5170d7")  # blue for e3sm family
        else:
            label_colors.append("#000000")  # black otherwise

    # --- Colormap/range selection (kept as-is) ---
    if stat in ["cor_xy"]:
        var_range = (0, 1.0)
        cmap_color = "viridis"
        cmap_bounds = np.linspace(0, 1, 21)
        cbar_label = "Pattern corr."
    elif stat in ["stdv_pc_ratio_to_obs"]:
        var_range = (0.5, 1.5)
        cmap_color = "jet"
        cmap_bounds = [r / 10 for r in range(5, 16, 1)]
        cbar_label = "PC σ ratio"
    elif stat in ["mae_xy"]:
        var_range = (-0.5, 0.5)
        cmap_color = "RdYlBu_r"
        cmap_bounds = np.linspace(-0.5, 0.5, 11)
        cbar_label = "Normalized MAE"
    elif stat in ["bias_xy"]:
        var_range = (-0.5, 0.5)
        cmap_color = "RdYlBu_r"
        cmap_bounds = np.linspace(-0.5, 0.5, 11)
        cbar_label = "Bias (norm.)"
    elif stat in ["rms_xy"]:
        var_range = (-0.5, 0.5)
        cmap_color = "RdYlBu_r"
        cmap_bounds = np.linspace(-0.5, 0.5, 11)
        cbar_label = "Normalized RMSE"
    else:
        var_range = (-0.5, 0.5)
        cmap_color = "RdYlBu_r"
        cmap_bounds = np.linspace(-0.5, 0.5, 11)
        cbar_label = f"{stat}"

    fig, ax, cbar = portrait_plot(
        data_all_nor,
        xaxis_labels=model_plot,
        yaxis_labels=var_plot,
        cbar_label=cbar_label,
        cbar_label_fontsize=fontsize * 0.95,
        cbar_tick_fontsize=fontsize * 0.95,
        box_as_square=box_as_square,
        vrange=var_range,
        figsize=figsize,
        cmap=cmap_color,
        cmap_bounds=cmap_bounds,
        cbar_kw={"extend": "both", "shrink": shrink, "pad": pad},
        missing_color=missing_color,
        legend_on=legend_on,
        legend_labels=legend_labels,
        legend_box_xy=legend_box_xy,
        legend_box_size=legend_box_size,
        legend_lw=legend_lw,
        legend_fontsize=legend_fontsize,
        logo_rect=logo_rect,
        logo_off=logo_off,
    )
    
    # --- Vertical separators now align with the reordering ---
    n_total = len(model_plot)  # <-- use reordered list length
    n_highlight = len(hi_non_means) + len(means)
    n_means = len(means)

    if n_highlight > 0:
        ax.axvline(x=n_total - n_highlight, color="k", linewidth=3)
    if 0 < n_means < n_highlight:
        ax.axvline(x=n_total - n_means, color="k", linewidth=3)

    # Use the reordered labels and colors
    ax.set_xticklabels(model_plot, rotation=45, va="bottom", ha="left")
    ax.set_yticklabels(var_plot,  rotation=0,  va="center", ha="right")
    for xtick, color in zip(ax.get_xticklabels(), label_colors):
        xtick.set_color(color)
        
    # 1) layout first
    fig.set_constrained_layout(False)   # guard, in case it's on
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    # 2) then place cbar + seasonal box
    realign_cbar_and_legend(
        fig, ax, cbar,
        cbar_side="right",
        auto_nudge=True,
        label_clearance_in=0.02,   # how much space beyond labels
        min_buffer_frac=0.0002,    # 0.0002 × 50 in ≈ 0.01"
        y_gap_above_cbar_in=0.005,
        top_guard=0.98,
    )

    # 3) center title over the data axes (exclude cbar/insets)
    ax_box = ax.get_position()
    x_center = 0.5 * (ax_box.x0 + ax_box.x1)
    # remove any existing suptitle to avoid double titles
    if getattr(fig, "_suptitle", None) is not None:
        fig._suptitle.remove()
        
    fig.text(x_center, 1.0, title, ha="center", va="top",
             fontsize=fontsize*1.1, fontweight="bold", transform=fig.transFigure)

    outdir = os.path.join(out_path, region)
    os.makedirs(outdir, exist_ok=True)
    outfile = f"{stat}_{region}_portrait_{group}.{fig_format}"
    fig.savefig(os.path.join(outdir, outfile), facecolor="w", bbox_inches="tight")
    plt.close(fig)
    
    return
                           
def parcoord_metric_plot(
    region,
    stat,
    metric,
    data_dict,
    stat_name,
    var_list,
    unit_list,
    run_list,
    ref_list,
    test_list,
    e3sm_list,
    out_path,
    fig_format,
    model_group=None,
    ref_group=None,
    mean_list=None,
    meanx_list=None,
    index_col="model",
    show_unit=False,
    base_fontsize=24,
    base_figsize=(50, 20),
    base_legend_lw=1.5,
    xcolors=None,
    style_cycle=None,
    color_map="tab20_r",
    identify_all_models=False,
    vertical_center=None,
    vertical_center_line=True,
    show_boxplot=False,
    show_violin=True,
    violin_colors=("lightgrey", "pink"),
    logo_rect=[0, 0, 0, 0],
    logo_off=True,
):
    """Function for parallel coordinate plots (revised: stable highlights & legend)."""
    # --- set up color and style cycles ---
    color_cycle = xcolors if xcolors is not None else [
        "#ff7f00", "#4daf4a", "#f781bf", "#a65628", "#984ea3", "#377eb8", "#dede00"
    ]
    line_styles = style_cycle if style_cycle is not None else ["solid", "dashed", "dashdot", "dotted"]

    if run_list is None or len(run_list) == 0:
        logger.warning("[ParCoord]: No models to plot (run_list empty); returning.")
        return

    # --- SMALL GUARD A: basic inputs present ---
    if not var_list:
        logger.warning("[ParCoord]: No variables to plot (var_list empty); returning.")
        return

    # --- Keep only existing, non-empty vars ---
    var_plot = sorted(v for v in var_list if (v in data_dict.columns) and data_dict[v].notna().any())
    if not var_plot:
        logger.warning(
            f"[ParCoord]: Nothing to plot for metric={metric}, region={region}, stat={stat}. "
            f"No valid variables found in metrics data (columns checked={len(var_list)})."
        )
        return

    # axis labels (align units to filtered vars)
    var_labels = var_plot.copy()
    if show_unit and unit_list is not None:
        unit_map = {v: u for v, u in zip(var_list, unit_list)}
        var_labels = [f"{v} [{unit_map.get(v, '')}]" if unit_map.get(v) else v for v in var_plot]

    # --- Build data matrix aligned to run_list if possible ---
    if index_col is not None and index_col in data_dict.columns:
        # De-dup by index to avoid multiple lines per label
        if data_dict[index_col].duplicated().any():
            logger.warning("[ParCoord]: Duplicated %s found; keeping first occurrence.", index_col)
        tmp = (
            data_dict
            .loc[~data_dict[index_col].duplicated(keep="first")]
            .set_index(index_col)
        )

        present = [m for m in run_list if m in tmp.index]
        missing = [m for m in run_list if m not in tmp.index]
        if missing:
            logger.warning("[ParCoord]: %d models from run_list missing in data and will be skipped: %s",
                           len(missing), missing)
        if not present:
            logger.warning("[ParCoord]: None of run_list found in data; returning.")
            return

        data_var = tmp.loc[present, var_plot].to_numpy()
        run_list = present  # keep labels aligned
    else:
        data_var = data_dict[var_plot].to_numpy()

    # Normalize possibly-None lists; de-dup while preserving order
    mean_list   = list(dict.fromkeys(mean_list or []))
    meanx_list  = list(dict.fromkeys(meanx_list or []))
    test_list  = list(dict.fromkeys(test_list or []))
    ref_list    = list(dict.fromkeys(ref_list or []))
    e3sm_list    = list(dict.fromkeys(e3sm_list or []))

    # --- Construct highlight model list (violin shading: group2 vs ref group) ---
    if model_group is None:
        model_group = "TEST"
    if ref_group is None:
        ref_group = "Reference"
        
    models_to_highlight1 = []
    for m in test_list:
        if m in run_list and m not in models_to_highlight1:
            models_to_highlight1.append(m)
    for m in mean_list:
        if m in run_list and m not in ref_list and m not in models_to_highlight1:
            models_to_highlight1.append(m)
    if not models_to_highlight1:
        logger.warning(f"[ParCoord]: No second model group found for model_group={model_group}")

    # --- Highlights (lines) with guaranteed alignment and stable mapping ---
    if len(mean_list) > 2:
        logger.warning(f"[ParCoord]: size of mean_list is supposed to be less than 2: {mean_list}")

    style_map = {}  # <-- INIT so we can fill it below
    combined = list(meanx_list) + [m for m in e3sm_list if m not in test_list and m not in meanx_list and m not in mean_list]
    for idx, m in enumerate(combined):
        style_map[m] = (color_cycle[idx % len(color_cycle)],
                        line_styles[idx % len(line_styles)])
        
    for m in mean_list:
        if m in ref_list:
            style_map[m] = ("#000000", "solid")  # ref mean in black
        elif m in test_list:
            style_map[m] = ("#e41a1c", "solid")  # test mean in red

    # --- Highlights (lines) with guaranteed alignment and stable mapping ---
    models_to_highlight2 = []
    hl_colors, hl_linestyles = [], []
    for m in style_map.keys():  # dict preserves insertion order (combined first, then means)
        if m in run_list and m not in models_to_highlight2:
            models_to_highlight2.append(m)
            c, ls = style_map[m]
            hl_colors.append(c)
            hl_linestyles.append(ls)
        
    # --- Reset violin when second group identical or no ref_list provided ---
    if sorted(models_to_highlight1) == sorted(run_list) or ref_list is None:
        logger.warning("[ParCoord]: second group identical to run_list; resetting highlights.")
        models_to_highlight1 = None
        violin_colors = ("lightgrey",)
        show_boxplot = True
        show_violin = False
    
    # === Figure scaling setup ===
    fscale = (len(var_labels) or 1) / 30.0
    fscale = max(0.6, min(fscale, 1.5))

    fontsize = base_fontsize
    figsize = (base_figsize[0] * fscale, base_figsize[1] * fscale)

    legend_ncol = max(1, int(7 * figsize[0] / 40.0))
    legend_position = (0.50, -0.14)

    xlabel = "Metric"
    ylabel = "{} ({})".format(stat_name, stat.upper())

    if "mean_climate" in [metric, region]:
        title = f"Model Performance of Mean Climatology ({stat.upper()}, {region.upper()})"
    elif "variability_modes" in [metric, region]:
        title = f"Model Performance of Modes Variability ({stat.upper()})"
    elif "enso" in [metric, region]:
        title = f"Model Performance of ENSO ({stat.upper()})"
    else:
        title = f"Model Performance ({stat.upper()}, {region.upper()})"

    if vertical_center is None:
        if stat in ["stdv_pc_ratio_to_obs"]:
            vertical_center = 1.0
        elif stat in ["cor_xy"]:
            vertical_center = 0.5
        elif stat in ["bias_xy"]:
            vertical_center = 0.0
        else:
            vertical_center = "median"

    fig, ax = parallel_coordinate_plot(
        data=data_var,
        metric_names=var_labels,
        model_names=run_list,
        model_names2=models_to_highlight1,
        group1_name=ref_group,
        group2_name=model_group,
        models_to_highlight=models_to_highlight2,
        models_to_highlight_colors=hl_colors,
        models_to_highlight_labels=models_to_highlight2,
        identify_all_models=identify_all_models,
        vertical_center=vertical_center,
        vertical_center_line=vertical_center_line,
        title="",
        figsize=figsize,
        colormap=color_map,
        show_boxplot=show_boxplot,
        show_violin=show_violin,
        violin_colors=violin_colors,
        legend_ncol=legend_ncol,
        legend_bbox_to_anchor=legend_position,
        legend_fontsize=fontsize * 0.85,
        xtick_labelsize=fontsize * 0.95,
        ytick_labelsize=fontsize * 0.95,
        logo_rect=logo_rect,
        logo_off=logo_off,
    )

    # apply neutral styling to violin/box elements when no second-group highlight
    if not models_to_highlight1:
        axes = np.atleast_1d(ax)
        for a in axes:
            for coll in getattr(a, "collections", []):
                if isinstance(coll, PolyCollection):
                    try:
                        n = len(coll.get_facecolors())
                        coll.set_facecolors([(0.75, 0.75, 0.75, 1.0)] * max(n, 1))
                    except Exception:
                        coll.set_facecolor((0.75, 0.75, 0.75, 1.0))
                    coll.set_edgecolor("black")
                    coll.set_linewidth(1.5)
            for p in getattr(a, "patches", []):
                if isinstance(p, PathPatch):
                    p.set_facecolor((0.80, 0.80, 0.80, 1.0))
                    p.set_edgecolor("black")
                    p.set_linewidth(1.5)
            for ln in getattr(a, "lines", []):
                if isinstance(ln, Line2D):
                    ln.set_linewidth(max(ln.get_linewidth(), 1.2))
                    ln.set_color("black")
        try:
            fig.canvas.draw_idle()
        except Exception:
            pass
        
    ax.set_xlabel(xlabel, fontsize=fontsize * 1.05)
    ax.set_ylabel(ylabel, fontsize=fontsize * 1.05)

    fig.suptitle(f"{title}", fontsize=fontsize * 1.05, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    outdir = os.path.join(out_path, region)
    os.makedirs(outdir, exist_ok=True)
    outfile = "{}_{}_parcoord_{}.{}".format(stat, region, metric, fig_format)
    fig.savefig(os.path.join(outdir, outfile), facecolor="w", bbox_inches="tight")
    plt.close(fig)
    return


