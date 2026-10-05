import os
import re
import glob
import json
import time
import numpy as np
import pandas as pd
import psutil
import matplotlib.pyplot as plt

import stat
import shutil
import datetime
import subprocess

import logging
import multiprocessing
import subprocess

import xarray as xr
import xcdat as xc
import pcmdi_metrics

from subprocess import Popen, PIPE
from itertools import chain
from typing import Dict, List, Tuple, Optional
from datetime import datetime
from copy import deepcopy

from collections import OrderedDict
from collections.abc import MutableMapping

from typing import Dict, List, Tuple, Optional, Literal, Sequence, Union, Iterable, Any
import tempfile
from pathlib import Path

from pcmdi_metrics.utils import sort_human, create_land_sea_mask

from pcmdi_metrics.io import (
    xcdat_open,
    base
)

from pcmdi_metrics.graphics import (
    Metrics,
    normalize_by_median,
    parallel_coordinate_plot,
    portrait_plot,
)

from pcmdi_metrics.enso.lib import (
    enso_portrait_plot
)

def count_child_processes(process=None):
    """
    Count the number of child processes for a given process.
    
    Parameters:
    - process (psutil.Process, optional): The process to check. If None, uses the current process.
    
    Returns:
    - int: Number of child processes.
    """
    if process is None:
        process = psutil.Process()
    
    children = process.children()
    return len(children)

def run_parallel_jobs(cmds: List[str], num_workers: int) -> List[Tuple[str, str, int]]:
    """
    Execute shell commands in parallel batches.

    Parameters:
    - cmds: List of command strings to run.
    - num_workers: Maximum number of subprocesses to run concurrently.

    Returns:
    - List of tuples: (stdout, stderr, return_code) for each command.
    """
    results = []
    procs = []

    for i, cmd in enumerate(cmds):
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
        procs.append((cmd, proc))

        # Run the batch if full or if it's the last command
        if len(procs) >= num_workers or i == len(cmds) - 1:
            print(f'Running {count_child_processes()} subprocesses...')
            for cmd, proc in procs:
                stdout, stderr = proc.communicate()
                return_code = proc.returncode

                if return_code != 0:
                    print(f"ERROR: Process failed: '{cmd}'\nError: {stderr.strip()}")
                    raise RuntimeError(f"Subprocess failed: {cmd}")

                results.append((stdout.strip(), stderr.strip(), return_code))

            time.sleep(0.25)  # Throttle before starting the next batch
            procs = []

    return results

def run_serial_jobs(cmds: List[str]) -> List[Tuple[str, str, int]]:
    """
    Execute shell commands one at a time (serially).

    Parameters:
    - cmds: List of command strings to run.

    Returns:
    - List of tuples: (stdout, stderr, return_code) for each command.
    """
    results = []

    for i, cmd in enumerate(cmds):
        print(f"Running [{i+1}/{len(cmds)}]: {cmd}")
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
        stdout, stderr = proc.communicate()
        return_code = proc.returncode

        if return_code != 0:
            print(f"ERROR: Process failed: '{cmd}'\nError: {stderr.strip()}")
            raise RuntimeError(f"Subprocess failed: {cmd}")

        results.append((stdout.strip(), stderr.strip(), return_code))

    return results

def setup_parallelization(num_workers=24):
    """Setup parallelization based on available workers."""
    available_workers = multiprocessing.cpu_count()
    logging.info(f"Available CPU cores: {available_workers}")
    
    if num_workers < 2 or num_workers > available_workers:
        multiprocessing_enabled = False
        logging.warning("Parallel processing disabled due to insufficient or excessive worker count.")
    else:
        multiprocessing_enabled = True
        logging.info(f"Parallel processing enabled with {num_workers} workers.")
    
    return multiprocessing_enabled

def shift_row_to_bottom(df, index_to_shift):
    """
    Moves the specified row to the bottom of the DataFrame and resets the index.

    Parameters:
        df (pd.DataFrame): The input DataFrame.
        index_to_shift (int): The index of the row to move to the bottom.

    Returns:
        pd.DataFrame: A new DataFrame with the row moved to the bottom and index reset.
    """
    if index_to_shift not in df.index:
        raise IndexError(f"Index {index_to_shift} not found in DataFrame.")

    df_top = df.drop(index=index_to_shift)
    df_bottom = df.loc[[index_to_shift]]

    new_df = pd.concat([df_top, df_bottom], ignore_index=True)
    return new_df

def derive_missing_variable(varin, path, model_id):
    """
    Derive variable with existing variables, preserving coordinates and attributes.

    Args:
        varin (str): Name of the derived variable (e.g., 'rstcre').
        path (str): Directory to look for/create the file.
        model_id (str): Identifier for constructing output filenames.
    """
    derived_var_map = {
        'rstcre': {'rsutcs': 1, 'rsut': -1},
        'rltcre': {'rlutcs': 1, 'rlut': -1},
    }

    if varin not in derived_var_map:
        return  # Nothing to derive

    var_dic = derived_var_map[varin]
    derived_data = None
    base_ds = None
    output_file = None


    for i, (src_var, scale) in enumerate(var_dic.items()):
        fpaths = sorted(glob.glob(os.path.join(path, f"*.{src_var}.*.nc")))
        if not fpaths:
            raise FileNotFoundError(f"No file found for source variable '{src_var}' in {path}")
        fpath = fpaths[0]
        ds = xcdat_open(fpath)
        data = ds[src_var] * scale

        if i == 0:
            base_ds = ds.copy(deep=True)
            derived_data = data.copy(deep=True)
            template = os.path.basename(fpath)
            output_file = os.path.join(path, template.replace(f".{src_var}.", f".{varin}."))
        else:
            derived_data = derived_data + data

    if base_ds is not None and derived_data is not None:
        derived_da = xr.DataArray(
            data=derived_data.data,
            coords=derived_data.coords,
            dims=derived_data.dims,
            attrs=derived_data.attrs
        )

        out_ds = base_ds.drop_vars(list(var_dic.keys()), errors="ignore")
        out_ds[varin] = derived_da

        # Optional: annotate global attributes
        out_ds.attrs.update({
            "derived_variable": varin,
            "derived_from": ", ".join(var_dic.keys()),
            "model_id": model_id
        })

        out_ds.to_netcdf(output_file)
        print(f"Derived variable '{varin}' written to {output_file}")

    return

def save_variable_regions(
        variables, 
        regions, 
        output_path="regions.json"
    ):
    """
    Maps each variable (simplified key) to a list of regions and saves to JSON.
    """
    region_map = OrderedDict()
    for var in variables:
        var_key = re.split(r"[_-]", var)[0] if "_" in var or "-" in var else var
        region_map[var_key] = regions

    with open(output_path, "w") as f:
        json.dump(region_map, f, sort_keys=False, indent=4, separators=(",", ": "))
    return region_map

def get_highlight_models(all_models, model_name):
    """
    Prioritize models containing 'e3sm' and then any additional specified models.

    Parameters:
        data_dict (dict): Dictionary with a 'model' key containing a list of model names.
        model_name (list): List of models to also highlight (after e3sm models).

    Returns:
        list: Ordered list of unique models to highlight.
    """
    highlight_model1 = []

    # First, collect all models that contain "e3sm" (case-insensitive)
    e3sm_models = [m for m in all_models if "e3sm" in m.lower()]

    # Then collect models in model_name that are not already in e3sm_models
    additional_models = [m for m in all_models if m in model_name and m not in e3sm_models]

    # Combine both lists
    highlight_model1 = e3sm_models + additional_models

    return highlight_model1

def generate_mean_clim_cmds(
        variables, 
        obs_dic, 
        case_id
    ):
    """
    Generates a list of shell commands for mean climate diagnostics.
    """
    commands = []
    for var in variables:
        var_key = re.split(r"[_-]", var)[0] if "_" in var or "-" in var else var
        if var_key in obs_dic:
            refset = obs_dic[var_key]["set"]
            cmd = " ".join([
                "mean_climate_driver.py",
                "-p parameterfile.py",
                "--vars", var,
                "-r", refset,
                "--case_id", case_id
            ])
            commands.append(cmd)
    return commands

def generate_varmode_cmds(
        modes,
        varOBS,
        reftyrs,
        reftyre,
        refname,
        refpath,
        case_id
    ):
    """Generates a list of command strings for variability modes processing."""

    # EOF mode overrides for specific variability modes (default is 1)
    eofn_map = {
        "NPO": 2,
        "NPGO": 2,
        "PSA1": 2,
        "PSA2": 3
    }

    commands = []

    for var_mode in modes:
        var_mode = var_mode.strip()
        # Use specified EOF number if in map, otherwise default to 1
        eofn = eofn_map.get(var_mode, 1)
        cmd = (
            f"variability_modes_driver.py -p parameterfile.py "
            f"--variability_mode {var_mode} "
            f"--eofn_mod {eofn} "
            f"--eofn_obs {eofn} "
            f"--varOBS {varOBS} "
            f"--osyear {reftyrs} "
            f"--oeyear {reftyre} "
            f"--reference_data_name {refname} "
            f"--reference_data_path {refpath} "
            f"--case_id {case_id}"
        )
        commands.append(cmd)

    return commands

def build_enso_obsvar_catalog(
        obs_dic: Dict, 
        variables: List[str], 
        output_file: str = "obs_catalogue.json"
    ) -> None:
    """
    Organize observational data for the ENSO driver based on the variable list.

    Parameters:
        obs_dic (dict): Dictionary mapping variable names to their observation sets and data files.
        variables (list): List of variable names to process.
        output_file (str): Output JSON file path to save the observation catalogue.
    """
    refr_dic = OrderedDict()

    for var in variables:
        vkey = re.split(r"[_-]", var)[0] if "_" in var or "-" in var else var

        if vkey not in obs_dic:
            raise KeyError(f"Variable key '{vkey}' not found in observation dictionary.")

        refset = obs_dic[vkey]['set']
        refname = obs_dic[vkey].get(refset)

        if not refname:
            raise KeyError(f"Reference name not found for variable '{vkey}' and set '{refset}'.")

        refr_dic.setdefault(refname, {})[vkey] = obs_dic[vkey][refname]

    with open(output_file, "w") as f:
        json.dump(refr_dic, f, indent=4, sort_keys=False, separators=(",", ": "))

    print(f"[INFO] Observation catalogue written to: {output_file}")

def build_enso_obsvar_landmask(
        obs_dic: Dict,
        variables: List[str],
        output_file: str = "obs_landmask.json",
        mask_dir: str = "fixed"
    ) -> None:
    """
    Organize observational land/sea mask mapping for ENSO diagnostics.

    Parameters:
        obs_dic (dict): Dictionary mapping variables to observation metadata.
        variables (list): List of variable names used in ENSO analysis.
        output_file (str): Path to output the landmask JSON.
        mask_dir (str): Directory prefix where the landmask files are located.
    """
    relf_dic = OrderedDict()

    for var in variables:
        vkey = re.split(r"[_-]", var)[0] if "_" in var or "-" in var else var

        if vkey not in obs_dic:
            raise KeyError(f"Variable key '{vkey}' not found in observation dictionary.")

        refset = obs_dic[vkey]['set']
        refname = obs_dic[vkey].get(refset)

        if not refname:
            raise KeyError(f"Reference name not found for variable '{vkey}' and set '{refset}'.")

        relf_dic.setdefault(refname, os.path.join(mask_dir, f"sftlf.{refname}.nc"))

    with open(output_file, "w") as f:
        json.dump(relf_dic, f, indent=4, sort_keys=False, separators=(",", ": "))

    print(f"[INFO] Landmask mapping written to: {output_file}")

def generate_enso_cmds(
        enso_groups_str, 
        case_id, 
        param_file="parameterfile.py", 
        driver_script="enso_driver.py"
    ):
    """
    Generate ENSO driver command-line strings for given metric groups.

    Parameters:
        enso_groups_str: Comma-separated list of ENSO metric groups.
        case_id: Case identifier.
        param_file: Parameter file used by the driver script.
        driver_script: ENSO driver script filename.

    Returns:
        cmds: List of shell command strings to run.
    """
    enso_groups = enso_groups_str.split(",")
    commands = [
        "{} -p {} --metricsCollection {} --case_id {}".format(driver_script, param_file, group, case_id)
        for group in enso_groups
    ]
    return commands

class ObservationLinker:
    def __init__(self, model_name, variables, obs_sets, ts_dir_ref_source, 
                 obstmp_dir, obs_alias_file, altobs_dic
        ):
        self.model_name = model_name
        self.variables = variables
        self.obs_sets = obs_sets
        self.ts_dir_ref_source = ts_dir_ref_source
        self.obstmp_dir = obstmp_dir
        self.obs_dic = json.load(open(obs_alias_file))
        self.altobs_dic = altobs_dic

    def _resolve_obs_file(self, varin, obsid):
        if varin not in self.obs_dic or obsid not in self.obs_dic[varin]:
            print(f"[Warning] No alias found for variable '{varin}' in obsid '{obsid}'")
            return None, None

        obsname = self.obs_dic[varin][obsid]
        obsstr = obsname.replace("_", "*").replace("-", "*") if "ceres_ebaf" in obsname else obsname
        pattern = os.path.join(self.ts_dir_ref_source, obsstr, f"{varin}_*.nc")
        fpaths = sorted(glob.glob(pattern))

        if fpaths and os.path.exists(fpaths[0]):
            return fpaths[0], varin

        # Try altobs mapping
        if varin in self.altobs_dic:
            alt_var = self.altobs_dic[varin]
            pattern_alt = os.path.join(self.ts_dir_ref_source, obsstr, f"{alt_var}_*.nc")
            fpaths = sorted(glob.glob(pattern_alt))
            if fpaths and os.path.exists(fpaths[0]):
                return fpaths[0], alt_var

        print(f"[Warning] Observation file not found for {varin} ({obsid})")
        return None, None

    def link_obs_data(self):
        for i, vv in enumerate(self.variables):
            varin = re.split(r"_|-", vv)[0] if "_" in vv or "-" in vv else vv
            if len(self.obs_sets) > 1 and len(self.obs_sets) == len(self.variables):
                obsid = self.obs_sets[i] 
            else:
                obsid = self.obs_sets[0]

            filepath, resolved_var = self._resolve_obs_file(varin, obsid)
            if filepath:
                template = os.path.basename(filepath)
                parts = template.replace(".nc", "").split("_")
                if len(parts) < 3:
                    print(f"[Error] Unexpected filename format: {template}")
                    continue
                yms, yme = parts[-2][:6], parts[-1][:6]
                obsname = self.obs_dic[varin][obsid].replace(".", "_")
                out = os.path.join(
                        self.obstmp_dir, 
                        f"{self.model_name.replace('%(model)', obsname)}.{varin}.{yms}-{yme}.nc"
                )

                if not os.path.exists(out):
                    os.makedirs(os.path.dirname(out), exist_ok=True)
                    if resolved_var == varin:
                        os.symlink(filepath, out)
                        print(f"[Info] Linked {resolved_var} → {out}")
                    else:
                        ds = xcdat_open(filepath)
                        ds = ds.rename({resolved_var: varin})
                        ds.to_netcdf(out)
                        print(f"[Info] Renamed and saved {resolved_var} as {varin} → {out}")
                else:
                    print(f"[Info] Skipping existing file: {out}")

    def derive_var(self, vout, var_dic):
        template = None
        out = None
        ds_out = None

        for i, (var, scale) in enumerate(var_dic.items()):
            fpaths = sorted(glob.glob(os.path.join(self.obstmp_dir, f"*.{var}.*.nc")))
            if not fpaths:
                print(f"[Warning] No file found for base variable '{var}' needed to derive '{vout}'")
                continue

            ds = xcdat_open(fpaths[0])
            if i == 0:
                template = os.path.basename(fpaths[0])
                out = os.path.join(self.obstmp_dir, template.replace(f".{var}.", f".{vout}."))
                shutil.copy(fpaths[0], out)
                ds_out = ds.rename_vars({var: vout})
                ds_out[vout] = ds_out[vout] * scale
            else:
                ds_other = xcdat_open(fpaths[0])
                ds_out[vout] = ds_out[vout] + ds_other[var] * scale

        if template and ds_out:
            ds_out.to_netcdf(out)
            print(f"[Info] Derived variable '{vout}' written to {out}")

    def process_derived_variables(self):
        for vv in self.variables:
            if vv in ['rltcre', 'rstcre']:
                fpaths = sorted(glob.glob(os.path.join(self.obstmp_dir, f"*{vv}_*.nc")))
                if not fpaths:
                    if vv == 'rstcre':
                        self.derive_var('rstcre', {'rsutcs': 1, 'rsut': -1})
                    elif vv == 'rltcre':
                        self.derive_var('rltcre', {'rlutcs': 1, 'rlut': -1})

class DataCatalogueBuilder:
    def __init__(self, test_path, test_set, ref_path, ref_set, variables, label, output_dir):
        self.test_path = test_path
        self.test_set = test_set
        self.ref_path = ref_path
        self.ref_set = ref_set
        self.variables = variables
        self.label = label
        self.output_dir = output_dir

        self.test_info = OrderedDict()
        self.ref_info = OrderedDict()

    def build_catalogues(self):
        for idx, var in enumerate(self.variables):
            varin = self._get_base_varname(var)
            test_files = sorted(glob.glob(os.path.join(self.test_path, f"*.{varin}.*.nc")))
            ref_files  = sorted(glob.glob(os.path.join(self.ref_path,  f"*.{varin}.*.nc")))

            if test_files and ref_files and os.path.exists(test_files[0]) and os.path.exists(ref_files[0]):
                for j, (fileset, info_dict, dataset, dataset_set) in enumerate([
                    (test_files[0], self.test_info, self.variables, self.test_set),
                    (ref_files[0],  self.ref_info,  self.variables, self.ref_set),
                ]):
                    metadata = self._extract_metadata(fileset, varin, var)
                    self._assign_metadata(info_dict, varin, var, dataset, dataset_set, idx, metadata)

        self._save_catalogue(self.test_path, self.test_info)
        self._save_catalogue(self.ref_path, self.ref_info)

        return self.test_info, self.ref_info

    def _get_base_varname(self, var):
        return re.split("_|-", var)[0] if ("_" in var or "-" in var) else var

    def _extract_metadata(self, filepath, varin, var):
        filename = os.path.basename(filepath)
        parts = filename.split(".")
        yymm_range = parts[6].split("-")
        return {
            "mip": parts[0],
            "exp": parts[1],
            "model": parts[2],
            "realization": parts[3],
            "tableID": parts[4],
            "yymms": yymm_range[0],
            "yymme": yymm_range[1],
            "var_in_file": varin,
            "var_name": var,
            "file_path": filepath,
            "template": filename
        }

    def _assign_metadata(self, target_dict, varin, var, dataset, dataset_names, idx, metadata):
        if varin not in target_dict:
            target_dict[varin] = {}
        kset = dataset_names[0] if len(dataset_names) != len(dataset) else dataset_names[idx]
        model = metadata["model"]

        target_dict[varin]['set'] = kset
        target_dict[varin][kset] = model
        target_dict[varin][model] = metadata

    def _save_catalogue(self, source_path, data_dict):
        filename = f"{source_path}_{self.label}_catalogue.json"
        filepath = os.path.join(self.output_dir, filename)
        with open(filepath, "w") as f:
            json.dump(data_dict, f, indent=4, sort_keys=False, separators=(",", ": "))

class LandSeaMaskGenerator:
    def __init__(self, test_path, ref_path, subsection, fixed_dir="fixed"):
        self.test_path = test_path
        self.ref_path = ref_path
        self.subsection = subsection
        self.fixed_dir = fixed_dir

    def run(self, enable_flag):
        if self._parse_flag(enable_flag):
            for group_path in [self.test_path, self.ref_path]:
                self._process_group(group_path)

    def _parse_flag(self, flag):
        return str(flag).lower() in ['true', 'y', 'yes']

    def _process_group(self, group):
        catalog_path = os.path.join(
            "pcmdi_diags",
            f"{group}_{self.subsection}_catalogue.json"
        )

        if not os.path.exists(catalog_path):
            print(f"Warning: Catalogue not found at {catalog_path}")
            return

        with open(catalog_path) as f:
            data_catalog = json.load(f)

        for var, meta in data_catalog.items():
            dataset = meta["set"]
            model = meta[dataset]
            input_file = meta[model]["file_path"]
            output_file = os.path.join(self.fixed_dir, f"sftlf.{model}.nc")

            if not os.path.exists(self.fixed_dir):
                os.makedirs(self.fixed_dir)

            if not os.path.exists(output_file):
                self._generate_mask(input_file, output_file, model)

    def _generate_mask(self, input_path, output_path, model_name):
        ds = xcdat_open(input_path, decode_times=True)
        ds = ds.bounds.add_missing_bounds()

        try:
            mask = create_land_sea_mask(ds, method="regionmask")
            print("Land mask estimated using regionmask method.")
        except Exception:
            mask = create_land_sea_mask(ds, method="pcmdi")
            print("Land mask estimated using PCMDI method.")

        mask = mask * 100.0
        mask.attrs.update({
            "long_name": "land_area_fraction",
            "units": "%",
            "id": "sftlf"
        })

        mask_ds = mask.to_dataset(name="sftlf").compute()
        mask_ds = mask_ds.bounds.add_missing_bounds()
        mask_ds = mask_ds.fillna(1.0e20)

        mask_ds.attrs.update({
            "model": model_name,
            "associated_files": input_path,
            "history": f"File processed: {datetime.now().strftime('%Y%m%d')}"
        })

        comp = dict(_FillValue=1.0e20, zlib=True, complevel=5)
        encoding = {v: comp for v in set(mask_ds.data_vars.keys()) | set(mask_ds.coords.keys())}

        mask_ds.to_netcdf(output_path, encoding=encoding)

        del ds, mask_ds, mask

class MeanClimateMetricsCollector:
    def __init__(self, regions, variables, fig_format, 
                 model_info, case_id, input_template, 
                 output_dir
        ):
        self.regions = regions
        self.variables = variables
        self.fig_format = fig_format
        self.mip, self.exp, self.model, self.relm = model_info
        self.case_id = case_id
        self.input_template = input_template
        self.output_dir = output_dir
        self.diag_metric = "mean_climate"
        self.seasons = ['DJF', 'MAM', 'JJA', 'SON', 'AC']
        self.model_name = f"{self.mip}.{self.exp}.{self.model}_{self.relm}"

    def collect(self):
        self._collect_figures()
        self._collect_metrics()
        self._collect_diags()

    def _collect_figures(self):
        fig_sets = OrderedDict()
        fig_sets['CLIM_patttern'] = ['graphics', '*']

        for fset, (fig_type, prefix) in fig_sets.items():
            for region in self.regions:
                for season in self.seasons:
                    for var in self.variables:
                        indir = self.input_template.replace('%(metric_type)', self.diag_metric)
                        indir = indir.replace('%(output_type)', fig_type)
                        search_path = os.path.join(
                                indir, var, f"{prefix}{region}_{season}*.{self.fig_format}"
                        )
                        fpaths = sorted(glob.glob(search_path))

                        for fpath in fpaths:
                            refname = os.path.basename(fpath).split("_")[0]
                            filname = f"{refname}_{region}_{season}.{self.fig_format}"
                            outpath = os.path.join(
                                self.output_dir.replace("%(group_type)", fset),
                                region, season
                            )
                            os.makedirs(outpath, exist_ok=True)
                            outfile = os.path.join(outpath, filname)
                            os.rename(fpath, outfile)

    def _collect_diags(self):
        inpath = self.input_template.replace('%(metric_type)', self.diag_metric)
        inpath = inpath.replace('%(output_type)', 'diagnostic_results')
        outpath = os.path.join(
                self.output_dir.replace('%(group_type)', 'metrics_data'), self.diag_metric
        )

        os.makedirs(outpath, exist_ok=True)
        fpaths = sorted(glob.glob(os.path.join(inpath, '*/*/*/*/*/*/*.nc')))

        for fpath in fpaths:
            filname = fpath.split("/")[-1]
            outfile = os.path.join(outpath, filname)
            os.rename(fpath, outfile)

    def _collect_metrics(self):
        inpath = self.input_template.replace('%(metric_type)', self.diag_metric)
        inpath = inpath.replace('%(output_type)', 'metrics_results')
        outpath = os.path.join(
                self.output_dir.replace('%(group_type)', 'metrics_data'), self.diag_metric
        )

        os.makedirs(outpath, exist_ok=True)
        fpaths = sorted(glob.glob(os.path.join(inpath, '*.json')))

        for fpath in fpaths:
            refname = os.path.basename(fpath).split("_")[:2]
            filname = f"{refname[0]}.{refname[1]}.{self.model_name}.{self.case_id}.json"
            outfile = os.path.join(outpath, filname)
            os.rename(fpath, outfile)

class VariabilityMetricsCollector:
    def __init__(self, modes, fig_format, mip, exp, model, relm,
                 case_id, input_dir, output_dir):
        self.modes = modes
        self.fig_format = fig_format
        self.mip = mip
        self.exp = exp
        self.model = model
        self.relm = relm
        self.case_id = case_id
        self.input_dir = input_dir.replace("%(metric_type)", "variability_modes")
        self.output_dir = output_dir
        self.model_name = f"{mip}.{exp}.{model}_{relm}"
        self.seasons = ['DJF', 'MAM', 'JJA', 'SON', 'yearly', 'monthly']
        self.fig_sets = OrderedDict({
            'MOV_eoftest': ['diagnostic_results', 'EG_Spec*'],
            'MOV_compose': ['graphics', '*compare_obs'],
            'MOV_telecon': ['graphics', '*teleconnection'],
            'MOV_pattern': ['graphics', '*']
        })

    def collect(self):
        self._collect_figures()
        self._collect_metrics()
        self._collect_diags() 

    def _collect_figures(self):
        for fig_set, (out_type, pattern_base) in self.fig_sets.items():
            for mode in self.modes:
                for season in self.seasons:
                    indir = self.input_dir.replace('%(output_type)', out_type)
                    template = (
                        f"{pattern_base}_{mode}_{season}*.{self.fig_format}"
                        if fig_set == "MOV_eoftest"
                        else f"{mode}_*_{season}_{pattern_base}.{self.fig_format}"
                    )
                    search_path = os.path.join(indir, mode, "*", template)
                    matched_files = sorted(glob.glob(search_path))

                    for fpath in matched_files:
                        filename = os.path.basename(fpath)
                        outfile = self._classify_output_name(fig_set, mode, season, filename)
                        outdir = os.path.join(
                            self.output_dir.replace("%(group_type)", "MOV_metric"),
                            fig_set, season
                        )
                        os.makedirs(outdir, exist_ok=True)
                        os.rename(fpath, os.path.join(outdir, outfile))

    def _classify_output_name(self, fig_set, mode, season, filename):
        suffix = "unknown"
        if "North_test" in filename:
            suffix = "EG_Spec"
        elif "_cbf_" in filename:
            suffix = "cbf"
        elif "EOF1" in filename:
            suffix = "eof1"
        elif "EOF2" in filename:
            suffix = "eof2"
        elif "EOF3" in filename:
            suffix = "eof3"
        return f"{fig_set}_{mode}_{season}_{suffix}.{self.fig_format}"

    def _collect_metrics(self):
        metrics_dir = self.input_dir.replace('%(output_type)', 'metrics_results')
        json_files = sorted(glob.glob(os.path.join(metrics_dir, '*/*/*.json')))

        for fpath in json_files:
            refmode = fpath.split("/")[-3]
            refname = fpath.split("/")[-2]
            reffile = fpath.split("/")[-1]

            eof_lookup = {"PSA1": "EOF2", "NPO": "EOF2", "NPGO": "EOF2", "PSA2": "EOF3"}
            refeof = eof_lookup.get(refmode, "EOF1")

            outdir = os.path.join(
                self.output_dir.replace('%(group_type)', 'metrics_data'),
                "variability_modes", refmode, refname
            )
            os.makedirs(outdir, exist_ok=True)

            base_name = f"var_mode_{refmode}.{refeof}.{self.model_name}.vs.{refname}.{self.case_id}"
            if 'diveDown' in reffile:
                outfile = os.path.join(outdir, f"{base_name}.diveDown.json")
            else:
                outfile = os.path.join(outdir, f"{base_name}.json")

            os.rename(fpath, outfile)

    def _collect_diags(self):
        diags_dir = self.input_dir.replace('%(output_type)', 'diagnostic_results')
        json_files = sorted(glob.glob(os.path.join(diags_dir, '*/*/*.nc')))

        for fpath in json_files:
            refmode = fpath.split("/")[-3]
            refname = fpath.split("/")[-2]
            reffile = fpath.split("/")[-1]

            eof_lookup = {"PSA1": "EOF2", "NPO": "EOF2", "NPGO": "EOF2", "PSA2": "EOF3"}
            refeof = eof_lookup.get(refmode, "EOF1")

            outdir = os.path.join(
                self.output_dir.replace('%(group_type)', 'metrics_data'),
                "variability_modes", refmode, refname
            )
            os.makedirs(outdir, exist_ok=True)

            outfile = os.path.join(outdir, reffile)

            os.rename(fpath, outfile)

class EnsoDiagnosticsCollector:
    def __init__(self, fig_format, refname, model_name_parts, case_id, input_dir, output_dir):
        self.fig_format = fig_format
        self.refname = refname
        self.mip, self.exp, self.model, self.relm = model_name_parts
        self.case_id = case_id
        self.model_name = f'{self.mip}.{self.exp}.{self.model}_{self.relm}'
        self.input_dir = input_dir.replace("%(metric_type)", "enso_metric")
        self.output_dir = output_dir
        self.diag_metric = "enso_metric"
        self.fig_sets = OrderedDict([("ENSO_metric", ['graphics', '*'])])

    def collect_figures(self, groups):
        for fset, (subdir, pattern) in self.fig_sets.items():
            for group in groups:
                fdir = self.input_dir.replace('%(output_type)', subdir)
                template = os.path.join(fdir, group, f"{pattern}.{self.fig_format}")
                fpaths = sorted(glob.glob(template))

                for fpath in fpaths:
                    tail = fpath.split("/")[-1].split(f"{self.model}_{self.relm}")[-1]
                    outpath = os.path.join(self.output_dir.replace("%(group_type)", fset), group)
                    os.makedirs(outpath, exist_ok=True)
                    outfile = f"{group}{tail}"
                    os.rename(fpath, os.path.join(outpath, outfile))

    def collect_metrics(self):
        inpath = self.input_dir.replace('%(output_type)', 'metrics_results')
        fpaths = sorted(glob.glob(os.path.join(inpath, '*/*.json')))

        for fpath in fpaths:
            refmode = fpath.split("/")[-2]
            reffile = fpath.split("/")[-1]
            outpath = os.path.join(
                self.output_dir.replace('%(group_type)', 'metrics_data'),
                self.diag_metric, refmode
            )
            os.makedirs(outpath, exist_ok=True)

            base_filename = f"{refmode}.{self.model_name}.vs.{self.refname}.{self.case_id}.json"
            outfile = base_filename.replace(".json", ".diveDown.json") if 'diveDown' in reffile else base_filename
            os.rename(fpath, os.path.join(outpath, outfile))

    def collect_diags(self):
        inpath = self.input_dir.replace('%(output_type)', 'diagnostic_results')
        fpaths = sorted(glob.glob(os.path.join(inpath, '*/*.nc')))

        for fpath in fpaths:
            refmode = fpath.split("/")[-2]
            reffile = fpath.split("/")[-1]
            outpath = os.path.join(
                self.output_dir.replace('%(group_type)', 'metrics_data'),
                self.diag_metric, refmode
            )
            os.makedirs(outpath, exist_ok=True)

            os.rename(fpath, os.path.join(outpath, reffile))

    def run(self, groups):
        self.collect_figures(groups)
        self.collect_metrics()
        self.collect_diags()

class SyntheticMetricsPlotter:
    def __init__(
        self,
        case_name,
        test_name,
        table_id,
        figure_format,
        figure_sets,
        metric_dict,
        save_data,
        base_test_input_path,
        results_dir=None,
        clim_vars=None,
        cmip_clim_dir=None,
        cmip_clim_set=None,
        mova_vars=None,
        movc_vars=None,
        cmip_movs_dir=None,
        cmip_movs_set=None,
        atm_modes=None,
        cpl_modes=None,
        enso_vars=None,
        cmip_enso_dir=None,
        cmip_enso_set=None,
        unit_check=True,
        badval_var_model=None,
        exclude_model_list=None,
        verbose=False,
        norm_method='default'
    ):
        self.case_name = case_name
        self.test_name = test_name
        self.table_id = table_id
        
        self.figure_format = figure_format
        self.figure_sets = figure_sets
        self.metric_dict = metric_dict
        
        self.save_data = save_data
        self.base_test_input_path = base_test_input_path
        self.results_dir = results_dir or "."
        
        self.clim_vars = clim_vars
        self.cmip_clim_dir = cmip_clim_dir
        self.cmip_clim_set = cmip_clim_set
        
        self.cmip_movs_dir = cmip_movs_dir
        self.cmip_movs_set = cmip_movs_set
        
        self.mova_vars = mova_vars
        self.movc_vars = movc_vars
        self.atm_modes = atm_modes.split(",") if atm_modes else []
        self.cpl_modes = cpl_modes.split(",") if cpl_modes else []
        
        self.enso_vars = enso_vars
        self.cmip_enso_dir = cmip_enso_dir
        self.cmip_enso_set = cmip_enso_set

        self.unit_check = unit_check
        self.badval_var_model = badval_var_model
        self.exclude_model_list = exclude_model_list
        self.verbose = verbose
        self.norm_method = norm_method
        
        self.parameter = self._initialize_parameter()

            
    def _initialize_parameter(self):
        parsed_test_names = []
        parsed_model_names = []

        for (raw_test,raw_case) in zip(self.test_name.split(","),self.case_name.split(",")):
            parts = raw_test.strip().split(".")
            if len(parts) != 4:
                raise ValueError(
                    f"Invalid test format '{raw_test}'. Expected 'a.b.c.d'"
                )

            # Construct strings
            test_id = f'{parts[0]}.{parts[1]}.{parts[2]}_{parts[3]}'
            parsed_test_names.append(test_id)

            parsed_model_names.append(raw_case)

        return OrderedDict({
            "save_data": self.save_data,
            "out_dir": os.path.join(self.results_dir, "ERROR_metric"),
            "test_name": parsed_test_names,
            "model_name": parsed_model_names,
            "tableID": [self.table_id],
        })

    def generate(self, metric_sets):
        print("Generating synthetic metrics plots ...")
        for metric in metric_sets:
            print(f"Processing metric: {metric}")
            self.parameter['test_path'] = self.base_test_input_path.replace('%(group_type)', metric)
            self.parameter['diag_vars'] = self.metric_dict[metric]

            if metric == "mean_climate":
                self._handle_mean_climate(metric)

            elif metric == "variability_modes":
                self._handle_variability_modes(metric)

            elif metric == "enso_metric":
                self._handle_enso_metric(metric)

    def _handle_mean_climate(self, metric):
        self.parameter.update({
            'cmip_path': self.cmip_clim_dir,
            'cmip_name': self.cmip_clim_set
        })
        
        # Instantiate the collector
        collector = ClimMetricsReader(
            self.parameter,
            unit_check=self.unit_check,
            badval_var_model=self.badval_var_model,       
            exclude_model_list=self.exclude_model_list,   
            verbose=self.verbose,
        )
        
        # Collect and merge metrics
        merge_lib = collector.collect()
        var_list = merge_lib.var_list
        var_unit_list = merge_lib.var_unit_list

        if self.clim_vars is not None:
            # Keep only variables requested in self.clim_vars
            # (and preserve order from self.clim_vars)
            var_list = [v for v in self.clim_vars if v in var_list]
            name_to_unit = dict(zip(merge_lib.var_list, merge_lib.var_unit_list))
            var_unit_list = [name_to_unit[v] for v in var_list]

        for stat, vars_ in self.metric_dict[metric].items():
            mean_climate_plot_driver(
                metric, stat,
                merge_lib.regions,
                self.parameter['model_name'],
                vars_,
                merge_lib.df_dict[stat],
                var_list,
                var_unit_list,
                self.parameter['save_data'],
                self.parameter['out_dir'],
                self.figure_format,
                self.norm_method
            )

    def _handle_variability_modes(self, metric):
        self.parameter.update({
            'cmip_path': self.cmip_movs_dir,
            'cmip_name': self.cmip_movs_set,
            'movs_mode': self.atm_modes + self.cpl_modes
        })

        reader = MoVsMetricsReader(self.parameter) 
        merge_lib, mode_season_list = reader.collect_metrics()
        
        for stat, vars_ in self.metric_dict[metric].items():
            variability_modes_plot_driver(
                metric, stat,
                self.parameter['model_name'],
                vars_,
                merge_lib[stat],
                mode_season_list,
                self.parameter['save_data'],
                self.parameter['out_dir'],
                self.figure_format
            )

    def _handle_enso_metric(self, metric):
        self.parameter.update({
            'cmip_path': self.cmip_enso_dir,
            'cmip_name': self.cmip_enso_set,
        })
        for stat in self.metric_dict[metric]:
            # Step 1: Load metrics JSON paths using the reader
            reader = EnsoMetricsReader(self.parameter,metric,stat)
            dict_json_path = reader.run()
            # Step 2: generate figures 
            enso_plot_driver(
                metric, stat,
                dict_json_path,
                self.parameter,
                self.figure_format
            )

def find_latest_file_list(
        path: str,
        file_pattern: str,
        var_pattern: str = r"\.(\w+)\.\d{8}\.nc$",
        time_pattern: str = r"\.(\d{8})\.nc$"
    ) -> List[str]:
    """
    Find the latest NetCDF file for each variable in the directory based on timestamps in filenames.

    Args:
        path (str): Directory to search.
        file_pattern (str): Regex to search file lists.
        var_pattern (str): Regex to extract variable name.
        time_pattern (str): Regex to extract date.

    Returns:
        List[str]: List of file paths, one for each variable (latest by timestamp).
    """
    latest_files = {}
    files = glob.glob(os.path.join(path, file_pattern))

    for f in files:
        fname = os.path.basename(f)
        var_match = re.search(var_pattern, fname)
        time_match = re.search(time_pattern, fname)

        if var_match and time_match:
            var = var_match.group(1)
            try:
                timestamp = datetime.strptime(time_match.group(1), "%Y%m%d")
            except ValueError:
                continue

            if var not in latest_files or timestamp > latest_files[var][0]:
                latest_files[var] = (timestamp, f)

    return [file for _, file in latest_files.values()]

class ClimMetricsMerger:
    def __init__(self, model_lib=None, cmip_lib=None, model_names=None):
        self.model_lib = model_lib or {}
        self.cmip_lib = cmip_lib or {}
        self.model_names = model_names or []
        self.merged_lib = None

    def merge(self):
        self._normalize_references()
        self._filter_regions()
        self._merge_and_standardize_units()
        self._highlight_and_sort_models()
        return self.merged_lib

    def _normalize_references(self):
        if hasattr(self.model_lib, 'references') and isinstance(self.model_lib.references, dict):
            self.model_lib.references = self._check_references(self.model_lib.references)
        if hasattr(self.cmip_lib, 'references') and isinstance(self.cmip_lib.references, dict):
            self.cmip_lib.references = self._check_references(self.cmip_lib.references)

    def _check_references(
        self,
        data_dict: MutableMapping[str, Optional[List[str]]],
        reference_alias: Optional[Dict[str, str]] = None
        ) -> MutableMapping[str, Optional[List[str]]]:
        if reference_alias is None:
            reference_alias = {
                'ceres_ebaf_toa_v4.1': 'ceres_ebaf_v4_1',
                'ceres_ebaf_toa_v4.0': 'ceres_ebaf_v4_0',
                'ceres_ebaf_toa_v2.8': 'ceres_ebaf_v2_8',
                'ceres_ebaf_surface_v4.1': 'ceres_ebaf_v4_1',
                'ceres_ebaf_surface_v4.0': 'ceres_ebaf_v4_0',
                'ceres_ebaf_surface_v2.8': 'ceres_ebaf_v2_8',
                'CERES-EBAF-4-1': 'ceres_ebaf_v4_1',
                'CERES-EBAF-4-0': 'ceres_ebaf_v4_0',
                'CERES-EBAF-2-8': 'ceres_ebaf_v2_8',
                'GPCP_v2.3': 'GPCP_v2_3',
                'GPCP_v2.2': 'GPCP_v2_2',
                'GPCP_v3.2': 'GPCP_v3_2',
                'GPCP-2-3': 'GPCP_v2_3',
                'GPCP-2-2': 'GPCP_v2_2',
                'GPCP-3-2': 'GPCP_v3_2',
                'NOAA_20C': 'NOAA-20C',
                'ERA-INT': 'ERA-Interim',
                'ERA-5': 'ERA5'
            }

        for key, values in data_dict.items():
            if isinstance(values, list):
                data_dict[key] = [reference_alias.get(val, val) for val in values]
            else:
                data_dict[key] = reference_alias.get(values, values)

        return data_dict

    def _filter_regions(self):
        self.model_lib, self.cmip_lib = self._check_regions(self.model_lib, self.cmip_lib)

    def _check_regions(self, data_lib, refr_lib):
        shared_regions = [region for region in data_lib.regions if region in refr_lib.regions]

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
                    if not df.empty and not df.isna().all().all()
                }
        return lib

    @staticmethod
    def _safe_merge_libs(lib1, lib2):
        """
        Merge two data libraries with nested dicts of DataFrames,
        gracefully handling missing or inconsistent keys, while
        avoiding FutureWarning due to all-NA/empty entries.
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

                    # Collect and clean valid DataFrames
                    valid_dfs = []
                    for df in (df1, df2):
                        if isinstance(df, pd.DataFrame) and not df.empty and not df.isna().all().all():
                            # Drop columns that are entirely NaN
                            df_clean = df.dropna(axis=1, how='all')
                            if not df_clean.empty and df_clean.shape[1] > 0:
                                valid_dfs.append(df_clean)

                    if valid_dfs:
                        merged_df = pd.concat(valid_dfs, ignore_index=True, sort=False)
                    else:
                        merged_df = pd.DataFrame()

                    merged.df_dict[stat][season][region] = merged_df

        return merged

    def _merge_and_standardize_units(self):
        # Prune empty or fully-NaN DataFrames from the model library
        cleaned_model_lib = self._prune_empty_dfs(self.model_lib)
        cleaned_model_lib = self._check_units(cleaned_model_lib)

        # Prune empty or fully-NaN DataFrames from the cmip library
        cleaned_cmip_lib = self._prune_empty_dfs(self.cmip_lib)
        cleaned_cmip_lib = self._check_units(cleaned_cmip_lib)

        # Safe merge with fallback for missing stats/seasons/regions
        self.merged_lib = self._safe_merge_libs(cleaned_cmip_lib, cleaned_model_lib)

        # Standardize units after merging
        self.merged_lib = self._check_units(self.merged_lib)

    def _check_units(self, data_lib, verbose=False):
        units_all = {
            "prw": "[kg m$^{-2}$]", "pr": "[mm d$^{-1}$]", "prsn": "[mm d$^{-1}$]",
            "prc": "[mm d$^{-1}$]", "hfls": "[W m$^{-2}$]", "hfss": "[W m$^{-2}$]",
            "clivi": "[kg $m^{-2}$]", "clwvi": "[kg $m^{-2}$]", "psl": "[Pa]",
            "rlds": "[W m$^{-2}$]", "rldscs": "[W $m^{-2}$]", "evspsbl": "[kg m$^{-2} s^{-1}$]",
            "rtmt": "[W m$^{-2}$]", "rsdt": "[W m$^{-2}$]", "rlus": "[W m$^{-2}$]",
            "rluscs": "[W m$^{-2}$]", "rlut": "[W m$^{-2}$]", "rlutcs": "[W m$^{-2}$]",
            "rsds": "[W m$^{-2}$]", "rsdscs": "[W m$^{-2}$]", "rstcre": "[W m$^{-2}$]",
            "rltcre": "[W m$^{-2}$]", "rsus": "[W m$^{-2}$]", "rsuscs": "[W m$^{-2}$]",
            "rsut": "[W m$^{-2}$]", "rsutcs": "[W m$^{-2}$]", "ts": "[K]",
            "tas": "[K]", "tauu": "[Pa]", "tauv": "[Pa]",
            "zg-500": "[m]", "ta-200": "[K]", "sfcWind": "[m s$^{-1}$]",
            "ta-850": "[K]", "ua-200": "[m s$^{-1}$]", "ua-850": "[m s$^{-1}$]",
            "va-200": "[m s$^{-1}$]", "va-850": "[m s$^{-1}$]", "uas": "[m s$^{-1}$]",
            "vas": "[m s$^{-1}$]", "tasmin": "[K]", "tasmax": "[K]", "clt": "[%]"
        }

        # Identify common variables and handle aliases like 'rt' or 'rmt'
        common_vars = [var for var in data_lib.var_list if var in units_all]
        if 'rtmt' not in common_vars and any(var in data_lib.var_list for var in ['rt', 'rmt']):
            common_vars.append('rtmt')

        # Collect units for these variables
        common_unts = [units_all[var] for var in common_vars if var in units_all]

        # Filter and correct reference list
        new_var_ref_dict = {}
        for var, ref in data_lib.var_ref_dict.items():
            if var in common_vars:
                new_var_ref_dict[var] = ref
            elif var in ['rt', 'rmt']:
                new_var_ref_dict['rtmt'] = ref
                if verbose:
                    print(f"Alias {var} mapped to 'rtmt' in references.")
        
        data_lib.var_ref_dict = self._check_references(new_var_ref_dict)

        # Clean DataFrames
        for stat, seasons in data_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in regions.items():
                    df = df.copy()
                    # Handle aliases
                    if 'rt' in df.columns:
                        df['rtmt'] = df['rt']
                    elif 'rmt' in df.columns:
                        df['rtmt'] = df['rmt']

                    # Drop irrelevant variables
                    drop_cols = [var for var in df.columns[3:] if var not in common_vars]
                    if drop_cols and verbose:
                        print(f"Dropping variables in {stat}/{season}/{region}: {drop_cols}")
                    df = df.drop(columns=drop_cols)
                    data_lib.df_dict[stat][season][region] = df

        data_lib.var_list = common_vars
        data_lib.var_unit_list = common_unts

        return data_lib

    def _highlight_and_sort_models(self):
        for stat, seasons in self.merged_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in regions.items():
                    df = pd.DataFrame(df)
                    highlight_models = get_highlight_models(df.get('model', []), self.model_names)
                    for model in highlight_models:
                        for idx in df[df["model"] == model].index:
                            df = shift_row_to_bottom(df, idx)
                    self.merged_lib.df_dict[stat][season][region] = df.fillna(np.nan)

class ClimMetricsReader:
    def __init__(
        self,
        parameter: Dict,
        unit_check: bool = True,
        verbose: bool = False,
        badval_var_model: Optional[Dict[str, List[str]]] = None,
        exclude_model_list: Optional[List[str]] = None,
    ):
        """
        Initialize the climate metrics collector.

        Args
        ----
        parameter : dict
            Required paths and identifiers.
        unit_check : bool
            Run unit consistency checks.
        verbose : bool
            Print extra progress info.
        badval_var_model : dict[str, list[str]] | None
            Mapping for known-bad variable values.
        exclude_model_list : list[str] | None
            Models to exclude from all DataFrames automatically.
        """
        self.parameter = parameter
        self.unit_check = unit_check
        self.verbose = verbose
        self.badval_var_model = badval_var_model
        self.exclude_model_list = exclude_model_list or []   # <- new

        self.cmip_lib = None
        self.all_lib = None
        self.all_names: List[str] = []

        self.var_pattern = re.compile(r"^([A-Za-z0-9\-]+)\.")
        self.time_pattern = re.compile(r"\.v(\d{8})\.json$")

        self._validate_parameter()

    # ---- public API ---------------------------------------------------------

    def collect(self):
        self._load_cmip_metrics()

        test_names = self.parameter["test_name"]
        model_names = self.parameter["model_name"]
        if len(test_names) != len(model_names):
            raise ValueError(
                f"'test_name' (n={len(test_names)}) and 'model_name' (n={len(model_names)}) must match length."
            )

        for test_name, model_name in zip(test_names, model_names):
            print(f"process metrics data for {model_name}: {test_name} ")
            model_lib = self._process_test_model(test_name, model_name)
            self.all_lib = model_lib.copy() if self.all_lib is None else self.all_lib.merge(model_lib)
            self.all_names.append(model_name)

        if self.verbose:
            print("Merging model metrics with CMIP reference metrics...")

        merger = ClimMetricsMerger(
            model_lib=self.all_lib,
            cmip_lib=self.cmip_lib,
            model_names=self.all_names,
        )
        return merger.merge()

    # ---- internals ----------------------------------------------------------

    def _validate_parameter(self) -> None:
        required = ["cmip_name", "cmip_path", "test_path", "test_name", "model_name"]
        missing = [k for k in required if k not in self.parameter]
        if missing:
            raise KeyError(f"Missing required parameter keys: {missing}")
        if not isinstance(self.parameter["test_name"], (list, tuple)) or not isinstance(
            self.parameter["model_name"], (list, tuple)
        ):
            raise TypeError("'test_name' and 'model_name' must be lists/tuples.")

    def _load_clim_metrics_from_files(
        self,
        file_paths: Sequence[str],
        *,
        badval_var_model: Optional[Dict[str, List[str]]] = None,
        exclude_model_list: Optional[List[str]] = None,
    ):
        """
        Load metrics -> apply bad-value nulling -> optional model exclusion.
        """
        if not file_paths:
            raise FileNotFoundError("No metric files provided to _load_clim_metrics_from_files().")

        lib = Metrics(file_paths)

        # 1) Null out known-bad values (uses per-instance default unless overridden)
        lib = self._check_badvals(
            lib,
            var_model=(badval_var_model if badval_var_model is not None else self.badval_var_model),
            verbose=self.verbose,
        )

        # 2) Drop selected models if specified (param wins; else use instance default)
        drop_models = exclude_model_list if exclude_model_list is not None else self.exclude_model_list
        if drop_models:
            lib = self.exclude_models(lib, drop_models)

        return lib

    def _load_cmip_metrics(self) -> None:
        cmip_id = self.parameter["cmip_name"]
        try:
            part0, part1, part2 = cmip_id.split(".")
        except ValueError:
            raise ValueError(
                f"parameter['cmip_name'] must look like 'MIP.activity.table', got: {cmip_id}"
            )

        cmip_dir = os.path.join(self.parameter["cmip_path"], part0, part1, part2)
        cmip_files = sorted(glob.glob(os.path.join(cmip_dir, f"*.{part2}.json")))
        if not cmip_files:
            raise FileNotFoundError(f"No CMIP metrics found in: {cmip_dir}")

        if self.verbose:
            print(f"Loading CMIP metrics from {len(cmip_files)} files in: {cmip_dir}")

        self.cmip_lib = self._load_clim_metrics_from_files(cmip_files)

    def _process_test_model(self, test_name: str, model_name: str):
        test_key = self._extract_test_key(test_name)
        test_path = self.parameter["test_path"].replace("%(model_name)", model_name)

        model_files = find_latest_file_list(
            path=test_path,
            file_pattern="*.v*.json",
            var_pattern=self.var_pattern,
            time_pattern=self.time_pattern,
        )
        if not model_files:
            raise FileNotFoundError(
                f"No mean climate metrics found for model '{model_name}' (searched in {test_path})"
            )

        if self.verbose:
            print(f"Reading metrics for model '{model_name}' from {len(model_files)} file(s)...")

        # Write remapped copies to a temp dir so Metrics reads the same structure your
        # original code produced—without touching the originals.
        tmpdir_ctx = tempfile.TemporaryDirectory(prefix="metrics_sidecars_")
        tmpdir = Path(tmpdir_ctx.name)

        sidecars = []
        for fp in model_files:
            try:
                with open(fp, "r") as f:
                    data = json.load(f)
                changed, new_data = self._remap_results_key(data, old_key=test_key, new_key="default")
                out_fp = tmpdir / (Path(fp).name.replace(".json", ".default.json") if changed else Path(fp).name)
                with open(out_fp, "w") as g:
                    json.dump(new_data, g, indent=2)
                sidecars.append(str(out_fp))
            except (OSError, json.JSONDecodeError) as e:
                print(f"Warning: Could not load {fp}: {e}")

        if not sidecars:
            raise FileNotFoundError(
                f"All discovered metric files failed to load for model '{model_name}'."
            )

        model_lib = self._load_clim_metrics_from_files(sidecars)

        # Standardize model names as before
        for stat, seasons in model_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, df in regions.items():
                    if not isinstance(df, pd.DataFrame):
                        df = pd.DataFrame(df)
                    if "model" in df.columns:
                        df.loc[:, "model"] = model_name
                    model_lib.df_dict[stat][season][region] = df

        # Keep the tempdir alive until the lib is fully used, or attach it to the object
        # (e.g., store tmpdir_ctx on self if you need it to persist longer)
        self._last_tmpdir_ctx = tmpdir_ctx  # prevent GC until reader is GC'd
        return model_lib

    def _extract_test_key(self, test_name: str) -> str:
        parts = test_name.split(".")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
        for p in parts:
            if p:
                return p
        return test_name

    @staticmethod
    def _remap_results_key(
        data: Dict, old_key: str, new_key: str = "default"
    ) -> Tuple[bool, Dict]:
        if not isinstance(data, dict):
            return False, data
        results = data.get("RESULTS")
        if not isinstance(results, dict) or old_key not in results:
            return False, data
        new_data = dict(data)
        new_results = dict(results)
        if new_key not in new_results:
            new_results[new_key] = new_results.pop(old_key)
        new_data["RESULTS"] = new_results
        return True, new_data

    # ---- your merged bad-value method --------------------------------------
    def _check_badvals(
        self,
        data_lib,
        var_model: Optional[Dict[str, List[str]]] = None,
        *,
        verbose: bool = False,
    ):
        """
        For each DataFrame in data_lib.df_dict[stat][season][region], set selected
        variable columns to NaN for specified models.

        Returns the same data_lib (mutated in place).
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
                            print(f"[WARN] Missing 'model' column for ({stat}, {season}, {region}); skipping.")
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    if not set(df["model"].unique()).intersection(target_models):
                        data_lib.df_dict[stat][season][region] = df
                        continue

                    cols_to_masks: Dict[str, np.ndarray] = {}

                    for model_name, cols in var_model.items():
                        model_mask = (df["model"] == model_name).to_numpy()
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
                        print(f"[INFO] ({stat}, {season}, {region}) nulled {total} values -> {detail}")

                    data_lib.df_dict[stat][season][region] = df

        return data_lib
    
    def exclude_models(self, data_lib, model_list: List[str]):
        """
        Exclude models listed in `model_list` from the DataFrames stored in data_lib.df_dict.
        """
        for stat, seasons in data_lib.df_dict.items():
            for season, regions in seasons.items():
                for region, table in regions.items():
                    df = pd.DataFrame(table)
                    mask = ~df.apply(lambda row: row.isin(model_list).any(), axis=1)
                    df = df.loc[mask].reset_index(drop=True)
                    data_lib.df_dict[stat][season][region] = df
        return data_lib

class MoVsMetricsReader:
    def __init__(self, parameter):
        self.parameter = parameter
        self.cmip_group, self.cmip_model, self.cmip_version = self.parameter['cmip_name'].split(".")
        self.movs_mode = parameter['movs_mode']
        self.var_pattern = re.compile(r"var_mode_(\w+)\.EOF\d+\..*\.json$")
        self.time_pattern = re.compile(r"\.v(\d{8})\.json$")

    def collect_metrics(self):
        cmip_files = self._get_cmip_files()
        if not cmip_files or not os.path.exists(cmip_files[0]):
            raise FileNotFoundError("ERROR: No Synthetic MoVs Metrics Data For CMIP, Aborting.")

        print("Found Synthetic MoVs Metrics Data For CMIP, Reading...")
        cmip_lib = self._load_movs_files(cmip_files)

        merge_lib = {}
        for stat, diag_vars in self.parameter['diag_vars'].items():
            merge_df, mode_season_list = self._movs_dict_to_df(cmip_lib, stat)

            for i, model_name in enumerate(self.parameter['model_name']):
                model_path = self.parameter['test_path'].replace("%(model_name)", model_name)
                model_files = find_latest_file_list(
                        path=f'{model_path}/*/*',
                        file_pattern="var_mode_*.json",
                        var_pattern=self.var_pattern,
                        time_pattern=self.time_pattern
                )
                if not model_files or not os.path.exists(model_files[0]):
                    raise FileNotFoundError(f"No Synthetic MoVs Metrics Data For {model_name}, Aborting.")

                print(f"Found Synthetic MoVs Metrics for {model_name}, Reading...")
                model_lib = self._load_movs_files(model_files)
                
                # Normalize model name key to match targets
                model_lib = {
                        mode: {model_name: next(iter(model_data.values()))}
                        for mode, model_data in model_lib.items()
                }

                # Convert dictionary to DataFrame
                model_df, _ = self._movs_dict_to_df(model_lib, stat)
                
                # Append to the merged DataFrame
                merge_df = pd.concat([merge_df, model_df], ignore_index=True)

            # Highlight and reorder models if applicable
            highlight_models = get_highlight_models(merge_df.get('model', []), self.parameter['model_name'])
            for model in merge_df["model"].tolist():
                if model in highlight_models:
                    for idx in merge_df[merge_df["model"] == model].index:
                        merge_df = shift_row_to_bottom(merge_df, idx)

            merge_lib[stat] = merge_df

        return merge_lib, mode_season_list

    def _get_cmip_files(self):
        return glob.glob(os.path.join(
            self.parameter['cmip_path'],
            self.cmip_group, self.cmip_model, self.cmip_version,
            "*/*/var_mode_*.json"
        ))

    def _load_movs_files(self, file_lists):
        json_lib = {}
        for mode in self.movs_mode:
            eof = {'PSA1': 'EOF2', 'NPO': 'EOF2', 'NPGO': 'EOF2', 'PSA2': 'EOF3'}.get(mode, 'EOF1')
            for json_file in file_lists:
                if mode in json_file and eof in json_file:
                    try:
                        with open(json_file, 'r') as fj:
                            data = json.load(fj)
                            json_lib[mode] = data.get('RESULTS', {})
                    except (FileNotFoundError, json.JSONDecodeError) as e:
                        print(f"Warning: Could not load {json_file}: {e}")
                    break
        return json_lib
    
    def _movs_dict_to_df(self, movs_dict, stat):
        models = sorted(movs_dict.get('NAM', {}).keys())
        df = pd.DataFrame({'model': models, 'num_runs': np.nan})
        mode_season_list = []

        for mode in self.movs_mode:
            seasons = (
                ['monthly'] if mode in ['PDO', 'NPGO'] else
                ['yearly'] if mode == 'AMO' else
                ['DJF', 'MAM', 'JJA', 'SON']
            )

            for season in seasons:
                col_name = f"{mode}_{season}"
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
                                run_stat = movs_dict[mode][model][run]['defaultReference'][mode][season]['cbf'][stat]
                                stat_values.append(run_stat)
                            except KeyError:
                                continue

                        if stat_values:
                            value = np.mean(stat_values)
                            num_runs = len(stat_values)

                    df.at[idx, col_name] = value
                    if np.isnan(df.at[idx, 'num_runs']):
                        df.at[idx, 'num_runs'] = num_runs
                    elif num_runs > 0:
                        df.at[idx, 'num_runs'] = max(df.at[idx, 'num_runs'], num_runs)

        return df, mode_season_list

def archive_data(region, stat, season, data_dict, model_name, var_names, var_units, outdir):
    """
    Archive processed data into a CSV file with variable units in column headers if available.

    Parameters:
        region (str): Region name.
        stat (str): Statistic type (e.g., mean, std).
        season (str): Season name.
        data_dict (dict or DataFrame): Data to archive.
        model_name (str): Model identifier.
        var_names (list): List of variable names.
        var_units (list): List of variable units (optional, same order as var_names).
        outdir (str): Directory to save the CSV file.
    """
    df = pd.DataFrame(data_dict)

    # Determine the index of the first variable column (assumes first 3 are metadata)
    metadata_cols = df.columns[:3].tolist()
    variable_cols = df.columns[3:]

    filtered_cols = []
    new_column_names = df.columns.tolist()

    for var in variable_cols:
        if var in var_names:
            filtered_cols.append(var)
            if var_units:
                idx = df.columns.get_loc(var)
                unit_label = var_units[var_names.index(var)]
                new_column_names[idx] = f"{var} ({unit_label})"

    # Subset dataframe and rename columns if units provided
    df = df[metadata_cols + filtered_cols]
    df.columns = new_column_names[:len(df.columns)]

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Construct and save the output filename
    outfile = f"{stat}_{region}_{season}_{model_name}.csv"
    df.to_csv(os.path.join(outdir, outfile), index=False)

    return

def drop_vars(data_dict, var_names, var_units=None):
    """
    Drop variables (columns) from data_dict where more than 90% of the values are NaN.

    Parameters:
        data_dict (pd.DataFrame): Data containing variable columns.
        var_names (list): List of variable names matching data_dict columns.
        var_units (list, optional): List of units for variables. Must match var_names in order.

    Returns:
        Tuple of (filtered_data_dict, updated_var_names, updated_var_units)
    """
    protected_columns = {'model', 'run', 'model_run', 'num_runs'}
    columns_to_drop = []

    for column in data_dict.columns:
        if column in protected_columns:
            continue
        nan_ratio = data_dict[column].isna().mean()
        if nan_ratio > 0.9:
            columns_to_drop.append(column)

    # Drop columns from DataFrame
    data_dict = data_dict.drop(columns=columns_to_drop)

    # Update var_names and var_units if applicable
    updated_var_names = [v for v in var_names if v not in columns_to_drop]
    updated_var_units = None
    if var_units is not None:
        # Keep units only for remaining variables
        name_to_unit = dict(zip(var_names, var_units))
        updated_var_units = [name_to_unit[v] for v in updated_var_names if v in name_to_unit]

    return data_dict, updated_var_names, updated_var_units

class EnsoMetricsReader:
    def __init__(self, parameter, metric, stat):
        self.parameter = parameter
        self.metric = metric
        self.stat = stat
        self.metric_dict = self.parameter['diag_vars'][stat]
        self.metrics_collections = self.metric_dict['collection']
        self.mips = [self.parameter['cmip_name'].split(".")[0]] + self.parameter['model_name']
        self.dict_json_path = {}
        
        self.var_pattern = re.compile(r"\.(\w+)\..*\.v(\d{8})\.json$")
        self.time_pattern = re.compile(r"\.v(\d{8})\.json$")

    def run(self):
        """Collect paths to ENSO metrics JSON files and return the mapping."""
        for mip in self.mips:
            self.dict_json_path[mip] = {}
            for metrics_collection in self.metrics_collections:
                if 'cmip' in mip:
                    self.dict_json_path[mip][metrics_collection] = self._get_cmip_json_path(mip, metrics_collection)
                else:
                    self.dict_json_path[mip][metrics_collection] = self._get_test_json_path(mip, metrics_collection)

            if len(self.dict_json_path[mip]) < 1:
                raise FileNotFoundError(f"No Synthetic ENSO Metrics Data for {mip}, aborting...")

        return self.dict_json_path

    def _get_cmip_json_path(self, mip, metrics_collection):
        path = os.path.join(
            self.parameter['cmip_path'],
            self.parameter['cmip_name'].split(".")[0],
            self.parameter['cmip_name'].split(".")[1],
            self.parameter['cmip_name'].split(".")[2],
            metrics_collection,
            f"{mip.lower()}_{self.parameter['cmip_name'].split('.')[1]}_{metrics_collection}_*.json"
        )
        matches = glob.glob(path)
        if not matches:
            raise FileNotFoundError(f"CMIP metrics file not found for {mip} and {metrics_collection}")
        return matches[0]

    def _get_test_json_path(self, mip, metrics_collection):
        for i, model_name in enumerate(self.parameter['model_name']):
            model_path = self.parameter['test_path'].replace("%(model_name)", model_name)
            model_files = find_latest_file_list(
                    path=f'{model_path}/{metrics_collection}',
                    file_pattern="*.json",
                    var_pattern=self.var_pattern,
                    time_pattern=self.time_pattern
            )
            print(f'{model_path}/{metrics_collection}')
            if not model_files or not os.path.exists(model_files[0]):
                raise FileNotFoundError(f"No Synthetic ENSO Metrics Data For {mip} {model_name}, Aborting.")
            
            for json_path in model_files:
                with open(json_path) as ff:
                    data_json = json.load(ff)

            old_key = list(data_json["RESULTS"]["model"].keys())[0]
            
            data_json["RESULTS"]["model"][mip] = data_json["RESULTS"]["model"].pop(old_key)
            
            with open(json_path, 'w', encoding='utf8') as ff:
                json.dump(data_json, ff, indent=4, separators=(",", ": "), sort_keys=True)

        return json_path

def enso_plot_driver(metric,stat,dict_json_path,parameter,fig_format):
    """
    Driver function to plot ENSO metrics based on specified type (e.g., portrait).
    """
    metric_dict = parameter['diag_vars'][stat]
    metrics_collections = metric_dict['collection']
    mips = [parameter['cmip_name'].split(".")[0]] + parameter['model_name']

    for mtype in metric_dict['type']:
        if mtype == "portrait":
            print(f"Processing Portrait Plots for {metric} {stat}...")

            list_project = mips
            list_obs = []  # fill in if observational references are needed
            outdir = os.path.join(parameter['out_dir'], metric)
            os.makedirs(outdir, exist_ok=True)

            outfile = f"{metric}_{stat}_portrait.{fig_format}"
            figure_name = os.path.join(outdir, outfile)

            fig, ref_info_dict = enso_portrait_plot(
                metrics_collections,
                list_project,
                list_obs,
                dict_json_path,
                figure_name=figure_name,
                reduced_set=True
            )

    return


def variability_modes_plot_driver(
        metric,stat,model_name,
        metric_dict,df_dict,
        mode_season_list,
        save_data,out_path,
        fig_format
    ):
    """Driver Function for the modes variability metrics plot"""
    season = "mon"
    if len(model_name) > 1:
        mout_name = model_name[0].split("_")[0]
    else:
        mout_name = model_name[0]

    for mtype in metric_dict['type']:
        if mtype == "portrait":
            print("Processing Portrait  Plots for {} {}....".format(metric,stat))
            if stat not in ["stdv_pc_ratio_to_obs"]:
                data_nor = normalize_by_median(
                        df_dict[mode_season_list].to_numpy().T, axis=1)
            else:
                data_nor = df_dict[mode_season_list].to_numpy().T
            if save_data:
                df_dict[mode_season_list] = data_nor.T
                outdir = os.path.join(out_path,metric)
                archive_data(metric,stat,season,df_dict,mout_name,mode_season_list,None,outdir)
            run_list = df_dict['model'].to_list() 
            stat_name = metric_dict['name']
            portrait_metric_plot(metric,stat,season,data_nor,
                                 stat_name,model_name,mode_season_list,
                                 run_list,out_path,fig_format)
        elif mtype == "parcoord":
            print("Processing Parallel Coordinate Plots for {} {}....".format(metric,stat))
            #drop data if all is NaNs
            data_dict,var_names,var_units = drop_vars(df_dict.copy(),mode_season_list.copy(),None)
            if save_data:
                outdir = os.path.join(out_path,metric)
                archive_data(metric,stat,season,data_dict,mout_name,mode_season_list,None,outdir)
            run_list = data_dict['model'].to_list() 
            stat_name = metric_dict['name']
            parcoord_metric_plot(metric,stat,season,data_dict,
                                 stat_name,model_name,var_names,var_units,
                                 run_list,out_path,fig_format)

    return 


def normalize_to_reference(
    df: pd.DataFrame,
    var_names,
    *,
    mode: str = "default",            # "default" (PCMDI) or "cmip_median"
    target: Optional[Union[str, List[str]]] = "E3SM", # regex/substring(s) to exclude from baseline
    model_col: str = "model",         # column name holding model IDs
    use_mad: bool = False,            # only used when mode=="cmip_median"
    case_insensitive: bool = False,   # regex matching flag
    baseline_mask: Optional[pd.Series] = None,   # advanced: True for rows to use in baseline
) -> pd.DataFrame:
    """
    Normalize selected metric columns and return a new DataFrame.

    - mode == "default": PCMDI normalize_by_median over ALL rows (exact call preserved)
    - mode == "cmip_median": compute baseline on rows NOT matching `target` (or per `baseline_mask`),
                             then apply to ALL rows; optionally use robust MAD z-scores.
    """
    if isinstance(var_names, str):
        var_names = [var_names]

    missing = [c for c in ([model_col] + list(var_names)) if c not in df.columns]
    if missing:
        raise KeyError(f"normalize_to_reference: missing required columns: {missing}")

    out = df.copy()

    if mode == "default":
        # ---- PCMDI default, bit-for-bit ----
        arr = df[var_names].to_numpy().T  # (N_vars, N_models)
        arr_norm = normalize_by_median(arr, axis=1)
        out.loc[:, var_names] = arr_norm.T
        return out

    if mode != "cmip_median":
        raise ValueError("mode must be 'default' or 'cmip_median'")

    # ---- CMIP-only baseline path ----
    # Compute which rows to use for the baseline
    if baseline_mask is not None:
        # True means "use in baseline"
        if baseline_mask.shape[0] != len(out):
            raise ValueError("baseline_mask length must match df length")
        use_rows = baseline_mask.astype(bool)
    else:
        if target is None:
            # If no target provided, fallback to using all rows (same as default baseline rows)
            use_rows = pd.Series(True, index=out.index)
        else:
            # Build a robust regex for one or many targets
            if isinstance(target, list):
                pat = "|".join(re.escape(t) for t in target)
            else:
                pat = target  # treat as a (possibly pre-escaped) regex or substring
            flags = 0 if not case_insensitive else re.IGNORECASE
            # rows NOT matching target are considered CMIP (baseline rows)
            matches = out[model_col].astype(str).str.contains(pat, regex=True, case=not case_insensitive)
            use_rows = ~matches

    cmip = out.loc[use_rows, var_names]
    if cmip.empty:
        raise ValueError(
            "normalize_to_reference: baseline subset is empty. "
            f"Check `target={target!r}` and `model_col='{model_col}'`."
        )

    med = cmip.median(axis=0, skipna=True)

    if use_mad:
        mad = (cmip - med).abs().median(axis=0, skipna=True)
        scale = (1.4826 * mad).replace(0.0, np.nan).fillna(1e-12)
        normed = (out[var_names] - med) / scale
    else:
        # (x - median) / median
        scale = med.replace(0.0, np.nan).fillna(1e-12)
        normed = (out[var_names] - med) / scale

    # Avoid infs if anything still slipped through
    normed = normed.where(np.isfinite(normed), np.nan)
    out.loc[:, var_names] = normed
    return out

def mean_climate_plot_driver(
    metric, stat, regions, model_name,
    metric_dict, df_dict,
    var_list, var_unit_list,
    save_data, out_path,
    fig_format,
    norm_method="default",
):
    """Driver Function for the mean climate metrics plot"""

    # Keep existing output-name logic
    mout_name = model_name[0].split("_")[0] if len(model_name) > 1 else model_name[0]

    # Precompute unit lookup to avoid repeated index() calls
    unit_map = {v: u for v, u in zip(var_list, var_unit_list)}

    for region in regions:
        # Skip early if region not requested
        if region not in metric_dict['region']:
            continue

        for mtype in metric_dict['type']:

            if mtype == "portrait":
                print(f"Processing Portrait  Plots for {metric} {region} {stat}....")

                var_names = sorted(var_list)
                var_units = [unit_map[v] for v in var_names]

                data_nor = {}
                # Build normalized arrays per season
                for season in metric_dict['season']:
                    # guard: skip if season/region missing
                    if season not in df_dict or region not in df_dict[season]:
                        continue

                    data_df = df_dict[season][region].copy()

                    # 1) Compute a (N_vars, N_models) array consistently
                    if stat == "cor_xy":
                        # Already in the correct order; transpose to (vars, models)
                        arr = data_df[var_names].to_numpy().T
                    else:
                        # Exclude ensemble (and E3SM) from baseline when not "default"
                        pattern = rf"^(?:{re.escape(mout_name)}|E3SM)"
                        df_norm = normalize_to_reference(
                            data_df, var_names,
                            mode=norm_method,
                            target=pattern,
                            model_col="model",
                        )
                        # Convert to (vars, models) array for plotting
                        arr = df_norm[var_names].to_numpy().T
                    data_nor[season] = arr  # always (N_vars, N_models)

                    # 2) Optional: save normalized data for this season
                    if save_data:
                        outdir = os.path.join(out_path, metric, region)
                        os.makedirs(outdir, exist_ok=True)
                        # Avoid KeyError if model_run is absent
                        outdic = data_df.drop(columns=["model_run"], errors="ignore").copy()
                        # Files are (N_models, N_vars) -> transpose back
                        outdic[var_names] = data_nor[season].T
                        archive_data(
                            region, stat, season, outdic, mout_name,
                            var_names, var_units, outdir
                        )

                # Choose a season that exists to get run_list
                run_list = None
                for season in metric_dict['season']:
                    if season in df_dict and region in df_dict[season]:
                        run_list = df_dict[season][region]['model'].to_list()
                        break
                if run_list is None:
                    continue

                stat_name = metric_dict['name']
                outdir = os.path.join(out_path, metric)
                os.makedirs(outdir, exist_ok=True)

                portrait_metric_plot(
                    region, stat, metric, data_nor,
                    stat_name, model_name, var_names,
                    run_list, outdir, fig_format
                )

            elif mtype == "parcoord":
                print(f"Processing Parallel Coordinate Plots for {metric} {region} {stat}....")

                for season in metric_dict['season']:
                    if season not in df_dict or region not in df_dict[season]:
                        continue

                    # Drop vars that are all-NaN; returns (data_df, kept_var_names, kept_var_units)
                    data_df, kept_var_names, kept_var_units = drop_vars(
                        df_dict[season][region].copy(),
                        var_list.copy(), var_unit_list.copy()
                    )

                    if data_df.empty or len(kept_var_names) == 0:
                        continue

                    if save_data:
                        outdir = os.path.join(out_path, metric, region)
                        os.makedirs(outdir, exist_ok=True)
                        outdic = data_df.drop(columns=["model_run"], errors="ignore").copy()
                        archive_data(
                            region, stat, season, outdic, mout_name,
                            kept_var_names, kept_var_units, outdir
                        )

                    run_list = data_df['model'].to_list()
                    stat_name = metric_dict['name']
                    outdir = os.path.join(out_path, metric)
                    os.makedirs(outdir, exist_ok=True)

                    parcoord_metric_plot(
                        region, stat, metric, data_df,
                        stat_name, model_name, kept_var_names, kept_var_units,
                        run_list, outdir, fig_format
                    )

    return

def append_mean_row(df: pd.DataFrame, source: pd.DataFrame, label: str) -> pd.DataFrame:
    if source.empty:
        return df
    num_cols = source.select_dtypes(include='number').columns
    mean_vals = source[num_cols].mean(numeric_only=True, skipna=True)
    row = {c: None for c in df.columns}
    row.update(mean_vals.to_dict())
    row['model'] = label
    return pd.concat([df, pd.DataFrame([row])], ignore_index=True)

def parcoord_metric_plot(
        region, stat, group, data_dict,
        stat_name, model_name,
        var_names, var_units,
        model_list, out_path,
        fig_format,
        vertical_center = "median",
        vertical_center_line = True,
        show_boxplot=False,
        show_violin=True,
        violin_colors=("lightgrey", "pink"),
        fontsize=18, 
        figsize = (40, 18),
        shrink = 0.8, 
        legend_box_xy = (1.08, 1.18),
        legend_box_size = 4,
        legend_lw = 1.5,
        identify_all_models = False,
        test_model_MMM = False,
        group1_name="CMIP",
        group2_name="E3SM",
        mean1_name=None,
        mean2_name=None,
        xcolors = None,
        color_map = "tab20_r",
        logo_rect=[0, 0, 0, 0],
        logo_off=True
    ):
    """ Function for parallel coordinate plots """
    legend_fontsize = fontsize * 0.8
    legend_ncol = int(7 * figsize[0] / 40.0)
    legend_posistion = (0.50, -0.14)
    
    if xcolors is None:
        xcolors = [
            "#e41a1c","#ff7f00","#4daf4a","#f781bf",
            "#a65628","#984ea3","#377eb8","#dede00"
        ]

    # --- Determine highlight models ---
    df_models = data_dict.get('model', [])
    highlight_model1 = get_highlight_models(df_models, model_name) or []

    # --- Work only on df (avoid mixing with data_dict) ---
    df = data_dict.reset_index(drop=True).copy()

    # Split CMIP vs E3SM
    highlight_set = set(highlight_model1)
    first_pos = next((i for i, m in enumerate(df['model'].tolist()) if m in highlight_set), None)
    if first_pos is None:
        cmip_slice = df
        test_slice = df.iloc[0:0]
    else:
        cmip_slice = df.iloc[:first_pos]
        test_slice = df.iloc[first_pos:]

    # --- Construct highlight lines & append means to df (not data_dict) ---
    highlight_model2 = []
    highlight_model2_label = []

    # Group mean names (defaults)
    mean1_name = mean1_name or f"{group1_name} (Mean)"
    mean2_name = mean2_name or f"{group2_name} (Mean)"

    if len(cmip_slice) > 0:
        df = append_mean_row(df, cmip_slice, mean1_name)
        highlight_model2.append(mean1_name)
        highlight_model2_label.append(mean1_name)

    subgroup_added = False
    # Either overall E3SM mean or version-specific subgroups (not both)
    if test_model_MMM and len(test_slice) > 0:
        df = append_mean_row(df, test_slice, mean2_name)
        highlight_model2.append(mean2_name)
        highlight_model2_label.append(mean2_name)
    elif len(test_slice) > 0:
        subgroup_patterns = OrderedDict({
            r'\bE3SM[-_ ]?1-1\b|\bv1(?:[._-]1\b)' : 'E3SMv1.1',
            r'\bE3SM[-_ ]?1-0\b|\bv1(?:[._-])'   : 'E3SMv1',
            r'\bE3SM[-_ ]?2-1\b|\bv2(?:[._-]1\b)' : 'E3SMv2.1',
            r'\bE3SM[-_ ]?2-0\b|\bv2(?:[._-])'   : 'E3SMv2',
            # includes v3- / v3. / v3_ and E3SM-3(-0)
            r'\bE3SM[-_ ]?3(?:-0)?\b|\bv3(?:[._-])' : 'E3SMv3',
            r'\bE3SM[-_ ]?4(?:-0)?\b|\bv4(?:[._-])' : 'E3SMv4',
        })

        remaining = test_slice.copy()
        test_subgroup_added = False

        for pattern, label in subgroup_patterns.items():
            mask = remaining['model'].astype(str).str.contains(
                pattern, regex=True, case=False, na=False
            )
            if mask.any():
                df = append_mean_row(df, remaining.loc[mask], label)
                highlight_model2.append(label)
                highlight_model2_label.append(label)
                test_subgroup_added = True
                # prevent rows from contributing to multiple subgroups
                remaining = remaining.loc[~mask]

        # Fallback if no version matched
        if not test_subgroup_added:
            df = append_mean_row(df, test_slice, mean2_name)
            highlight_model2.append(mean2_name)
            highlight_model2_label.append(mean2_name)

    # Colors for highlight lines (use actual mean labels)
    lncolors = []
    for i, label in enumerate(highlight_model2):
        if label == mean1_name:
            lncolors.append("#000000")     # CMIP mean
        elif label == mean2_name:
            lncolors.append("#5b5b5b")     # E3SM mean
        else:
            lncolors.append(xcolors[i % len(xcolors)])

    # Preserve user-provided variable order; filter to existing columns
    var_name1 = [v for v in var_names if v in df.columns]

    # Align labels with units
    var_labels = []
    for v in var_name1:
        idx = var_names.index(v)
        var_labels.append(f"{var_names[idx]}\n{var_units[idx]}" if var_units is not None else var_names[idx])

    # Build data array from df (after appends)
    data_var = df[var_name1].to_numpy()

    # Recompute model list from df to include appended means
    model_list_final = df['model'].astype(str).tolist()

    xlabel = "Metric"
    ylabel = f"{stat_name} ({stat.upper()})"
    if "mean_climate" in [group, region]:
        title = f"Model Performance of Annual Climatology ({stat.upper()}, {region.upper()})"
    elif "variability_modes" in [group, region]:
        title = f"Model Performance of Modes Variability ({stat.upper()})"
    elif "enso" in [group, region]:
        title = f"Model Performance of ENSO ({stat.upper()})"
    else:
        title = f"Model Performance ({stat.upper()})"

    fig, ax = parallel_coordinate_plot(
        data_var,
        var_labels,
        model_list_final,               # <-- updated
        model_names2=highlight_model1,
        group1_name=group1_name,
        group2_name=group2_name,
        models_to_highlight=highlight_model2,
        models_to_highlight_colors=lncolors,
        models_to_highlight_labels=highlight_model2_label,
        identify_all_models=identify_all_models,
        vertical_center=vertical_center,
        vertical_center_line=vertical_center_line,
        title=title,
        figsize=figsize,
        colormap=color_map,
        show_boxplot=show_boxplot,
        show_violin=show_violin,
        violin_colors=violin_colors,
        legend_ncol=legend_ncol,
        legend_bbox_to_anchor=legend_posistion,
        legend_fontsize=fontsize * 0.85,
        xtick_labelsize=fontsize * 0.95,
        ytick_labelsize=fontsize * 0.95,
        logo_rect=logo_rect,
        logo_off=logo_off
    )

    ax.set_xlabel(xlabel, fontsize=fontsize * 1.1)
    ax.set_ylabel(ylabel, fontsize=fontsize * 1.1)
    ax.set_title(title, fontsize=fontsize * 1.1)

    outdir = os.path.join(out_path, region)
    os.makedirs(outdir, exist_ok=True)
    
    outfile = f"{stat}_{region}_parcoord_{group}.{fig_format}"
    
    fig.savefig(os.path.join(outdir, outfile), facecolor="w", bbox_inches="tight")
    
    plt.close(fig)
    
    return

def portrait_metric_plot(
        region, stat, group, data_dict,
        stat_name, model_name,
        var_list, model_list,
        out_path, fig_format
    ):
    # process figure
    fontsize = 20
    add_vertical_line = True
    figsize = (40, 18)
    legend_box_xy = (1.08, 1.18)
    legend_box_size = 4
    legend_lw = 1.5
    shrink = 0.8 * len(var_list) / 30.0
    legend_fontsize = fontsize * 0.8

    # Prepare data array and legend controls
    if group == "mean_climate":
        # Expect seasonal dict: djf/mam/jja/son
        required_keys = ("djf", "mam", "jja", "son")
        if not all(k in data_dict for k in required_keys):
            raise KeyError(f"mean_climate group expects keys {required_keys}, got {tuple(data_dict.keys())}")
        data_all_nor = np.stack(
            [data_dict["djf"], data_dict["mam"], data_dict["jja"], data_dict["son"]]
        )
        legend_on = True
        legend_labels = ["DJF", "MAM", "JJA", "SON"]
    else:
        # Assume already a 2D array-like or 3D with extra dimension(s) handled by portrait_plot
        data_all_nor = data_dict
        legend_on = False
        legend_labels = []

    highlight_models = get_highlight_models(model_list, model_name)

    # Label colors: highlight target models, then any "e3sm", then default black
    label_colors = []
    for model in model_list:
        if model in model_name:
            label_colors.append("#FC5A50")
        elif "e3sm" in model.lower():
            label_colors.append("#5170d7")
        else:
            label_colors.append("#000000")

    # Colormap/range per stat
    if stat in ["cor_xy"]:
        var_range = (0, 1.0)
        cmap_color = "viridis"
        cmap_bounds = np.linspace(0, 1, 21)
    elif stat in ["stdv_pc_ratio_to_obs"]:
        var_range = (0.5, 1.5)
        cmap_color = "jet"
        # fine-grained 0.5 → 1.5 step 0.1
        cmap_bounds = [r / 10 for r in range(5, 16, 1)]
    else:
        var_range = (-0.5, 0.5)
        cmap_color = "RdYlBu_r"
        cmap_bounds = np.linspace(-0.5, 0.5, 11)

    fig, ax, cbar = portrait_plot(
        data_all_nor,
        xaxis_labels=model_list,
        yaxis_labels=var_list,
        cbar_label=stat_name,
        cbar_label_fontsize=fontsize * 1.0,
        cbar_tick_fontsize=fontsize,
        box_as_square=True,
        vrange=var_range,
        figsize=figsize,
        cmap=cmap_color,
        cmap_bounds=cmap_bounds,
        cbar_kw={"extend": "both", "shrink": shrink},
        missing_color="white",
        legend_on=legend_on,
        legend_labels=legend_labels,
        legend_box_xy=legend_box_xy,
        legend_box_size=legend_box_size,
        legend_lw=legend_lw,
        legend_fontsize=legend_fontsize,
        logo_rect=[0, 0, 0, 0],
        logo_off=True
    )

    # Optional vertical divider between non-highlight and highlight models
    if add_vertical_line and highlight_models:
        # Assumes highlight models are appended at the end of xaxis_labels
        xpos = len(model_list) - len(highlight_models)
        if 0 < xpos < len(model_list):
            ax.axvline(x=xpos, color="k", linewidth=3)

    # Ticks & colors
    ax.set_xticklabels(model_list, rotation=45, va="bottom", ha="left")
    ax.set_yticklabels(var_list, rotation=0, va="center", ha="right")
    for xtick, color in zip(ax.get_xticklabels(), label_colors):
        xtick.set_color(color)

    # Save figure
    outdir = os.path.join(out_path, region)
    os.makedirs(outdir, exist_ok=True)
    outfile = f"{stat}_{region}_portrait_{group}.{fig_format}"
    fig.savefig(os.path.join(outdir, outfile), facecolor="w", bbox_inches="tight")
    plt.close(fig)
    return


