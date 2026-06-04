from __future__ import annotations
import os 
import re 
import json
from pathlib import Path
from typing import Dict, Any, Iterable, List, Optional, Tuple

class CLIMJSONMerger:
    """
    Merge mean-climate per-model JSONs into consolidated files.

    Input modes:
      A) member_dirs -> each dir contains files directly:
         <...>/pcmdi_diags/model_vs_obs/metrics_data/mean_climate/*.json
      B) fallback single tree:
         pmprdir/metrics_results/mean_climate/{mip}/{exp}/{case_id}/*.json

    Output (per var):
      {out_path}/metrics_results/mean_climate/{mip}/{exp}/{case_id}/{var}.{mip}.{exp}.{case_id}.json

    Filenames like:
      {var}.{grid}.{mip}.{exp}.{member_or_model}.{case_id}.json
      e.g. zg-500.2.5x2.5.e3sm.historical.v3-LR_0321.v20251015.json
    """

    def __init__(
        self,
        mips: Optional[List[str]] = None,
        exps: Optional[List[str]] = None,
        case_id: str = "v20230202",
        pmprdir: str = "/p/user_pub/pmp/pmp_results/pmp_v1.1.2",
        *,
        model_rename: Optional[Tuple[str, str, bool]] = None,
        strict: bool = False,
        verbose: bool = True,
        dry_run: bool = False,
        out_path: str = "./",
        member_dirs: Optional[List[str]] = None,
    ):
        self.mips = mips or ["cmip6"]
        self.exps = exps or ["historical"]
        self.case_id = case_id
        self.pmprdir = Path(pmprdir)
        self.strict = strict
        self.verbose = verbose
        self.dry_run = dry_run
        self.out_path = Path(out_path)
        self.member_dirs = [Path(p) for p in (member_dirs or [])]
        self.model_rename = model_rename

        # Fallback input root
        self.in_root_pattern = (
            self.pmprdir / "metrics_results" / "mean_climate" / "{mip}" / "{exp}" / self.case_id
        )

        # Output root (kept simple as requested)
        self.out_root_pattern = self.out_path / "{mip}" / "{exp}" / self.case_id

    # ==========================================================
    # Public API
    # ==========================================================
    def merge_all(self) -> None:
        """Discover variables and merge per MIP/EXP/VAR."""
        for mip in self.mips:
            for exp in self.exps:
                try:
                    variables = self._discover_variables(mip, exp)
                    self._log(f"[CLIM] {mip}/{exp} variables: {sorted(variables)}")
                    out_root = Path(str(self.out_root_pattern).format(mip=mip, exp=exp))
                    out_root.mkdir(parents=True, exist_ok=True)
                    for var in sorted(variables):
                        final_name = f"{var}.{mip}.{exp}.{self.case_id}.json"
                        final_path = out_root / final_name
                        self._merge_one(mip, exp, var, final_path)
                except Exception as err:
                    self._log(f"[CLIM][ERROR] {mip}/{exp} -> {err}")
                    if self.strict:
                        raise

    # ==========================================================
    # Core logic
    # ==========================================================
    def _merge_one(self, mip: str, exp: str, var: str, outfile: Path) -> None:
        json_files: List[Path] = []
        var_glob = f"{var}.*.{mip}.{exp}.*.{self.case_id}.json"

        if self.member_dirs:
            for root in self.member_dirs:
                if root.is_dir():
                    json_files.extend(sorted(root.glob(var_glob)))
            in_desc = f"[members x{len(self.member_dirs)}] var={var}"
        else:
            in_dir = Path(str(self.in_root_pattern).format(mip=mip, exp=exp))
            json_files = sorted(in_dir.glob(var_glob))
            in_desc = str(in_dir)

        if self.verbose:
            self._log(f"[CLIM] Search in: {in_desc}")
            self._log(f"[CLIM] {len(json_files)} files matched for var={var}")

        jsons = self._filter_individual_jsons(json_files)
        if not jsons:
            self._log(f"[CLIM][skip] No individual JSONs to merge for var={var}")
            return

        merged: Dict[str, Any] = {}
        for j, path in enumerate(jsons):
            if self.verbose and (j < 3 or j == len(jsons) - 1):
                self._log(f"  [{j+1}/{len(jsons)}] {path}")
            with path.open("r") as f:
                d = json.load(f)
            if j == 0:
                merged = d.copy()
            else:
                self._dict_merge(merged, d)

        if self.verbose:
            self._log(f"[CLIM] → Writing: {outfile}")
        outfile.parent.mkdir(parents=True, exist_ok=True)
        if not self.dry_run:
            with outfile.open("w") as fp:
                json.dump(merged, fp, sort_keys=True, indent=4)
                
        # Optional: rename model key inside RESULTS after merge
        if self.model_rename:
            old_key, new_key, with_ensemble = self.model_rename
            self.rename_results_model(
                outfile,
                old_key,
                new_key,
                suffix_with_ensemble = with_ensemble,
            )
            
    # ==========================================================
    # Rename RESULTS["model"] key (post-merge)
    # ==========================================================
    def rename_results_model(
        self,
        outfile: Path,  # was untyped
        old_key: str,
        new_key: str,
        *,
        replace_inside_strings: bool = False,
        overwrite_branch: bool = True,
        suffix_with_ensemble: bool = False,
    ) -> None:
        """
        Rename the top-level model key under RESULTS in a merged JSON file.

        Modes
        -----
        1) Simple rename (default):
             RESULTS[old_key] -> RESULTS[new_key]

        2) Fan-out by ensemble (suffix_with_ensemble=True):
             For each ensemble id 'ens' under RESULTS[old_key][obs][ens],
             create RESULTS[f"{new_key}-{ens}"][obs][ens] = original data.

        Options
        -------
        replace_inside_strings : also replace `old_key`->`new_key` within strings in the moved branch
        overwrite_branch       : allow overwriting an existing RESULTS[new_key] (or suffixed keys)
        """
        fpath = Path(outfile)

        if not fpath.exists():
            self._log(f"[CLIM][WARN] merged file not found: {fpath}")
            return

        with fpath.open("r") as fp:
            data = json.load(fp)

        results = data.get("RESULTS")
        if not isinstance(results, dict):
            self._log(f"[CLIM][WARN] no RESULTS in {fpath}")
            return

        if old_key not in results:
            self._log(f"[CLIM][WARN] RESULTS[{old_key}] not found in {fpath}")
            return

        # -------- Mode 1: simple rename --------
        if not suffix_with_ensemble:
            if (new_key in results) and not overwrite_branch:
                self._log(f"[CLIM][skip] RESULTS[{new_key}] exists (set overwrite_branch=True): {fpath}")
                return

            branch = results.pop(old_key)
            if replace_inside_strings:
                branch = self._replace_token_in_obj(branch, old=old_key, new=new_key)

            results[new_key] = branch
            data["RESULTS"] = results

        # -------- Mode 2: fan out by ensemble --------
        else:
            src = results.pop(old_key)
            if not isinstance(src, dict):
                self._log(f"[CLIM][WARN] RESULTS[{old_key}] not a dict in {fpath}")
                return
            # defensive checks and fan-out
            for obs, ens_map in (src.items() if isinstance(src, dict) else []):
                if not isinstance(ens_map, dict):
                    # unexpected, just skip this obs entry
                    continue
                for ens, payload in ens_map.items():
                    model_name = f"{new_key}-{ens}"
                    if (model_name in results) and not overwrite_branch:
                        # do not overwrite this model_name branch
                        self._log(f"[CLIM][skip] RESULTS[{model_name}] exists (set overwrite_branch=True): {fpath}")
                        continue

                    # ensure nesting RESULTS[model_name][obs][ens]
                    mbranch: Dict[str, Any] = results.setdefault(model_name, {})
                    obs_branch: Dict[str, Any] = mbranch.setdefault(obs, {})
                    # clone payload; optionally replace tokens inside
                    new_payload = payload
                    if replace_inside_strings:
                        new_payload = self._replace_token_in_obj(payload, old=old_key, new=model_name)

                    obs_branch[ens] = new_payload

            data["RESULTS"] = results

        # write back
        if not self.dry_run:
            with fpath.open("w") as fp:
                json.dump(data, fp, sort_keys=True, indent=4)

        action = (
            f"RESULTS[{old_key}] → RESULTS[{new_key}]"
            if not suffix_with_ensemble
            else f"RESULTS[{old_key}] → RESULTS[{new_key}-<ens>]"
        )
        self._log(f"[CLIM][write] {fpath} ({action})")
        
    def _replace_token_in_obj(self, obj: Any, old: str, new: str) -> Any:
        """
        Recursively replace `old` with `new` in all string values,
        list elements, and dictionary keys.
        """
        if isinstance(obj, dict):
            out: Dict[Any, Any] = {}
            for k, v in obj.items():
                nk = k.replace(old, new) if isinstance(k, str) else k
                out[nk] = self._replace_token_in_obj(v, old, new)
            return out

        if isinstance(obj, list):
            return [self._replace_token_in_obj(x, old, new) for x in obj]

        if isinstance(obj, str):
            return obj.replace(old, new)

        return obj

    # ==========================================================
    # Helpers
    # ==========================================================
    def _discover_variables(self, mip: str, exp: str) -> List[str]:
        """Collect unique leading tokens (var) from *.{mip}.{exp}.*.{case_id}.json."""
        vars_found = set()
        token_glob = f"*.{mip}.{exp}.*.{self.case_id}.json"

        if self.member_dirs:
            for root in self.member_dirs:
                if not root.is_dir():
                    self._log(f"[CLIM][WARN] missing dir: {root}")
                    continue
                for p in root.glob(token_glob):
                    var = p.name.split(".", 1)[0]
                    if var:
                        vars_found.add(var)
            return sorted(vars_found)

        in_dir = Path(str(self.in_root_pattern).format(mip=mip, exp=exp))
        if not in_dir.is_dir():
            self._log(f"[CLIM][WARN] input dir not found: {in_dir}")
            return []

        for p in in_dir.glob(token_glob):
            var = p.name.split(".", 1)[0]
            if var:
                vars_found.add(var)
        return sorted(vars_found)

    @staticmethod
    def _filter_individual_jsons(paths: Iterable[Path]) -> List[Path]:
        """Drop pre-merged/diveDown files (stem contains 'allModels'/'allRuns')."""
        kept: List[Path] = []
        for p in paths:
            parts = p.stem.split("_")
            if "allModels" in parts or "allRuns" in parts:
                continue
            kept.append(p)
        return kept

    @staticmethod
    def _dict_merge(dct: Dict[str, Any], merge_dct: Dict[str, Any]) -> None:
        """Recursive dict merge (in-place)."""
        for k, v in merge_dct.items():
            if isinstance(v, dict) and isinstance(dct.get(k), dict):
                CLIMJSONMerger._dict_merge(dct[k], v)
            else:
                dct[k] = v

    def _log(self, *args) -> None:
        if self.verbose:
            print(*args)


class MOVSJSONMerger:
    """
    Merge member-level variability-mode JSONs into consolidated files.

    Input (member mode):
      Each entry in `member_dirs` points at .../pcmdi_diags/model_vs_obs/metrics_data
      Files are located at:
        {member_dir}/variability_modes/{mode}/{obs}/
          var_mode_{mode}.{eof}.{mip}.{exp}.{member}.vs.{obs}.{case_id}.json

    Output:
      {out_path}/{mip}/{exp}/{case_id}/{mode}/{obs}/
        var_mode_{mode}.{eof}.{mip}.{exp}.allModels_allRuns.{syear}-{eyear}.json
    """

    def __init__(
        self,
        mips: Optional[List[str]] = None,
        exps: Optional[List[str]] = None,
        case_id: str = "v20230202",
        pmprdir: str = "/p/user_pub/pmp/pmp_results/pmp_v1.1.2",
        *,
        model_rename: Optional[Tuple[str, str, bool]] = None,
        movs_obses : Optional[List[str]] = None,
        movs_modes : Optional[List[str]] = None,
        syear: int = 1900,
        eyear: int = 2005,
        strict: bool = False,
        verbose: bool = True,
        dry_run: bool = False,
        out_path: str = "./",
        member_dirs: Optional[List[str]] = None,
    ):
        self.mips = mips or ["cmip6"]
        self.exps = exps or ["historical"]
        self.case_id = case_id
        self.pmprdir = pmprdir

        self.model_rename = model_rename
        self.member_dirs = [Path(p) for p in (member_dirs or [])]

        self.movs_obses = movs_obses or []   # keep your explicit selection
        self.movs_modes = movs_modes or []   # keep your explicit selection
        self.period = f"{syear}-{eyear}"

        self.strict = strict
        self.verbose = verbose
        self.dry_run = dry_run
        self.out_path = Path(out_path)

        # var_mode_{mode}.{eof}.{mip}.{exp}.{member}.vs.{obs}.{case_id}.json
        self.collect_glob = "var_mode_{mode}.*.{mip}.{exp}.*.vs.{obs}.{case_id}.json"

    # =========================
    # Public API
    # =========================
    def merge_all(self) -> None:
        """Discover (mode, obs, mip, exp, eof) combos from member_dirs and merge each."""
        combos = self._discover_combinations(
            self.movs_modes,
            self.movs_obses,
            self.period,
            self.case_id,
            self.member_dirs,
        )

        if self.verbose:
            self._log(f"[MOVS] discovered {len(combos)} (mode,obs,mip,exp,eof) combos")
        for (mode, obs, mip, exp, eof) in sorted(combos):
            try:
                final_dir = (self.out_path / mip / exp / self.case_id / mode / obs)  # FIX: Path
                final_dir.mkdir(parents=True, exist_ok=True)
                final_name = f"var_mode_{mode}.{eof}.{mip}.{exp}.allModels_allRuns.{self.period}.json"
                final_path = final_dir / final_name
                self._merge_one(mode, obs, mip, exp, eof, final_path)
            except Exception as err:
                self._log(f"[ERROR] {mip}/{exp}/{mode}/{obs} -> {err}")
                if self.strict:
                    raise

    # =========================
    # Core logic
    # =========================
    # ==========================================================
    # Rename RESULTS["model"] key (post-merge)
    # ==========================================================
    def rename_results_model(
        self,
        outfile: Path,  # was untyped
        old_key: str,
        new_key: str,
        *,
        replace_inside_strings: bool = False,
        overwrite_branch: bool = True,
        suffix_with_ensemble: bool = False,
    ) -> None:
        """
        Rename the top-level model key under RESULTS in a merged JSON file.

        Modes
        -----
        1) Simple rename (default):
             RESULTS[old_key] -> RESULTS[new_key]

        2) Fan-out by ensemble (suffix_with_ensemble=True):
             For each ensemble id 'ens' under RESULTS[old_key][obs][ens],
             create RESULTS[f"{new_key}-{ens}"][obs][ens] = original data.

        Options
        -------
        replace_inside_strings : also replace `old_key`->`new_key` within strings in the moved branch
        overwrite_branch       : allow overwriting an existing RESULTS[new_key] (or suffixed keys)
        """
        fpath = Path(outfile)

        if not fpath.exists():
            self._log(f"[CLIM][WARN] merged file not found: {fpath}")
            return

        with fpath.open("r") as fp:
            data = json.load(fp)

        results = data.get("RESULTS")
        if not isinstance(results, dict):
            self._log(f"[CLIM][WARN] no RESULTS in {fpath}")
            return

        if old_key not in results:
            self._log(f"[CLIM][WARN] RESULTS[{old_key}] not found in {fpath}")
            return

        # -------- Mode 1: simple rename --------
        if not suffix_with_ensemble:
            if (new_key in results) and not overwrite_branch:
                self._log(f"[CLIM][skip] RESULTS[{new_key}] exists (set overwrite_branch=True): {fpath}")
                return

            branch = results.pop(old_key)
            if replace_inside_strings:
                branch = self._replace_token_in_obj(branch, old=old_key, new=new_key)

            results[new_key] = branch
            data["RESULTS"] = results

        # -------- Mode 2: fan out by ensemble --------
        else:
            src = results.pop(old_key)
            if not isinstance(src, dict):
                self._log(f"[MOVS][WARN] RESULTS[{old_key}] not a dict in {fpath}")
                return
            # defensive checks and fan-out
            for ens, payload in (src.items() if isinstance(src, dict) else []):
                model_name = f"{new_key}-{ens}"
                if (model_name in results) and not overwrite_branch:
                    # do not overwrite this model_name branch
                    self._log(f"[MOVS][skip] RESULTS[{model_name}] exists (set overwrite_branch=True): {fpath}")
                    continue

                # ensure nesting RESULTS[model_name][obs][ens]
                mbranch: Dict[str, Any] = results.setdefault(model_name, {})
                # clone payload; optionally replace tokens inside
                new_payload = payload
                if replace_inside_strings:
                    new_payload = self._replace_token_in_obj(payload, old=old_key, new=model_name)

                mbranch[ens] = new_payload

            data["RESULTS"] = results

        # write back
        if not self.dry_run:
            with fpath.open("w") as fp:
                json.dump(data, fp, sort_keys=True, indent=4)

        action = (
            f"RESULTS[{old_key}] → RESULTS[{new_key}]"
            if not suffix_with_ensemble
            else f"RESULTS[{old_key}] → RESULTS[{new_key}-<ens>]"
        )
        self._log(f"[MOVS][write] {fpath} ({action})")
        
    def _replace_token_in_obj(self, obj: Any, old: str, new: str) -> Any:
        """
        Recursively replace `old` with `new` in all string values,
        list elements, and dictionary keys.
        """
        if isinstance(obj, dict):
            out: Dict[Any, Any] = {}
            for k, v in obj.items():
                nk = k.replace(old, new) if isinstance(k, str) else k
                out[nk] = self._replace_token_in_obj(v, old, new)
            return out

        if isinstance(obj, list):
            return [self._replace_token_in_obj(x, old, new) for x in obj]

        if isinstance(obj, str):
            return obj.replace(old, new)

        return obj
    
    def _merge_one(self, mode: str, obs: str, mip: str, exp: str, eof: str, outfile: Path) -> None:
        # Collect candidate files from all member roots
        json_files: List[Path] = []
        pattern = self.collect_glob.format(mode=mode, mip=mip, exp=exp, obs=obs, case_id=self.case_id)

        for root in self.member_dirs:
            search_dir = root / mode / obs  # FIX: include variability_modes
            json_files.extend(sorted(search_dir.glob(pattern)))

        if not json_files:
            self._log(f"[skip] No files for {mip}/{exp}/{mode}/{obs} (pattern: {pattern})")
            return

        files = self._filter_individual_jsons(json_files)
        if not files:
            self._log(f"[skip] Only merged/diveDown files found for {mip}/{exp}/{mode}/{obs}")
            return

        if self.verbose:
            self._log(f"[MOVS] Merging {len(files)} files: {mip}/{exp}/{mode}/{obs} ({eof})")
            for j, p in enumerate(files[:3], start=1):
                self._log(f"  [{j}/{len(files)}] {p}")

        merged: Dict[str, Any] = {}
        for idx, path in enumerate(files, start=1):
            with path.open("r") as f:
                d = json.load(f)
            if idx == 1:
                merged = d.copy()
            else:
                self._dict_merge(merged, d)

        if self.verbose:
            self._log(f"→ Writing: {outfile}")
        if not self.dry_run:
            with outfile.open("w") as fp:   # FIX: outfile is Path
                json.dump(merged, fp, sort_keys=True, indent=4)
                
        # Optional: rename model key inside RESULTS after merge
        if self.model_rename:
            old_key, new_key, with_ensemble = self.model_rename
            self.rename_results_model(
                outfile,
                old_key,
                new_key,
                suffix_with_ensemble = with_ensemble,
            )
            
    # =========================
    # Discovery & helpers
    # =========================
    def _discover_combinations(
        self,
        movs_modes: List[str],
        movs_obses: List[str],
        period: str,
        case_id: str,
        member_dirs: List[Path],
    ) -> List[tuple]:
        """
        Scan member_dirs to discover unique (mode, obs, mip, exp, eof) tuples.
        Filename:
          var_mode_{mode}.{eof}.{mip}.{exp}.{member}.vs.{obs}.{case_id}.json
        """
        combos = set()
        for root in member_dirs:
            if not root.is_dir():
                self._log(f"[WARN] missing dir: {base}")
                continue
                
            for mode, obs in zip(movs_modes,movs_obses):
                mode_dir = root / mode / obs
                if not mode_dir.is_dir():
                    self._log(f"[WARN] missing dir: {mode_dir}")
                    continue

                for p in mode_dir.glob(f"var_mode_{mode}.*.*.*.*.vs.{obs}.{case_id}.json"):  # FIX: Path.glob
                    parts = p.name.split(".")
                    if len(parts) < 9:
                        continue
                    try:
                        eof = parts[1]
                        mip = parts[2]
                        exp = parts[3]
                        case = parts[7]
                        if case != self.case_id:
                            continue
                        if self.mips and (mip not in self.mips):
                            continue
                        if self.exps and (exp not in self.exps):
                            continue
                        combos.add((mode, obs, mip, exp, eof))
                    except Exception:
                        continue
        return list(combos)

    @staticmethod
    def _filter_individual_jsons(paths: Iterable[Path]) -> List[Path]:
        """Drop already-merged files (member token equals 'allModels_allRuns')."""
        kept: List[Path] = []
        for p in paths:
            parts = p.name.split(".")
            if len(parts) >= 5 and parts[4] == "allModels_allRuns":
                continue
            kept.append(p)
        return kept

    @staticmethod
    def _dict_merge(dct: Dict[str, Any], merge_dct: Dict[str, Any]) -> None:
        """Recursive dict merge (in-place)."""
        for k, v in merge_dct.items():
            if isinstance(v, dict) and isinstance(dct.get(k), dict):
                MOVSJSONMerger._dict_merge(dct[k], v)
            else:
                dct[k] = v

    def _log(self, *args) -> None:
        if self.verbose:
            print(*args)

class ENSOJSONMerger:
    """
    Merge per-model ENSO metric JSONs into a single 'allModels_allRuns' JSON.

    Input filename shape:
      {mc}.{mip}.{exp}.{member}.vs.{obs}.{case_id}.json
    Example:
      ENSO_perf.e3sm.historical.v3-LR_0321.vs.ERSST.v20251015.json
    """

    def __init__(
        self,
        mips: Optional[List[str]] = None,
        exps: Optional[List[str]] = None,
        case_id: str = "v20230202",
        pmprdir: str = "/p/user_pub/pmp/pmp_results/pmp_v1.1.2",
        *,
        model_rename: Optional[Tuple[str, str, bool]] = None,
        collections: Optional[List[str]] = None,
        observations: Optional[List[str]] = None,
        strict: bool = False,
        verbose: bool = True,
        dry_run: bool = False,
        skip_tokens: Optional[Iterable[str]] = None,
        out_path: str = "./",
        member_dirs: Optional[List[str]] = None,
    ):
        self.mips = mips or ["cmip6"]
        self.exps = exps or ["historical"]
        self.case_id = case_id

        self.pmprdir = Path(pmprdir)
        self.out_path = Path(out_path)
        self.member_dirs = [Path(p) for p in (member_dirs or [])]

        self.model_rename = model_rename
        self.collections = collections or []
        self.observations = observations or []

        self.strict = strict
        self.verbose = verbose
        self.dry_run = dry_run

        self.skip_tokens: Tuple[str, ...] = tuple(skip_tokens) if skip_tokens is not None else (
            "diveDown",
            "allModels",
            "allRuns",
        )

        # {mc}.{mip}.{exp}.{member}.vs.{obs}.{case_id}.json
        self.enso_glob = "{mc}.{mip}.{exp}.*.vs.{obs}.{case_id}.json"

    # -------------------- helpers --------------------
    @staticmethod
    def _filter_individual_jsons(paths: Iterable[Path]) -> List[Path]:
        """Drop already-merged files (member token equals 'allModels_allRuns')."""
        kept: List[Path] = []
        for p in paths:
            parts = p.name.split(".")
            if len(parts) >= 5 and parts[3] == "allModels_allRuns":  # mc mip exp [member]
                continue
            kept.append(p)
        return kept

    @staticmethod
    def _dict_merge(dct: Dict[str, Any], merge_dct: Dict[str, Any]) -> None:
        """Recursive dict merge (in-place)."""
        for k, v in merge_dct.items():
            if isinstance(v, dict) and isinstance(dct.get(k), dict):
                ENSOJSONMerger._dict_merge(dct[k], v)
            else:
                dct[k] = v

    def _log(self, *args) -> None:
        if self.verbose:
            print(*args)

    def rename_results_model(
        self,
        outfile: Path,
        new_key: str,
        *,
        replace_inside_strings: bool = True,
        overwrite_branch: bool = True,
        suffix_with_ensemble: bool = False,
    ) -> None:
        """
        Rename model branch in ENSO merged JSON (no OBS layer).

        Supports two layouts:
          A) data["RESULTS"][<model>] = { ens: payload, ... }
          B) data["RESULTS"]["model"][<model>] = { ens: payload, ... }

        If suffix_with_ensemble=True, fan out into RESULTS[f"{new_key}-{ens}"] (same layout).
        """
        fpath = Path(outfile)
        if not fpath.exists():
            self._log(f"[ENSO][WARN] merged file not found: {fpath}")
            return

        with fpath.open("r") as fp:
            data = json.load(fp)

        results_root = data.get("RESULTS")
        if not isinstance(results_root, dict):
            self._log(f"[ENSO][WARN] no RESULTS in {fpath}")
            return

        # Where are model branches stored?
        if isinstance(results_root.get("model"), dict):
            container_key = "model"          # write back to RESULTS["model"]
            container = results_root["model"]
        else:
            container_key = "RESULTS"        # write back to RESULTS
            container = results_root

        if not isinstance(container, dict):
            self._log(f"[ENSO][WARN] model container is not a dict in {fpath}")
            return

        model_keys = list(container.keys())
        if len(model_keys) == 0:
            self._log(f"[ENSO][WARN] no model keys under {container_key} in {fpath}")
            return
        if len(model_keys) > 1 and not suffix_with_ensemble:
            self._log(f"[ENSO][WARN] multiple model keys {model_keys}; "
                      f"ambiguous which to rename → '{new_key}'")
            return

        keys_to_process = model_keys if suffix_with_ensemble else [model_keys[0]]

        for old_key in list(keys_to_process):
            if old_key not in container:
                self._log(f"[ENSO][WARN] {container_key}[{old_key}] not found")
                continue

            src = container.pop(old_key)
            if not isinstance(src, dict):
                self._log(f"[ENSO][WARN] {container_key}[{old_key}] not a dict")
                continue

            if not suffix_with_ensemble:
                # Simple rename: RESULTS[old_key] -> RESULTS[new_key]
                if (new_key in container) and not overwrite_branch:
                    self._log(f"[ENSO][skip] {container_key}[{new_key}] exists "
                              f"(set overwrite_branch=True)")
                    container[old_key] = src  # restore
                    continue
                branch = (self._replace_token_in_obj(src, old=old_key, new=new_key)
                          if replace_inside_strings else src)
                container[new_key] = branch
            else:
                # Fan-out: src is {ens: payload}
                for ens, payload in src.items():
                    model_name = f"{new_key}-{ens}"
                    if (model_name in container) and not overwrite_branch:
                        self._log(f"[ENSO][skip] {container_key}[{model_name}] exists "
                                  f"(set overwrite_branch=True)")
                        continue
                    new_payload = (self._replace_token_in_obj(payload, old=old_key, new=model_name)
                                   if replace_inside_strings else payload)
                    # Keep same “no-obs” layout: {model_name: {ens: payload}}
                    mbranch = container.setdefault(model_name, {})
                    mbranch[ens] = new_payload

        # Write back in same layout
        if container_key == "model":
            results_root["model"] = container
            data["RESULTS"] = results_root
        else:
            data["RESULTS"] = container

        if not self.dry_run:
            with fpath.open("w") as fp:
                json.dump(data, fp, sort_keys=True, indent=4)

        action = "simple rename" if not suffix_with_ensemble else "fan-out by ensemble"
        self._log(f"[ENSO][write] {fpath} ({action} → '{new_key}')")

    def _replace_token_in_obj(self, obj: Any, old: str, new: str) -> Any:
        """Recursively replace tokens in strings, list elements, and dict keys."""
        if isinstance(obj, dict):
            out: Dict[Any, Any] = {}
            for k, v in obj.items():
                nk = k.replace(old, new) if isinstance(k, str) else k
                out[nk] = self._replace_token_in_obj(v, old, new)
            return out
        if isinstance(obj, list):
            return [self._replace_token_in_obj(x, old, new) for x in obj]
        if isinstance(obj, str):
            return obj.replace(old, new)
        return obj

    # -------------------- merging --------------------
    def merge_one(self, mip: str, exp: str, mc: str, obs: str, case_id: str, out_file: Path) -> Path:
        """
        Merge all per-model JSONs for one (mip, exp, mc, obs, case_id) into out_file.
        """
        pattern = self.enso_glob.format(mc=mc, mip=mip, exp=exp, obs=obs, case_id=case_id)

        # Try both {member_dir}/{mc} and {member_dir}/{mc}/{obs}
        json_files: List[Path] = []
        for root in self.member_dirs:
            for search_dir in (root / mc, root / mc / obs):
                if search_dir.is_dir():
                    json_files.extend(sorted(search_dir.glob(pattern)))

        if not json_files:
            self._log(f"[ENSO][skip] No files for {mip}/{exp}/{mc}/{obs} (pattern: {pattern})")
            return out_file

        files = [p for p in self._filter_individual_jsons(json_files)
                 if not any(tok in p.stem for tok in self.skip_tokens)]
        if not files:
            self._log(f"[ENSO][skip] Only merged/diveDown files found for {mip}/{exp}/{mc}/{obs}")
            return out_file

        if self.verbose:
            self._log(f"[ENSO] Merging {len(files)} files: {mip}/{exp}/{mc}/{obs}")
            for j, p in enumerate(files[:3], start=1):
                self._log(f"  [{j}/{len(files)}] {p}")

        merged: Dict[str, Any] = {}
        for idx, path in enumerate(files, start=1):
            with path.open("r") as f:
                d = json.load(f)
            if idx == 1:
                merged = d.copy()
            else:
                self._dict_merge(merged, d)

        out_file.parent.mkdir(parents=True, exist_ok=True)
        if not self.dry_run:
            with out_file.open("w") as fp:
                json.dump(merged, fp, sort_keys=True, indent=4)

        # Optional in-file model rename
        if self.model_rename:
            _, new_key, with_ens = self.model_rename
            self.rename_results_model(
                out_file, new_key, suffix_with_ensemble=with_ens
            )

        return out_file

    def merge_all(self) -> None:
        """Process all (mip, exp, mc, obs) pairs and write merged JSONs."""
        if len(self.collections) != len(self.observations):
            raise ValueError("collections and observations must be the same length (paired).")

        for mip in self.mips:
            for exp in self.exps:
                for mc, obs in zip(self.collections, self.observations):
                    self._log(f"[INFO] mip={mip} exp={exp} mc={mc} obs={obs} case_id={self.case_id}")
                    try:
                        final_dir = self.out_path / mip / exp / self.case_id / mc
                        out_name = f"{mip}_{exp}_{mc}_{obs}_{self.case_id}_allModels_allRuns.json"
                        out_file = final_dir / out_name
                        self.merge_one(mip, exp, mc, obs, self.case_id, out_file)
                    except Exception as err:
                        self._log(f"[ERROR] {mip}/{exp}/{mc}/{self.case_id} -> {err}")
                        if self.strict:
                            raise
