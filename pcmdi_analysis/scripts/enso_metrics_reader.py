import glob
import json
import os
import re

from logger import _setup_child_logger
from utils import find_latest_file_list

logger = _setup_child_logger(__name__)


class EnsoMetricsReader:
    def __init__(self, parameter, metric_dict=None, mips=None, collections=None):
        # keep extra args for backwards-compat; unused by design
        self.ref_path = parameter["ref_path"]
        self.ref_name = parameter["ref_name"]
        self.test_path = parameter["test_path"]
        self.diag_vars = parameter["diag_vars"]

        self.var_pattern = re.compile(r"\.(\w+)\..*\.v(\d{8})\.json$")
        self.time_pattern = re.compile(r"\.v(\d{8})\.json$")
        self.test_combined = parameter.get("test_combined", False)
        
        if self.test_combined: 
            self.test_name = parameter["test_name"]
        else:
            self.mips = parameter["test_mip"]
            self.tests = parameter["test_name"]
            self.caseids = parameter["test_id"]
        
    def _get_combined_json_path(self, data_path, mip, mod, cid, metrics_collection):
        path = os.path.join(
            data_path,
            mod,
            mip,
            cid,
            metrics_collection,
            f"{mod.lower()}_{mip}_{metrics_collection}_*{cid}*.json",
        )
        matches = glob.glob(path)
        if not matches:
            raise FileNotFoundError(
                f"Combined metrics file not found for {mod} {mip} {cid} [{metrics_collection}]"
            )
        return matches[0]

    def _get_test_json_path(self, data_path, mip, mod, cid, metrics_collection):
        model_path = data_path.replace("put_model_here", mod)
        dir_path = f"{model_path}/{metrics_collection}"
        model_files = find_latest_file_list(
            path=dir_path,
            file_pattern=f"*{mip}*{cid}.json",
            var_pattern=self.var_pattern,
            time_pattern=self.time_pattern,
        )
        logger.info(dir_path)
        if not model_files or not os.path.exists(model_files[0]):
            raise FileNotFoundError(
                f"No Synthetic ENSO Metrics Data For {mip} {mod}, Aborting."
            )

        # Normalize model key inside JSON(s) to the actual model name
        for json_path in model_files:
            try:
                with open(json_path) as ff:
                    data_json = json.load(ff)
            except json.JSONDecodeError:
                logger.warning(f"[enso] Skipping unreadable JSON: {json_path}")
                continue

            results = data_json.get("RESULTS", {})
            models = results.get("model", {})
            if not isinstance(models, dict) or not models:
                logger.warning(f"[enso] Missing RESULTS.model in {json_path}; skipping rename.")
                continue

            old_key = next(iter(models.keys()))
            if old_key != mod:
                models[mod] = models.pop(old_key)
                # write back only if we actually changed the key
                with open(json_path, "w", encoding="utf8") as ff:
                    json.dump(
                        data_json, ff, indent=4, separators=(",", ": "), sort_keys=True
                    )

        # return the newest file selected by find_latest_file_list
        return model_files[0]

    def run(self, stat):
        # Validate diag_vars
        if not isinstance(self.diag_vars, dict):
            logger.error(
                f"[enso] parameter['diag_vars'] must be a dict; got {type(self.diag_vars).__name__}."
            )
            return {}

        metric_dict = self.diag_vars.get(stat, {})
        if not metric_dict:
            logger.warning(f"[enso] No variables configured for stat='{stat}'. Skipping.")
            return {}

        # --- Collections (optional config) ---
        enso_collections = metric_dict.get("collection", [])
        if not isinstance(enso_collections, (list, tuple)):
            logger.warning(
                f"[enso] 'collection' should be list/tuple; got {type(enso_collections).__name__}. Using empty list."
            )
            enso_collections = []

        logger.debug(f"[enso] stat='{stat}', collections={list(enso_collections)}")

        # --- Collect paths to ENSO metrics JSON files and return the mapping. ---
        dict_json_path = {}
        
        # Reference model (optional)
        ref_mips = []
        if isinstance(self.ref_name, str) and self.ref_name.strip():
            parts = self.ref_name.split(".")
            if len(parts) >= 3:
                ref_mod, ref_mip = parts[0], parts[1]
                ref_cid = ".".join(parts[2:])
                ref_mips = [(ref_mip, self.ref_name, ref_cid)]
                
            dict_json_path.setdefault(self.ref_name, {})
            for metrics_collection in enso_collections:
                dict_json_path[self.ref_name][metrics_collection] = self._get_combined_json_path(
                    self.ref_path, ref_mip, ref_mod, ref_cid, metrics_collection
                )
            if enso_collections and len(dict_json_path[self.ref_name]) < 1:
                raise FileNotFoundError(
                    f"No Synthetic ENSO Metrics Data for reference {ref_mod} {ref_mip} {ref_cid}."
                )
                 
        # TEST model 
        test_mips = []
        if self.test_combined and isinstance(self.test_name, str) and self.test_name.strip(): 
            parts = self.test_name.split(".")
            if len(parts) >= 3:
                test_mod, test_mip = parts[0], parts[1]
                test_cid = ".".join(parts[2:])
                test_mips = [(test_mip, self.test_name, test_cid)]
                
            dict_json_path.setdefault(self.test_name, {})
            for metrics_collection in enso_collections:
                dict_json_path[self.test_name][metrics_collection] = self._get_combined_json_path(
                    self.test_path, test_mip, test_mod, test_cid, metrics_collection
                )
            if enso_collections and len(dict_json_path[self.test_name]) < 1:
                raise FileNotFoundError(
                    f"No Synthetic ENSO Metrics Data for test {test_mod} {test_mip} {test_cid}."
                )
                    
        else:
            # Pair up test configurations (ensure aligned lengths to avoid silent truncation)
            if not (len(self.mips) == len(self.tests) == len(self.caseids)):
                raise ValueError(
                    "[enso] test_mip, test_name, and test_id must be the same length "
                    f"(got {len(self.mips)}, {len(self.tests)}, {len(self.caseids)})"
                )
            test_mips = list(zip(self.mips, self.tests, self.caseids))
            
            for mip, mod, cid in test_mips:
                dict_json_path.setdefault(mod, {})
                for metrics_collection in enso_collections:
                    dict_json_path[mod][metrics_collection] = self._get_test_json_path(
                        self.test_path, mip, mod, cid, metrics_collection
                    )
                if enso_collections and len(dict_json_path[mod]) < 1:
                    raise FileNotFoundError(
                        f"No Synthetic ENSO Metrics Data for test {mod} {mip} {cid}."
                    )
                            
        if ref_mips: 
            enso_mips = ref_mips + test_mips
        else:
            enso_mips = list(test_mips)

        return enso_mips, metric_dict, dict_json_path
