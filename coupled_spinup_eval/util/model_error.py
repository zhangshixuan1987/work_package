"""e3sm_diags metrics analysis: drift of PCC/RMSE/bias over 50-yr windows and
cross-experiment comparison (SpinupMetricAnalyzer).

Options:
  bias_ratio_metric : "Mean_Bias" -> (bias - ref)/|ref|  (9_model_error_compare_analysis)
                      "Bias"      -> bias / ref          (8_model_error_drift_analysis)
"""
import os
import re
import glob
import string
import json
import numpy as np
import xarray as xr
import cftime
import fnmatch
import pandas as pd
import seaborn as sns
from collections import OrderedDict
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, FixedLocator
from matplotlib.colors import TwoSlopeNorm, BoundaryNorm,ListedColormap
import matplotlib as mpl
from matplotlib.patches import Polygon, Patch, Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable


class SpinupMetricAnalyzer:
    """
    Pipeline:
      1) extract_metrics()  -> read CSVs, assemble pandas tables (no ingest normalization)
      2) build_dataset()    -> create xr.Dataset with raw and normalized fields
      3) to_netcdf(path)    -> save the diagnostic file
      4) load_dataset(path) -> load the saved file (later sessions)
      5) plot_heatmap(...) / plot_heatmap_quadseason(...) / plot_heatmap_experiments_period(...)

    Heatmap orientation: columns=X=variables, rows=Y=time periods

    IMPORTANT:
    - Windows whose START year > picontrol_switch_year are read from piControl (offset).
    - We only build windows up to spinup_years (no automatic +500 extension).

    Normalization:
    - Legacy: single reference period label (e.g., "2451-2500")
    - New: window mean, e.g., reference_window=(2001, 2500) uses mean across all segments
           whose START is within [2001, 2500].
    - External: call set_external_reference_vectors(...) to use a precomputed baseline (e.g., 500-yr piControl)
      for per-variable reference in ratio/relative metrics.
    """

    # ----------------------- init -----------------------
    def __init__(
        self,
        base_dir,
        spinup_name,
        picontrol_name=None,
        subdir=None,
        season="ANN",
        variables=None,
        metrics=None,
        normalize_across_rows=True,     # build ratio/relative fields
        reference_period="0001-0050",
        reference_window=None,          # (start, end) -> window-mean reference
        spinup_years=2000,
        picontrol_switch_year=None,     # defaults to spinup_years
        bias_ratio_metric="Mean_Bias",  # "Mean_Bias": (bias-ref)/|ref| ; "Bias": bias/ref (drifting_error)
    ):
        self.bias_ratio_metric = bias_ratio_metric
        self.base_dir = base_dir
        self.spinup_name = spinup_name
        self.picontrol_name = picontrol_name
        self.subdir_pattern = subdir or (
            "{exp}/e3sm_diags/atm_monthly_180x360_aave/model_vs_obs_{period}/"
            "viewer/table-data/{season}_metrics_table.csv"
        )
        self.season = season
        self.normalize_across_rows = normalize_across_rows
        self.reference_period = reference_period
        self.reference_window = reference_window
        self.spinup_years = int(spinup_years)
        self.picontrol_switch_year = (
            int(picontrol_switch_year) if picontrol_switch_year is not None else int(spinup_years)
        )

        # content
        self.variables = variables or []    # full row labels initially; later base vars
        self.units = {}                     # {var: unit}
        self.metrics = metrics or ["Mean_Bias", "RMSE", "Correlation"]
        self.metric_short = {"Mean_Bias": "Bias", "RMSE": "RMSE", "Correlation": "PCC"}

        # filled by extract_metrics()
        self.tables = {m: pd.DataFrame() for m in self.metrics}  # rows=variables, cols=period labels
        self.obs_std = {}  # per-base-variable observational std (filled after ingest)

        # filled by build_dataset()/load_dataset()
        self.ds = None

        # ingest bookkeeping
        self.obs_std_raw = {}          # {full_row_label: std}
        self._once = set()

        # optional plotting catalog
        self._plot_catalog = None           # normalized {var: {region: [sources...]}}
        self._plot_order_bases = None       # ordered list of base variables to preserve order

        # external reference map: {"Mean_Bias": {var: val}, "RMSE": {...}, "Correlation": {...}}
        self._external_ref = None

    # ----------------------- external reference -----------------------
    def set_external_reference_vectors(self, ref_map: dict):
        """
        Provide per-metric, per-variable reference values (e.g., 500-yr piControl means).
        Example:
            {
              "Mean_Bias": {"PRECT": 0.12, ...},
              "RMSE": {"PRECT": 1.8, ...},
              "Correlation": {"PRECT": 0.67, ...},
              "RMSE_sigma": {"PRECT": 0.03, ...},   # optional, see piControl_sigma_vector()
            }
        These are used by ratio/relative-normalized fields (Bias/RMSE *_ratio_ref, *_rel_ref);
        "RMSE_sigma" (piControl internal variability of RMSE) enables RMSE_z_ref.
        """
        self._external_ref = {}
        for k in ["Mean_Bias", "RMSE", "Correlation", "RMSE_sigma"]:
            self._external_ref[k] = {str(v): float(val) for v, val in (ref_map.get(k, {}) or {}).items()}

    def _external_ref_vector(self, metric_name: str, variables: list):
        """
        Build a per-variable reference vector using self._external_ref if available.
        Missing variables get NaN (leading to NaNs in normalized outputs for those columns).
        """
        if not self._external_ref or metric_name not in self._external_ref:
            return None
        mp = self._external_ref[metric_name]
        out = np.array([mp.get(str(v), np.nan) for v in variables], dtype=float)
        return out

    # ----------------------- helpers -----------------------
    def _parse_period(self, label: str):
        """'0001-0050' -> (1, 50) for numeric sorting."""
        a, b = label.split("-")
        return int(a), int(b)

    def _period_label(self, start: int, end: int):
        return f"{start:04d}-{end:04d}"

    def _get_csv_paths(self):
        """
        Build list of (display_label, csv_path) for all 50-year windows from 1..spinup_years.
        Windows whose START > picontrol_switch_year come from piControl (offset indexing).
        """
        paths = []
        last_end = self.spinup_years

        for start_year in range(1, last_end + 1, 50):
            end_year = start_year + 49
            if end_year > last_end:
                break

            label = self._period_label(start_year, end_year)

            if start_year <= self.picontrol_switch_year:
                exp_name = self.spinup_name
                source_period = label
            else:
                if self.picontrol_name is None:
                    break
                exp_name = self.picontrol_name
                start_raw = start_year - self.picontrol_switch_year
                end_raw   = end_year   - self.picontrol_switch_year
                source_period = f"{int(start_raw):04d}-{int(end_raw):04d}"

            csv_path = os.path.join(
                self.base_dir,
                self.subdir_pattern.format(exp=exp_name, period=source_period, season=self.season),
            )
            if os.path.isfile(csv_path):
                paths.append((label, csv_path))
            else:
                print(f"[WARNING] Missing file for {label}: {csv_path}")

        paths.sort(key=lambda t: self._parse_period(t[0])[0])

        if paths:
            first_spin_label = next((p for p, _ in paths if self._parse_period(p)[0] <= self.picontrol_switch_year), None)
            first_ctrl_label = next((p for p, _ in paths if self._parse_period(p)[0] >  self.picontrol_switch_year), None)
            print(f"[INFO] windows: first spin-up={first_spin_label}, first piControl={first_ctrl_label}, last={paths[-1][0]}")
        return paths

    def _sort_period_columns(self, df_wide: pd.DataFrame):
        cols = sorted([c for c in df_wide.columns], key=lambda c: self._parse_period(c)[0])
        return df_wide[cols]

    def _split_var_obs(self, var: str):
        """
        Return (base, obs) for labels like:
          - 'PRECT (GPCP_v2.3)'
          - 'PRECT global GPCP_v2.3'
          - 'U_850mb global ERA5'
          - 'CLDHGH_CAL global CALIPSOCOSP'
        Fallback: obs=None and base=first token.
        """
        s = str(var).strip()
        m = re.match(r"^(?P<base>.+?)\s*\((?P<obs>[^)]+)\)\s*$", s)
        if m: return m.group("base").strip(), m.group("obs").strip()
        m = re.match(r"^(?P<base>[^ ]+?)\s+(global|land|ocean)\s+(?P<obs>.+)$", s, flags=re.IGNORECASE)
        if m: return m.group("base").strip(), m.group("obs").strip()
        m = re.match(r"^(?P<base>.+?)\s+(global|land|ocean)\s+(?P<obs>.+)$", s, flags=re.IGNORECASE)
        if m: return m.group("base").strip(), m.group("obs").strip()
        parts = s.split()
        return (parts[0], None) if len(parts) >= 1 else (s, None)

    def _parse_triplet_from_row(self, row_label: str):
        """
        Parse a CSV 'Variables' row into (var, region, source).
        Accepts:
          - 'VAR global ERA5' / 'VAR land ERA5' / 'VAR ocean ERA5'
          - 'VAR (ERA5)'          -> region='global'
          - 'VAR'                 -> region='global', source=None
        """
        s = str(row_label).strip()
        m = re.match(r"^(?P<var>.+?)\s*\((?P<src>[^)]+)\)\s*$", s)
        if m:
            return m.group("var").strip(), "global", m.group("src").strip()
        m = re.match(r"^(?P<var>[^ ]+)\s+(?P<region>global|land|ocean)\s+(?P<src>.+)$", s, flags=re.IGNORECASE)
        if m:
            return m.group("var").strip(), m.group("region").lower().strip(), m.group("src").strip()
        parts = s.split()
        var = parts[0] if parts else s
        return var, "global", None

    # ----------------------- plotting catalog -----------------------
    def use_plot_catalog(self, path: str):
        """
        Provide a JSON file defining {var: {region: [sources...]}} (or wrapper with 'regions').
        Preserves the order of base variables in the file for plotting.
        """
        with open(path, "r", encoding="utf-8") as f:
            raw = json.load(f)

        catalog = {}
        order = []
        for var, val in raw.items():  # dict preserves insertion order in Py3.7+
            regmap = val["regions"] if (isinstance(val, dict) and "regions" in val and isinstance(val["regions"], dict)) else val
            norm = {}
            for region, srcs in regmap.items():
                if srcs is None:
                    norm[region] = [None]
                elif isinstance(srcs, list):
                    clean = [s for s in srcs if s] or [None]
                    norm[region] = list(dict.fromkeys(clean))  # dedupe, preserve order
                else:
                    norm[region] = [str(srcs)]
            catalog[str(var)] = norm
            order.append(str(var))

        self._plot_catalog = catalog
        self._plot_order_bases = order

    def _allowed_triplets(self):
        """
        Build a set of allowed (var, region, source_or_None) from self._plot_catalog.
        """
        if not self._plot_catalog:
            return None
        allowed = set()
        for var, regmap in self._plot_catalog.items():
            for region, srcs in regmap.items():
                region_norm = (region or "global").lower()
                for src in (srcs or [None]):
                    src_norm = src if (src is None or str(src).strip() == "") else str(src).strip()
                    allowed.add((str(var), region_norm, src_norm if src_norm else None))
        return allowed

    # ----------------------- Stage 1: gather -----------------------
    def extract_metrics(self, select_vars=None, exact_match=True):
        """
        Read all CSVs into self.tables[metric] (rows=variables, cols=period labels).
        If a plot catalog was provided via use_plot_catalog(), rows are filtered to that set,
        and the display order of base variables will follow the catalog.
        """
        selected = set(select_vars) if select_vars else None

        for period, path in self._get_csv_paths():
            df = pd.read_csv(
                path,
                na_values=[999.999, 9.99999e2, 999.99, 9.99e2, -999, -999.0, -999.99, "nan", "NaN"],
            )
            df.columns = [str(c).strip().replace("\ufeff","") for c in df.columns]
            if "Variables" not in df.columns:
                for alt in ("Variable", "Var", "Name"):
                    if alt in df.columns:
                        df = df.rename(columns={alt: "Variables"})
                        break
            if "Variables" not in df.columns:
                print(f"[WARNING] No 'Variables' column in {path}; skipping.")
                continue

            df = df.replace([999.999, 9.99999e2, -999, -999.0], np.nan)

            # Variable selection first
            all_vars = df["Variables"].astype(str)
            if selected:
                if exact_match:
                    mask = all_vars.isin(selected)
                else:
                    mask = all_vars.apply(lambda v: any(s.lower() in v.lower() for s in selected))
                df = df.loc[mask]

            if df.empty:
                print(f"[INFO] No variables kept for {period} (selection={len(selected) if selected else 'all'})")
                continue

            # Initial variable list + units snapshot
            if not self.variables:
                self.variables = df["Variables"].tolist()
                if len(df.columns) > 1:
                    self.units = dict(zip(df["Variables"], df.iloc[:, 1]))

            df = df.set_index("Variables")

            # record Ref._STD if present (per full row label)
            if "Ref._STD" in df.columns:
                for var_full, val in df["Ref._STD"].items():
                    try:
                        v = float(val)
                    except Exception:
                        v = np.nan
                    if np.isfinite(v) and v > 0:
                        self.obs_std_raw[var_full] = v

            # Append metric columns
            for metric in self.metrics:
                if metric not in df.columns:
                    continue
                if self.tables[metric].empty:
                    self.tables[metric] = pd.DataFrame(index=df.index)
                series = pd.to_numeric(df[metric], errors="coerce")
                self.tables[metric][period] = series.reindex(df.index)

        for m in self.metrics:
            if not self.tables[m].empty:
                self.tables[m] = self._sort_period_columns(self.tables[m])

        # If a catalog is provided, filter rows to triplets and set base-var display order
        allowed = self._allowed_triplets()
        if allowed is not None:
            matched_rows = set()
            matched_bases = []
            base_seen = set()

            for metric, df in self.tables.items():
                if df is None or df.empty:
                    continue
                keep_rows = []
                for row in df.index:
                    var, region, src = self._parse_triplet_from_row(row)
                    trip = (var, region, src if src else None)
                    if trip in allowed:
                        keep_rows.append(row)
                        base, _ = self._split_var_obs(row)
                        matched_rows.add(row)
                        if base not in base_seen:
                            base_seen.add(base)
                            matched_bases.append(base)

                self.tables[metric] = df.loc[keep_rows]

            # honor catalog base order if present
            if self._plot_order_bases:
                ordered = [b for b in self._plot_order_bases if b in matched_bases]
            else:
                ordered = matched_bases
            self.variables = ordered

            print(f"[INFO] Catalog filter: matched rows={len(matched_rows)}, matched base vars={len(self.variables)}")

        # ---- Fold Ref._STD values to base-level obs std (median over multiple obs) ----
        from collections import defaultdict
        acc = defaultdict(list)
        for full, val in self.obs_std_raw.items():
            base, _ = self._split_var_obs(full)
            if np.isfinite(val) and val > 0:
                acc[base].append(float(val))
        self.obs_std = {b: float(np.nanmedian(vs)) for b, vs in acc.items() if len(vs)}

        print(f"[INFO] Loaded metrics for {len(self.variables)} variables "
              f"({len(select_vars) if select_vars else 'all'} selected).")

    # ----------------------- variable collapsing -----------------------
    def set_obs_priority(self, obs_list):
        """Optional: set preferred observation order for collapsing."""
        self.obs_priority = list(obs_list)

    def _best_coverage(self, df, rows):
        counts = {r: df.loc[r].notna().sum() for r in rows if r in df.index}
        return max(counts, key=counts.get) if counts else rows[0]

    def _best_by_metric(
        self, df, rows, metric_name, strategy,
        use_ref, refp, prefer_zeros_as_invalid_for_pcc=True,
    ):
        if use_ref and refp in df.columns:
            col = df[refp]
            vals = {}
            for r in rows:
                v = col.loc[r] if r in col.index else np.nan
                if metric_name == "Correlation":
                    vals[r] = v if (pd.notna(v) and np.isfinite(v) and abs(v) <= 1.0) else np.nan
                else:
                    vals[r] = v
            finite = {r: v for r, v in vals.items() if pd.notna(v) and np.isfinite(v)}
            if not finite:
                return None
            if strategy == "min_abs":
                return min(finite, key=lambda r: abs(finite[r]))
            if strategy == "min":
                return min(finite, key=lambda r: finite[r])
            if strategy == "max":
                if metric_name == "Correlation" and prefer_zeros_as_invalid_for_pcc:
                    finite = {r: v for r, v in finite.items() if v != 0}
                    if not finite:
                        return None
                return max(finite, key=lambda r: finite[r])
            return None
        else:
            vals = {}
            for r in rows:
                if r not in df.index:
                    vals[r] = np.nan; continue
                s = df.loc[r]
                if metric_name == "Correlation":
                    s = s.where(s.abs() <= 1)
                    if prefer_zeros_as_invalid_for_pcc:
                        s = s.replace(0, np.nan)
                if strategy == "min_abs":
                    vals[r] = s.abs().mean(skipna=True)
                elif strategy in ("min","max"):
                    vals[r] = s.mean(skipna=True)
            finite = {r: v for r, v in vals.items() if pd.notna(v) and np.isfinite(v)}
            if not finite:
                return None
            return (min if strategy in ("min_abs","min") else max)(finite, key=lambda r: finite[r])

    def collapse_variable_groups_per_metric(
        self,
        strategy_map=None,
        use_ref=True,
        ref_period=None,
        prefer_zeros_as_invalid_for_pcc=True,
    ):
        """
        Collapse multiple obs per base variable to ONE row PER METRIC.
        Bias picks min |bias|, RMSE picks min rmse, PCC picks max corr (at ref if available).
        The kept row (obs) per metric and base variable is stored in self.selected_rows.
        """
        if all(tbl.empty for tbl in self.tables.values()):
            print("[WARN] collapse_variable_groups_per_metric(): no tables to collapse.")
            return

        strategy_map = strategy_map or {"Mean_Bias": "min_abs", "RMSE": "min", "Correlation": "max"}
        refp = ref_period or self.reference_period

        by_base = {}
        for r in list(self.tables.values())[0].index:
            base, _ = self._split_var_obs(r)
            by_base.setdefault(base, []).append(r)

        new_tables = {}
        self.selected_rows = {}
        for metric, df in self.tables.items():
            if df.empty:
                new_tables[metric] = df
                continue

            strat = strategy_map.get(metric, "best_coverage")
            out_rows = []
            self.selected_rows[metric] = {}
            for base, rows in by_base.items():
                keep = self._best_by_metric(
                    df=df, rows=rows, metric_name=metric, strategy=strat,
                    use_ref=use_ref, refp=refp, prefer_zeros_as_invalid_for_pcc=prefer_zeros_as_invalid_for_pcc
                )
                if keep is None and strat == "priority" and hasattr(self, "obs_priority"):
                    # honor explicit priority if provided
                    for pref in self.obs_priority:
                        for r in rows:
                            _, obs = self._split_var_obs(r)
                            if (obs or "").lower() == pref.lower():
                                keep = r; break
                        if keep: break
                if keep is None:
                    keep = self._best_coverage(df, rows) or rows[0]

                self.selected_rows[metric][base] = keep
                s = df.loc[keep].copy() if keep in df.index else pd.Series(index=df.columns, dtype=float)
                s.name = base
                out_rows.append(s)

            new_df = pd.DataFrame(out_rows).reindex(columns=df.columns)
            new_tables[metric] = new_df

        self.tables = new_tables
        # preserve catalog order if available, else first appearance order
        if self._plot_order_bases:
            self.variables = [b for b in self._plot_order_bases if b in list(self.tables[self.metrics[0]].index)]
        else:
            base_order, seen = [], set()
            for r in list(self.tables[self.metrics[0]].index):
                b, _ = self._split_var_obs(r)
                if b not in seen:
                    base_order.append(b); seen.add(b)
            self.variables = base_order

        print(f"[INFO] Collapsed per metric: kept {len(self.variables)} base variables with per-metric best obs.")

    # ----------------------- arrays & normalization -----------------------
    def as_array(self, all_periods, variables, m):
        """Return (n_periods, n_variables) for metric m."""
        if self.tables[m].empty:
            return np.full((len(all_periods), len(variables)), np.nan)
        df = self.tables[m].reindex(columns=all_periods).reindex(index=variables)
        return df.T.to_numpy(dtype=float)

    def _robust_sigma(self, arr):
        """
        Robust spread across periods (per variable).
        Uses max(MAD*1.4826, IQR/1.349). For columns with no finite data, returns eps.
        """
        eps = 1e-3
        x = np.asarray(arr, dtype=float)  # (period, variable)
        if x.ndim != 2 or x.size == 0:
            return np.full(x.shape[1] if x.ndim == 2 else 1, eps, dtype=float)

        col_has_data = np.isfinite(x).any(axis=0)
        sigma = np.full(x.shape[1], eps, dtype=float)
        if col_has_data.any():
            x_ok = x[:, col_has_data]
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                med = np.nanmedian(x_ok, axis=0)
                mad = np.nanmedian(np.abs(x_ok - med), axis=0)
                s_mad = 1.4826 * mad
                q75 = np.nanpercentile(x_ok, 75, axis=0)
                q25 = np.nanpercentile(x_ok, 25, axis=0)
                s_iqr = (q75 - q25) / 1.349
                s = np.nanmax(np.stack([s_mad, s_iqr], axis=0), axis=0)
            s = np.where((~np.isfinite(s)) | (s <= 0), eps, s)
            sigma[col_has_data] = s
        return sigma

    @staticmethod
    def piControl_sigma_vector(df, detrend=True, method="std", floor_frac=0.001, min_windows=5):
        """
        Internal-variability σ of a metric across piControl windows, per variable.

        df: rows = period labels ("0001-0050", ...), columns = variables (e.g. ds["RMSE"].to_pandas()).
        Per variable:
          1) optionally remove a least-squares linear trend vs window-centre year (residual drift
             would otherwise inflate σ);
          2) σ of the residuals, with the degrees of freedom of the fit removed:
               method="std"    sample std (ddof = 1 + trend); most efficient for ~10 windows
               method="robust" max(1.4826*MAD, IQR/1.349), scaled by sqrt((n-1)/(n-1-k));
                               outlier resistant but noisy for ~10 windows (can underestimate σ)
               method="max"    max of the two (never smaller than either)
          3) floor σ at floor_frac * |mean| so constant metrics (σ = 0) do not blow up.
        Variables with fewer than min_windows finite windows get NaN.
        Returns {var: σ}.
        """
        centers = np.array([np.mean([int(t) for t in str(p).split("-")]) for p in df.index], float)
        out = {}
        for v in df.columns:
            y = df[v].to_numpy(dtype=float)
            ok = np.isfinite(y)
            n = int(ok.sum())
            if n < min_windows:
                out[str(v)] = np.nan
                continue
            x, y = centers[ok], y[ok]
            k = 0
            if detrend:
                slope, icpt = np.polyfit(x - x.mean(), y, 1)
                r = y - (slope * (x - x.mean()) + icpt)
                k = 1
            else:
                r = y - np.median(y)
            s_std = float(np.std(r, ddof=1 + k))
            mad = 1.4826 * np.median(np.abs(r - np.median(r)))
            iqr = (np.percentile(r, 75) - np.percentile(r, 25)) / 1.349
            s_rob = max(mad, iqr) * np.sqrt((n - 1) / (n - 1 - k))
            s = {"std": s_std, "robust": s_rob, "max": max(s_std, s_rob)}[method]
            out[str(v)] = float(max(s, floor_frac * abs(np.mean(y))))
        return out

    def _ref_vector(self, arr, all_periods):
        """
        Reference vector for normalization:
          - If self.reference_window is set: mean across segments whose start in [w0, w1]
          - Else if self.reference_period is in all_periods: that row
          - Else: first row
        """
        x = np.asarray(arr, dtype=float)
        starts = np.array([self._parse_period(p)[0] for p in all_periods], dtype=int)
        if self.reference_window is not None:
            w0, w1 = map(int, self.reference_window)
            idx = np.where((starts >= w0) & (starts <= w1))[0]
            if idx.size > 0:
                return np.nanmean(x[idx, :], axis=0)
        if (self.reference_period is not None) and (self.reference_period in all_periods):
            return x[all_periods.index(self.reference_period), :].astype(float)
        return x[0, :].astype(float)

    def ratio_to_ref(self, arr, metric_name, all_periods):
        """
        Normalized ratio vs reference.
        If external reference vectors are set, use them; otherwise use _ref_vector().
        RMSE -> percent ((X - Xref)/Xref * 100); Bias -> fraction of |Xref|; PCC -> ratio to ref.
        """
        out = np.full_like(arr, np.nan, dtype=float)
        if not self.normalize_across_rows:
            return out

        variables = self.variables  # current variable order
        ext_vec = self._external_ref_vector(metric_name, variables)
        ref_vals = ext_vec if ext_vec is not None else self._ref_vector(arr, all_periods)

        eps = 1e-6
        mask = np.isfinite(ref_vals) & (np.abs(ref_vals) > eps)
        if not np.any(mask):
            return out
        if metric_name == "RMSE":
            out[:, mask] = (arr[:, mask] - ref_vals[mask]) / ref_vals[mask] * 100.0
        elif metric_name == "Mean_Bias":
            out[:, mask] = (arr[:, mask] - ref_vals[mask]) / np.abs(ref_vals[mask])
        else:
            out[:, mask] = arr[:, mask] / ref_vals[mask]
        return out

    def rel_change_to_ref(self, arr, all_periods):
        """
        Relative change vs reference: (X_t - X_ref) / X_ref.
        Uses external ref if available; else falls back to window/single-row ref.
        """
        out = np.full_like(arr, np.nan, dtype=float)
        if not self.normalize_across_rows:
            return out

        variables = self.variables
        ext_vec = self._external_ref_vector("RMSE", variables)  # metric name not used here; any key works
        ref_vals = ext_vec if ext_vec is not None else self._ref_vector(arr, all_periods)

        mask = np.isfinite(ref_vals) & (ref_vals != 0)
        if not np.any(mask):
            return out
        out[:, mask] = (arr[:, mask] - ref_vals[mask]) / ref_vals[mask]
        return out

    def z_to_ref(self, arr, metric_name="RMSE"):
        """
        (X - X_ref) / σ_ref with the external reference mean and σ (e.g. 500-yr piControl,
        σ from piControl_sigma_vector()). All NaN unless both external vectors are set.
        """
        out = np.full_like(arr, np.nan, dtype=float)
        ref_vals = self._external_ref_vector(metric_name, self.variables)
        sig_vals = self._external_ref_vector(f"{metric_name}_sigma", self.variables)
        if ref_vals is None or sig_vals is None:
            return out
        mask = np.isfinite(ref_vals) & np.isfinite(sig_vals) & (sig_vals > 0)
        out[:, mask] = (arr[:, mask] - ref_vals[mask]) / sig_vals[mask]
        return out

    def _delta_sigma(self, arr, periods, ref_period=None, is_bias=False):
        """
        Δσ = (X(t) - X_ref) / σ_temporal  (if is_bias: use |X| - |X_ref|)
        Uses window-mean reference if self.reference_window is set; otherwise ref_period/single-row.
        (External ref is NOT applied here; Δσ measures change vs internal time reference.)
        """
        x = np.asarray(arr, dtype=float)
        sigma = self._robust_sigma(x)
        denom = np.where((~np.isfinite(sigma)) | (sigma <= 0), 1e-3, sigma)

        if self.reference_window is not None:
            ref_vec = self._ref_vector(x, periods)
        else:
            if (ref_period is not None) and (ref_period in periods):
                ref_vec = x[periods.index(ref_period), :]
            else:
                ref_vec = self._ref_vector(x, periods)

        num = (np.abs(x) - np.abs(ref_vec)) if is_bias else (x - ref_vec)
        return num / denom

    # ----------------------- Stage 2: build diagnostic dataset -----------------------
    def build_dataset(self):
        """
        Create an xr.Dataset with:
          - Bias, RMSE, PCC (raw vs obs)
          - Bias_sigma_obs, RMSE_sigma_obs (absolute level in obs-σ units)
          - Bias_delta_obs, RMSE_delta_obs (change from reference in obs-σ units)
          - Bias_ratio_ref, RMSE_ratio_ref (global/external reference)
          - Bias_rel_ref, RMSE_rel_ref (zero-centered relative change vs ref, in %)
          - Bias_delta_sigma, RMSE_delta_sigma (change from reference in temporal-σ units)
        """
        if all(tbl.empty for tbl in self.tables.values()):
            raise RuntimeError("[build_dataset] No metric tables. Did you run extract_metrics()?")

        all_periods = sorted(
            set().union(*[set(t.columns) for t in self.tables.values() if not t.empty]),
            key=lambda c: self._parse_period(c)[0]
        )
        if not all_periods:
            raise RuntimeError("[build_dataset] No periods found after gathering metrics.")

        source = ["spinup" if self._parse_period(p)[0] <= self.picontrol_switch_year else "piControl"
                  for p in all_periods]
        variables = self.variables[:]

        p_start = [self._parse_period(p)[0] for p in all_periods]
        p_end   = [self._parse_period(p)[1] for p in all_periods]

        bias = self.as_array(all_periods, variables, "Mean_Bias")
        rmse = self.as_array(all_periods, variables, "RMSE")
        pcc  = self.as_array(all_periods, variables, "Correlation")

        # info logs for all-nan columns
        self._log_allnan_columns_once("RMSE", rmse, variables)
        self._log_allnan_columns_once("Bias", bias, variables)
        self._log_allnan_columns_once("PCC",  pcc,  variables)

        # ratios, rel changes (honor external reference if set)
        bias_ratio = self.ratio_to_ref(bias, self.bias_ratio_metric, all_periods)
        rmse_ratio = self.ratio_to_ref(rmse, "RMSE", all_periods)

        # obs-normalize (fallback 1e-3)
        obs_sigma_vec = np.array([self.obs_std.get(v, 1e-3) for v in variables], dtype=float)
        obs_sigma_vec = np.where((~np.isfinite(obs_sigma_vec)) | (obs_sigma_vec <= 0), 1e-3, obs_sigma_vec)
        bias_sigma_obs = bias / obs_sigma_vec
        rmse_sigma_obs = rmse / obs_sigma_vec

        # deltas in obs-σ units (internal reference)
        if self.reference_window is not None:
            ref_sigma_vec_bias = self._ref_vector(bias_sigma_obs, all_periods)
            ref_sigma_vec_rmse = self._ref_vector(rmse_sigma_obs, all_periods)
            bias_delta_obs = bias_sigma_obs - ref_sigma_vec_bias
            rmse_delta_obs = rmse_sigma_obs - ref_sigma_vec_rmse
        else:
            ref_idx = all_periods.index(self.reference_period) if self.reference_period in all_periods else 0
            bias_delta_obs = bias_sigma_obs - bias_sigma_obs[ref_idx, :]
            rmse_delta_obs = rmse_sigma_obs - rmse_sigma_obs[ref_idx, :]
            bias_delta_obs[ref_idx, :] = 0.0
            rmse_delta_obs[ref_idx, :] = 0.0

        # temporal-σ deltas
        bias_delta_sigma = self._delta_sigma(bias, all_periods, ref_period=self.reference_period, is_bias=True)
        rmse_delta_sigma = self._delta_sigma(rmse, all_periods, ref_period=self.reference_period, is_bias=False)

        # zero-centered relative change vs ref (%)
        bias_rel_ref = self.rel_change_to_ref(bias, all_periods) * 100.0
        rmse_rel_ref = self.rel_change_to_ref(rmse, all_periods) * 100.0

        # change vs external ref in units of its internal-variability σ
        rmse_z_ref = self.z_to_ref(rmse, "RMSE")

        ds = xr.Dataset(
            data_vars={
                # raw
                "Bias":              (("period","variable"), bias),
                "RMSE":              (("period","variable"), rmse),
                "PCC":               (("period","variable"), pcc),

                # obs-anchored absolute and deltas
                "Bias_sigma_obs":    (("period","variable"), bias_sigma_obs),
                "RMSE_sigma_obs":    (("period","variable"), rmse_sigma_obs),
                "Bias_delta_obs":    (("period","variable"), bias_delta_obs),
                "RMSE_delta_obs":    (("period","variable"), rmse_delta_obs),

                # ratio to global/external ref (note units below)
                "Bias_ratio_ref":    (("period","variable"), bias_ratio),
                "RMSE_ratio_ref":    (("period","variable"), rmse_ratio),

                # relative change vs ref (%)
                "Bias_rel_ref":      (("period","variable"), bias_rel_ref),
                "RMSE_rel_ref":      (("period","variable"), rmse_rel_ref),

                # change vs external ref in units of its σ (needs RMSE_sigma reference)
                "RMSE_z_ref":        (("period","variable"), rmse_z_ref),

                # temporal-σ deltas (significance of drift)
                "Bias_delta_sigma":  (("period","variable"), bias_delta_sigma),
                "RMSE_delta_sigma":  (("period","variable"), rmse_delta_sigma),
            },
            coords={
                "period":        ("period", all_periods),
                "period_start":  ("period", p_start),
                "period_end":    ("period", p_end),
                "variable":      ("variable", variables),
                "source":        ("period", source),
            },
            attrs={
                "has_rel_ref": 1,
                "spinup_name": self.spinup_name,
                "picontrol_name": self.picontrol_name or "",
                "season": self.season,
                "reference_mode": "window" if getattr(self, "reference_window", None) else "single",
                "reference_period": str(self.reference_period) if self.reference_period is not None else "",
                "reference_window": (
                    f"{int(self.reference_window[0])}-{int(self.reference_window[1])}"
                    if getattr(self, "reference_window", None) else ""
                ),
                "normalize_ratio_to_ref": json.dumps(bool(self.normalize_across_rows)),
                "units_json": json.dumps(self.units, ensure_ascii=False),
                "metric_short_json": json.dumps(self.metric_short),
                "picontrol_switch_year": int(self.picontrol_switch_year),
                "spinup_years": int(self.spinup_years),
                # Clarify units for normalized fields:
                "rmse_ratio_ref_units": "%",                 # ((X - Xref)/Xref)*100
                "rmse_z_ref_units": "sigma_ref",             # (X - Xref)/σref
                "bias_ratio_ref_units": "fraction_of_abs_ref" # (X - Xref)/|Xref|
            }
        )
        self.ds = ds
        return ds

    # ----------------------- IO -----------------------
    def to_netcdf(
        self, path,
        mode="w",
        engine="netcdf4",
        overwrite=True,
        atomic=True,
        compress=True,
        dtype=float
    ):
        """Safely write NetCDF."""
        if self.ds is None:
            raise RuntimeError("No dataset to save. Call build_dataset() first.")
        out_dir = os.path.dirname(path) or "."
        os.makedirs(out_dir, exist_ok=True)

        enc = None
        if compress:
            enc = {name: {"zlib": True, "complevel": 4, "dtype": dtype}
                   for name in self.ds.data_vars}

        try:
            self.ds.close()  # harmless for in-memory
        except Exception:
            pass

        target = path
        tmp = f"{path}.tmp" if atomic else path

        if not atomic and overwrite and os.path.exists(target):
            try:
                os.remove(target)
            except PermissionError as e:
                raise PermissionError(f"[to_netcdf] Cannot remove existing file: {target} ({e})")

        try:
            self.ds.to_netcdf(tmp, mode=mode, engine=engine, encoding=enc)
        except PermissionError as e:
            raise PermissionError(
                f"[to_netcdf] Permission denied when writing '{tmp}'. "
                f"Check directory/file ACLs and group ownership.\n{e}"
            )

        if atomic:
            try:
                os.replace(tmp, target)  # atomic on POSIX
            except PermissionError as e:
                raise PermissionError(
                    f"[to_netcdf] Cannot replace '{target}' with temp file. "
                    f"Check write perms on directory.\n{e}"
                )
        print(f"[INFO] Wrote diagnostic NetCDF: {target}")

    def load_dataset(self, path, engine="netcdf4", chunks=None, cache=False):
        """Load the diagnostic dataset. Use cache=False to avoid lingering open handles."""
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
        try:
            if self.ds is not None:
                self.ds.close()
        except Exception:
            pass
        self.ds = xr.open_dataset(path, engine=engine, chunks=chunks, cache=cache)
        try:
            self.units = json.loads(self.ds.attrs.get("units_json", "{}"))
            self.metric_short = json.loads(self.ds.attrs.get("metric_short_json", "{}"))
        except Exception:
            pass
        print(f"[INFO] Loaded diagnostic NetCDF: {path}")

    def close(self):
        """Explicitly close the currently loaded dataset to release file handles."""
        try:
            if isinstance(self.ds, xr.Dataset):
                self.ds.close()
            print("[INFO] Closed diagnostic dataset handle.")
        except Exception:
            pass

    def summary(self):
        """Print diagnostic info about the loaded dataset."""
        if self.ds is None:
            print("[WARN] No dataset loaded.")
            return
        ref_mode = self.ds.attrs.get("reference_mode", "")
        ref_per  = self.ds.attrs.get("reference_period", "")
        ref_win  = self.ds.attrs.get("reference_window", "")
        print("\n[=== Dataset Diagnostic Summary ===]")
        print(f"Experiment: {self.spinup_name} vs {self.picontrol_name}")
        print(f"Season:     {self.season}")
        if ref_mode == "window" and ref_win:
            print(f"Reference:  window-mean over {ref_win}")
        elif ref_per:
            print(f"Reference:  single period {ref_per}")
        else:
            print("Reference:  (not set)")
        print(f"Periods:    {len(self.ds['period'])} ({self.ds['period'].values[0]} → {self.ds['period'].values[-1]})")
        print(f"Variables:  {len(self.ds['variable'])}")
        print(f"Metrics:    {[v for v in self.ds.data_vars]}")
        print("---------------------------------------")
        for key in ["Bias", "RMSE", "PCC"]:
            if key in self.ds:
                arr = self.ds[key].values
                finite = np.isfinite(arr)
                ffrac = finite.mean()
                vmin = np.nanmin(arr) if ffrac > 0 else np.nan
                vmax = np.nanmax(arr) if ffrac > 0 else np.nan
                print(f"{key:15s} min={vmin:8.3f}  max={vmax:8.3f}  finite={ffrac:5.2f}")
        print("=======================================\n")

    # ----------------------- selection & plotting -----------------------
    def list_variables(self, verbose: bool = True):
        """Return a list of all variables currently in the dataset."""
        if self.ds is None:
            raise RuntimeError("No dataset loaded. Call build_dataset() or load_dataset() first.")
        vars_ = list(self.ds["variable"].values)
        if verbose:
            print(f"[INFO] Found {len(vars_)} variables:")
            for v in vars_:
                print(f"  - {v}")
        return vars_
    
    def _field_kind_and_label(self, metric: str, field: str):
        """
        Return ('kind', colorbar_label, unit_hint)
        kind ∈ {'pcc','ratio_pct','ratio_frac','delta_sigma','sigma_obs','raw'}
        """
        if metric == "Correlation" or field == "PCC":
            return "pcc", "PCC", None
        if field.endswith("_rel_ref"):
            return "ratio_pct", "Δ / ref (%)", "%"
        if field.endswith("_ratio_ref"):
            return ("ratio_pct", "Δ / ref (%)", "%") if metric == "RMSE" else ("ratio_frac", "Δ / |ref|", None)
        if field.endswith("_z_ref"):
            return "delta_sigma", "Δ / σ_ref", None
        if field.endswith("_delta_sigma") or field.endswith("_delta_obs"):
            return "delta_sigma", "Δ / σ", None
        if field.endswith("_sigma_obs"):
            return "sigma_obs", f"{metric} / σ_obs", None
        return "raw", self.metric_short.get(metric, metric), None

    def select_variables(self, keep):
        keep = [v for v in keep if v in list(self.ds["variable"].values)]
        self.ds = self.ds.sel(variable=keep)

    # -------- Period selection helpers (new) --------
    def select_periods(self, keep_periods):
        """
        Keep only selected period labels in the gathered tables (pre-build).
        Call this after extract_metrics() and before build_dataset().
        """
        if not isinstance(keep_periods, (list, tuple, set)):
            keep_periods = [keep_periods]
        keep = [str(p) for p in keep_periods]
        for m, df in self.tables.items():
            if not df.empty:
                cols = [c for c in df.columns if c in keep]
                if cols:
                    self.tables[m] = df[cols]
                else:
                    self.tables[m] = pd.DataFrame(index=df.index)
        # adjust spinup_years bookkeeping (optional)
        try:
            starts = [self._parse_period(p)[0] for p in keep if isinstance(p, str) and "-" in p]
            if starts:
                self.spinup_years = max(starts) + 49
        except Exception:
            pass

    def select_periods_on_ds(self, keep_periods):
        """Keep only selected period labels on the built dataset (post-build)."""
        if self.ds is None:
            raise RuntimeError("No dataset loaded/built. Call build_dataset() or load_dataset() first.")
        if not isinstance(keep_periods, (list, tuple, set)):
            keep_periods = [keep_periods]
        keep = [p for p in keep_periods if p in list(self.ds["period"].values)]
        if not keep:
            raise ValueError(f"None of the requested periods are present: {keep_periods}")
        self.ds = self.ds.sel(period=keep)

    def _select_field(self, metric: str, mode: str = "auto") -> str:
        if metric == "Correlation":
            return "PCC"
        if mode == "delta_obs":
            return "Bias_delta_obs" if metric == "Mean_Bias" else "RMSE_delta_obs"
        if mode == "sigma_obs":
            return "Bias_sigma_obs" if metric == "Mean_Bias" else "RMSE_sigma_obs"
        if mode == "raw":
            return "Bias" if metric == "Mean_Bias" else "RMSE"
        if mode == "delta":
            return "Bias_delta_sigma" if metric == "Mean_Bias" else "RMSE_delta_sigma"
        if mode == "ratio":
            return "Bias_ratio_ref" if metric == "Mean_Bias" else "RMSE_ratio_ref"
        if mode == "rel":
            return "Bias_rel_ref" if metric == "Mean_Bias" else "RMSE_rel_ref"
        if mode == "z" and metric == "RMSE":
            return "RMSE_z_ref"

        pref = {
            "Mean_Bias": ["Bias_rel_ref","Bias_delta_obs","Bias_delta_sigma","Bias_ratio_ref","Bias_sigma_obs","Bias"],
            "RMSE":      ["RMSE_rel_ref","RMSE_delta_obs","RMSE_delta_sigma","RMSE_ratio_ref","RMSE_sigma_obs","RMSE"],
        }
        for name in pref.get(metric, []):
            if self.ds is not None and name in self.ds:
                return name
        return "PCC"

    def _log_allnan_columns_once(self, name, arr, variables):
        key = ("allnan", name)
        if key in self._once: return
        col_has_data = np.isfinite(arr).any(axis=0)
        bad = [variables[i] for i, ok in enumerate(col_has_data) if not ok]
        if bad:
            head = ", ".join(map(str, bad[:12])) + ("..." if len(bad) > 12 else "")
            print(f"[INFO] {name}: {len(bad)} variables have no finite values across periods: {head}")
        self._once.add(key)

    # ---------- Color utilities ----------
    def _finite_or_default(self, arr, default=0.0):
        """Return finite subset of arr; if empty, return [default]."""
        a = np.asarray(arr, dtype=float)
        a = a[np.isfinite(a)]
        if a.size == 0:
            return np.array([default], dtype=float)
        return a

    def _finite_extent(self, vals, pct=98, fallback=1.0):
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return fallback
        return max(np.nanpercentile(np.abs(vals), pct), 1e-6)

    def _make_colormap_and_norm(
        self, metric: str, field: str, arr_all: np.ndarray,
        discrete: bool = True, bins=None, vmin=None, vmax=None,
        pcc_range: str = "pos"  # 'pos'=[0,1], 'sym'=[-1,1]
    ):
        finite = self._finite_or_default(arr_all, default=0.0)
        kind, cbar_label, _ = self._field_kind_and_label(metric, field)

        # ==========================================================
        # --- PCC (sequential, typically Viridis) ------------------
        # ==========================================================
        if kind == "pcc":
            base_cmap = mpl.colormaps.get_cmap("viridis")
            if discrete:
                edges = np.asarray(bins, float) if bins is not None else (
                    np.arange(0.0, 1.0001, 0.05)
                    if pcc_range == "pos" else np.arange(-1.0, 1.0001, 0.1)
                )
                if edges.size < 3 or not np.all(np.diff(edges) > 0):
                    edges = np.arange(0.0, 1.0001, 0.05)

                N = len(edges) - 1
                # sample N+2 colors, use middle N for bins
                cols_full = np.asarray(base_cmap(np.linspace(0, 1, N + 2)))
                bin_cols = cols_full[1:-1]
                cmap = ListedColormap(bin_cols, name="pcc_bins")
                cmap.set_under(cols_full[0])
                cmap.set_over(cols_full[-1])

                norm = BoundaryNorm(edges, N, clip=False)
            else:
                if pcc_range == "sym":
                    lo, hi = -1.0, 1.0
                else:
                    lo = max(0.0, np.nanpercentile(finite, 2))
                    hi = min(1.0, np.nanpercentile(finite, 98))
                cmap = base_cmap
                norm = mpl.colors.Normalize(vmin=lo, vmax=hi)
                edges = None

            cmap = cmap.copy()
            cmap.set_bad((0.85, 0.85, 0.85, 1.0))
            return cmap, norm, cbar_label, (edges if discrete else None)

        # ==========================================================
        # --- Diverging (ratio %, frac, Δ/σ) ------------------------
        # ==========================================================
        if kind in {"ratio_pct", "ratio_frac", "delta_sigma"}:
            base_cmap = mpl.colormaps.get_cmap("RdBu_r")
            if discrete:
                if bins is None:
                    edges = (np.array(
                        [-20, -16, -12, -8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8, 12, 16, 20],
                        float
                    ) if kind == "ratio_pct" else np.array(
                        [-3, -2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 3],
                        float
                    ))
                else:
                    edges = np.asarray(bins, float)
                    if edges.size < 3 or not np.all(np.diff(edges) > 0):
                        edges = np.array(
                            [-10, -8, -6, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 6, 8, 10], float
                        )

                N = len(edges) - 1
                cols_full = np.asarray(base_cmap(np.linspace(0, 1, N + 2)))
                bin_cols = cols_full[1:-1]

                tmp_cmap = ListedColormap(bin_cols, name="div_bins_tmp")
                tmp_cmap = self._whiten_bins_around_zero(tmp_cmap, edges, n_each_side=1)
                bin_cols = np.asarray(tmp_cmap.colors)

                cmap = ListedColormap(bin_cols, name="div_bins")
                cmap.set_under(cols_full[0])
                cmap.set_over(cols_full[-1])

                norm = BoundaryNorm(edges, N, clip=False)
            else:
                spread = self._finite_extent(finite, pct=98, fallback=10.0)
                vvmin = -(vmin if vmin is not None else spread)
                vvmax = +(vmax if vmax is not None else spread)
                cmap = base_cmap
                norm = TwoSlopeNorm(vmin=vvmin, vcenter=0.0, vmax=vvmax)
                edges = None

            cmap = cmap.copy()
            cmap.set_bad((0.85, 0.85, 0.85, 1.0))
            return cmap, norm, cbar_label, (edges if discrete else None)

        # ==========================================================
        # --- Raw / sigma_obs (sequential, Blues) -------------------
        # ==========================================================
        base_cmap = mpl.colormaps.get_cmap("Blues")
        if discrete:
            edges = np.asarray(bins, float) if bins is not None else None
            if edges is None or edges.size < 3 or not np.all(np.diff(edges) > 0):
                lo = np.nanpercentile(finite, 2) if np.isfinite(finite).any() else 0.0
                hi = np.nanpercentile(finite, 98) if np.isfinite(finite).any() else 1.0
                edges = np.linspace(lo, hi, 12)

            N = len(edges) - 1
            cols_full = np.asarray(base_cmap(np.linspace(0, 1, N + 2)))
            bin_cols = cols_full[1:-1]
            cmap = ListedColormap(bin_cols, name="seq_bins")
            cmap.set_under(cols_full[0])
            cmap.set_over(cols_full[-1])
            norm = BoundaryNorm(edges, N, clip=False)
        else:
            lo = np.nanpercentile(finite, 2) if np.isfinite(finite).any() else 0.0
            hi = np.nanpercentile(finite, 98) if np.isfinite(finite).any() else 1.0
            if vmin is not None:
                lo = vmin
            if vmax is not None:
                hi = vmax
            cmap = base_cmap
            norm = mpl.colors.Normalize(vmin=lo, vmax=hi)
            edges = None

        cmap = cmap.copy()
        cmap.set_bad((0.85, 0.85, 0.85, 1.0))
        return cmap, norm, cbar_label, (edges if discrete else None)

    def _whiten_bins_around_zero(self, cmap, edges, n_each_side=1):
        """
        Make the bins nearest 0 (one on each side) white for a discrete, diverging map.
        """
        edges = np.asarray(edges, dtype=float)
        if edges.ndim != 1 or edges.size < 3 or not np.all(np.diff(edges) > 0):
            N = max(2, int(getattr(cmap, "N", 256)))
            cols = cmap(np.linspace(0, 1, N))
            return ListedColormap(cols, name=getattr(cmap, "name", "listed"))

        N = edges.size - 1
        cols = np.array(cmap(np.linspace(0, 1, N)))
        cols = cols[:, :4]  # RGBA

        if not (edges[0] < 0.0 < edges[-1]):
            return ListedColormap(cols, name=getattr(cmap, "name", "listed"))

        k0 = int(np.searchsorted(edges, 0.0, side="right") - 1)
        white = np.array([1.0, 1.0, 1.0, 1.0])

        for i in range(n_each_side):
            k_pos = k0 + i
            if 0 <= k_pos < N:
                cols[k_pos] = white
        for i in range(n_each_side):
            k_neg = k0 - 1 - i
            if 0 <= k_neg < N:
                cols[k_neg] = white

        return ListedColormap(cols, name=getattr(cmap, "name", "listed"))

    # ---------- Rendering helper ----------
    def render_panel(
        self, *, E, V, common_vars, seasons, 
        experiments, explabels, 
        ax, mats, sm, title_panel: str,
        add_xlabel: bool, ylabel: str,
        draw_grid: bool, x_step: int, 
        fontz: int, hatch_nan: bool
    ):
        def tri_top(j, i):    return [(j, i), (j+1, i), (j+0.5, i+0.5)]   # DJF
        def tri_right(j, i):  return [(j+1, i), (j+1, i+1), (j+0.5, i+0.5)]# MAM
        def tri_bottom(j, i): return [(j, i+1), (j+1, i+1), (j+0.5, i+0.5)]# JJA
        def tri_left(j, i):   return [(j, i), (j, i+1), (j+0.5, i+0.5)]    # SON
        tri_for_season = {"DJF": tri_top, "MAM": tri_right, "JJA": tri_bottom, "SON": tri_left}

        edgecolor = "grey"; lw = 0.001
        for i in range(E):
            for j in range(V):
                for s in seasons:
                    val = mats[s][i, j]
                    if not np.isfinite(val):
                        if hatch_nan:
                            ax.add_patch(Polygon(tri_for_season[s](j, i), closed=True,
                                                 facecolor=(0.95, 0.95, 0.95, 1.0),
                                                 edgecolor=edgecolor, linewidth=lw, hatch="//",
                                                 rasterized=True))   # PDF hatch tiles do not scale in LaTeX
                        continue
                    ax.add_patch(Polygon(tri_for_season[s](j, i), closed=True,
                                         facecolor=sm.to_rgba(val),
                                         edgecolor=edgecolor, linewidth=lw))
        ax.set_xlim(0, V); ax.set_ylim(E, 0)

        # X ticks
        ax.set_xticks(np.arange(V) + 0.5)
        if x_step > 1:
            ax.set_xticks(np.arange(0, V, x_step) + 0.5)
            xticklabels = [str(common_vars[k]) for k in range(0, V, x_step)]
        else:
            xticklabels = [str(v) for v in common_vars]
        ax.set_xticklabels(xticklabels, rotation=90, fontsize=fontz * 0.90)

        # Y ticks
        ax.set_yticks(np.arange(E) + 0.5)
        ax.set_yticklabels([str(e) for e in explabels], fontsize=fontz * 0.95)

        if ylabel:
            ax.set_ylabel(ylabel, fontsize=fontz)
        if not add_xlabel:
            # leave space; shared xlabel will be added on the figure
            ax.set_xlabel("")

        if draw_grid:
            for j in range(V + 1):
                ax.plot([j, j], [0, E], color="grey", lw=0.01, zorder=0)
            for i in range(E + 1):
                ax.plot([0, V], [i, i], color="grey", lw=0.01, zorder=0)

        if title_panel:
            ax.set_title(title_panel, fontsize=fontz * 1.02, pad=6, loc="left")

    # --- helper: choose RMSE field by mode (add once if you don't already have it) ---
    def _field_for_rmse_from_mode(self, rmse_mode: str) -> str:
        return {
            "ratio":     "RMSE_ratio_ref",
            "rel":       "RMSE_rel_ref",
            "z":         "RMSE_z_ref",
            "delta":     "RMSE_delta_sigma",
            "sigma_obs": "RMSE_sigma_obs",
            "raw":       "RMSE",
        }.get(str(rmse_mode or "ratio").lower(), "RMSE_ratio_ref")


    # --- helper: draw one season diamond legend (DJF/MAM/JJA/SON) ---
    def _draw_season_box(self, ax, fontz: int = 12):
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_aspect("equal"); ax.axis("off")

        # frame
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="0.25", linewidth=1.1))

        tri_defs = {
            "DJF": [(0.0, 1.0), (1.0, 1.0), (0.5, 0.5)],
            "MAM": [(1.0, 0.0), (1.0, 1.0), (0.5, 0.5)],
            "JJA": [(0.0, 0.0), (1.0, 0.0), (0.5, 0.5)],
            "SON": [(0.0, 0.0), (0.0, 1.0), (0.5, 0.5)],
        }
        for s in ["DJF", "MAM", "JJA", "SON"]:
            ax.add_patch(Polygon(tri_defs[s], closed=True,
                                 facecolor="white", edgecolor="0.4", lw=1.0))
        label_pos = {"DJF": (0.50, 0.80), "MAM": (0.78, 0.50),
                     "JJA": (0.50, 0.20), "SON": (0.22, 0.50)}
        for s, (tx, ty) in label_pos.items():
            ax.text(tx, ty, s, ha="center", va="center",
                    fontsize=fontz * 0.90, weight="normal", color="black")


    # --- helper: create two perfectly aligned colorbars and a single season box ---
    def _add_aligned_colorbars(
        self, fig, ax_top, ax_bot,
        sm_top, sm_bot,
        label_top: str, label_bot: str,
        *, edges_top=None, edges_bot=None, ticks_bot=None, fontz: int = 12, show_season_box: bool = True
    ):
        """Create two perfectly aligned colorbars (with optional discrete boundaries) and one season box."""
        pos1 = ax_top.get_position(fig)   # Bbox in figure coords
        pos2 = ax_bot.get_position(fig)

        gap   = 0.012      # space between axes and colorbar (figure fraction)
        width = 0.016      # colorbar width (figure fraction)

        x = max(pos1.x1, pos2.x1) + gap            # common x for both
        h1 = min(0.60, pos1.height * 0.95)         # keep a small margin inside each subplot height
        h2 = min(0.60, pos2.height * 0.95)

        # vertically center each cbar within its subplot height
        y1 = pos1.y0 + 0.5 * (pos1.height - h1)
        y2 = pos2.y0 + 0.5 * (pos2.height - h2)

        # colorbar axes
        cax1 = fig.add_axes([x, y1, width, h1])
        cax2 = fig.add_axes([x, y2, width, h2])

        # top colorbar (PCC)
        if edges_top is not None:
            cb1 = fig.colorbar(sm_top, cax=cax1, boundaries=edges_top, extend="both")
            # tidy ticks for lots of bins
            if len(edges_top) > 16:
                step = max(1, len(edges_top)//12)
                idx = list(range(0, len(edges_top), step))
                if idx[-1] != len(edges_top)-1:
                    idx[-1] = len(edges_top)-1
                cb1.set_ticks([edges_top[k] for k in idx])
            else:
                cb1.set_ticks(edges_top)
            cb1.set_ticklabels([f"{v:g}" for v in cb1.get_ticks()])
        else:
            cb1 = fig.colorbar(sm_top, cax=cax1, extend="both")

        # bottom colorbar (RMSE)
        if edges_bot is not None:
            cb2 = fig.colorbar(sm_bot, cax=cax2, boundaries=edges_bot, extend="both")
            if len(edges_bot) > 16:
                step = max(1, len(edges_bot)//12)
                idx = list(range(0, len(edges_bot), step))
                if idx[-1] != len(edges_bot)-1:
                    idx[-1] = len(edges_bot)-1
                cb2.set_ticks([edges_bot[k] for k in idx])
            else:
                cb2.set_ticks(edges_bot)
            cb2.set_ticklabels([f"{v:g}" for v in cb2.get_ticks()])
        else:
            cb2 = fig.colorbar(sm_bot, cax=cax2, extend="both")
        if ticks_bot is not None:
            cb2.set_ticks(list(ticks_bot))
            cb2.set_ticklabels([f"{v:g}" for v in ticks_bot])

        cb1.ax.tick_params(labelsize=fontz*0.90, length=3)
        cb2.ax.tick_params(labelsize=fontz*0.90, length=3)
        if label_top:
            cb1.set_label(label_top, fontsize=fontz, labelpad=5)
        if label_bot:
            cb2.set_label(label_bot, fontsize=fontz, labelpad=5)

        # draw ONE season box above the top colorbar
        if show_season_box:
            # make it noticeably bigger and lower
            box_w = width * 4.5        # wider box (was 3.6)
            box_h = h1 * 1.05          # taller box (was 0.85)
            box_gap = 0.055            # more horizontal space

            # lower vertical placement (move closer to middle of panels)
            total_h = (y1 + h1) - y2   # total vertical span
            box_y = y2 + 0.5 * total_h - (box_h * 0.55)  # lower by ~half its height

            # position horizontally
            box_x = x + width + box_gap
            box_ax = fig.add_axes([box_x, box_y, box_w, box_h])
            self._draw_season_box(box_ax, fontz=fontz * 0.85)  # larger text for clarity

    def plot_heatmap_experiments_quadseason_dual(
        self,
        analyzers_by_exp_by_season: dict,
        experiments: list,
        explabels: list,
        period_label: str,
        fig_path,
        *,
        # PCC panel config
        pcc_bins=None,
        pcc_discrete=True,
        pcc_label="PCC (unitless)",
        # RMSE panel config
        rmse_mode="ratio",
        rmse_bins=None,
        rmse_discrete=True,
        rmse_label="Norm. RMSE (%)",
        rmse_ticks=None,           # explicit RMSE colorbar ticks (default: matplotlib's choice)
        # shared figure/config
        title: str = "",
        xlabel: str = "",          # shared at bottom (we put it on the bottom axes)
        ylabel_left: str = "",     # y-label for both rows (left)
        fontz: int = 12,
        figsize: tuple[float, float] | None = None,
        x_step: int = 1,
        show_season_box: bool = True,
        draw_grid: bool = False,
        hatch_nan: bool = False,
    ):
        seasons = ["DJF", "MAM", "JJA", "SON"]
        field_pcc  = "PCC"
        field_rmse = self._field_for_rmse_from_mode(rmse_mode)

        # --- common variables across all (exp × season) ---
        var_orders = []
        for e in experiments:
            row_vars = []
            for s in seasons:
                an = analyzers_by_exp_by_season[e][s]
                for fld in (field_pcc, field_rmse):
                    if fld not in an.ds:
                        raise RuntimeError(f"[{e}:{s}] field '{fld}' missing")
                row_vars.append(list(map(str, an.ds["variable"].values)))
            base = row_vars[0]
            common = [v for v in base if all(v in rv for rv in row_vars[1:])]
            var_orders.append(common)

        base = var_orders[0]
        common_vars = [v for v in base if all(v in vo for vo in var_orders[1:])]
        if not common_vars:
            raise RuntimeError("No common variables across experiments & seasons.")

        # --- gather matrices ---
        def _gather(field):
            mats, all_vals = {}, []
            for s in seasons:
                rows = []
                for e in experiments:
                    df = analyzers_by_exp_by_season[e][s].ds[field].to_pandas()
                    if period_label not in df.index:
                        raise RuntimeError(f"[{e}:{s}] period '{period_label}' not found for '{field}'")
                    rows.append(df.loc[period_label, common_vars].to_numpy(dtype=float))
                M = np.vstack(rows)  # (E,V)
                mats[s] = M
                all_vals.append(M[np.isfinite(M)])
            arr_all = np.concatenate(all_vals) if all_vals and all_vals[0].size else np.array([0.0])
            return mats, arr_all

        mats_pcc,  vals_pcc  = _gather(field_pcc)
        mats_rmse, vals_rmse = _gather(field_rmse)

        # --- colormaps/norms ---
        cmap_pcc,  norm_pcc,  cbl_pcc,  edges_pcc  = self._make_colormap_and_norm(
            metric="Correlation", field=field_pcc,  arr_all=vals_pcc,
            discrete=pcc_discrete, bins=pcc_bins
        )
        cmap_rmse, norm_rmse, cbl_rmse, edges_rmse = self._make_colormap_and_norm(
            metric="RMSE",       field=field_rmse, arr_all=vals_rmse,
            discrete=rmse_discrete, bins=rmse_bins
        )
        if pcc_label  is not None: cbl_pcc  = pcc_label
        if rmse_label is not None: cbl_rmse = rmse_label

        sm_pcc  = mpl.cm.ScalarMappable(norm=norm_pcc,  cmap=cmap_pcc);  sm_pcc.set_array(vals_pcc)
        sm_rmse = mpl.cm.ScalarMappable(norm=norm_rmse, cmap=cmap_rmse); sm_rmse.set_array(vals_rmse)

        E, V = next(iter(mats_pcc.values())).shape
        if figsize is None:
            figsize = (max(15, 0.70 * V), max(7.5, 2.0 + 1.6 * E))

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, constrained_layout=False, sharex=True)

        # --- draw panels (top: no xlabel; bottom: yes) ---
        self.render_panel(
            E=E, V=V, common_vars=common_vars, seasons=seasons, 
            experiments=experiments,explabels=explabels,
            ax=ax1, mats=mats_pcc, sm=sm_pcc,
            title_panel="(a) Spatial Pattern Correlation Coefficient (PCC)",
            add_xlabel=False, ylabel=ylabel_left,
            draw_grid=draw_grid, x_step=x_step, fontz=fontz, hatch_nan=hatch_nan
        )

        self.render_panel(
            E=E, V=V, common_vars=common_vars, seasons=seasons, 
            experiments=experiments,explabels=explabels,
            ax=ax2, mats=mats_rmse, sm=sm_rmse,
            title_panel="(b) Normalized Spatial Root-Mean-Square Error (RMSE)",
            add_xlabel=True, ylabel=ylabel_left,
            draw_grid=draw_grid, x_step=x_step, fontz=fontz, hatch_nan=hatch_nan
        )

        # --- consistent layout for stacked panels + long x labels ---
        fig.subplots_adjust(
            left=0.09, right=0.90,
            hspace=0.22,
            bottom=0.28 if V > 28 else (0.24 if V > 22 else 0.18),
            top=0.90
        )

        # --- colorbars (aligned) + single season box ---
        self._add_aligned_colorbars(
            fig, ax1, ax2,
            sm_top=sm_pcc, sm_bot=sm_rmse,
            label_top=cbl_pcc, label_bot=cbl_rmse, ticks_bot=rmse_ticks,
            fontz=fontz, show_season_box=show_season_box,
        )
        if title:
            fig.suptitle(title, fontsize=fontz * 1.05, fontweight="bold", y=0.995)
        # --- save ---
        out_dir = os.path.dirname(os.fspath(fig_path)) or "."
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig(fig_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
        plt.show(); plt.close()
        print(f"[INFO] Saved figure: {fig_path}")

    def plot_heatmap_regions_pcc_rmse(
        self,
        analyzers_by_region: dict,
        region_titles: dict,
        experiments: list,
        explabels: list,
        period_label: str,
        fig_path,
        *,
        pcc_bins=None,
        pcc_discrete=True,
        pcc_label="PCC (unitless)",
        rmse_mode="ratio",
        rmse_bins=None,
        rmse_discrete=True,
        rmse_label="Norm. RMSE (%)",
        title: str = "",
        fontz: int = 12,
        figsize: tuple[float, float] | None = None,
        x_step: int = 1,
        show_season_box: bool = True,
        draw_grid: bool = False,
        hatch_nan: bool = False,
        panel_wspace: float = 0.10,
        panel_hspace: float = 0.30,
    ):
        """Regions x {PCC, normalized RMSE} grid: one row per region, PCC in the left column,
        RMSE in the right column (same quad-season cells as plot_heatmap_experiments_quadseason_dual).

        analyzers_by_region : {region: {exp: {season: SpinupMetricAnalyzer with .ds loaded}}}
        region_titles       : {region: title}; rows follow this dict's order.
        Each column has ONE color scale shared by all regions so rows are directly comparable.
        Variables shown are those present for every region, experiment and season.
        """
        seasons = ["DJF", "MAM", "JJA", "SON"]
        regions = [r for r in region_titles if r in analyzers_by_region]
        if not regions:
            raise RuntimeError("No regions to plot.")
        field_pcc = "PCC"
        field_rmse = self._field_for_rmse_from_mode(rmse_mode)

        # --- variables common to every region x experiment x season ---
        common_vars = None
        for r in regions:
            for e in experiments:
                for s in seasons:
                    ds = analyzers_by_region[r][e][s].ds
                    for fld in (field_pcc, field_rmse):
                        if fld not in ds:
                            raise RuntimeError(f"[{r}:{e}:{s}] field '{fld}' missing")
                    vs = [str(v) for v in ds["variable"].values]
                    common_vars = vs if common_vars is None else [v for v in common_vars if v in vs]
        if not common_vars:
            raise RuntimeError("No variables common to all regions, experiments and seasons.")

        def _gather(region, field):
            mats = {}
            for s in seasons:
                rows = []
                for e in experiments:
                    df = analyzers_by_region[region][e][s].ds[field].to_pandas()
                    if period_label not in df.index:
                        raise RuntimeError(f"[{region}:{e}:{s}] period '{period_label}' not found for '{field}'")
                    rows.append(df.loc[period_label, common_vars].to_numpy(dtype=float))
                mats[s] = np.vstack(rows)          # (E, V)
            return mats

        mats_pcc = {r: _gather(r, field_pcc) for r in regions}
        mats_rmse = {r: _gather(r, field_rmse) for r in regions}

        def _all_vals(mats_by_region):
            vals = [m[np.isfinite(m)] for mats in mats_by_region.values() for m in mats.values()]
            vals = [v for v in vals if v.size]
            return np.concatenate(vals) if vals else np.array([0.0])

        vals_pcc, vals_rmse = _all_vals(mats_pcc), _all_vals(mats_rmse)

        # --- one color scale per column, shared by all regions ---
        cmap_pcc, norm_pcc, cbl_pcc, edges_pcc = self._make_colormap_and_norm(
            metric="Correlation", field=field_pcc, arr_all=vals_pcc, discrete=pcc_discrete, bins=pcc_bins)
        cmap_rmse, norm_rmse, cbl_rmse, edges_rmse = self._make_colormap_and_norm(
            metric="RMSE", field=field_rmse, arr_all=vals_rmse, discrete=rmse_discrete, bins=rmse_bins)
        if pcc_label is not None:
            cbl_pcc = pcc_label
        if rmse_label is not None:
            cbl_rmse = rmse_label
        sm_pcc = mpl.cm.ScalarMappable(norm=norm_pcc, cmap=cmap_pcc); sm_pcc.set_array(vals_pcc)
        sm_rmse = mpl.cm.ScalarMappable(norm=norm_rmse, cmap=cmap_rmse); sm_rmse.set_array(vals_rmse)

        E, V, R = len(experiments), len(common_vars), len(regions)
        if figsize is None:
            figsize = (2 * max(12.0, 0.42 * V) + 2.0, R * (1.4 + 0.75 * E) + 3.0)
        fig, axes = plt.subplots(R, 2, figsize=figsize, sharex=True, squeeze=False)

        letters = string.ascii_lowercase
        for i, r in enumerate(regions):
            last = i == R - 1
            for j, (mats, sm, name) in enumerate(((mats_pcc[r], sm_pcc, "PCC"),
                                                  (mats_rmse[r], sm_rmse, "Norm. RMSE"))):
                ax = axes[i, j]
                self.render_panel(
                    E=E, V=V, common_vars=common_vars, seasons=seasons,
                    experiments=experiments, explabels=explabels,
                    ax=ax, mats=mats, sm=sm,
                    title_panel=f"({letters[2 * i + j]}) {region_titles[r]}: {name}",
                    add_xlabel=last, ylabel="",
                    draw_grid=draw_grid, x_step=x_step, fontz=fontz, hatch_nan=hatch_nan,
                )
                if not last:
                    ax.tick_params(labelbottom=False)
                if j == 1:
                    ax.tick_params(labelleft=False)

        bottom = 0.30 if V > 28 else (0.26 if V > 22 else 0.20)
        fig.subplots_adjust(left=0.07, right=0.98, top=0.92 if title else 0.96,
                            bottom=bottom, wspace=panel_wspace, hspace=panel_hspace)

        # --- horizontal colorbar under each column; season key between them ---
        cb_h, cb_gap = 0.018, 0.075
        for j, (sm, edges, label) in enumerate(((sm_pcc, edges_pcc, cbl_pcc),
                                                (sm_rmse, edges_rmse, cbl_rmse))):
            pos = axes[-1, j].get_position(fig)
            cax = fig.add_axes([pos.x0 + 0.08 * pos.width, max(0.01, bottom - cb_gap - cb_h - 0.13),
                                0.84 * pos.width, cb_h])
            if edges is not None:
                cb = fig.colorbar(sm, cax=cax, orientation="horizontal", boundaries=edges, extend="both")
                ticks = list(edges) if len(edges) <= 16 else list(edges[::max(1, len(edges) // 12)])
                cb.set_ticks(ticks)
                cb.set_ticklabels([f"{v:g}" for v in ticks])
            else:
                cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend="both")
            cb.ax.tick_params(labelsize=fontz * 0.90, length=3)
            cb.set_label(label, fontsize=fontz, labelpad=4)

        if show_season_box:
            p0, p1 = axes[-1, 0].get_position(fig), axes[-1, 1].get_position(fig)
            bw = 0.045
            bh = bw * figsize[0] / figsize[1]           # square on the page
            bx = 0.5 * (p0.x1 + p1.x0) - bw / 2
            by = max(0.005, bottom - cb_gap - cb_h - 0.13 - bh * 0.35)
            self._draw_season_box(fig.add_axes([bx, by, bw, bh]), fontz=fontz * 0.75)

        if title:
            fig.suptitle(title, fontsize=fontz * 1.05, fontweight="bold", y=0.995)
        out_dir = os.path.dirname(os.fspath(fig_path)) or "."
        os.makedirs(out_dir, exist_ok=True)
        plt.savefig(fig_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
        plt.show(); plt.close()
        print(f"[INFO] Saved figure: {fig_path}")

    # ---- from the error_drifting variant ----
    def _setup_colorbar_and_season_box(
        self,
        fig,
        ax,
        sm,
        edges=None,
        cbar_label=None,
        fontz: int = 10,
        *,
        show_season_box: bool = True,
        bar_width_frac: float = 0.016,   # ~1.6% of figure width
        bar_height_frac: float = None,   # None => auto from figure height
        bar_x_gap_frac: float = 0.012,   # gap between ax and colorbar (fig fraction)
        bar_y_bias: float = 0.52,        # 0..1: vertical anchor; >0.5 lowers the bar slightly
        season_box_w_mult: float = 2.0,  # box width = w_mult * bar width
        season_box_h_mult: float = 0.40, # box height = h_mult * bar height
        season_box_gap_frac: float = 0.010,  # gap above bar (figure fraction)
    ):
        """
        Build a colorbar whose size & position are in *figure* coordinates (so it
        won't be stretched by layout engines), and optionally add the season box
        above it—also in figure coords. Returns (cbar, cax, box_ax).
        """
        fig.canvas.draw_idle()
        ax_pos = ax.get_position(fig)  # Bbox in figure coords (0..1)

        # ---- Choose bar size in figure coords ----
        # width (unchanged)
        bw = float(np.clip(bar_width_frac, 0.008, 0.04))

        # height: auto if None, otherwise clamp   <-- CHANGE THIS BLOCK
        if bar_height_frac is None:
            # was: min(0.30, ax_pos.height * 1.0)
            bh = min(0.65, ax_pos.height * 1.15)   # allow a bit longer auto bar
        else:
            bh = float(bar_height_frac)

        # was: bh = float(np.clip(bh, 0.12, 0.40))
        bh = float(np.clip(bh, 0.12, 0.60))         # raise the max to 0.60

        # position: to the right of the main axes, centered vertically then biased down a bit
        x0 = ax_pos.x1 + float(bar_x_gap_frac)
        y_center = ax_pos.y0 + ax_pos.height * float(bar_y_bias)  # bias > 0.5 means slightly lower
        y0 = float(np.clip(y_center - bh / 2.0, 0.02, 0.98 - bh))

        # ---- Create colorbar axis in figure coords (decoupled from ax height) ----
        cax = fig.add_axes([x0, y0, bw, bh])

        # ---- Draw colorbar ----
        if edges is not None:
            cbar = fig.colorbar(sm, cax=cax, boundaries=edges, extend="both")
            if len(edges) > 16:
                step = max(1, len(edges)//12)
                idx = list(range(0, len(edges), step))
                if idx[-1] != len(edges)-1:
                    idx[-1] = len(edges)-1
                cbar.set_ticks([edges[k] for k in idx])
            else:
                cbar.set_ticks(edges)
            cbar.set_ticklabels([f"{v:g}" for v in cbar.get_ticks()])
        else:
            cbar = fig.colorbar(sm, cax=cax, extend="both")

        cbar.ax.tick_params(labelsize=fontz * 0.90, length=3)
        if cbar_label:
            cbar.set_label(cbar_label, fontsize=fontz, labelpad=5)

        # ---- Season box above the bar (also in figure coords) ----
        box_ax = None
        if show_season_box:
            gap = float(season_box_gap_frac)
            bx = x0 - 0.10 * bw                            # tiny left shift for aesthetics
            by = float(np.clip(y0 + bh + gap, 0.0, 0.98))
            bw_box = float(bw * season_box_w_mult)
            bh_box = float(bh * season_box_h_mult)

            box_ax = fig.add_axes([bx, by, bw_box, bh_box])
            box_ax.set_xlim(0, 1); box_ax.set_ylim(0, 1)
            box_ax.set_aspect("equal"); box_ax.axis("off")

            # Frame
            box_ax.add_patch(
                Rectangle((0, 0), 1, 1, facecolor="white",
                          edgecolor="0.25", linewidth=1.1)
            )
            # Triangles
            tri_defs = {
                "DJF": [(0.0, 1.0), (1.0, 1.0), (0.5, 0.5)],
                "MAM": [(1.0, 0.0), (1.0, 1.0), (0.5, 0.5)],
                "JJA": [(0.0, 0.0), (1.0, 0.0), (0.5, 0.5)],
                "SON": [(0.0, 0.0), (0.0, 1.0), (0.5, 0.5)],
            }
            for s in ["DJF", "MAM", "JJA", "SON"]:
                box_ax.add_patch(
                    Polygon(tri_defs[s], closed=True,
                            facecolor="white", edgecolor="0.4", lw=1.0)
                )
            label_pos = {"DJF": (0.50, 0.80), "MAM": (0.78, 0.50),
                         "JJA": (0.50, 0.20), "SON": (0.22, 0.50)}
            for s, (tx, ty) in label_pos.items():
                box_ax.text(tx, ty, s, ha="center", va="center",
                            fontsize=fontz * 0.90, weight="normal", color="black")

        return cbar, cax, box_ax

    # ---- from the error_drifting variant ----
    def plot_heatmap_quadseason(
        self,
        analyzers_by_season: dict,      # {"DJF": an_djf, "MAM": an_mam, "JJA": an_jja, "SON": an_son}
        metric: str,
        fig_path,
        mode: str | None = "auto",
        fontz: int = 10,
        title: str = "",
        xlabel: str = "Variable",
        ylabel: str = "Time Period",
        cblabel: str = None,
        x_step: int = 1,
        y_step: int = 1,
        separator_year: int | None = None,
        separator_color: str = "black",
        separator_style: str = "--",
        separator_width: float = 1.0,
        discrete: bool = True,
        bins=None,
        vmin=None,
        vmax=None,
        figsize: tuple[float, float] | None = None,
        transpose: bool = False,          # time on x (True) vs y (False)
        show_season_box: bool = True,     # mini season legend (square)
    ):
        """
        Quad-season heatmap. Each cell is split into triangles:
          DJF=top, MAM=right, JJA=bottom, SON=left.
        If transpose=False (default): variables on x, time on y.
        If transpose=True:            time on x, variables on y.
        """

        if isinstance(fig_path, (tuple, list)):
            fig_path = fig_path[0]
        fig_path = os.fspath(fig_path)

        season_order = ["DJF", "MAM", "JJA", "SON"]
        missing = [s for s in season_order if s not in analyzers_by_season]
        if missing:
            raise ValueError(f"Missing analyzers for seasons: {', '.join(missing)}")

        field = self._select_field(metric, mode=mode if isinstance(mode, str) else "auto")

        # --- Collect dataframes per season ---
        dfs = {}
        for s in season_order:
            an_s = analyzers_by_season[s]
            if an_s.ds is None or field not in an_s.ds:
                raise RuntimeError(f"[{s}] No dataset or field '{field}' not present.")
            df = an_s.ds[field].to_pandas().replace([999.999, 9.99999e2, -999, -999.0], np.nan)
            dfs[s] = df

        # Intersect rows/columns across seasons (rows=periods, cols=variables)
        common_rows, common_cols = None, None
        for df in dfs.values():
            common_rows = df.index if common_rows is None else common_rows.intersection(df.index)
            common_cols = df.columns if common_cols is None else common_cols.intersection(df.columns)
        if len(common_rows) == 0 or len(common_cols) == 0:
            raise RuntimeError("No overlapping periods/variables across seasons.")
        for s in season_order:
            dfs[s] = dfs[s].loc[common_rows, common_cols]

        # Values for scaling
        arr_all = np.concatenate([np.ravel(dfs[s].to_numpy(dtype=float)) for s in season_order])
        arr_all = arr_all[np.isfinite(arr_all)]

        # Colormap/norm (edges for discrete)
        cmap, norm, cbar_label, edges = self._make_colormap_and_norm(
            metric=metric, field=field, arr_all=arr_all,
            discrete=discrete, bins=bins, vmin=vmin, vmax=vmax,
        )
        if cblabel is not None:
            cbar_label = cblabel

        # === minimal tweak: feed real data to ScalarMappable to avoid stray "bad" color in cbar
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array(arr_all if arr_all.size else np.array([0.0]))

        # --- Figure size based on orientation ---
        n_time = len(common_rows)  # periods
        n_vars = len(common_cols)  # variables
        # === minimal tweak: remove pre-use of fig_w/fig_h
        if figsize is None:
            figsize = (max(14, 0.35 * (n_time if transpose else n_vars)),
                       8 if transpose else 10)
        fig_w, fig_h = figsize

        fig, ax = plt.subplots(figsize=figsize, constrained_layout=False)

        # --- Triangle geometry (rotate mapping if transposed) ---
        if not transpose:
            # x=j (variable), y=i (time)
            tri_map = {
                "DJF": lambda j,i: [(j, i), (j+1, i), (j+0.5, i+0.5)],
                "MAM": lambda j,i: [(j+1, i), (j+1, i+1), (j+0.5, i+0.5)],
                "JJA": lambda j,i: [(j, i+1), (j+1, i+1), (j+0.5, i+0.5)],
                "SON": lambda j,i: [(j, i), (j, i+1), (j+0.5, i+0.5)],
            }
            for i, row in enumerate(common_rows):      # time
                for j, col in enumerate(common_cols):  # variable
                    for s in season_order:
                        val = dfs[s].loc[row, col]
                        poly = Polygon(
                            tri_map[s](j, i), closed=True, facecolor=sm.to_rgba(val),
                            edgecolor="black", linewidth=0.01
                        )
                        ax.add_patch(poly)

            # Axis setup (diagnostic layout)
            ax.set_xlim(0, n_vars)
            ax.set_ylim(n_time, 0)
            ax.set_xticks(np.arange(n_vars) + 0.5)
            ax.set_xticklabels([str(c) for c in common_cols], rotation=90, fontsize=fontz*0.95)
            ax.set_yticks(np.arange(n_time) + 0.5)
            ax.set_yticklabels([str(r) for r in common_rows], fontsize=fontz*0.95)
            ax.set_xlabel(xlabel or "", fontsize=fontz)
            ax.set_ylabel(ylabel or "", fontsize=fontz)

        else:
            # x=i (time), y=j (variable)
            tri_map = {
                "DJF": lambda i,j: [(i, j), (i, j+1), (i+0.5, j+0.5)],
                "MAM": lambda i,j: [(i, j+1), (i+1, j+1), (i+0.5, j+0.5)],
                "JJA": lambda i,j: [(i+1, j), (i+1, j+1), (i+0.5, j+0.5)],
                "SON": lambda i,j: [(i, j), (i+1, j), (i+0.5, j+0.5)],
            }
            for i, row in enumerate(common_rows):      # time
                for j, col in enumerate(common_cols):  # variable
                    for s in season_order:
                        val = dfs[s].loc[row, col]
                        poly = Polygon(
                            tri_map[s](i, j), closed=True, facecolor=sm.to_rgba(val),
                            edgecolor='grey', linewidth=0.01
                        )
                        ax.add_patch(poly)

            # Axis setup (publication layout)
            ax.set_xlim(0, n_time)
            ax.set_ylim(n_vars, 0)
            ax.set_xticks(np.arange(n_time) + 0.5)
            ax.set_xticklabels([str(r) for r in common_rows], rotation=45, ha="right", fontsize=fontz*0.95)
            ax.set_yticks(np.arange(n_vars) + 0.5)
            ax.set_yticklabels([str(c) for c in common_cols], rotation=0, fontsize=fontz*0.95)
            ax.set_xlabel(xlabel or "", fontsize=fontz)
            ax.set_ylabel(ylabel or "", fontsize=fontz)

        if x_step > 1:
            ticks = ax.get_xticks()
            labels = [lbl.get_text() for lbl in ax.get_xticklabels()]
            ax.set_xticks(ticks[::x_step])
            ax.set_xticklabels(labels[::x_step], fontsize=fontz*0.95)

        if y_step > 1:
            ticks = ax.get_yticks()
            labels = [lbl.get_text() for lbl in ax.get_yticklabels()]
            ax.set_yticks(ticks[::y_step])
            ax.set_yticklabels(labels[::y_step], fontsize=fontz*0.95)

        ax.set_title(title, fontsize=fontz*1.1, loc="left")
        
        cbar, cax, box_ax = self._setup_colorbar_and_season_box(
            fig, ax, sm, edges=edges,
            cbar_label=(cblabel or cbar_label),
            fontz=fontz,
            show_season_box=show_season_box,
            bar_width_frac=0.016,      # thicker (was 0.016)
            bar_height_frac=3.0,      # longer (was 0.95)
            bar_x_gap_frac=-0.02,      # same small gap
            bar_y_bias=0.45,
            season_box_w_mult=4.5,
            season_box_h_mult=0.5,
            season_box_gap_frac=-0.05,
        )
        
        # ===================== Spinup vs piControl separator =====================
        if separator_year is None and hasattr(self, "picontrol_switch_year"):
            separator_year = self.picontrol_switch_year

        if separator_year is not None:
            an0 = analyzers_by_season[season_order[0]]
            full_periods = list(an0.ds["period"].values)
            full_starts  = list(an0.ds["period_start"].values)
            sep_index = None
            for r_i, r_label in enumerate(common_rows):
                try:
                    full_idx = full_periods.index(r_label)
                    if full_starts[full_idx] > separator_year:
                        sep_index = r_i
                        break
                except ValueError:
                    continue
            if sep_index is not None:
                if transpose:
                    ax.axvline(sep_index, color=separator_color,
                               linestyle=separator_style, linewidth=separator_width)
                else:
                    ax.axhline(sep_index, color=separator_color,
                               linestyle=separator_style, linewidth=separator_width)

        fig.subplots_adjust(left=0.12, right=0.85, bottom=0.10, top=0.95)

        # Save
        out_dir = os.path.dirname(fig_path) or "."
        os.makedirs(out_dir, exist_ok=True)        
        plt.savefig(fig_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
        plt.show()
        plt.close()
        print(f"[INFO] Saved figure: {fig_path}")
