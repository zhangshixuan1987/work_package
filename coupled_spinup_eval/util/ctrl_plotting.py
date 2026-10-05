"""Figures for the ctrl analyses (one experiment: spin-up followed by its piControl).

  MPASTimeDiagnosticsPlotter (shared base):
    OceanDriftPlotter  (14_ocn_ts_ctrl_analysis)
    OHCDriftPlotter    (15_ohc_ts_ctrl_analysis)
  AtmosphereBalancePlotter (12_atm_flux_ctrl_analysis)
  CoupledBudgetPlotter     (13_ocn_conserve_ctrl_analysis)
"""
import os, re, glob, warnings, string
import numpy as np
import xarray as xr
import cftime
import contextlib
import pandas as pd
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, FixedLocator
import os
import glob
import string
import warnings


class MPASTimeDiagnosticsPlotter:
    """Shared trend/decay plotting for annual tidy time series (from ctrl_processing builders). Use a subclass."""

    def __init__(self, data_var, var_name=None, var_label=None, var_unit=None):
        """
        Parameters
        ----------
        data_var : str
            The actual NetCDF variable name in the tidy file (e.g., MPAS SST var).
        var_label : str, optional
            Pretty label for titles/axis (e.g., "Sea Surface Temperature").
        var_unit : str, optional
            Unit string for y-axis (e.g., "°C" or "Sv").
        """
        self.data_var  = data_var
        self.var_name  = var_name 
        self.var_label = var_label or data_var
        self.var_unit  = var_unit

    def apply_centered_filter(self, y, window=11):
        w = max(1, int(window))
        if w % 2 == 0: w += 1
        pad = w // 2
        y_arr = np.asarray(y, dtype=float)
        isfin = np.isfinite(y_arr).astype(float)
        y_filled = np.where(np.isfinite(y_arr), y_arr, 0.0)
        k = np.ones(w, dtype=float)
        num = np.convolve(y_filled, k, mode="same")
        den = np.convolve(isfin,   k, mode="same")
        out = np.full_like(y_arr, np.nan, dtype=float)
        g = den > 0
        out[g] = num[g] / den[g]
        out[:pad] = np.nan; out[-pad:] = np.nan
        return out

    def sliding_window_trend(self, y, x, window):
        n = len(x)
        w = max(2, int(window))
        if w % 2 == 0: w += 1
        half = w // 2
        out = np.full(n, np.nan, dtype=float)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        for i in range(n):
            i0 = max(0, i - half); i1 = min(n, i + half + 1)
            xx = x[i0:i1]; yy = y[i0:i1]
            m = np.isfinite(xx) & np.isfinite(yy)
            if m.sum() >= 3:
                out[i] = np.polyfit(xx[m], yy[m], 1)[0]  # per year
        return out

    def _exp_decay_fit(self, X, Y, fit_range=None, use_filtered_for_fit=True,
                       default_window=400, tau_min=5, tau_max=5000, n_tau=200,
                       min_points=20):
        """
        Robust single-exponential fit using grid search over tau.
        Returns dict with tau, A, y_inf, fit_range, and yhat over all X.
        """
        X = np.asarray(X, float); Y = np.asarray(Y, float)
        if fit_range is None:
            x0 = float(np.nanmin(X))
            x1 = x0 + default_window
        else:
            x0, x1 = map(float, fit_range)

        mfit = np.isfinite(X) & np.isfinite(Y) & (X >= x0) & (X <= x1)
        if mfit.sum() < min_points:
            return None

        # de-noise for fitting, but keep original Y for plotting
        Yfit = Y.copy()
        if use_filtered_for_fit:
            Yfit = self.apply_centered_filter(Yfit, window=11)
        mfit &= np.isfinite(Yfit)
        if mfit.sum() < min_points:
            return None

        t = X[mfit] - x0  # start fit at t=0 to improve conditioning
        y = Yfit[mfit]

        # candidate taus (log-spaced)
        taus = np.exp(np.linspace(np.log(tau_min), np.log(tau_max), n_tau))
        best = {"rmse": np.inf, "tau": np.nan, "A": np.nan, "y_inf": np.nan}

        for tau in taus:
            e = np.exp(-t / tau)
            # Solve [1, e] * [y_inf, A] = y  → least squares
            G = np.column_stack([np.ones_like(e), e])
            coef, *_ = np.linalg.lstsq(G, y, rcond=None)
            y_inf, A = float(coef[0]), float(coef[1])
            yhat = y_inf + A * e
            rmse = float(np.sqrt(np.nanmean((y - yhat) ** 2)))
            if rmse < best["rmse"]:
                best.update(rmse=rmse, tau=float(tau), A=A, y_inf=y_inf)

        # build prediction over all X (use original axis)
        e_all = np.exp(-(X - x0) / best["tau"])
        yhat_all = best["y_inf"] + best["A"] * e_all
        best["yhat"] = yhat_all
        best["fit_range"] = (x0, x1)
        return best

    def summarize_trends_by_window(self, x, trend_vals, box_window=50):
        x = np.asarray(x, dtype=float)
        if len(x) == 0: return []
        start = int(np.nanmin(x)); end = int(np.nanmax(x))
        bw = max(1, int(box_window))
        edges = list(range(start, end + 1, bw))
        if edges[-1] < end: edges.append(end)
        return [{"start_year": edges[i], "end_year": edges[i+1]} for i in range(len(edges)-1)]

    def _find_plot_csvs(self, data_out_dir, fig_base):
        ts_csv  = os.path.join(data_out_dir, f"{fig_base}_timeseries.csv")
        box_csv = os.path.join(data_out_dir, f"{fig_base}_trend_boxes.csv")
        ok = os.path.exists(ts_csv) and os.path.exists(box_csv)
        return ts_csv, box_csv, ok

    def _prepare_from_csv(self, ts_csv, box_window, window):
        df = pd.read_csv(ts_csv)
        x = df["year"].to_numpy()
        y = df["value"].to_numpy()
        want = f"sliding_trend_{int(window)}yr"
        if want in df.columns:
            trend = df[want].to_numpy()
        else:
            trend = self.sliding_window_trend(y, x, window) * 10.0
        stats = self.summarize_trends_by_window(x, trend, box_window=box_window)
        for s in stats:
            s["values"] = trend[(x >= s["start_year"]) & (x <= s["end_year"]) & np.isfinite(trend)]
        return {"x": x, "y": y, "trend": trend, "stats": stats}

    def _save_csvs(self, trend_dict, data_out_dir, fig_base, show_filter, filter_window, window):
        os.makedirs(data_out_dir, exist_ok=True)
        x = trend_dict["x"]; y = trend_dict["y"]; trend = trend_dict["trend"]
        out = {"year": x.astype(int, copy=False), "value": y}
        if show_filter:
            out[f"filtered_{int(filter_window)}yr"] = self.apply_centered_filter(y, window=filter_window)
        out[f"sliding_trend_{int(window)}yr"] = trend
        pd.DataFrame(out).to_csv(os.path.join(data_out_dir, f"{fig_base}_timeseries.csv"), index=False)

        rows = []
        for s in trend_dict["stats"]:
            vals = np.asarray(s["values"], dtype=float)
            if vals.size == 0 or not np.isfinite(vals).any():
                rows.append({"start_year": int(s["start_year"]), "end_year": int(s["end_year"]),
                             "mean": np.nan, "p10": np.nan, "p90": np.nan, "count": 0})
            else:
                rows.append({"start_year": int(s["start_year"]), "end_year": int(s["end_year"]),
                             "mean": float(np.nanmean(vals)),
                             "p10": float(np.nanpercentile(vals, 10)),
                             "p90": float(np.nanpercentile(vals, 90)),
                             "count": int(np.isfinite(vals).sum())})
        pd.DataFrame(rows).to_csv(os.path.join(data_out_dir, f"{fig_base}_trend_boxes.csv"), index=False)

    def plot_from_tidy(
        self,
        tidy_nc,
        out_pdf,
        fig_dir=".",
        experiment="MPAS_spinup",
        window=10,
        box_window=50,
        filter_window=11,
        syear=None,
        eyear=None,
        show_filter=True,
        nxlab=250,
        spinup_markers=None,             # e.g., [2000] or (2000, 2500)
        rgn_label=None,
        show_decay=False,
        decay_fit_range=None,     # e.g., (0, 350) or (0, spinup_markers[0])
        decay_tail_years=200,
        use_filtered_for_fit=True,
        default_window=400,
        tau_min=5,
        tau_max=5000,
        n_tau=200,
        legend_one_row=True,      # both panels legend in one row
        data_out_dir=None,         # optional CSV cache dir
        fontz=16,
        figsize = (30,10),
        nrows = 1,
        ncols = 2,
        fig_id = 0,
        pad_inches=0, 
        dpi=300,
        # --- add these new bounds ---
        ymin=None,
        ymax=None,
        dmin=None,
        dmax=None,
        fig=None,
        axes=None,          # expected shape: (nrows, 2) or list-like
        row_idx=0,          # which row to draw into if axes provided
        do_save=True,       # save only when you want (combined figure: False inside loop)
        do_show=True,       # show only at end for combined
        do_close=True,      # close only at end for combined
        out_path_override=None,  # optional: save combined to a different filename
        show_trend_range_box=True,
        trend_mode = "full_period",
        trend_fit_range=None,        # e.g., (syear, eyear) or (spinup_markers[0], eyear)
        use_filtered_for_trend=True, # use yf if available
        trend_x=None,
        trend_y=None,
    ):
        """
        Reads a *single* tidy NetCDF (combined spin-up + piControl).
        Produces a 2-panel figure: (1) annual + filtered series, (2) sliding trend + box summaries.
        """
        ds = xr.open_dataset(tidy_nc)

        # locate series and years
        if self.data_var not in ds.data_vars:
            ds.close()
            raise KeyError(f"'{self.data_var}' not found in tidy file. Found: {list(ds.data_vars)}")

        years = ds["year"].values if "year" in ds.coords else (
                ds["time"].dt.year.values if "time" in ds.coords else np.arange(ds.dims.get("time", len(ds[self.data_var])))
        )
        y = ds[self.data_var].values

        # sensible x-range defaults
        if syear is None: syear = int(np.nanmin(years))
        if eyear is None: eyear = int(np.nanmax(years))

        # prefer precomputed trend over recompute
        pre_trend = f"sliding_trend_{int(window)}yr"
        pre_trend_dec = f"{pre_trend}_per_decade"
        if pre_trend_dec in ds:
            trend = ds[pre_trend_dec].values
        elif pre_trend in ds:
            trend = ds[pre_trend].values * 10.0
        else:
            trend = self.sliding_window_trend(y, years, window) * 10.0

        # prefer precomputed filtered series if present
        pre_filt = f"filtered_{int(filter_window)}yr"
        yf = ds[pre_filt].values if (show_filter and pre_filt in ds) else (
             self.apply_centered_filter(y, window=filter_window) if show_filter else None
        )

        # CSV cache reuse/save
        tidy_mtime = os.path.getmtime(tidy_nc) if os.path.exists(tidy_nc) else -1
        fig_base = os.path.splitext(os.path.basename(out_pdf))[0]
        if data_out_dir:
            ts_csv, box_csv, ok = self._find_plot_csvs(data_out_dir, fig_base)
            can_reuse = ok and (os.path.getmtime(ts_csv) >= tidy_mtime) and (os.path.getmtime(box_csv) >= tidy_mtime)
        else:
            can_reuse = False

        if can_reuse:
            rebuilt = self._prepare_from_csv(ts_csv, box_window=box_window, window=window)
            X = rebuilt["x"]; Y = rebuilt["y"]; trend_use = rebuilt["trend"]; stats = rebuilt["stats"]
        else:
            X = years; Y = y; trend_use = trend
            stats = self.summarize_trends_by_window(X, trend_use, box_window=box_window)
            for s in stats:
                s["values"] = trend_use[(X >= s["start_year"]) & (X <= s["end_year"]) & np.isfinite(trend_use)]
                
        # ---- plotting
        created_fig = False
        if axes is None:
            fig, axes_local = plt.subplots(1, 2, figsize=figsize)
            ax, ax2 = axes_local[0], axes_local[1]
            created_fig = True
        else:
            # axes expected as (nrows, 2) from plt.subplots(nrows, 2, ...)
            ax  = axes[row_idx, 0]
            ax2 = axes[row_idx, 1]
            
        # ---- plotting        
        color = {"line": "#000000", "weak": "#7f7f7f", "box": "#87ceeb", "mean": "#d62728"}
        panel_id = string.ascii_lowercase[fig_id]
        
        # Left panel: time series
        ax.plot(X, Y, color=color["line"], lw=1.5, alpha=0.35, label="Annual mean")
        
        if show_filter and yf is not None:
            ax.plot(X, yf, color=color["line"], lw=1.6, alpha=1.0, label=f"{filter_window}-yr filtered")
            
        # --- overlay exponential transient fit (panel a1) ---
        if show_decay:
            fit = self._exp_decay_fit(
                X, Y,
                fit_range=decay_fit_range,          # e.g., (syear, spinup_markers[0])
                use_filtered_for_fit=use_filtered_for_fit,          # denoise for stability
                default_window=default_window,                 # if no fit_range
                tau_min=tau_min, tau_max=tau_max, n_tau=n_tau  # search space
            )
            if fit is not None and np.isfinite(fit["tau"]):
                tau = fit["tau"]; yhat = fit["yhat"]; x0, x1 = fit["fit_range"]
                ax.plot(X, yhat, color="#d62728", lw=2.0, alpha=0.9,
                        label=f"Exp fit (τ≈{tau:.0f} yr)")
                #ax.axvspan(x0, x1, color="#d62728", alpha=0.06, lw=0)

        if ymin is None or ymax is None: 
            ymin, ymax = np.nanmin(Y), np.nanmax(Y)
            if np.isfinite(ymin) and np.isfinite(ymax) and ymax > ymin:
                margin = 0.05 * (ymax - ymin)
                ax.set_ylim(ymin - margin, ymax + margin)
        else:
            ax.set_ylim(ymin, ymax)
            
        ax.set_xlim(syear - box_window/5, eyear + box_window/5)
        
        # tick labels & axes
        xticks = np.arange(max(syear, (syear // nxlab) * nxlab), eyear + 1, nxlab)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(int(v)) for v in xticks], fontsize=fontz * 0.95)
        ax.tick_params(axis="y", labelsize=fontz * 0.95)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f" if (np.nanmax(np.abs(Y)) < 10) else "%.1f"))

        unit_lbl = f"{self.var_unit}" if self.var_unit else ""
        ax.set_title(f"({panel_id}1) {self.var_label} ({rgn_label})", fontsize=fontz*1.05, loc="left")
        ax.set_ylabel(f"{unit_lbl}", fontsize=fontz * 0.95)
        ax.set_xlabel("Model Year", fontsize=fontz * 0.95)

        ax.axhline(0, color=color["weak"], ls="--", lw=1.2, alpha=0.7)
        if spinup_markers is not None:
            for s in (spinup_markers if isinstance(spinup_markers, (list, tuple)) else [spinup_markers]):
                ax.axvline(s, color=color["weak"], ls="--", lw=1.2, alpha=0.8)

        # one-line legend
        leg = ax.legend(loc="upper right", fontsize=fontz * 0.95, frameon=True,
                        ncol=3 if legend_one_row else 1, columnspacing=1.2, handlelength=2.0)
        
        # Ensure single row if only one/two entries
        if legend_one_row and len(leg.legend_handles) > 0:
            leg._ncols = len(leg.legend_handles)

        # Right panel: sliding trend + boxes
        ax2.plot(X, trend_use, color=color["line"], lw=1.5, alpha=0.35, label=f"Sliding Trend ({window} yr)")
        mids = [(s["start_year"] + s["end_year"]) // 2 for s in stats]
        box_vals = [s["values"] for s in stats]
        ax2.boxplot(
            box_vals,
            positions=mids,
            widths=max(10, box_window * 0.8),
            patch_artist=True,
            showfliers=False,
            whis=[10, 90],
            boxprops=dict(facecolor=color["box"], edgecolor="black", linewidth=1.0),
            medianprops=dict(color="black", linewidth=1.0),
            whiskerprops=dict(color="black", linewidth=1.0),
            capprops=dict(color="black", linewidth=1.0),
        )
        
        means = [np.nanmean(b) if (len(b) and np.isfinite(b).any()) else np.nan for b in box_vals]
        ax2.scatter(mids, means, color=color["mean"], s=22, zorder=3, label=f"Mean Trend ({box_window} yr)")
        
        # ---- summary box: long-term trend OR range of window-mean trends ----
        if show_trend_range_box:

            # Default text position
            if trend_x is None or trend_y is None:
                trend_x, trend_y = 0.5, 0.12   # looks nicer near bottom center usually

            # Choose mode
            if trend_mode == "range":
                marr = np.asarray(means[2:], dtype=float)
                marr = marr[np.isfinite(marr)]
                if marr.size > 0:
                    mmin = np.nanmin(marr)
                    mmax = np.nanmax(marr)
                    txt = (f"Mean trend: "
                           f"[{mmin:.2f},{mmax:.2f}] {unit_lbl} decade$^{{-1}}$")
                else:
                    txt = "Mean trend: n/a"

            else:
                # ---- LONG-TERM trend (single value) ----
                # pick series for fitting
                y_fit = yf if (use_filtered_for_trend and show_filter and (yf is not None)) else Y

                # choose fit range
                if trend_fit_range is None:
                    x0, x1 = float(syear), float(eyear)
                else:
                    x0, x1 = map(float, trend_fit_range)

                m = np.isfinite(X) & np.isfinite(y_fit) & (X >= x0) & (X <= x1)
                if m.sum() >= 3:
                    slope_per_year = np.polyfit(X[m], y_fit[m], 1)[0]
                    slope_per_dec  = 10.0 * slope_per_year
                    txt = (f"Long-term Trend ({int(x0)}–{int(x1)}yr): "
                           f"{slope_per_dec:.5f} {unit_lbl} decade$^{{-1}}$")
                else:
                    txt = "Long-term trend: n/a"

            ax2.text(
                trend_x, trend_y, txt,
                transform=ax2.transAxes,
                ha="center", va="center",
                fontsize=fontz * 0.9,
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.35",
                    facecolor="white",
                    edgecolor="red",
                    linewidth=1.5,
                    alpha=0.95,
                ),
                zorder=10,
            )


        if dmin is None or dmax is None:
            dmin, dmax = np.nanmin(trend_use), np.nanmax(trend_use)
            if np.isfinite(dmin) and np.isfinite(dmax) and dmax > dmin:
                dmargin = 0.05 * (dmax - dmin)
                ax2.set_ylim(dmin - dmargin, dmax + dmargin)
        else:
            ax2.set_ylim(dmin, dmax)
            
        ax2.set_xlim(syear - box_window/5, eyear + box_window/5)
        ax2.set_xticks(xticks)
        ax2.set_xticklabels([str(int(v)) for v in xticks], fontsize=fontz * 0.95)
        ax2.tick_params(axis="y", labelsize=fontz * 0.95)
        ax2.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

        ax2.set_title(f"({panel_id}2) {self.var_label} Trend ({rgn_label})", fontsize=fontz*1.05, loc="left")
        ax2.set_ylabel(fr"{unit_lbl} decade$^{-1}$", fontsize=fontz * 0.95)
        ax2.set_xlabel("Model Year", fontsize=fontz * 0.95)

        ax2.axhline(0, color=color["weak"], ls="--", lw=1.2, alpha=0.7)
        if spinup_markers is not None:
            for s in (spinup_markers if isinstance(spinup_markers, (list, tuple)) else [spinup_markers]):
                ax2.axvline(s, color=color["weak"], ls="--", lw=1.2, alpha=0.7)

        # one-line legend
        h, l = ax2.get_legend_handles_labels()
        ax2.legend(h, l, loc="upper right", fontsize=fontz * 0.95, frameon=True,
                   ncol=len(l) if legend_one_row else 1, columnspacing=1.2, handlelength=2.0)

        fig.subplots_adjust(
            left=0.045, right=0.995, bottom=0.17, top=0.91,
            wspace=0.15, hspace=0.35
        )
        # save tightly to out_path
        out_path = out_path_override if out_path_override else (
            os.path.join(fig_dir, out_pdf) if os.path.isdir(fig_dir) else out_pdf
        )

        if do_save:
            fig.savefig(out_path, bbox_inches="tight", pad_inches=pad_inches, dpi=dpi)
            print(f"[PLOT] Saved: {out_path}")

        if do_show:
            plt.show()

        # Only close if we created this fig, or if caller explicitly requests
        if do_close and (created_fig or axes is not None):
            plt.close(fig)
            
        # --- Hide x-axis labels for combined multi-row figures (except bottom row) ---
        if axes is not None:
            if row_idx < (nrows - 1):
                # left panel
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=False)
        
                # right panel
                ax2.set_xlabel("")
                ax2.tick_params(axis="x", labelbottom=False)
        
        # save CSVs if we didn’t reuse cache
        if data_out_dir and not can_reuse:
            self._save_csvs({"x": X, "y": Y, "trend": trend_use, "stats": stats},
                            data_out_dir, fig_base, show_filter, filter_window, window)

        ds.close()
        return out_path


class OceanDriftPlotter(MPASTimeDiagnosticsPlotter):
    """Ocean / sea-ice drift panels: annual series, filtered, sliding trends, e-folding fit."""


class OHCDriftPlotter(MPASTimeDiagnosticsPlotter):
    """OHC drift panels with implied heat-flux (W m-2) right axis."""

    def plot_from_tidy(
        self,
        tidy_nc,
        out_pdf,
        fig_dir=".",
        experiment="MPAS_spinup",
        window=10,
        box_window=50,
        filter_window=11,
        syear=None,
        eyear=None,
        show_filter=True,
        nxlab=250,
        spinup_markers=None,
        rgn_label=None,
        show_decay=False,
        decay_fit_range=None,
        decay_tail_years=200,
        use_filtered_for_fit=True,
        default_window=400,
        tau_min=5,
        tau_max=5000,
        n_tau=200,
        legend_one_row=True,
        data_out_dir=None,
        fontz=16,
        figsize=(30,10),
        nrows=1,
        ncols=2,
        fig_id=0,
        pad_inches=0,
        dpi=300,
        # --- add these new bounds ---
        ymin=None,
        ymax=None,
        dmin=None,
        dmax=None,
        # --- NEW: equivalent flux slope lines (W/m^2) ---
        flux_lines=None,            # e.g. [-0.2,-0.1,-0.05,0.05,0.1,0.2]
        show_flux_axis=True,
        flux_area="ocean",          # "ocean" or "earth"
        flux_area_m2=None,          # optional override
        flux_anchor_year=None,      # t0 (default: syear)
        flux_anchor_value=0.0,      # ΔOHC0 at t0
        flux_ref_year=None,         # tref (default: eyear)
        flux_axis_color="#1f77b4",  # match your example blue
        flux_axis_fmt="%.3f",
        flux_style=":",
        flux_lw=2,
        flux_alpha=1.0,
        sec_per_year = 365.0 * 24.0 * 3600.0, 
        # --- NEW: allow external figure/axes for multi-panel layout ---
        fig=None,
        axes=None,          # expected shape: (nrows, 2) or list-like
        row_idx=0,          # which row to draw into if axes provided
        do_save=True,       # save only when you want (combined figure: False inside loop)
        do_show=True,       # show only at end for combined
        do_close=True,      # close only at end for combined
        out_path_override=None,  # optional: save combined to a different filename
        show_trend_range_box=True,
        trend_mode = "full_period",
        trend_fit_range=None,        # e.g., (syear, eyear) or (spinup_markers[0], eyear)
        use_filtered_for_trend=True, # use yf if available
        trend_x=None,
        trend_y=None,
    ):

        """
        Reads a *single* tidy NetCDF (combined spin-up + piControl).
        Produces a 2-panel figure: (1) annual + filtered series (+ optional exp fit),
        (2) sliding trend + box summaries.
        """
        
        ds = xr.open_dataset(tidy_nc)
        # locate series and years
        if self.data_var not in ds.data_vars:
            ds.close()
            raise KeyError(f"'{self.data_var}' not found in tidy file. Found: {list(ds.data_vars)}")

        years = ds["year"].values if "year" in ds.coords else (
            ds["time"].dt.year.values if "time" in ds.coords
            else np.arange(ds.dims.get("time", len(ds[self.data_var])))
        )
        y = ds[self.data_var].values

        # sensible x-range defaults
        if syear is None: syear = int(np.nanmin(years))
        if eyear is None: eyear = int(np.nanmax(years))

        # precomputed trend (per decade) if present
        pre_trend = f"sliding_trend_{int(window)}yr"
        pre_trend_dec = f"{pre_trend}_per_decade"
        if pre_trend_dec in ds:
            trend = ds[pre_trend_dec].values
        elif pre_trend in ds:
            trend = ds[pre_trend].values * 10.0
        else:
            trend = self.sliding_window_trend(y, years, window) * 10.0

        # precomputed filtered series if present
        pre_filt = f"filtered_{int(filter_window)}yr"
        yf = ds[pre_filt].values if (show_filter and pre_filt in ds) else (
             self.apply_centered_filter(y, window=filter_window) if show_filter else None
        )

        # Units fallback from tidy file if not provided
        unit_lbl = self.var_unit
        if not unit_lbl:
            try:
                unit_lbl = str(ds[self.data_var].attrs.get("units", "") or "")
            except Exception:
                unit_lbl = ""

        # CSV cache reuse/save
        tidy_mtime = os.path.getmtime(tidy_nc) if os.path.exists(tidy_nc) else -1
        fig_base = os.path.splitext(os.path.basename(out_pdf))[0]
        if data_out_dir:
            ts_csv, box_csv, ok = self._find_plot_csvs(data_out_dir, fig_base)
            can_reuse = ok and (os.path.getmtime(ts_csv) >= tidy_mtime) and (os.path.getmtime(box_csv) >= tidy_mtime)
        else:
            can_reuse = False

        if can_reuse:
            rebuilt = self._prepare_from_csv(ts_csv, box_window=box_window, window=window)
            X = rebuilt["x"]; Y = rebuilt["y"]; trend_use = rebuilt["trend"]; stats = rebuilt["stats"]
        else:
            X = years; Y = y; trend_use = trend
            stats = self.summarize_trends_by_window(X, trend_use, box_window=box_window)
            for s in stats:
                s["values"] = trend_use[(X >= s["start_year"]) & (X <= s["end_year"]) & np.isfinite(trend_use)]

        # ---- plotting
        created_fig = False
        if axes is None:
            fig, axes_local = plt.subplots(1, 2, figsize=figsize)
            ax, ax2 = axes_local[0], axes_local[1]
            created_fig = True
        else:
            # axes expected as (nrows, 2) from plt.subplots(nrows, 2, ...)
            ax  = axes[row_idx, 0]
            ax2 = axes[row_idx, 1]

        color = {"line": "#000000", "weak": "#7f7f7f", "box": "#87ceeb", "mean": "#d62728"}
        panel_id = string.ascii_lowercase[fig_id]
        
        # Left panel: time series
        ax.plot(X, Y, color=color["line"], lw=1.5, alpha=0.35, label="Annual mean")
        
        if ymin is None or ymax is None: 
            ymin, ymax = np.nanmin(Y), np.nanmax(Y)
            if np.isfinite(ymin) and np.isfinite(ymax) and ymax > ymin:
                margin = 0.05 * (ymax - ymin)
                ax.set_ylim(ymin - margin, ymax + margin)
        else:
            ax.set_ylim(ymin, ymax)
            
        if show_filter and yf is not None:
            ax.plot(X, yf, color=color["line"], lw=1.6, alpha=1.0, label=f"{filter_window}-yr filtered")
            
        # --- Equivalent flux right axis + constant-flux slope lines (W/m^2) ---
        if show_flux_axis:
            # area
            if flux_area_m2 is not None:
                A = float(flux_area_m2)
            else:
                fa = str(flux_area).lower()
                A = 3.61e14 if fa.startswith("ocean") else 5.10e14
        
            t0   = float(syear if flux_anchor_year is None else flux_anchor_year)
            y0   = float(flux_anchor_value)
            tref = float(eyear if flux_ref_year is None else flux_ref_year)
            dt   = max(1.0, tref - t0)  # years (constant scaling, matches your old logic)
        
            # right axis: flux that would produce (y - y0) over dt years
            def ohc_to_flux(y_left):
                return ((np.asarray(y_left) - y0) * 1.0e22) / (A * sec_per_year * dt)
        
            def flux_to_ohc(F):
                return (np.asarray(F) * A * sec_per_year * dt) / 1.0e22 + y0
        
            axr = ax.secondary_yaxis('right', functions=(ohc_to_flux, flux_to_ohc))
            axr.set_ylabel(r"W/m$^2$", fontsize=fontz * 0.95, color=flux_axis_color)
            axr.tick_params(axis="y", labelsize=fontz * 0.95, colors=flux_axis_color)
            axr.spines["right"].set_color(flux_axis_color)
            axr.yaxis.set_major_formatter(FormatStrFormatter(flux_axis_fmt))
        
            # slope lines: use right-axis ticks (W/m^2) as guides
            fig = ax.figure
            fig.canvas.draw()  # ensure ticks exist
            Fticks = [t for t in axr.get_yticks() if np.isfinite(t)]
        
            for F in Fticks:
                # slope in (1e22 J)/yr for constant flux F
                slope = (float(F) * A * sec_per_year) / 1.0e22
                yline = y0 + slope * (X - t0)
                ax.plot(X, yline, color=flux_axis_color, ls=flux_style, lw=flux_lw,
                        alpha=flux_alpha, zorder=20)


        # Optional exponential transient fit
        if show_decay:
            fit = self._exp_decay_fit(
                X, Y,
                fit_range=decay_fit_range,
                use_filtered_for_fit=use_filtered_for_fit,
                default_window=default_window,
                tau_min=tau_min, tau_max=tau_max, n_tau=n_tau
            )
            if fit is not None and np.isfinite(fit["tau"]):
                tau = fit["tau"]; yhat = fit["yhat"]; x0, x1 = fit["fit_range"]
                ax.plot(X, yhat, color="#d62728", lw=2.0, alpha=0.9,
                        label=f"Exp fit (τ≈{tau:.0f} yr)")
                # Optionally highlight used fit window:
                # ax.axvspan(x0, x1, color="#d62728", alpha=0.06, lw=0)
            
        ax.set_xlim(syear - box_window/5, eyear + box_window/5)

        # ticks & labels
        xticks = np.arange(max(syear, (syear // nxlab) * nxlab), eyear + 1, nxlab)
        ax.set_xticks(xticks)
        ax.set_xticklabels([str(int(v)) for v in xticks], fontsize=fontz * 0.95)
        ax.tick_params(axis="y", labelsize=fontz * 0.95)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f" if (np.nanmax(np.abs(Y)) < 10) else "%.1f"))

        ax.set_title(f"({panel_id}1) {self.var_label} ({rgn_label})", fontsize=fontz*1.05, loc="left")
        ax.set_ylabel(f"{unit_lbl}", fontsize=fontz * 0.95)
        ax.set_xlabel("Model Year", fontsize=fontz * 0.95)

        ax.axhline(0, color=color["weak"], ls="--", lw=1.2, alpha=0.7)
        if spinup_markers is not None:
            for s in (spinup_markers if isinstance(spinup_markers, (list, tuple)) else [spinup_markers]):
                ax.axvline(s, color=color["weak"], ls="--", lw=1.2, alpha=0.8)

        # one-line legend
        leg = ax.legend(loc="lower left", fontsize=fontz * 0.95, frameon=True,
                        ncol=3 if legend_one_row else 1, columnspacing=1.2, handlelength=2.0)
        if legend_one_row and len(leg.legend_handles) > 0:
            leg._ncols = len(leg.legend_handles)

        # Right panel: sliding trend + boxes
        ax2.plot(X, trend_use, color=color["line"], lw=1.5, alpha=0.35, label=f"Sliding Trend ({window} yr)")
        mids = [(s["start_year"] + s["end_year"]) // 2 for s in stats]
        box_vals = [s["values"] for s in stats]
        ax2.boxplot(
            box_vals,
            positions=mids,
            widths=max(10, box_window * 0.8),
            patch_artist=True,
            showfliers=False,
            whis=[10, 90],
            boxprops=dict(facecolor=color["box"], edgecolor="black", linewidth=1.0),
            medianprops=dict(color="black", linewidth=1.0),
            whiskerprops=dict(color="black", linewidth=1.0),
            capprops=dict(color="black", linewidth=1.0),
        )
        means = [np.nanmean(b) if (len(b) and np.isfinite(b).any()) else np.nan for b in box_vals]
        ax2.scatter(mids, means, color=color["mean"], s=22, zorder=3, label=f"Mean Trend ({box_window} yr)")
        
        # ---- summary box: long-term trend OR range of window-mean trends ----
        if show_trend_range_box:
            
            # Default text position
            if trend_x is None or trend_y is None:
                trend_x, trend_y = 0.5, 0.12   # looks nicer near bottom center usually

            # Choose mode
            if trend_mode == "range":
                marr = np.asarray(means[2:], dtype=float)
                marr = marr[np.isfinite(marr)]
                if marr.size > 0:
                    mmin = np.nanmin(marr)
                    mmax = np.nanmax(marr)
                    txt = (f"Mean trend: "
                           f"[{mmin:.2f},{mmax:.2f}]x{unit_lbl} decade$^{{-1}}$")
                else:
                    txt = "Mean trend: n/a"

            else:
                # ---- LONG-TERM trend (single value) ----
                # pick series for fitting
                y_fit = yf if (use_filtered_for_trend and show_filter and (yf is not None)) else Y

                # choose fit range
                if trend_fit_range is None:
                    x0, x1 = float(syear), float(eyear)
                else:
                    x0, x1 = map(float, trend_fit_range)

                m = np.isfinite(X) & np.isfinite(y_fit) & (X >= x0) & (X <= x1)
                if m.sum() >= 3:
                    slope_per_year = np.polyfit(X[m], y_fit[m], 1)[0]
                    slope_per_dec  = 10.0 * slope_per_year
                    txt = (f"Long-term Trend ({int(x0)}–{int(x1)}yr): "
                           f"{slope_per_dec:.5f}x{unit_lbl} decade$^{{-1}}$")
                else:
                    txt = "Long-term trend: n/a"

            ax2.text(
                trend_x, trend_y, txt,
                transform=ax2.transAxes,
                ha="center", va="center",
                fontsize=fontz * 0.9,
                color="black",
                bbox=dict(
                    boxstyle="round,pad=0.35",
                    facecolor="white",
                    edgecolor="red",
                    linewidth=1.5,
                    alpha=0.95,
                ),
                zorder=10,
            )

        if dmin is None or dmax is None:
            dmin, dmax = np.nanmin(trend_use), np.nanmax(trend_use)
            if np.isfinite(dmin) and np.isfinite(dmax) and dmax > dmin:
                dmargin = 0.05 * (dmax - dmin)
                ax2.set_ylim(dmin - dmargin, dmax + dmargin)
        else:
            ax2.set_ylim(dmin, dmax)
    
        ax2.set_xlim(syear - box_window/5, eyear + box_window/5)
        ax2.set_xticks(xticks)
        ax2.set_xticklabels([str(int(v)) for v in xticks], fontsize=fontz * 0.95)
        ax2.tick_params(axis="y", labelsize=fontz * 0.95)
        ax2.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

        ax2.set_title(f"({panel_id}2) {self.var_label} Trend ({rgn_label})", fontsize=fontz*1.05, loc="left")
        ax2.set_ylabel(fr"{unit_lbl} decade$^{-1}$", fontsize=fontz * 0.95)
        ax2.set_xlabel("Model Year", fontsize=fontz * 0.95)

        ax2.axhline(0, color=color["weak"], ls="--", lw=1.2, alpha=0.7)
        if spinup_markers is not None:
            for s in (spinup_markers if isinstance(spinup_markers, (list, tuple)) else [spinup_markers]):
                ax2.axvline(s, color=color["weak"], ls="--", lw=1.2, alpha=0.7)

        # one-line legend
        h, l = ax2.get_legend_handles_labels()
        ax2.legend(h, l, loc="upper right", fontsize=fontz * 0.95, frameon=True,
                   ncol=len(l) if legend_one_row else 1, columnspacing=1.2, handlelength=2.0)

        fig.subplots_adjust(left=0.047, right=0.995, bottom=0.17, top=0.91, wspace=0.30, hspace=0.35)

        # save tightly to out_path
        out_path = out_path_override if out_path_override else (
            os.path.join(fig_dir, out_pdf) if os.path.isdir(fig_dir) else out_pdf
        )

        if do_save:
            fig.savefig(out_path, bbox_inches="tight", pad_inches=pad_inches, dpi=dpi)
            print(f"[PLOT] Saved: {out_path}")

        if do_show:
            plt.show()

        # Only close if we created this fig, or if caller explicitly requests
        if do_close and (created_fig or axes is not None):
            plt.close(fig)
            
        # --- Hide x-axis labels for combined multi-row figures (except bottom row) ---
        if axes is not None:
            if row_idx < (nrows - 1):
                # left panel
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=False)
        
                # right panel
                ax2.set_xlabel("")
                ax2.tick_params(axis="x", labelbottom=False)
        
        # save CSVs if we didn’t reuse cache
        if data_out_dir and not can_reuse:
            self._save_csvs({"x": X, "y": Y, "trend": trend_use, "stats": stats},
                            data_out_dir, fig_base, show_filter, filter_window, window)
            
        ds.close()
        return out_path


class AtmosphereBalancePlotter:
    def __init__(
        self, 
        combined_tidy, 
        region_map=None,
        var_dict=None,
    ):
        self.combined_tidy = combined_tidy
        self.region_map = self._get_region(region_map)
        self.var_dict = self._get_var(var_dict)
        #make sure unit in certain quantities are consistent with plot target
        self.default_units = {
            "TAS": "degC", 
            "TREFHT": "degC", 
            "TS": "degC",
            "MASSFLUX": "mm yr-1",
        }
        
    def _get_region(self, region_map):
        if region_map is None:
            region_map = {0: "Global", 1: "Northern Hemisphere", 2: "Southern Hemisphere"}
        return region_map
        
    def _get_var(self, var_dict):
        if  var_dict is None:
            var_dict = {
                "net_toa_flux_restom":             {"name": "RESTOM",   'label': "RESTOM",  'unit': 'W m$^{-2}$'},
                "net_sfc_flux_ressurf":            {"name": "RESSURF",  'label': "RESSURF", 'unit': 'W m$^{-2}$'},
                "global_surface_skin_temperature": {"name": "TS",       'label': "TS",      'unit': '$^{o}$C'},
                "global_surface_air_temperature":  {"name": "TAS",      'label': "TAS",     'unit': '$^{o}$C'},
                "toa_radiation":                   {"name": "RESTOA",   'label': "RESTOA",  'unit': 'W m$^{-2}$'},
                "net_atm_energy_imbalance":        {"name": "RESNET",   'label': "RESNET",  'unit': 'W m$^{-2}$'},
                "net_atm_water_imbalance":         {"name": "MASSFLUX", 'label': "E - P",   'unit': 'mm yr$^{-1}$'},
            }
        return var_dict
        
    def _normalize_unit_label(self, u):
        if u is None:
            return None
        u = u.strip().lower().replace("/", " ").replace("−", "-")
        aliases = {
            "kg m^-2 s^-1": "kg m-2 s-1",
            "kg/m2/s": "kg m-2 s-1",
            "kelvin": "k",
            "°c": "degc",
            "c": "degc",
            "mm/yr": "mm yr-1",
            "mm per year": "mm yr-1",
        }
        return aliases.get(u, u)

    def _check_and_warn_units(self, varname, units_map, y=None, where=None):
        """
        Check-only: compare current units (from units_map) to expected defaults.
        Emits warnings; does NOT convert data or mutate units_map.
        """
        expected = self._normalize_unit_label(self.default_units.get(varname))
        current  = self._normalize_unit_label((units_map or {}).get(varname))

        loc = f" [{where}]" if where else ""

        if current is None:
            warnings.warn(
                f"Unit missing for '{varname}'{loc}. Expected '{expected}'. "
                f"Data will be used as-is; set df.attrs['units']['{varname}'] to silence.",
                stacklevel=2
            )
            # Optional heuristic for temps
            if varname in {"TAS","TREFHT","TS"} and y is not None and np.nanmin(y) > 100:
                warnings.warn(
                    f"'{varname}'{loc} values look like Kelvin (min>100) but unit is unknown. "
                    f"Expected '{expected}'.",
                    stacklevel=2
                )
            return

        if expected and current != expected:
            warnings.warn(
                f"Unit mismatch for '{varname}'{loc}: got '{current}', expected '{expected}'. "
                f"Data will be used as-is. If already converted, update df.attrs['units']['{varname}']='{expected}'.",
                stacklevel=2
            )

        # Light sanity check for temperatures only
        if varname in {"TAS","TREFHT","TS"} and y is not None:
            ymin = np.nanmin(y)
            if current == "degc" and ymin > 100:
                warnings.warn(
                    f"'{varname}'{loc} labeled 'degC' but values (min={ymin:.2f}) look like Kelvin.",
                    stacklevel=2
                )
            if current in {"k"} and ymin < 100:
                warnings.warn(
                    f"'{varname}'{loc} labeled 'K' but values (min={ymin:.2f}) look like Celsius.",
                    stacklevel=2
                )
                
    def sliding_window_trend(self, y_vals, x_vals=None, win=10):
        y_vals = np.asarray(y_vals)
        trends = np.full_like(y_vals, np.nan, dtype=np.float64)
        x_vals = np.arange(len(y_vals)) if x_vals is None else np.asarray(x_vals)
        half = win // 2
        for i in range(half, len(y_vals) - half):
            y_window = y_vals[i - half:i + half]
            x_window = x_vals[i - half:i + half]
            if np.any(np.isnan(y_window)) or np.any(np.isnan(x_window)):
                continue
            slope, *_ = linregress(x_window, y_window)
            trends[i] = slope
        return trends

    def apply_centered_filter(self, y, window=11):
        y = np.asarray(y)
        filtered = np.full_like(y, np.nan, dtype=np.float64)
        half = window // 2
        for i in range(half, len(y) - half):
            window_vals = y[i - half:i + half + 1]
            if np.all(np.isfinite(window_vals)):
                filtered[i] = np.mean(window_vals)
        return filtered

    def summarize_trends_by_window(self, years, trend_vals, box_window=50, verbose=True):
        trend_years = np.array(years)
        valid_mask = np.isfinite(trend_vals)
        trend_years = trend_years[valid_mask]
        trend_vals = trend_vals[valid_mask]
        stats = []
        if len(trend_years) == 0:
            return stats
        start_year = int(trend_years[0]) // box_window * box_window
        end_year = int(trend_years[-1]) + 1
        bin_edges = np.arange(start_year, end_year + box_window, box_window)
        for i in range(len(bin_edges) - 1):
            yr_start, yr_end = bin_edges[i], bin_edges[i + 1]
            mask = (trend_years >= yr_start) & (trend_years < yr_end)
            bin_data = trend_vals[mask]
            if bin_data.size > 0:
                stats_dict = {
                    "start_year": yr_start,
                    "end_year": yr_end - 1,
                    "mean": float(np.mean(bin_data)),
                    "std": float(np.std(bin_data)),
                    "median": float(np.median(bin_data)),
                    "min": float(np.min(bin_data)),
                    "max": float(np.max(bin_data)),
                    "n": int(bin_data.size),
                }
                stats.append(stats_dict)
                if verbose:
                    print(f"[{yr_start}–{yr_end - 1}] n={stats_dict['n']}, mean={stats_dict['mean']:.3e}")
        return stats
        
    def _prepare_trend_data(self, df, varnames, region_map, window, box_window):
        units_map = {**self.default_units, **df.attrs.get("units", {})}

        trend_all = {}
        for varname in varnames:
            trend_all[varname] = {}
            for rgn, name in region_map.items():
                dfr = df[df["rgn"] == rgn].set_index("year")
                if varname not in dfr:
                    continue

                x = dfr.index.values
                y = dfr[varname].values

                # Check only; warns if something is off, does NOT change y
                self._check_and_warn_units(varname, units_map, y=y, where=f"region={name}")

                trend_vals = self.sliding_window_trend(y, x, window) * 10
                stats = self.summarize_trends_by_window(
                    x, trend_vals, box_window=box_window, verbose=False
                )
                for s in stats:
                    s["values"] = trend_vals[(x >= s["start_year"]) & (x <= s["end_year"]) & np.isfinite(trend_vals)]

                trend_all[varname][name] = {"x": x, "y": y, "trend": trend_vals, "stats": stats}
        return trend_all

    def fit_best_exponential(self, X, Y, fit_range=None, use_filtered_for_fit=True,
                       default_window=400, tau_min=5, tau_max=5000, n_tau=200,
                       min_points=20):
        """
        Robust single-exponential fit using grid search over tau.
        Returns dict with tau, A, y_inf, fit_range, and yhat over all X.
        """
        X = np.asarray(X, float); Y = np.asarray(Y, float)
        if fit_range is None:
            x0 = float(np.nanmin(X))
            x1 = x0 + default_window
        else:
            x0, x1 = map(float, fit_range)

        mfit = np.isfinite(X) & np.isfinite(Y) & (X >= x0) & (X <= x1)
        if mfit.sum() < min_points:
            return None

        # de-noise for fitting, but keep original Y for plotting
        Yfit = Y.copy()
        if use_filtered_for_fit:
            Yfit = self.apply_centered_filter(Yfit, window=11)
        mfit &= np.isfinite(Yfit)
        if mfit.sum() < min_points:
            return None

        t = X[mfit] - x0  # start fit at t=0 to improve conditioning
        y = Yfit[mfit]

        # candidate taus (log-spaced)
        taus = np.exp(np.linspace(np.log(tau_min), np.log(tau_max), n_tau))
        best = {"rmse": np.inf, "tau": np.nan, "A": np.nan, "y_inf": np.nan}

        for tau in taus:
            e = np.exp(-t / tau)
            # Solve [1, e] * [y_inf, A] = y  → least squares
            G = np.column_stack([np.ones_like(e), e])
            coef, *_ = np.linalg.lstsq(G, y, rcond=None)
            y_inf, A = float(coef[0]), float(coef[1])
            yhat = y_inf + A * e
            rmse = float(np.sqrt(np.nanmean((y - yhat) ** 2)))
            if rmse < best["rmse"]:
                best.update(rmse=rmse, tau=float(tau), A=A, y_inf=y_inf)

        # build prediction over all X (use original axis)
        e_all = np.exp(-(X - x0) / best["tau"])
        yhat_all = best["y_inf"] + best["A"] * e_all
        best["yhat"] = yhat_all
        best["fit_range"] = (x0, x1)
        return best

    def plot_time_series_panels(
            self, varnames, fig_name,
            text_position_dict=None, 
            region_labels=None,
            window=10, box_window=50, 
            syear=0, eyear=2500,
            spinup_period=[220, 2000], 
            filter_window=11, 
            init_window=0,
            show_mean_box=True,
            vars_for_mean_box=None, 
            box_x=0.02,
            box_y=0.10, # top-left in axes coords
            show_filter=True, 
            show_decay=True,
            decay_fit_range=None, # e.g., (0, 350) or (0, spinup_markers[0])
            decay_tail_years=200,
            use_filtered_for_fit=True,
            show_trend_range_box=True,
            trend_region="Global",
            trend_mode="full_period",
            trend_fit_range=None, # e.g., (syear, eyear) or (spinup_markers[0], eyear)
            use_filtered_for_trend=True, # use yf if available
            plot_dict=None, 
            fig_id=0,
            figsize=(30, 10),
            nxlab=2, 
            fontz=14,
            dpi=600, 
            bbox_inches="tight", 
            pad_inches=0.05,
            wspace=0.15, 
            hspace=0.35
        ):
        
        ds = xr.open_dataset(self.combined_tidy)
        df = ds.to_dataframe().reset_index()

        region_map = region_labels or self.region_map
        trend_all = self._prepare_trend_data(
            df, varnames, region_map, window, box_window
        )

        color_dict = {
            'black': '#000000', 'blue': '#1f77b4', 'red': '#d62728',
            'skyblue': '#87ceeb', 'gray': '#7f7f7f'
        }

        style_dict = {
            'raw': ['-', 1.5, 0.3, 1],
            'filter': ['-', 1.5, 1.0, 2],
            'mean': ['o', 40, 0, 4],
            'decay': ['--', 1.5, 1.0, 1],
        }

        nrows, ncols = len(varnames), 2
        fig1, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1:
            axes = np.expand_dims(axes, axis=0)

        xticks = np.arange(max(syear, (syear // nxlab) * nxlab), eyear + 1, nxlab)
        xlabels = [str(y) for y in xticks]

        for row, varname in enumerate(varnames):
            panel_id = string.ascii_lowercase[row]
            label_text = varname 
            var_text = varname
            if plot_dict and varname in plot_dict:
                var_text = f"{plot_dict[varname]['label']}"
                label_text = f"{plot_dict[varname]['unit']}"
                show_decay = plot_dict[varname]["decay"]
                
            # === Column 1: Time series (NH, SH, Global) ===
            ax = axes[row, 0]
            decay_labels = []
            ymin, ymax = np.inf, -np.inf
            for rgn, color, alpha in zip(["Northern Hemisphere", "Southern Hemisphere"],
                                         [color_dict["blue"], color_dict["red"]],
                                         [0.4, 0.4]):
                if rgn not in trend_all[varname]:
                    continue
                x = trend_all[varname][rgn]["x"]
                y = trend_all[varname][rgn]["y"]
            
                ymin = min(ymin, np.nanmin(y))
                ymax = max(ymax, np.nanmax(y))
                ax.plot(x, y, label=rgn, color=color, linestyle=style_dict['raw'][0],
                        linewidth=style_dict['raw'][1], alpha=style_dict['raw'][2],
                        zorder=style_dict['raw'][3])

                if show_filter:
                    yf = self.apply_centered_filter(y, window=filter_window)
                    ax.plot(x, yf, color=color, linestyle=style_dict['filter'][0],
                            linewidth=style_dict['filter'][1], alpha=style_dict['filter'][2],
                            zorder=style_dict['filter'][3], label=f"{rgn} ({filter_window}yr MA)")

                if show_decay:
                    x_valid = x[init_window:]
                    y_valid = y[init_window:]
                    fit_result = self.fit_best_exponential(
                            x_valid, 
                            y_valid, 
                            fit_range=decay_fit_range,          
                            use_filtered_for_fit=use_filtered_for_fit
                    )
                    tau = fit_result["tau"]; yhat = fit_result["yhat"]; x0, x1 = fit_result["fit_range"]
                    decay_labels.append(f"τ≈{tau:.1f}yr" if not np.isnan(tau) else "N/A")
                    ax.plot(x_valid, yhat, color=color, linestyle=style_dict['decay'][0],
                            linewidth=style_dict['decay'][1], alpha=style_dict['decay'][2],
                            zorder=style_dict['decay'][3])

            # === Global ===
            if "Global" in trend_all[varname]:
                x = trend_all[varname]["Global"]["x"]
                y = trend_all[varname]["Global"]["y"]
                global_x = x
                global_y = y
                global_yf = None
                ymin = min(ymin, np.nanmin(y))
                ymax = max(ymax, np.nanmax(y))
                ax.plot(x, y, label="Global", color=color_dict["black"], alpha=0.2, linewidth=2.2, zorder=10)

                if show_filter:
                    yf = self.apply_centered_filter(y, window=filter_window)
                    global_yf = yf
                    ax.plot(x, yf, color=color_dict["black"], linestyle=style_dict['filter'][0],
                            linewidth=style_dict['filter'][1], alpha=style_dict['filter'][2],
                            zorder=style_dict['filter'][3], label=f"Global ({filter_window}yr MA)")

                if show_decay:
                    x_valid = x[init_window:]
                    y_valid = y[init_window:]
                    fit_result = self.fit_best_exponential(
                            x_valid, 
                            y_valid, 
                            fit_range=decay_fit_range,          
                            use_filtered_for_fit=use_filtered_for_fit
                    )
                    tau = fit_result["tau"]; yhat = fit_result["yhat"]; x0, x1 = fit_result["fit_range"]
                    decay_labels.append(f"τ≈{tau:.1f}yr" if not np.isnan(tau) else "N/A")
                    tau = fit_result["tau"]; yhat = fit_result["yhat"]; x0, x1 = fit_result["fit_range"]
                    ax.plot(x_valid, yhat, color=color_dict["black"], linestyle=style_dict['decay'][0],
                            linewidth=style_dict['decay'][1], alpha=style_dict['decay'][2],
                            zorder=style_dict['decay'][3])

            if plot_dict and varname in plot_dict:
                ymin = plot_dict[varname]["min"]
                ymax = plot_dict[varname]["max"]
                ax.set_ylim(ymin, ymax)
                ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            else:
                # === After all plotting ===
                if np.isfinite(ymin) and np.isfinite(ymax):
                    y_margin = 0.05 * (ymax - ymin)
                    ax.set_ylim(ymin - y_margin, ymax + y_margin)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

            ax.set_xlim(syear - box_window/5, eyear + box_window/5)
            ax.set_xticks(xticks)
            ax.set_xticklabels(xlabels, fontsize=fontz * 0.95)
            
            ax.set_title(f"({panel_id}1) {var_text} (mean)", fontsize=fontz*1.05, loc='left')
            ax.set_ylabel(label_text, fontsize=fontz * 0.95)
            ax.set_xlabel("Model Year", fontsize=fontz * 0.95)
            
            ax.tick_params(axis='y', labelsize=fontz * 0.95)
            ax.minorticks_on()
            ax.xaxis.set_minor_locator(FixedLocator([(xticks[i] + xticks[i+1]) / 2 for i in range(len(xticks)-1)]))
            ax.axhline(0, color=color_dict["gray"], linestyle='--', linewidth=1.5, alpha=0.7)
            for period in spinup_period:
                ax.axvline(period, color=color_dict["gray"], linestyle='--', linewidth=1.5, alpha=0.8)

            region_legend = [
                Line2D([0], [0], color=color_dict["black"], lw=2.2, alpha=1.0, label='Global'),
                Line2D([0], [0], color=color_dict["blue"], lw=1.5, alpha=1.0, label='NH'),
                Line2D([0], [0], color=color_dict["red"], lw=1.5, alpha=1.0, label='SH')
            ]
            legend1 = ax.legend(handles=region_legend, loc="upper center", ncol=3, fontsize=fontz * 0.95,
                frameon=True, bbox_to_anchor=(0.5, 1.02) if varname == "TAS" else None)
            legend1.get_frame().set_alpha(0.6)
            ax.add_artist(legend1)

            if show_decay and decay_labels:
                decay_colors = [color_dict["black"], color_dict["blue"], color_dict["red"]]
                exp_line = [
                    Line2D([0], [0], 
                           color=decay_colors[i], 
                           linestyle=style_dict['decay'][0],
                           linewidth=style_dict['decay'][1], 
                           alpha=style_dict['decay'][2], 
                           label=lbl)
                    for i, lbl in enumerate(decay_labels)
                ]
                ax.legend(handles=exp_line, loc="lower center", 
                          ncol=3, fontsize=fontz * 0.95,
                          frameon=True, bbox_to_anchor=(0.5, 0.02)
                         )

            if show_mean_box:
                # --- Mean ± std over long-term window for mean time-series (Global, NH, SH) ---
                if varname in vars_for_mean_box:
                    # averaging window: use trend_fit_range if provided
                    if trend_fit_range is None:
                        x0, x1 = float(syear), float(eyear)
                    else:
                        x0, x1 = map(float, trend_fit_range)

                    region_order = ["Global", "Northern Hemisphere", "Southern Hemisphere"]
                    region_color_names = {
                        "Global": "black",
                        "Northern Hemisphere": "blue",
                        "Southern Hemisphere": "red",
                    }

                    # list of (color_name, mean_val, std_val)
                    mean_entries = []

                    for rname in region_order:
                        if rname not in trend_all[varname]:
                            continue

                        if rname == "Global" and "Global" in trend_all[varname]:
                            # reuse already-defined global series
                            x_mean = global_x
                            if show_filter and (global_yf is not None):
                                y_mean_series = global_yf
                            else:
                                y_mean_series = global_y
                        else:
                            # NH / SH: get directly from trend_all and optionally filter
                            x_mean = trend_all[varname][rname]["x"]
                            y_raw  = trend_all[varname][rname]["y"]
                            if show_filter:
                                y_mean_series = self.apply_centered_filter(
                                    y_raw, window=filter_window
                                )
                            else:
                                y_mean_series = y_raw

                        m_mean = (
                            np.isfinite(x_mean)
                            & np.isfinite(y_mean_series)
                            & (x_mean >= x0)
                            & (x_mean <= x1)
                        )

                        if m_mean.sum() > 0:
                            mean_val = float(np.nanmean(y_mean_series[m_mean]))
                            std_val  = float(np.nanstd(y_mean_series[m_mean]))
                            cname = region_color_names.get(rname, "black")
                            mean_entries.append((cname, mean_val, std_val))

                    if mean_entries:
                        # base unit for display
                        if plot_dict and varname in plot_dict:
                            base_unit = plot_dict[varname]["unit"]
                        else:
                            base_unit = label_text

                        # --- draw one boxed label (black) ---
                        base_label = f"Mean ± Std ({int(x0)}–{int(x1)}yr):"
                        ax.text(
                            box_x, box_y, base_label,
                            transform=ax.transAxes,
                            ha="left", va="top",
                            fontsize=fontz * 0.9,
                            color="black",
                            bbox=dict(
                                boxstyle="round,pad=0.35",
                                facecolor="white",
                                edgecolor="black",
                                linewidth=1.5,
                                alpha=0.95,
                            ),
                            zorder=10,
                        )

                        # --- colored mean±std values + unit ---
                        x_start = box_x + 0.46  # tweak spacing if needed
                        dx      = 0.14          # horizontal spacing between regions

                        for i, (cname, mean_val, std_val) in enumerate(mean_entries):
                            text_val = f"{mean_val:.2f}±{std_val:.2f}"
                            ax.text(
                                x_start + i * dx, box_y,
                                text_val,
                                transform=ax.transAxes,
                                ha="left", va="top",
                                fontsize=fontz * 0.9,
                                color=cname,
                                zorder=11,
                            )

                        # unit at the end in black
                        ax.text(
                            x_start + len(mean_entries) * dx + 0.02, box_y,
                            base_unit,
                            transform=ax.transAxes,
                            ha="left", va="top",
                            fontsize=fontz * 0.9,
                            color="black",
                            zorder=11,
                        )
                        
            # === Column 2: Global trend boxplot + overlay ===
            ax2 = axes[row, 1]
            if trend_region not in trend_all[varname]:
                ax2.set_visible(False)
                continue
            trend_labels = []
            trend_x = trend_all[varname][trend_region]["x"]
            trend_y = trend_all[varname][trend_region]["trend"]
            trend_labels.append(f"Sliding Trend ({window}yr)")
            stats = trend_all[varname][trend_region]["stats"]
            mid_years = [(s["start_year"] + s["end_year"]) // 2 for s in stats]
            box_data = [s["values"] for s in stats]
            ax2.boxplot(
                box_data,
                positions=mid_years,
                widths=40,
                patch_artist=True,
                showfliers=False,
                whis=[10, 90],
                boxprops=dict(facecolor="skyblue", edgecolor="black", linewidth=1.0, alpha=1.0),
                medianprops=dict(color="black", linewidth=1.0),
                whiskerprops=dict(color="black", linewidth=1.0),
                capprops=dict(color="black", linewidth=1.0)
            )

            ax2.plot(
                trend_x, trend_y,
                color=color_dict["black"],
                linestyle='-',
                linewidth=1.5,
                alpha=0.3,
                zorder=1
            )
            
            if show_decay:
                x_valid = trend_x[init_window:]
                y_valid = trend_y[init_window:]
                fit_result = self.fit_best_exponential(
                        x_valid, 
                        y_valid, 
                        fit_range=decay_fit_range,          
                        use_filtered_for_fit=use_filtered_for_fit
                )
                tau = fit_result["tau"]; yhat = fit_result["yhat"]; x0, x1 = fit_result["fit_range"]
                    
            # Overlay mean trend per box
            means = [np.nanmean(trends) for trends in box_data]
            trend_labels.append(f"Mean Trend ({box_window}yr)")
            ax2.scatter(
                mid_years, means,
                color=color_dict['red'],
                marker=style_dict['mean'][0],
                s=style_dict['mean'][1],
                linewidth=style_dict['mean'][2],
                zorder=style_dict['mean'][3],
            )
            
            # Format axis
            if plot_dict and varname in plot_dict:
                tmin = plot_dict[varname]["dmin"]
                tmax = plot_dict[varname]["dmax"]
                ax2.set_ylim(tmin, tmax)
                ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            else:
                tmin, tmax = np.nanmin(trend_y), np.nanmax(trend_y)
                if np.isfinite(tmin) and np.isfinite(tmax):
                    margin = 0.05 * (tmax - tmin)
                    ax2.set_ylim(tmin - margin, tmax + margin)
                    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    
            ax2.set_xlim(syear - box_window/5, eyear + box_window/5)
            ax2.set_xticks(xticks)
            ax2.set_xticklabels([str(y) for y in xticks], fontsize=fontz * 0.95)
            ax2.set_title(f"({panel_id}2) {var_text} Trend (Global)", fontsize=fontz*1.05, loc='left')
            ax2.set_ylabel(f"{label_text} decade$^{{-1}}$", fontsize=fontz * 0.95)
            ax2.set_xlabel("Model Year", fontsize=fontz * 0.95)
            ax2.tick_params(axis='y', labelsize=fontz * 0.95)
            ax2.axhline(0, color=color_dict["gray"], linestyle="--", linewidth=1.5, alpha=0.7)
            for period in spinup_period:
                ax2.axvline(period, color=color_dict["gray"], linestyle="--", linewidth=1.5, alpha=0.7)

            # Legend
            legend_handles = [
                Line2D([0], [0], color=color_dict["black"], linestyle=style_dict['filter'][0],
                       linewidth=style_dict['filter'][1], alpha=style_dict['filter'][2], label=trend_labels[0]),
                Line2D([0], [0], marker=style_dict['mean'][0], color=color_dict["red"],
                       markersize=style_dict['mean'][1] * 0.3, linewidth=0, label=trend_labels[-1])
            ]
            
            ax2.legend(
                handles=legend_handles,
                loc="upper center",
                fontsize=fontz * 0.95,
                frameon=True,
                labelspacing=0.2,
                columnspacing=0.5,
                handletextpad=0.2,
                ncol=3
            )
            
            # === Filter trends between year 350 and 2000 ===
            if show_trend_range_box:
                # Default text position
                text_x, text_y = 0.02, 0.98
                if text_position_dict is not None: 
                    pos = text_position_dict.get(varname, {'x': 0.02, 'y': 0.98})
                    text_x, text_y = pos['x'], pos['y']
                                    
                if plot_dict and varname in plot_dict:
                    trend_unit = plot_dict[varname]["unit"] + " decade$^{{-1}}$"
                else:
                    trend_unit = f"{label_text} decade$^{{-1}}$"
                    
                # Choose mode
                if trend_mode == "range":
                    valid_inds = [(y >= 350) and (y <= 2000) for y in mid_years]
                    filtered_means = [np.nanmean(box_data[i]) for i, valid in enumerate(valid_inds) if valid and np.isfinite(box_data[i]).any()]
                    if filtered_means:
                        mean_min = np.nanmin(filtered_means)
                        mean_max = np.nanmax(filtered_means)
                        txt = f"Mean trend: [{mean_min:.2f},{mean_max:.2f}] {trend_unit}"
                    else:
                        txt = f"Mean trend: n/a"
                else:
                    # ---- LONG-TERM trend (single value) ----
                    # pick series for fitting
                    x = global_x
                    y = global_y
                    yf = global_yf
                    y_fit = yf if (use_filtered_for_trend and show_filter and (yf is not None)) else y
                    # choose fit range
                    if trend_fit_range is None:
                        x0, x1 = float(syear), float(eyear)
                    else:
                        x0, x1 = map(float, trend_fit_range)
                        
                    m = np.isfinite(x) & np.isfinite(y_fit) & (x >= x0) & (x <= x1)
                    if m.sum() >= 3:
                        slope_per_year = np.polyfit(x[m], y_fit[m], 1)[0]
                        slope_per_dec  = 10.0 * slope_per_year
                        txt = (f"Long-term Trend ({int(x0)}–{int(x1)}yr): "
                               f"{slope_per_dec:.5f} {trend_unit}")
                    else:
                        txt = "Long-term trend: n/a"
                            
                ax2.text(
                    text_x, text_y, txt,
                    transform=ax2.transAxes,
                    ha="center", va="center",
                    fontsize=fontz * 0.9,
                    color="black",
                    bbox=dict(
                        boxstyle="round,pad=0.35",
                        facecolor="white",
                        edgecolor="red",
                        linewidth=1.5,
                        alpha=0.95,
                    ),
                    zorder=10,
                )
                
        fig1.subplots_adjust(wspace=wspace, hspace=hspace)
        fig1.savefig(fig_name, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches)
        print(f"[PLOT] Saved: {fig_name}")
        plt.show()

        return 

    def plot_hemisphere_panels(
            self, varnames, fig_name,        
            text_position_dict=None,
            region_labels=None,
            window=10, box_window=50, 
            syear=0, eyear=2500,
            spinup_period=[220, 2000], 
            filter_window=11, 
            init_window=0,
            show_filter=False, 
            show_decay=False,
            decay_fit_range=None,    
            decay_tail_years=200,
            use_filtered_for_fit=True,
            show_trend_range_box=True,
            trend_region=None,
            trend_mode="full_period",
            trend_fit_range=None,         
            use_filtered_for_trend=None,  
            plot_dict=None, 
            fig_id=0,
            figsize=(30, 10),
            nxlab=2, 
            fontz=14,
            dpi=600, 
            bbox_inches="tight", 
            pad_inches=0.05,
            wspace=0.15, 
            hspace=0.35
        ):
        
        ds = xr.open_dataset(self.combined_tidy)
        df = ds.to_dataframe().reset_index()

        region_map = region_labels or self.region_map
        trend_all = self._prepare_trend_data(
            df, varnames, region_map, window, box_window
        )
        if trend_region is None:
            region_order = ['Northern Hemisphere', 'Southern Hemisphere']
        elif isinstance(trend_region, str):
            region_order = [trend_region]
        else:
            region_order = list(trend_region)
            
        color_dict = {
            'black': '#000000', 'blue': '#1f77b4', 'red': '#d62728',
            'skyblue': '#87ceeb', 'gray': '#7f7f7f'
        }
        
        style_dict = {
            'filter': ['-', 1.5, 0.3, 2],
            'mean': ['o', 40, 0, 4],
            'decay': ['--', 1.5, 1.0, 2],
        }
        # === NH & SH Trend Panel ===
        nrows, ncols = len(varnames), 2
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows == 1:
            axes = np.expand_dims(axes, axis=0)
        
        for row, varname in enumerate(varnames):
            panel_id = string.ascii_lowercase[row]
            if plot_dict and varname in plot_dict:
                var_text = f"{plot_dict[varname]['label']}"
                label_text = f"{plot_dict[varname]['unit']}"
                show_decay = plot_dict[varname]["decay"]

            for col, rname in enumerate(region_order):
                ax = axes[row, col]
                if rname not in trend_all[varname]:
                    ax.set_visible(False)
                    continue
                    
                trend_labels = [f"Sliding Trend ({window}yr)"]
                trend_x = trend_all[varname][rname]["x"]
                trend_y = trend_all[varname][rname]["trend"]  
                ax.plot(
                    trend_x, trend_y,
                    color=color_dict["black"],
                    linestyle=style_dict['filter'][0],
                    linewidth=style_dict['filter'][1],
                    alpha=style_dict['filter'][2],
                    zorder=style_dict['filter'][3],
                    label=trend_labels[0]
                )

                if show_decay:
                    x_valid = trend_x[init_window:]
                    y_valid = trend_y[init_window:]
                    fit_result = self.fit_best_exponential(
                            x_valid, 
                            y_valid, 
                            fit_range=decay_fit_range,          
                            use_filtered_for_fit=use_filtered_for_fit
                    )
                    tau = fit_result["tau"]; yhat = fit_result["yhat"]; x0, x1 = fit_result["fit_range"]
                    
                # Boxplot data
                stats = trend_all[varname][rname]["stats"]
                mid_years = [(s["start_year"] + s["end_year"]) // 2 for s in stats]
                box_data = [s["values"] for s in stats]
                ax.boxplot(
                    box_data,
                    positions=mid_years,
                    widths=40,
                    patch_artist=True,
                    showfliers=False,
                    whis=[10, 90],
                    boxprops=dict(facecolor="skyblue", edgecolor="black", linewidth=1.0, alpha=1.0),
                    medianprops=dict(color="black", linewidth=1.0),
                    whiskerprops=dict(color="black", linewidth=1.0),
                    capprops=dict(color="black", linewidth=1.0)
                )

                # Overlay mean trend per box
                means = [np.nanmean(trends) for trends in box_data]
                trend_labels.append(f"Mean Trend ({box_window}yr)")
                ax.scatter(
                    mid_years, means,
                    color=color_dict['red'],
                    marker=style_dict['mean'][0],
                    s=style_dict['mean'][1],
                    linewidth=style_dict['mean'][2],
                    zorder=style_dict['mean'][3],
                    label=trend_labels[-1]
                )
                
                # Format axis
                if plot_dict and varname in plot_dict:
                    tmin = plot_dict[varname]["dmin"]
                    tmax = plot_dict[varname]["dmax"]
                    ax.set_ylim(tmin, tmax)
                    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                else:
                    tmin, tmax = np.nanmin(trend_y), np.nanmax(trend_y)
                    if np.isfinite(tmin) and np.isfinite(tmax):
                        margin = 0.05 * (tmax - tmin)
                        ax.set_ylim(tmin - margin, tmax + margin)
                        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                        
                ax.set_xlim(syear - box_window, eyear + box_window)
                xticks = np.arange(max(syear, (syear // nxlab) * nxlab), eyear + 1, nxlab)
                ax.set_xticks(xticks)
                ax.set_xticklabels([str(y) for y in xticks], fontsize=fontz * 0.95)
                ax.set_xlabel("Model Year", fontsize=fontz * 0.95)
                ax.set_ylabel(f"{label_text} decade$^{{-1}}$", fontsize=fontz * 0.95)
                ax.tick_params(axis='y', labelsize=fontz * 0.95)
                ax.axhline(0, color=color_dict["black"], linestyle="--", linewidth=1.5, alpha=0.7)
                for period in spinup_period:
                    ax.axvline(period, color=color_dict["gray"], linestyle="--", linewidth=1.5, alpha=0.7)
                ax.set_title(f"({panel_id}{col + 3}) {var_text} Trends ({rname})", fontsize=fontz*1.05, loc='left')
                
                # Legend
                legend_handles = [
                    Line2D([0], [0], color=color_dict["black"], linestyle=style_dict['filter'][0],
                           linewidth=style_dict['filter'][1], alpha=style_dict['filter'][2], label=trend_labels[0]),
                    Line2D([0], [0], marker=style_dict['mean'][0], color=color_dict["red"],
                           markersize=style_dict['mean'][1] * 0.3, linewidth=0, label=trend_labels[-1])
                ]

                ax.legend( handles=legend_handles,loc="upper center", fontsize=fontz * 0.95, frameon=True, 
                           labelspacing=0.2, columnspacing=0.5, handletextpad=0.2, ncol=3 )
                
                # Add annotation box with mean trend range
                if show_trend_range_box:
                    # Default text position
                    text_x, text_y = 0.02, 0.98
                    if text_position_dict is not None: 
                        pos = text_position_dict.get(varname, {'x': 0.02, 'y': 0.98})
                        text_x, text_y = pos['x'], pos['y']
                    
                    if plot_dict and varname in plot_dict:
                        trend_unit = plot_dict[varname]["unit"] + " decade$^{{-1}}$"
                    else:
                        trend_unit = f"{label_text} decade$^{{-1}}$"
                        
                    # Choose mode
                    if trend_mode == "range":
                        valid_inds = [(y >= 350) and (y <= 2000) for y in mid_years]
                        filtered_means = [np.nanmean(box_data[i]) for i, valid in enumerate(valid_inds) if valid]
                        if filtered_means:
                            txt = f"Mean trend: [{np.min(filtered_means):.2f}, {np.max(filtered_means):.2f}] {trend_unit}"
                        else:
                            txt = f"Mean trend: n/a"
                    else:
                        # ---- LONG-TERM trend (single value) ----
                        x = trend_all[varname][rname]["x"]
                        y = trend_all[varname][rname]["y"]
                        if use_filtered_for_trend:
                            yf = self.apply_centered_filter(y, window=filter_window)

                        # pick series for fitting
                        y_fit = yf if (use_filtered_for_trend and show_filter and (yf is not None)) else y
                        # choose fit range
                        if trend_fit_range is None:
                            x0, x1 = float(syear), float(eyear)
                        else:
                            x0, x1 = map(float, trend_fit_range)
                            
                        m = np.isfinite(x) & np.isfinite(y_fit) & (x >= x0) & (x <= x1)
                        if m.sum() >= 3:
                            slope_per_year = np.polyfit(x[m], y_fit[m], 1)[0]
                            slope_per_dec  = 10.0 * slope_per_year
                            txt = (f"Long-term Trend ({int(x0)}–{int(x1)}yr): "
                                   f"{slope_per_dec:.5f} {trend_unit}")
                        else:
                            txt = "Long-term trend: n/a"
                            
                    ax.text(
                        text_x, text_y, txt,
                        transform=ax.transAxes,
                        ha="center", va="center",
                        fontsize=fontz * 0.9,
                        color="black",
                        bbox=dict(
                            boxstyle="round,pad=0.35",
                            facecolor="white",
                            edgecolor="red",
                            linewidth=1.5,
                            alpha=0.95,
                        ),
                        zorder=10,
                    )
        
        fig.subplots_adjust(wspace=wspace, hspace=hspace)
        fig.savefig(fig_name, dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches)
        print(f"[PLOT] Saved: {fig_name}")
        plt.show()
        return


class CoupledBudgetPlotter:
    """
    Plotter for Carolyn’s comments:
      (1) P-E (or E-P) consistency with ocean volume/SSH drift
      (2) integrated ocean salt conservation (plus FW-equiv implied by salt drift)

    Notes on sign conventions (global-mean):
      - If plotting E-P: positive means freshwater LOSS from the ocean (SSH tends to fall)
        Compare against:  -(1/Ao) dV/dt  (freshwater loss positive)
      - If plotting P-E: positive means freshwater GAIN by the ocean (SSH tends to rise)
        Compare against:  +(1/Ao) dV/dt  (freshwater gain positive)

    Update (2026-02):
      - Combine the old (b) and (d) into a single integrated panel that compares:
          ΔSSH, ∫(flux)dt, and h_FW^salt
      - Switch figure layout to 1x3 panels: (a) rates, (b) integrated+salt, (c) salt conservation
    """

    def __init__(self, xr_engine="netcdf4"):
        self.xr_engine = xr_engine

    @staticmethod
    def _require(path, label):
        if not os.path.exists(path):
            raise FileNotFoundError(f"[ERROR] Missing {label}: {path}")
        return path

    @staticmethod
    def _load_tidy_series(nc, var_name, year_name="year", engine=None):
        """
        Load a tidy variable and return (years, vals).

        Behavior:
          - If 'rgn' dimension exists, select rgn=0 (GLOBAL).
          - Otherwise return full array.
        """
        ds = xr.open_dataset(nc, decode_times=False, engine=engine)
        try:
            # --- year axis ---
            if year_name in ds:
                years = np.asarray(ds[year_name].values)
            elif year_name in ds.coords:
                years = np.asarray(ds.coords[year_name].values)
            else:
                raise KeyError(f"'{year_name}' not found in {nc}")

            # --- variable ---
            if var_name not in ds:
                raise KeyError(f"'{var_name}' not in {nc}. Found vars: {list(ds.data_vars)}")

            da = ds[var_name]

            # --- if regional dimension exists, select GLOBAL (rgn=0) ---
            if "rgn" in da.dims:
                da = da.isel(rgn=0)

            vals = np.asarray(da.values, dtype=float)
        finally:
            ds.close()

        return years.astype(int, copy=False), vals

    @staticmethod
    def _match_years(y_ref, y_other):
        """
        Intersect on years to compare two time series.
        Returns (common_years, idx_in_ref, idx_in_other).
        """
        y_ref = np.asarray(y_ref, dtype=int)
        y_other = np.asarray(y_other, dtype=int)

        # Ensure sorted years for searchsorted
        s_ref = np.argsort(y_ref)
        s_oth = np.argsort(y_other)
        y_ref_s = y_ref[s_ref]
        y_oth_s = y_other[s_oth]

        common = np.intersect1d(y_ref_s, y_oth_s)
        i_ref = s_ref[np.searchsorted(y_ref_s, common)]
        i_oth = s_oth[np.searchsorted(y_oth_s, common)]
        return common, i_ref, i_oth

    @staticmethod
    def _trend_per_year(x_years, y):
        """
        Simple linear trend (units per year) ignoring NaNs.
        """
        x = np.asarray(x_years, dtype=float)
        y = np.asarray(y, dtype=float)
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 3:
            return np.nan
        coef = np.polyfit(x[m], y[m], 1)
        return float(coef[0])

    def plot_global_fw_and_salt_budget(
        self,
        atm_tidy_nc,
        ocn_conserve_nc,
        out_path,
        flux_var,
        flux_units="mm/yr",
        flux_form="P-E",
        flux_data_form="E-P",
        filter_window=None,
        trend_fit_range=None,
        figsize=None,                 # <-- now optional
        layout="1x3",                 # <-- NEW: "1x3" or "3x1"
        syear=0,
        eyear=2500,
        dpi=150,
        title=None,
        do_show=False,
        flux_ylim=None,
        fluxint_ylim=None,
        saltfrac_ylim=None,
        add_zero_lines=True,
        fontz=12,
    ):
        """
        Makes a 1x3 figure:
          (a) flux (P-E or E-P) (mm/yr) vs volume-based SSH rate term (mm/yr)
          (b) Integrated comparison (mm): ΔSSH, ∫(flux)dt, and salt FW-equivalent + constants/equation
          (c) Ocean salt mass change frac (×1e-6)
        """

        flux_form = str(flux_form).strip()
        flux_data_form = str(flux_data_form).strip()
        if flux_data_form not in ("E-P", "P-E") or flux_form not in ("E-P", "P-E"):
            raise ValueError("flux_form and flux_data_form must be 'E-P' or 'P-E'")

        # --- load ocean combined diagnostics ---
        ds_o = xr.open_dataset(ocn_conserve_nc, decode_times=False, engine=self.xr_engine)
        try:
            y_o = np.asarray(ds_o["year"].values, dtype=int)

            ssh_mm   = np.asarray(ds_o["ocean_volume_change_equiv_mm"].values, dtype=float)
            ssh_rate = np.asarray(ds_o["ocean_volume_tendency_equiv_mm_per_yr"].values, dtype=float)

            salt_frac  = np.asarray(ds_o["ocean_salt_mass_change_frac"].values, dtype=float)
            fw_salt_mm = np.asarray(ds_o["ocean_salt_freshwater_equiv_mm"].values, dtype=float)

            Ao     = float(ds_o.attrs.get("ocean_area_m2", np.nan))
            rho0   = float(ds_o.attrs.get("rho0_kg_m3_for_salt_mass", np.nan))
            S0     = float(ds_o.attrs.get("S0_kgkg_for_fw_equiv", np.nan))
            S0_src = str(ds_o.attrs.get("S0_source_for_fw_equiv", ""))

            salt_mass_kg = None
            if "ocean_salt_mass_kg" in ds_o:
                salt_mass_kg = np.asarray(ds_o["ocean_salt_mass_kg"].values, dtype=float)

            iS0_attr = ds_o.attrs.get("baseline_year_index_salt", None)
            try:
                iS0_attr = int(iS0_attr) if iS0_attr is not None else None
            except Exception:
                iS0_attr = None
        finally:
            ds_o.close()

        # --- load atmosphere flux tidy (rgn=0 if present) ---
        y_a, flux = self._load_tidy_series(
            atm_tidy_nc, flux_var, year_name="year", engine=self.xr_engine
        )

        # --- optional: use filtered series if present (also rgn=0 if present) ---
        if filter_window is not None and filter_window > 1:
            ds = xr.open_dataset(atm_tidy_nc, decode_times=False, engine=self.xr_engine)
            try:
                key = f"filtered_{int(filter_window)}yr"
                if key in ds:
                    da_f = ds[key]
                    if "rgn" in da_f.dims:
                        da_f = da_f.isel(rgn=0)
                    flux = np.asarray(da_f.values, dtype=float)
            finally:
                ds.close()

        # --- convert data to requested plotting form ---
        if flux_data_form != flux_form:
            flux = -flux

        flux_label = r"$P-E$" if flux_form == "P-E" else r"$E-P$"
        flux_int_label = r"$\int (P-E)\,dt$" if flux_form == "P-E" else r"$\int (E-P)\,dt$"

        # --- intersect years ---
        y_common, ia, io = self._match_years(y_a, y_o)
        flux_c       = np.asarray(flux[ia], dtype=float)
        ssh_rate_c   = np.asarray(ssh_rate[io], dtype=float)
        ssh_mm_c     = np.asarray(ssh_mm[io], dtype=float)
        salt_frac_c  = np.asarray(salt_frac[io], dtype=float)
        fw_salt_mm_c = np.asarray(fw_salt_mm[io], dtype=float)

        # --- Msalt0-related vars (ONLY add/compute; does not affect plotting) ---
        Msalt0 = np.nan
        salt_mass_trend_kg_yr = np.nan
        salt_mass_kg_c = None
        if salt_mass_kg is not None:
            salt_mass_kg_c = np.asarray(salt_mass_kg[io], dtype=float)

            # baseline index: prefer the one saved in attrs; otherwise first finite in the intersected series
            if (iS0_attr is not None) and (0 <= iS0_attr < len(salt_mass_kg)):
                Msalt0 = float(salt_mass_kg[iS0_attr])
            else:
                m0 = np.isfinite(salt_mass_kg_c)
                if m0.any():
                    Msalt0 = float(salt_mass_kg_c[int(np.argmax(m0))])

        # --- optional plotting window ---
        if (syear is not None) or (eyear is not None):
            lo = -np.inf if syear is None else int(syear)
            hi =  np.inf if eyear is None else int(eyear)
            wplot = (y_common >= lo) & (y_common <= hi)
            y_common = y_common[wplot]
            flux_c = flux_c[wplot]
            ssh_rate_c = ssh_rate_c[wplot]
            ssh_mm_c = ssh_mm_c[wplot]
            salt_frac_c = salt_frac_c[wplot]
            fw_salt_mm_c = fw_salt_mm_c[wplot]
            if salt_mass_kg_c is not None:
                salt_mass_kg_c = salt_mass_kg_c[wplot]

        # --- fit mask ---
        if trend_fit_range is not None:
            y0, y1 = trend_fit_range
            fit_mask = (y_common >= int(y0)) & (y_common <= int(y1))
        else:
            fit_mask = np.isfinite(y_common)

        # --- Msalt trend (kg/yr) now that fit_mask exists ---
        if salt_mass_kg_c is not None:
            salt_mass_trend_kg_yr = self._trend_per_year(
                y_common[fit_mask],
                salt_mass_kg_c[fit_mask],
            )

        # --- ocean sign convention to match chosen flux_form in rate panel ---
        if flux_form == "P-E":
            ssh_rate_match = ssh_rate_c
            fw_salt_mm_c = -fw_salt_mm_c
            ssh_rate_label = r"$\frac{1}{A_o}\,\frac{dV}{dt}$"
            ssh_note = r"$\frac{1}{A_o}\,\frac{dV}{dt} \equiv \frac{d\,\mathrm{SSH}}{dt}$"
        else:
            ssh_rate_match = -ssh_rate_c
            fw_salt_mm_c = fw_salt_mm_c 
            ssh_rate_label = r"$-\frac{1}{A_o}\,\frac{dV}{dt}$"
            ssh_note = r"$-\frac{1}{A_o}\,\frac{dV}{dt} \equiv -\frac{d\,\mathrm{SSH}}{dt}$"

        flux_tr = self._trend_per_year(y_common[fit_mask], flux_c[fit_mask])
        ssh_tr  = self._trend_per_year(y_common[fit_mask], ssh_rate_match[fit_mask])

        # --- integrate flux (mm/yr) to mm using actual year spacing ---
        flux_int = np.full_like(flux_c, np.nan, dtype=float)
        m = np.isfinite(flux_c) & np.isfinite(y_common)
        if m.any():
            i0 = int(np.argmax(m))
            flux_int[i0] = 0.0
            for i in range(i0 + 1, len(flux_c)):
                if (np.isfinite(flux_c[i]) and np.isfinite(flux_int[i - 1]) and
                    np.isfinite(y_common[i]) and np.isfinite(y_common[i - 1])):
                    dy = float(y_common[i] - y_common[i - 1])
                    flux_int[i] = flux_int[i - 1] + flux_c[i] * dy

        # --- SSH change; shift to 0 at first finite overlap year ---
        ssh_int = ssh_mm_c.copy()
        m2 = np.isfinite(ssh_int)
        if m2.any():
            j0 = int(np.argmax(m2))
            ssh_int = ssh_int - ssh_int[j0]

        # Direct comparison in integrated panel:
        if flux_form == "E-P":
            ssh_int_plot = -ssh_int
            ssh_int_label = r"$-\frac{1}{A_o}\,\frac{dV}{dt} \equiv -\frac{d\,\mathrm{SSH}}{dt}$"
        else:
            ssh_int_plot = ssh_int
            ssh_int_label = r"$\frac{1}{A_o}\,\frac{dV}{dt} \equiv \frac{d\,\mathrm{SSH}}{dt}$"

        # --- plot (layout) ---
        if layout == "1x3":
            nrows, ncols = 1, 3
            default_figsize = (18, 5)
            ax_order = (0, 1, 2)
        elif layout == "3x1":
            nrows, ncols = 3, 1
            default_figsize = (6, 14)
            ax_order = (0, 1, 2)
        else:
            raise ValueError("layout must be '1x3' or '3x1'")

        if figsize is None:
            figsize = default_figsize

        fig, ax = plt.subplots(nrows, ncols, figsize=figsize, dpi=dpi)

        # normalize axes to 1D list
        ax = np.atleast_1d(ax).ravel()
        ax0, ax1, ax2 = ax[ax_order[0]], ax[ax_order[1]], ax[ax_order[2]]
        ax0.set_xlim(syear-1, eyear+1)
        ax1.set_xlim(syear-1, eyear+1)
        ax2.set_xlim(syear-1, eyear+1)
        
        # (a) rates
        ax0.plot(y_common, ssh_rate_match, label=ssh_rate_label, color="tab:blue")
        ax0.plot(y_common, flux_c, label=flux_label, color="tab:red")
        ax0.set_title("(a) Global Freshwater Flux vs Ocean Volume Tendency", fontsize=fontz, loc="left")
        ax0.set_ylabel(r"$\mathrm{mm\ yr^{-1}}$", fontsize=fontz * 0.95)
        ax0.set_xlabel("Year", fontsize=fontz * 0.95)
        ax0.tick_params(axis="both", which="major", labelsize=fontz * 0.9)
        if flux_ylim is not None:
            ax0.set_ylim(flux_ylim)
        if add_zero_lines:
            ax0.axhline(0.0, lw=0.8, linestyle=":", color="black", alpha=0.8, zorder=0)
        ax0.grid(True, alpha=0.3)
        leg = ax0.legend(
            loc="best",
            frameon=True,
            framealpha=0.5,
            fontsize=fontz * 0.85,
        )
        leg.get_frame().set_linewidth(0.6)
        
        if np.any(fit_mask):
            yfit_min = int(np.nanmin(y_common[fit_mask]))
            yfit_max = int(np.nanmax(y_common[fit_mask]))
        else:
            yfit_min = int(np.nanmin(y_common))
            yfit_max = int(np.nanmax(y_common))

        txt = (
            rf"Trend({flux_form}, {yfit_min}-{yfit_max}) = {flux_tr*10:.3e} mm yr$^{{-1}}$ decade$^{{-1}}$" "\n"
            rf"Trend({ssh_rate_label}, {yfit_min}-{yfit_max}) = {ssh_tr*10:.3e} mm yr$^{{-1}}$ decade$^{{-1}}$"
        )
        ax0.text(0.02, 0.02, txt, transform=ax0.transAxes,
                 va="bottom", ha="left", fontsize=fontz * 0.85)
        ax0.text(0.02, 0.85, ssh_note, transform=ax0.transAxes,
                 va="bottom", ha="left", fontsize=fontz * 0.85)

        # (b) integrated combined: ΔSSH, ∫flux dt, h_FW^salt
        eqn = (
            r"$h_{\mathrm{FW}}^{\mathrm{salt}}(t)"
            r" = -\frac{\Delta M_{\mathrm{salt}}(t)}{\rho_0\,S_0\,A_o}$"
        )
        ax1.plot(y_common, ssh_int_plot, label=ssh_int_label, color="tab:blue")
        ax1.plot(y_common, flux_int, label=flux_int_label, color="tab:red")
        ax1.plot(y_common, fw_salt_mm_c, label=rf"{eqn}", color="black")
        ax1.set_title("(b) Integrated Freshwater Budget", fontsize=fontz, loc="left")
        ax1.set_ylabel(r"$\mathrm{mm}$", fontsize=fontz * 0.95)
        ax1.set_xlabel("Year", fontsize=fontz * 0.95)
        ax1.tick_params(axis="both", which="major", labelsize=fontz * 0.9)
        if fluxint_ylim is not None:
            ax1.set_ylim(fluxint_ylim)
        if add_zero_lines:
            ax1.axhline(0.0, lw=0.8, linestyle=":", color="black", alpha=0.8, zorder=0)
        ax1.grid(True, alpha=0.3)
        leg = ax1.legend(
            loc="best",
            frameon=True,
            framealpha=0.5,          # transparency
            fontsize=fontz * 0.85,
        )
        leg.get_frame().set_linewidth(0.6)

        const = (
            fr"$A_o={Ao:.3g}\,\mathrm{{m}}^2,\ "
            fr"\rho_0={rho0:.4g}\,\mathrm{{kg\,m}}^{{-3}},\ "
            fr"S_0={S0:.4g}\,\mathrm{{PSU}}$"
        )
        ax1.text(
            0.35, 0.90, const,
            transform=ax1.transAxes,
            ha="left", va="bottom",
            fontsize=fontz * 0.85,
        )

        # (c) salt conservation
        ax2.plot(
            y_common, salt_frac_c * 1e6,
            label=r"$\frac{\Delta M_{\mathrm{salt}}}{M_{\mathrm{salt},0}}$",
            color="black",
        )
        ax2.set_title("(c) Global Ocean Salt Conservation", fontsize=fontz, loc="left")
        ax2.set_ylabel(r"$\mathrm{fraction\ (\times10^{-6})}$", fontsize=fontz * 0.95)
        ax2.set_xlabel("Year", fontsize=fontz * 0.95)
        ax2.tick_params(axis="both", which="major", labelsize=fontz * 0.9)
        if saltfrac_ylim is not None:
            ax2.set_ylim(saltfrac_ylim)
        if add_zero_lines:
            ax2.axhline(0.0, lw=0.8, linestyle=":", color="black", alpha=0.8, zorder=0)
        ax2.grid(True, alpha=0.3)
        leg = ax2.legend(
            loc="upper left",
            frameon=True,
            framealpha=0.5,
            fontsize=fontz * 0.85,
        )
        leg.get_frame().set_linewidth(0.6)
        
        if np.isfinite(salt_frac_c).any():
            txt_c = (
                fr"$M_{{\mathrm{{salt}},0}} = {Msalt0:.3e}\ \mathrm{{kg}}$" "\n"
                fr"$\frac{{d\,M_{{\mathrm{{salt}}}}}}{{dt}} = {salt_mass_trend_kg_yr:.3e}\ \mathrm{{kg\,yr^{{-1}}}}$"
            )
            ax2.text(
                0.65, 0.25,
                txt_c,
                transform=ax2.transAxes,
                va="top", ha="left",
                fontsize=fontz * 0.85,
            )

        if title:
            fig.suptitle(title, y=1.02)

        fig.tight_layout()
        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")
        if do_show:
            plt.show()
        plt.close(fig)

        return out_path
