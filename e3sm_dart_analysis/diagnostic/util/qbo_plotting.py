"""Figures for the QBO sections written by ``util.qbo_diagnostics``.

``QboPlotter`` draws, for one variable (``ua`` or ``ta``) and the files
``qbo_<source>_<var>_<period>.nc``:

* time-pressure sections of a reference and all experiments (RMSE/PCC vs the
  reference annotated) and their differences from the reference;
* amplitude profiles;
* RMSE/PCC skill bars vs one reference, or averaged over several references
  (block-bootstrap 95% CI or reference half-range as error bars);
* persistence of the initial vertical profile and valid-time pattern
  correlation with the reference vs lead month;
* vertical coherence: corr_t(X(p), X(p0)).

Sections can be plotted raw (``field="qbo"``) or as anomalies
(``field="qbo_anom"``) from each source's own monthly climatology
(``own_clim``), the reference climatology (``ref_clim``, comparable across
experiments), or the first month (``first_month``). Temperature defaults to
``ref_clim`` anomalies.

``plot_t_shear_corr_profiles`` and ``plot_t_shear_scatter`` diagnose the
thermal-wind coupling of the QBO: corr_t(T', dU/dln p) by level and scatter
in a layer.
"""
from __future__ import annotations

import math
import string
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from configs.hind_experiment_config import PERIOD, STYLE_MAP
from util.qbo_diagnostics import QBO_VARIABLES, qbo_section_path

Field = Literal["qbo", "qbo_anom"]
AnomMode = Literal["own_clim", "ref_clim", "first_month"]
AmpMethod = Literal["p2p", "rms", "rms_detrend"]
_AMP_VARS = {"p2p": "amp_p2p", "rms": "amp_rms", "rms_detrend": "amp_rms_detrend"}
_EMPHASIS = ("DARTOP", "ERA5OP", "CTRL")


@dataclass
class QboPlotParameter:
    main_title: str
    variable: str = "ua"
    period: str = PERIOD
    default_field: Field = "qbo"
    default_anom_mode: AnomMode = "ref_clim"

    @property
    def var_label(self) -> str:
        return QBO_VARIABLES[self.variable]["label"]

    @property
    def var_units(self) -> str:
        return QBO_VARIABLES[self.variable]["units"]

    @classmethod
    def for_variable(cls, variable: str, **kwargs) -> "QboPlotParameter":
        meta = QBO_VARIABLES[variable]
        kwargs.setdefault("main_title", f"QBO {meta['long_name']} (5°S–5°N)")
        if variable == "ta":
            kwargs.setdefault("default_field", "qbo_anom")
        return cls(variable=variable, **kwargs)


@dataclass
class SavedQboDiag:
    name: str
    qbo: xr.DataArray        # (time, plev)
    amplitude: xr.DataArray  # (plev,)


def _line_kwargs(name: str, style_map: dict | None, ref: str | None = None) -> dict[str, Any]:
    kw = dict((style_map or STYLE_MAP).get(name, {}))
    kw.setdefault("linestyle", "-")
    kw.setdefault("linewidth", 2.2 if name in (ref, *_EMPHASIS) else 1.6)
    return kw


def _yyyymm(tvals) -> np.ndarray:
    return np.array([int(t.year) * 100 + int(t.month) if hasattr(t, "year")
                     else int(np.datetime_as_string(t, unit="M").replace("-", "")) for t in tvals])


def _pcorr(a: np.ndarray, b: np.ndarray, min_n: int = 3) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(m) < min_n:
        return np.nan
    a, b = a[m] - a[m].mean(), b[m] - b[m].mean()
    den = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / den) if den > 0 else np.nan


class QboPlotter:
    """Plot saved QBO sections of one variable (see module docstring)."""

    def __init__(self, parameter: QboPlotParameter, qbo_dir):
        self.parameter = parameter
        self.qbo_dir = Path(qbo_dir)

    # ------------------------------------------------------------------ I/O
    def load(self, source: str, amplitude_method: AmpMethod = "rms") -> SavedQboDiag:
        path = qbo_section_path(self.qbo_dir, source, self.parameter.variable, self.parameter.period)
        if not path.is_file():
            raise FileNotFoundError(f"Missing QBO section: {path}")
        ds = xr.open_dataset(path)
        qbo = ds["qbo"].transpose("time", "plev").load()
        amp = ds[_AMP_VARS[amplitude_method]].squeeze().load().rename("amplitude")
        ds.close()
        return SavedQboDiag(name=source, qbo=qbo, amplitude=amp)

    def load_set(self, ref: str, sources, amplitude_method: AmpMethod = "rms",
                 skip_bad: bool = True) -> tuple[SavedQboDiag, list[SavedQboDiag]]:
        """Reference and the sources that exist and share its months and levels."""
        ref_d = self.load(ref, amplitude_method)
        tests = []
        for s in sources:
            try:
                d = self.load(s, amplitude_method)
                self.check_same_grid(d.qbo, ref_d.qbo)
                tests.append(d)
            except Exception as e:
                if not skip_bad:
                    raise
                print(f"[SKIP] {s}: {e}")
        return ref_d, tests

    @staticmethod
    def check_same_grid(a: xr.DataArray, b: xr.DataArray) -> None:
        if not np.array_equal(_yyyymm(a.time.values), _yyyymm(b.time.values)):
            raise ValueError("months differ from the reference")
        if not np.array_equal(a.plev.values, b.plev.values):
            raise ValueError("pressure levels differ from the reference")

    # ------------------------------------------------------------ fields
    @staticmethod
    def _remove_monthly_clim(x: xr.DataArray, clim_from: xr.DataArray | None = None) -> xr.DataArray:
        clim = (clim_from if clim_from is not None else x).groupby("time.month").mean("time")
        return (x.groupby("time.month") - clim).drop_vars("month", errors="ignore")

    def to_field(self, qbo: xr.DataArray, ref: xr.DataArray, field: Field | None = None,
                 anom_mode: AnomMode | None = None) -> xr.DataArray:
        field = field or self.parameter.default_field
        anom_mode = anom_mode or self.parameter.default_anom_mode
        if field == "qbo":
            return qbo
        if anom_mode == "own_clim":
            return self._remove_monthly_clim(qbo)
        if anom_mode == "first_month":
            return qbo - qbo.isel(time=0)
        return self._remove_monthly_clim(qbo, clim_from=ref)

    def _field_tag(self, field, anom_mode) -> str:
        field = field or self.parameter.default_field
        anom_mode = anom_mode or self.parameter.default_anom_mode
        return "raw" if field == "qbo" else f"anom_{anom_mode}"

    # ------------------------------------------------------------ metrics
    @staticmethod
    def rmse_and_pcc(test: xr.DataArray, ref: xr.DataArray) -> tuple[float, float]:
        """RMSE and Pearson correlation over all (time, plev) points."""
        a, b = np.asarray(test.values).ravel(), np.asarray(ref.values).ravel()
        m = np.isfinite(a) & np.isfinite(b)
        a, b = a[m], b[m]
        rmse = float(np.sqrt(np.mean((a - b) ** 2))) if a.size else np.nan
        pcc = float(np.corrcoef(a, b)[0, 1]) if a.size >= 2 else np.nan
        return rmse, pcc

    @classmethod
    def block_bootstrap_skill(cls, test: xr.DataArray, refs: list[xr.DataArray], *, n_boot: int = 500,
                              block_len: int = 6, seed: int = 0):
        """Reference-averaged RMSE/PCC and their 95% block-bootstrap CI (blocks of months)."""
        nt = test.sizes["time"]
        rng = np.random.default_rng(seed)
        rmse0 = float(np.nanmean([cls.rmse_and_pcc(test, r)[0] for r in refs]))
        pcc0 = float(np.nanmean([cls.rmse_and_pcc(test, r)[1] for r in refs]))
        starts = np.arange(0, nt, block_len)
        nblocks = int(np.ceil(nt / block_len))
        boot = np.empty((n_boot, 2))
        for b in range(n_boot):
            idx = np.concatenate([np.arange(s, min(s + block_len, nt)) for s in rng.choice(starts, nblocks)])[:nt]
            tb = test.isel(time=idx)
            stats = [cls.rmse_and_pcc(tb, r.isel(time=idx)) for r in refs]
            boot[b] = np.nanmean(stats, axis=0)
        lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
        return rmse0, pcc0, (lo[0], hi[0]), (lo[1], hi[1])

    # ------------------------------------------------------------ helpers
    @staticmethod
    def _grid(n: int, nrows: int | None, ncols: int | None) -> tuple[int, int]:
        if nrows is None and ncols is None:
            ncols = 3
        if nrows is None:
            return int(math.ceil(n / ncols)), int(ncols)
        if ncols is None:
            return int(nrows), int(math.ceil(n / nrows))
        if nrows * ncols < n:
            raise ValueError(f"Grid {nrows}x{ncols} too small for {n} panels")
        return int(nrows), int(ncols)

    @staticmethod
    def _time_ticks(ax, tvals, every: int, fontsize: int) -> None:
        idx = np.arange(0, len(tvals), every)
        ax.set_xticks(idx)
        ax.set_xticklabels([f"{t.year:04d}-{t.month:02d}" if hasattr(t, "year")
                            else np.datetime_as_string(t, unit="M") for t in np.asarray(tvals)[idx]],
                           fontsize=fontsize, rotation=30, ha="right")

    def _default_levels(self, diff: bool):
        if self.parameter.variable == "ta":
            return (np.arange(-4, 4.1, 0.5), [-4, -2, -1, 0, 1, 2, 4]) if diff else \
                   (np.arange(-6, 6.1, 0.5), [-6, -3, -1, 1, 3, 6])
        return (np.arange(-20, 21, 2), [-20, -10, -5, 0, 5, 10, 20]) if diff else \
               (np.arange(-50, 51, 5), [-50, -25, -5, 5, 25, 50])

    def _section_grid(self, panels, *, levels, ticks, cbar_label, nrows, ncols, w_per_col, h_per_row,
                      font_scale, time_every, cmap, wspace, hspace, annotate_ref=None, fig_id=0):
        n = len(panels)
        nrows, ncols = self._grid(n, nrows, ncols)
        fs = dict(title=round(font_scale), label=round(font_scale * 0.95), tick=round(font_scale * 0.95),
                  anno=round(font_scale * 0.8))
        fig, axes = plt.subplots(nrows, ncols, figsize=(w_per_col * ncols, h_per_row * nrows), squeeze=False)
        cs = None
        for i, (title, da) in enumerate(panels):
            r, c = divmod(i, ncols)
            ax = axes[r, c]
            p = da.plev.values
            cs = ax.contourf(np.arange(da.sizes["time"]), p, da.transpose("plev", "time").values,
                             levels=levels, cmap=cmap, extend="both")
            ax.set_title(f"({string.ascii_lowercase[fig_id + i]}) {title}", fontsize=fs["title"], loc="left")
            ax.set_title(self.parameter.main_title, fontsize=fs["title"] * 0.95, loc="right")
            ax.set_yscale("log")
            ax.set_ylim(100, 1)
            if c == 0:
                ax.set_ylabel("hPa", fontsize=fs["label"])
                ax.tick_params(axis="y", labelsize=fs["tick"])
            else:
                ax.tick_params(axis="y", labelleft=False)
            if r == nrows - 1 or i + ncols >= n:
                self._time_ticks(ax, da.time.values, time_every, fs["tick"])
                ax.set_xlabel("Year-Month", fontsize=fs["label"])
            else:
                ax.set_xticks([])
            if annotate_ref is not None and i > 0:
                rmse, pcc = self.rmse_and_pcc(da, annotate_ref)
                ax.text(0.68, 0.05, f"RMSE: {rmse:.2f} {self.parameter.var_units}\nPCC: {pcc:.2f}",
                        transform=ax.transAxes, fontsize=fs["anno"],
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.7", alpha=0.85))
        for j in range(n, nrows * ncols):
            axes.flat[j].axis("off")
        fig.subplots_adjust(left=0.06, right=0.88, bottom=0.08, top=0.90, wspace=wspace, hspace=hspace)
        cbar = fig.colorbar(cs, cax=fig.add_axes([0.90, 0.18, 0.02, 0.64]), ticks=ticks)
        cbar.ax.tick_params(labelsize=fs["tick"])
        cbar.set_label(cbar_label, fontsize=fs["label"])
        return fig

    # ------------------------------------------------------------ figures
    def plot_timepressure(self, ref: str, sources, *, field: Field | None = None, anom_mode: AnomMode | None = None,
                          levels=None, ticks=None, nrows=None, ncols=3, w_per_col=5.2, h_per_row=3.2,
                          font_scale=12, time_every=6, cmap="RdBu_r", wspace=0.25, hspace=0.35, fig_id=0,
                          skip_bad=True) -> plt.Figure:
        """Reference + every source; RMSE/PCC vs the reference in each source panel."""
        ref_d, tests = self.load_set(ref, sources, skip_bad=skip_bad)
        ref_f = self.to_field(ref_d.qbo, ref_d.qbo, field, anom_mode)
        panels = [(ref, ref_f)] + [(t.name, self.to_field(t.qbo, ref_d.qbo, field, anom_mode)) for t in tests]
        d_levels, d_ticks = self._default_levels(diff=False)
        p = self.parameter
        return self._section_grid(panels, levels=levels if levels is not None else d_levels,
                                  ticks=ticks or d_ticks, cbar_label=f"{p.var_label} ({p.var_units})",
                                  nrows=nrows, ncols=ncols, w_per_col=w_per_col, h_per_row=h_per_row,
                                  font_scale=font_scale, time_every=time_every, cmap=cmap, wspace=wspace,
                                  hspace=hspace, annotate_ref=ref_f, fig_id=fig_id)

    def plot_timepressure_diffs(self, ref: str, sources, *, field: Field | None = None,
                                anom_mode: AnomMode | None = None, levels=None, ticks=None, nrows=None, ncols=3,
                                w_per_col=5.2, h_per_row=3.2, font_scale=12, time_every=6, cmap="RdBu_r",
                                wspace=0.25, hspace=0.35, fig_id=0, skip_bad=True) -> plt.Figure:
        """Source − reference sections."""
        ref_d, tests = self.load_set(ref, sources, skip_bad=skip_bad)
        ref_f = self.to_field(ref_d.qbo, ref_d.qbo, field, anom_mode)
        panels = [(f"{t.name} − {ref}", self.to_field(t.qbo, ref_d.qbo, field, anom_mode) - ref_f) for t in tests]
        if not panels:
            raise ValueError("No sources to difference")
        d_levels, d_ticks = self._default_levels(diff=True)
        p = self.parameter
        return self._section_grid(panels, levels=levels if levels is not None else d_levels,
                                  ticks=ticks or d_ticks, cbar_label=f"Δ{p.var_label} ({p.var_units})",
                                  nrows=nrows, ncols=ncols, w_per_col=w_per_col, h_per_row=h_per_row,
                                  font_scale=font_scale, time_every=time_every, cmap=cmap, wspace=wspace,
                                  hspace=hspace, fig_id=fig_id)

    def plot_amplitude(self, ref: str, sources, *, amplitude_method: AmpMethod = "rms", style_map=None,
                       figsize=(6.8, 7.8), font_scale=12, xlim=None, ylim=(100, 1), legend_outside=True,
                       fig_id=0, skip_bad=True) -> plt.Figure:
        """Amplitude profiles of every source with the reference on top."""
        ref_d, tests = self.load_set(ref, sources, amplitude_method, skip_bad)
        fig, ax = plt.subplots(figsize=figsize)
        for d in [*tests, ref_d]:
            ax.plot(d.amplitude.values, d.amplitude.plev.values, label=d.name,
                    **_line_kwargs(d.name, style_map, ref))
        p = self.parameter
        ax.set_title(f"({string.ascii_lowercase[fig_id]}) {p.main_title}", fontsize=font_scale * 1.05, loc="left")
        ax.set_xlabel(f"{amplitude_method.upper()} amplitude ({p.var_units})", fontsize=font_scale)
        ax.set_ylabel("hPa", fontsize=font_scale)
        ax.set_yscale("log")
        ax.set_ylim(*ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)
        ax.grid(True, which="both", alpha=0.25)
        ax.tick_params(labelsize=font_scale * 0.95)
        if legend_outside:
            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=font_scale * 0.9, frameon=False)
        else:
            ax.legend(loc="best", fontsize=font_scale * 0.9, frameon=False)
        return fig

    def skill_table(self, refs, sources, *, n_boot=500, block_len=6, seed=0, skip_bad=True):
        """RMSE/PCC of every source averaged over ``refs`` (list) with CI and reference half-range."""
        refs = [refs] if isinstance(refs, str) else list(refs)
        ref_ds = [self.load(r) for r in refs]
        rows = []
        for s in sources:
            try:
                t = self.load(s)
                for r in ref_ds:
                    self.check_same_grid(t.qbo, r.qbo)
            except Exception as e:
                if not skip_bad:
                    raise
                print(f"[SKIP] {s}: {e}")
                continue
            per = np.array([self.rmse_and_pcc(t.qbo, r.qbo) for r in ref_ds])
            rmse, pcc, rmse_ci, pcc_ci = self.block_bootstrap_skill(
                t.qbo, [r.qbo for r in ref_ds], n_boot=n_boot, block_len=block_len, seed=seed)
            rows.append(dict(source=s, rmse=rmse, pcc=pcc, rmse_lo=rmse_ci[0], rmse_hi=rmse_ci[1],
                             pcc_lo=pcc_ci[0], pcc_hi=pcc_ci[1],
                             rmse_ref_halfrange=0.5 * np.ptp(per[:, 0]), pcc_ref_halfrange=0.5 * np.ptp(per[:, 1])))
        import pandas as pd
        return pd.DataFrame(rows)

    def plot_skill_bars(self, refs, sources, *, errorbar: Literal["bootstrap", "ref_spread"] = "bootstrap",
                        sort_by: Literal["input", "rmse", "pcc"] = "rmse", baseline: str = "CTRL",
                        n_boot=500, block_len=6, seed=0, figsize=(11.0, 4.6), font_scale=12, nrows=1, ncols=2,
                        rmse_ylim=(0, 20), pcc_ylim=(-0.8, 0.8), fig_id=0, skip_bad=True, table=None):
        """RMSE and PCC bars vs one reference or averaged over several.

        Error bars: 95% block-bootstrap CI (``bootstrap``) or half the range across
        references (``ref_spread``, needs >= 2 references). RMSE bars are annotated
        with the % change vs ``baseline``; with two references their own
        disagreement is shown as a floor. ``table``: a precomputed ``skill_table``
        for the same ``refs`` (``sources`` and the bootstrap options are then unused).
        Returns (figure, table).
        """
        refs = [refs] if isinstance(refs, str) else list(refs)
        tab = (table.copy() if table is not None else
               self.skill_table(refs, sources, n_boot=n_boot, block_len=block_len, seed=seed, skip_bad=skip_bad))
        if tab.empty:
            raise ValueError("No sources to score")
        if sort_by == "rmse":
            tab = tab.sort_values("rmse")
        elif sort_by == "pcc":
            tab = tab.sort_values("pcc", ascending=False)
        if errorbar == "ref_spread":
            if len(refs) < 2:
                raise ValueError("errorbar='ref_spread' needs at least two references")
            rmse_err, pcc_err = tab.rmse_ref_halfrange.values, tab.pcc_ref_halfrange.values
        else:
            rmse_err = np.vstack([tab.rmse - tab.rmse_lo, tab.rmse_hi - tab.rmse])
            pcc_err = np.vstack([tab.pcc - tab.pcc_lo, tab.pcc_hi - tab.pcc])
        x = np.arange(len(tab))
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        ax, ax2 = axes.flat[0], axes.flat[1]
        tick_fs = font_scale * 0.95
        ref_txt = refs[0] if len(refs) == 1 else f"avg over {', '.join(refs)}"
        ax.bar(x, tab.rmse, yerr=rmse_err, color="0.6", capsize=3)
        ax.set_title(f"({string.ascii_lowercase[fig_id]}) RMSE vs {ref_txt}", fontsize=font_scale * 1.05, loc="left")
        ax.set_ylabel(f"RMSE ({self.parameter.var_units})", fontsize=font_scale)
        if rmse_ylim:
            ax.set_ylim(*rmse_ylim)
        if baseline in tab.source.values:
            b = float(tab.loc[tab.source == baseline, "rmse"].iloc[0])
            dy = 0.06 * np.diff(ax.get_ylim())[0]
            for xi, v in zip(x, tab.rmse):
                ax.text(xi, v + dy, f"{100 * (v - b) / b:+.0f}%", ha="center", va="bottom", fontsize=tick_fs * 0.85)
        ax2.bar(x, tab.pcc, yerr=pcc_err, color="0.6", capsize=3)
        ax2.axhline(0.0, linewidth=0.8, alpha=0.6)
        ax2.set_title(f"({string.ascii_lowercase[fig_id + 1]}) PCC vs {ref_txt}", fontsize=font_scale * 1.05, loc="left")
        ax2.set_ylabel("PCC", fontsize=font_scale)
        if pcc_ylim:
            ax2.set_ylim(*pcc_ylim)
        if len(refs) == 2:
            r0, r1 = (self.load(r).qbo for r in refs)
            rmse_b, pcc_b = self.rmse_and_pcc(r0, r1)
            ax.axhline(rmse_b, color="k", linestyle="--", linewidth=0.8)
            ax.text(0.02, 0.96, f"RMSE({refs[0]},{refs[1]}) = {rmse_b:.2f}", transform=ax.transAxes,
                    va="top", fontsize=tick_fs * 0.9)
            ax2.text(0.02, 0.04, f"PCC({refs[0]},{refs[1]}) = {pcc_b:.2f}", transform=ax2.transAxes,
                     fontsize=tick_fs * 0.9)
        for a in (ax, ax2):
            a.set_xticks(x)
            a.set_xticklabels(tab.source, rotation=45, ha="right", fontsize=tick_fs)
            a.tick_params(axis="y", labelsize=tick_fs)
            a.grid(True, axis="y", alpha=0.25)
        err_txt = "95% block-bootstrap CI" if errorbar == "bootstrap" else "half-range across references"
        fig.text(0.99, 0.01, f"Error bars: {err_txt}", ha="right", fontsize=tick_fs * 0.8)
        fig.tight_layout()
        return fig, tab

    def plot_lagcorr_persistence(self, ref: str, sources, *, init_months=3, max_lag=18, style_map=None,
                                 figsize=(9.5, 8.0), font_scale=12, nrows=2, ncols=1, ylim=(-1.0, 1.2),
                                 fig_id=0, skip_bad=True) -> plt.Figure:
        """(a) corr_p(initial profile, profile at lag L); (b) corr_p(source, reference) at lag L."""
        ref_d, tests = self.load_set(ref, sources, skip_bad=skip_bad)
        ref_arr = ref_d.qbo.values
        lmax = min(max_lag, ref_arr.shape[0] - 1) if max_lag is not None else ref_arr.shape[0] - 1
        lags = np.arange(lmax + 1)
        series = [ref_d, *tests]
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, squeeze=False)
        ax1, ax2 = axes.flat[0], axes.flat[1]
        for d in series:
            arr = d.qbo.values
            x0 = np.nanmean(arr[:init_months], axis=0)
            kw = _line_kwargs(d.name, style_map, ref)
            ax1.plot(lags, [_pcorr(x0, arr[lag]) for lag in lags], label=d.name, **kw)
            ax2.plot(lags, [_pcorr(arr[lag], ref_arr[lag]) for lag in lags], label=d.name, **kw)
        for k, (ax, title) in enumerate(((ax1, f"Persistence vs initial {init_months}-month profile"),
                                         (ax2, f"Correlation with {ref}"))):
            ax.axhline(0.0, linewidth=0.8, alpha=0.6)
            ax.set_xlim(0, lags[-1])
            ax.set_ylim(*ylim)
            ax.set_title(f"({string.ascii_lowercase[fig_id + k]}) {title}", fontsize=font_scale * 1.05, loc="left")
            ax.set_ylabel("Correlation over levels", fontsize=font_scale)
            ax.tick_params(labelsize=font_scale * 0.95)
            ax.grid(True, alpha=0.25)
        ax2.set_xlabel("Lag (months)", fontsize=font_scale)
        if ncols > 1:
            ax1.set_xlabel("Lag (months)", fontsize=font_scale)
        h, lab = ax1.get_legend_handles_labels()
        fig.legend(h, lab, loc="center left", bbox_to_anchor=(0.88, 0.5), fontsize=font_scale * 0.9, frameon=False)
        fig.subplots_adjust(left=0.09, right=0.86, bottom=0.10, top=0.92, hspace=0.25, wspace=0.28)
        return fig

    def plot_vertical_coherence(self, ref: str, sources, *, ref_level_hpa=30.0, field: Field | None = None,
                                anom_mode: AnomMode | None = None, style_map=None, figsize=(6.8, 7.8),
                                font_scale=12, ylim=(100, 1), fig_id=0, skip_bad=True) -> plt.Figure:
        """r(p) = corr_t(X(t, p), X(t, p0)) for the reference and every source."""
        ref_d, tests = self.load_set(ref, sources, skip_bad=skip_bad)
        fig, ax = plt.subplots(figsize=figsize)
        for d in [*tests, ref_d]:
            x = self.to_field(d.qbo, ref_d.qbo, field, anom_mode)
            i0 = int(np.argmin(np.abs(x.plev.values - ref_level_hpa)))
            x0 = x.isel(plev=i0).values
            r = [_pcorr(x.isel(plev=i).values, x0, min_n=2) for i in range(x.sizes["plev"])]
            ax.plot(r, x.plev.values, label=d.name, **_line_kwargs(d.name, style_map, ref))
        ax.axvline(0.0, linewidth=0.8, alpha=0.5)
        ax.axhline(ref_level_hpa, linewidth=0.8, alpha=0.5, linestyle="--")
        ax.set_title(f"({string.ascii_lowercase[fig_id]}) {self.parameter.main_title}: coherence with "
                     f"{ref_level_hpa:g} hPa", fontsize=font_scale * 1.05, loc="left")
        ax.set_xlabel("Correlation over time", fontsize=font_scale)
        ax.set_ylabel("hPa", fontsize=font_scale)
        ax.set_xlim(-1, 1)
        ax.set_yscale("log")
        ax.set_ylim(*ylim)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=font_scale * 0.9, frameon=False)
        return fig


# ---------------------------------------------------------------------------
# Temperature - wind-shear coupling
# ---------------------------------------------------------------------------
def shear_dlnp(u: xr.DataArray) -> xr.DataArray:
    """dU/dln(p) by centred differences (one-sided at the end levels)."""
    lnp = np.log(u.plev.values.astype(float))
    U = u.transpose("time", "plev").values
    S = np.gradient(U, lnp, axis=1, edge_order=1)
    return xr.DataArray(S, dims=("time", "plev"), coords={"time": u.time, "plev": u.plev}, name="dU_dlnp")


def _t_and_shear(qbo_dir, source: str, period: str, ref_t: xr.DataArray | None):
    t = QboPlotter(QboPlotParameter.for_variable("ta", period=period), qbo_dir)
    u = QboPlotter(QboPlotParameter.for_variable("ua", period=period), qbo_dir)
    T = t.load(source).qbo
    U = u.load(source).qbo
    QboPlotter.check_same_grid(T, U)
    T_anom = t.to_field(T, ref_t if ref_t is not None else T, "qbo_anom",
                        "ref_clim" if ref_t is not None else "own_clim")
    return T_anom, shear_dlnp(U)


def t_shear_corr_profile(qbo_dir, source: str, period: str = PERIOD, ref_t: xr.DataArray | None = None):
    """corr_t(T'(p), dU/dln p (p)) of one source."""
    T, S = _t_and_shear(qbo_dir, source, period, ref_t)
    r = [_pcorr(T.isel(plev=i).values, S.isel(plev=i).values, min_n=2) for i in range(T.sizes["plev"])]
    return xr.DataArray(r, dims="plev", coords={"plev": T.plev}, name="corr_T_shear")


def plot_t_shear_corr_profiles(qbo_dir, ref: str, sources, *, period: str = PERIOD, band_hpa=(10.0, 70.0),
                               style_map=None, figsize=(6.8, 7.8), font_scale=12, xlim=(-1, 1), ylim=(100, 1),
                               fig_id=0, skip_bad=True) -> plt.Figure:
    """Profiles of corr[T', dU/dln p]; T' are anomalies from the reference monthly climatology.

    The legend gives each profile's mean over ``band_hpa`` (shaded).
    """
    ref_t = QboPlotter(QboPlotParameter.for_variable("ta", period=period), qbo_dir).load(ref).qbo
    lo, hi = min(band_hpa), max(band_hpa)
    fig, ax = plt.subplots(figsize=figsize)
    for s in [*sources, ref]:
        try:
            r = t_shear_corr_profile(qbo_dir, s, period, ref_t)
        except Exception as e:
            if not skip_bad:
                raise
            print(f"[SKIP] {s}: {e}")
            continue
        score = float(r.sel(plev=slice(lo, hi)).mean())
        ax.plot(r.values, r.plev.values, label=f"{s} (mean r={score:.2f})", **_line_kwargs(s, style_map, ref))
    ax.axvline(0.0, color="k", linewidth=0.3)
    ax.axhspan(lo, hi, alpha=0.12)
    ax.set_title(fr"({string.ascii_lowercase[fig_id]}) Corr[$T^\prime$, $\partial U/\partial \ln p$] "
                 f"(mean {lo:g}–{hi:g} hPa)", fontsize=font_scale * 1.05, loc="left")
    ax.set_xlabel("Correlation over time", fontsize=font_scale)
    ax.set_ylabel("hPa", fontsize=font_scale)
    ax.set_xlim(*xlim)
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=font_scale * 0.85, frameon=False)
    return fig


def plot_t_shear_scatter(qbo_dir, ref: str, sources, *, period: str = PERIOD, layer_hpa=(20.0, 50.0),
                         ncols=3, figsize=None, font_scale=12, cloud_color="0.6", mean_color="tab:red",
                         fig_id=0, skip_bad=True) -> plt.Figure:
    """T' vs dU/dln p in a layer: grey cloud of (month, level) points and the layer-mean
    series (one point per month), whose correlation and slope are annotated (the
    honest sample size is the number of months)."""
    ref_t = QboPlotter(QboPlotParameter.for_variable("ta", period=period), qbo_dir).load(ref).qbo
    names = [ref, *sources]
    nrows = int(math.ceil(len(names) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize or (5.2 * ncols, 4.2 * nrows), squeeze=False)
    lo, hi = min(layer_hpa), max(layer_hpa)
    for i, s in enumerate(names):
        ax = axes.flat[i]
        try:
            T, S = _t_and_shear(qbo_dir, s, period, ref_t)
        except Exception as e:
            if not skip_bad:
                raise
            print(f"[SKIP] {s}: {e}")
            ax.axis("off")
            continue
        Tl, Sl = T.sel(plev=slice(lo, hi)), S.sel(plev=slice(lo, hi))
        ax.scatter(Sl.values.ravel(), Tl.values.ravel(), s=12, alpha=0.45, color=cloud_color)
        xm, ym = Sl.mean("plev").values, Tl.mean("plev").values
        ax.scatter(xm, ym, s=28, color=mean_color, edgecolor="k", linewidth=0.4, zorder=3)
        if xm.size >= 2:
            a, b = np.polyfit(xm, ym, 1)
            xx = np.linspace(xm.min(), xm.max(), 50)
            ax.plot(xx, a * xx + b, color=mean_color, linewidth=2)
            ax.text(0.04, 0.04, f"r={_pcorr(xm, ym, 2):.2f}\nslope={a:.2g}", transform=ax.transAxes,
                    fontsize=font_scale * 0.9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
        ax.set_title(f"({string.ascii_lowercase[fig_id + i]}) {s}", fontsize=font_scale, loc="left")
        ax.set_title(f"{lo:g}–{hi:g} hPa", fontsize=font_scale, loc="right")
        ax.set_xlabel(r"$\partial U/\partial \ln p$ (m s$^{-1}$)", fontsize=font_scale)
        if i % ncols == 0:
            ax.set_ylabel("T′ (K)", fontsize=font_scale)
        ax.grid(True, alpha=0.2)
    for j in range(len(names), nrows * ncols):
        axes.flat[j].axis("off")
    fig.tight_layout()
    return fig
