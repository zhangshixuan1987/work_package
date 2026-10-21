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
* vertical coherence: corr_t(X(p), X(p0));
* memory and timing (``qbo_timing_table``, ``plot_qbo_timing``): months of
  memory of the initial QBO phase, onset months of westerly/easterly winds at
  each level and their error, and the descent rate of the westerly shear zone.
  These separate the initialized memory from the model's own (too fast) QBO,
  which a single RMSE/PCC over the record mixes.

Sections can be plotted raw (``field="qbo"``) or as anomalies
(``field="qbo_anom"``) from each source's own monthly climatology
(``own_clim``), the reference climatology (``ref_clim``, comparable across
experiments), or the first month (``first_month``). Temperature defaults to
``ref_clim`` anomalies.

``plot_hovmoller_ua_ta`` draws the U and T sections in one figure (one row each)
against a single reference; ``plot_amplitude_ua_ta``, ``plot_skill_bars_ua_ta``,
``plot_persistence_ua_ta`` and ``plot_coherence_ua_ta`` do the same for the
other metrics (T scored as anomalies, see ``QboPlotter.to_field``).

``plot_t_shear_corr_profiles`` and ``plot_t_shear_scatter`` diagnose the
thermal-wind coupling of the QBO: corr_t(T', dU/dln p) by level and scatter
in a layer.
"""
from __future__ import annotations

import math
import string
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


_PRESSURE_TICKS = (1, 2, 5, 10, 20, 30, 50, 70, 100)


def _pressure_axis(ax, ylim=(100, 1), labelsize=None) -> None:
    """Log pressure axis, surface down, labelled in plain hPa (1, 2, 5, 10, ...) instead of 10^n."""
    ticks = [t for t in _PRESSURE_TICKS if min(ylim) <= t <= max(ylim)]
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.yaxis.set_minor_locator(plt.NullLocator())
    if labelsize is not None:
        ax.tick_params(axis="y", labelsize=labelsize)


def _errorbar_text(errorbar: str) -> str:
    return "95% block-bootstrap CI" if errorbar == "bootstrap" else "half-range across references"


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
        starts = np.arange(max(1, nt - block_len + 1))   # moving (overlapping) blocks
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
            _pressure_axis(ax, labelsize=fs["tick"])
            if c == 0:
                ax.set_ylabel("Pressure (hPa)", fontsize=fs["label"])
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
                       fig_id=0, skip_bad=True, ax=None, legend=True) -> plt.Figure:
        """Amplitude profiles of every source with the reference on top (into ``ax`` if given)."""
        ref_d, tests = self.load_set(ref, sources, amplitude_method, skip_bad)
        fig, ax = plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax)
        for d in [*tests, ref_d]:
            ax.plot(d.amplitude.values, d.amplitude.plev.values, label=d.name,
                    **_line_kwargs(d.name, style_map, ref))
        p = self.parameter
        ax.set_title(f"({string.ascii_lowercase[fig_id]}) {p.main_title}", fontsize=font_scale * 1.05, loc="left")
        ax.set_xlabel(f"{amplitude_method.upper()} amplitude ({p.var_units})", fontsize=font_scale)
        ax.set_ylabel("Pressure (hPa)", fontsize=font_scale)
        _pressure_axis(ax, ylim)
        if xlim is not None:
            ax.set_xlim(*xlim)
        ax.grid(True, which="both", alpha=0.25)
        ax.tick_params(labelsize=font_scale * 0.95)
        if legend:
            ax.legend(**(dict(loc="center left", bbox_to_anchor=(1.02, 0.5)) if legend_outside else dict(loc="best")),
                      fontsize=font_scale * 0.9, frameon=False)
        return fig

    def skill_table(self, refs, sources, *, n_boot=500, block_len=6, seed=0, skip_bad=True):
        """RMSE/PCC of every source averaged over ``refs`` (list) with CI and reference half-range.

        Scored on the plotted field (``to_field``; anomalies take the climatology of ``refs[0]``).
        """
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
            tf, rfs = self.to_field(t.qbo, ref_ds[0].qbo), [self.to_field(r.qbo, ref_ds[0].qbo) for r in ref_ds]
            per = np.array([self.rmse_and_pcc(tf, r) for r in rfs])
            rmse, pcc, rmse_ci, pcc_ci = self.block_bootstrap_skill(
                tf, rfs, n_boot=n_boot, block_len=block_len, seed=seed)
            rows.append(dict(source=s, rmse=rmse, pcc=pcc, rmse_lo=rmse_ci[0], rmse_hi=rmse_ci[1],
                             pcc_lo=pcc_ci[0], pcc_hi=pcc_ci[1],
                             rmse_ref_halfrange=0.5 * np.ptp(per[:, 0]), pcc_ref_halfrange=0.5 * np.ptp(per[:, 1])))
        return pd.DataFrame(rows)

    def plot_skill_bars(self, refs, sources, *, errorbar: Literal["bootstrap", "ref_spread"] = "bootstrap",
                        sort_by: Literal["input", "rmse", "pcc"] = "rmse", baseline: str = "CTRL",
                        n_boot=500, block_len=6, seed=0, figsize=(11.0, 4.6), font_scale=12, nrows=1, ncols=2,
                        rmse_ylim=(0, 20), pcc_ylim=(-0.8, 0.8), fig_id=0, skip_bad=True, table=None, axes=None):
        """RMSE and PCC bars vs one reference or averaged over several.

        Error bars: 95% block-bootstrap CI (``bootstrap``) or half the range across
        references (``ref_spread``, needs >= 2 references). RMSE bars are annotated
        with the % change vs ``baseline``; with two references their own
        disagreement is shown as a floor. ``table``: a precomputed ``skill_table``
        for the same ``refs`` (``sources`` and the bootstrap options are then unused).
        ``axes``: an (RMSE, PCC) pair to draw into, without the footnote and layout.
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
        own = axes is None
        if own:
            fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
        ax, ax2 = np.ravel(axes)[:2]
        fig = ax.figure
        tick_fs = font_scale * 0.95
        var = self.parameter.var_label
        ref_txt = refs[0] if len(refs) == 1 else f"avg over {', '.join(refs)}"
        ax.bar(x, tab.rmse, yerr=rmse_err, color="0.6", capsize=3)
        ax.set_title(f"({string.ascii_lowercase[fig_id]}) {var} RMSE vs {ref_txt}", fontsize=font_scale * 1.05, loc="left")
        ax.set_ylabel(f"RMSE ({self.parameter.var_units})", fontsize=font_scale)
        if rmse_ylim:
            ax.set_ylim(*rmse_ylim)
        else:  # room above the bars for the % labels and the reference-floor note
            ax.set_ylim(0, 1.35 * float(np.nanmax(tab.rmse.values + np.atleast_2d(rmse_err)[-1])))
        if baseline in tab.source.values:
            b = float(tab.loc[tab.source == baseline, "rmse"].iloc[0])
            dy = 0.06 * np.diff(ax.get_ylim())[0]
            for xi, v in zip(x, tab.rmse):
                ax.text(xi, v + dy, f"{100 * (v - b) / b:+.0f}%", ha="center", va="bottom", fontsize=tick_fs * 0.85)
        ax2.bar(x, tab.pcc, yerr=pcc_err, color="0.6", capsize=3)
        ax2.axhline(0.0, linewidth=0.8, alpha=0.6)
        ax2.set_title(f"({string.ascii_lowercase[fig_id + 1]}) {var} PCC vs {ref_txt}", fontsize=font_scale * 1.05, loc="left")
        ax2.set_ylabel("PCC", fontsize=font_scale)
        if pcc_ylim:
            ax2.set_ylim(*pcc_ylim)
        if len(refs) == 2:
            r0, r1 = (self.load(r).qbo for r in refs)
            rmse_b, pcc_b = self.rmse_and_pcc(self.to_field(r1, r0), self.to_field(r0, r0))
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
        if own:
            fig.text(0.99, 0.01, f"Error bars: {_errorbar_text(errorbar)}", ha="right", fontsize=tick_fs * 0.8)
            fig.tight_layout()
        return fig, tab

    def plot_lagcorr_persistence(self, ref: str, sources, *, init_months=3, max_lag=18, style_map=None,
                                 figsize=(9.5, 8.0), font_scale=12, nrows=2, ncols=1, ylim=(-1.0, 1.2),
                                 fig_id=0, skip_bad=True, axes=None) -> plt.Figure:
        """(a) corr_p(initial profile, profile at lag L); (b) corr_p(source, reference) at lag L.

        Profiles of the plotted field (``to_field``: T anomalies, so the mean vertical structure
        does not dominate the correlation). ``axes``: a pair to draw into, without legend and layout.
        """
        ref_d, tests = self.load_set(ref, sources, skip_bad=skip_bad)
        field = {d.name: self.to_field(d.qbo, ref_d.qbo).transpose("time", "plev").values for d in [ref_d, *tests]}
        ref_arr = field[ref_d.name]
        lmax = min(max_lag, ref_arr.shape[0] - 1) if max_lag is not None else ref_arr.shape[0] - 1
        lags = np.arange(lmax + 1)
        series = [ref_d, *tests]
        own = axes is None
        if own:
            fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, squeeze=False)
        ax1, ax2 = np.ravel(axes)[:2]
        fig = ax1.figure
        var = self.parameter.var_label
        for d in series:
            arr = field[d.name]
            x0 = np.nanmean(arr[:init_months], axis=0)
            kw = _line_kwargs(d.name, style_map, ref)
            ax1.plot(lags, [_pcorr(x0, arr[lag]) for lag in lags], label=d.name, **kw)
            ax2.plot(lags, [_pcorr(arr[lag], ref_arr[lag]) for lag in lags], label=d.name, **kw)
        for k, (ax, title) in enumerate(((ax1, f"{var} persistence vs initial {init_months}-month profile"),
                                         (ax2, f"{var} correlation with {ref}"))):
            ax.axhline(0.0, linewidth=0.8, alpha=0.6)
            ax.set_xlim(0, lags[-1])
            ax.set_ylim(*ylim)
            ax.set_title(f"({string.ascii_lowercase[fig_id + k]}) {title}", fontsize=font_scale * 1.05, loc="left")
            ax.set_ylabel("Correlation over levels", fontsize=font_scale)
            ax.tick_params(labelsize=font_scale * 0.95)
            ax.grid(True, alpha=0.25)
        ax2.set_xlabel("Lag (months)", fontsize=font_scale)
        if not own:
            return fig
        if ncols > 1:
            ax1.set_xlabel("Lag (months)", fontsize=font_scale)
        h, lab = ax1.get_legend_handles_labels()
        fig.legend(h, lab, loc="center left", bbox_to_anchor=(0.88, 0.5), fontsize=font_scale * 0.9, frameon=False)
        fig.subplots_adjust(left=0.09, right=0.86, bottom=0.10, top=0.92, hspace=0.25, wspace=0.28)
        return fig

    def plot_vertical_coherence(self, ref: str, sources, *, ref_level_hpa=30.0, field: Field | None = None,
                                anom_mode: AnomMode | None = None, style_map=None, figsize=(6.8, 7.8),
                                font_scale=12, ylim=(100, 1), fig_id=0, skip_bad=True, ax=None,
                                legend=True) -> plt.Figure:
        """r(p) = corr_t(X(t, p), X(t, p0)) for the reference and every source (into ``ax`` if given)."""
        ref_d, tests = self.load_set(ref, sources, skip_bad=skip_bad)
        fig, ax = plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax)
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
        ax.set_ylabel("Pressure (hPa)", fontsize=font_scale)
        ax.set_xlim(-1, 1)
        _pressure_axis(ax, ylim)
        ax.grid(True, which="both", alpha=0.25)
        ax.tick_params(labelsize=font_scale * 0.95)
        if legend:
            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=font_scale * 0.9, frameon=False)
        return fig


# ---------------------------------------------------------------------------
# U and T in one figure
# ---------------------------------------------------------------------------
def ua_ta_plotters(qbo_dir, *, period: str = PERIOD, ta_field: Field = "qbo_anom",
                   ta_anom_mode: AnomMode = "ref_clim") -> dict[str, QboPlotter]:
    """{"ua": raw U plotter, "ta": T plotter drawing ``ta_field``/``ta_anom_mode``}."""
    return {v: QboPlotter(QboPlotParameter.for_variable(v, period=period, default_field=f,
                                                        default_anom_mode=ta_anom_mode), qbo_dir)
            for v, f in (("ua", "qbo"), ("ta", ta_field))}


def _legend_right(fig, ax, font_scale) -> None:
    h, lab = ax.get_legend_handles_labels()
    fig.legend(h, lab, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=font_scale * 0.9, frameon=False)


def plot_amplitude_ua_ta(qbo_dir, ref: str, sources, *, amplitude_method: AmpMethod = "rms",
                         xlim=((0, 30), None), figsize=(12.0, 7.0), font_scale=12, **kw) -> plt.Figure:
    """(a) U and (b) T amplitude profiles (``kw``: ``ua_ta_plotters`` options, ``skip_bad``)."""
    skip_bad = kw.pop("skip_bad", True)
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    for k, (pl, ax, lim) in enumerate(zip(ua_ta_plotters(qbo_dir, **kw).values(), axes, xlim)):
        pl.plot_amplitude(ref, sources, amplitude_method=amplitude_method, xlim=lim, font_scale=font_scale,
                          fig_id=k, skip_bad=skip_bad, ax=ax, legend=False)
    axes[1].set_ylabel("")
    fig.tight_layout()
    _legend_right(fig, axes[0], font_scale)
    return fig


def plot_coherence_ua_ta(qbo_dir, ref: str, sources, *, ref_level_hpa=30.0, figsize=(12.0, 7.0), font_scale=12,
                         **kw) -> plt.Figure:
    """Vertical coherence with ``ref_level_hpa`` of (a) U and (b) T."""
    skip_bad = kw.pop("skip_bad", True)
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    for k, (pl, ax) in enumerate(zip(ua_ta_plotters(qbo_dir, **kw).values(), axes)):
        pl.plot_vertical_coherence(ref, sources, ref_level_hpa=ref_level_hpa, font_scale=font_scale, fig_id=k,
                                   skip_bad=skip_bad, ax=ax, legend=False)
        ax.set_title(f"({string.ascii_lowercase[k]}) {pl.parameter.var_label}: coherence with "
                     f"{ref_level_hpa:g} hPa", fontsize=font_scale * 1.05, loc="left")
    axes[1].set_ylabel("")
    fig.tight_layout()
    _legend_right(fig, axes[0], font_scale)
    return fig


def plot_skill_bars_ua_ta(qbo_dir, refs, sources, *, tables: Mapping[str, pd.DataFrame] | None = None,
                          errorbar: Literal["bootstrap", "ref_spread"] = "bootstrap", n_boot=500, block_len=6,
                          rmse_ylim=((0, 20), None), pcc_ylim=(-0.8, 0.8), figsize=(11.0, 8.4), font_scale=12,
                          **kw) -> plt.Figure:
    """RMSE (left) and PCC (right) bars of U (top) and T (bottom), sources ordered by U RMSE.

    ``tables``: precomputed ``skill_table`` per variable for the same ``refs``.
    """
    skip_bad = kw.pop("skip_bad", True)
    plotters = ua_ta_plotters(qbo_dir, **kw)
    tabs = {v: tables[v] if tables is not None else
            pl.skill_table(refs, sources, n_boot=n_boot, block_len=block_len, skip_bad=skip_bad)
            for v, pl in plotters.items()}
    order = list(tabs["ua"].sort_values("rmse").source)
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    for r, (v, pl) in enumerate(plotters.items()):
        tab = tabs[v].set_index("source").loc[[s for s in order if s in set(tabs[v].source)]].reset_index()
        pl.plot_skill_bars(refs, sources, table=tab, sort_by="input", errorbar=errorbar, rmse_ylim=rmse_ylim[r],
                           pcc_ylim=pcc_ylim, font_scale=font_scale, fig_id=2 * r, axes=axes[r])
    for ax in axes[0]:
        ax.tick_params(axis="x", labelbottom=False)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.text(0.99, 0.005, f"Error bars: {_errorbar_text(errorbar)}", ha="right", fontsize=font_scale * 0.75)
    return fig


def plot_persistence_ua_ta(qbo_dir, ref: str, sources, *, init_months=3, max_lag=18, figsize=(14.0, 8.0),
                           font_scale=12, **kw) -> plt.Figure:
    """Persistence of the initial profile (top) and correlation with ``ref`` (bottom) for U (left) and T (right)."""
    skip_bad = kw.pop("skip_bad", True)
    fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=True, sharey=True)
    for c, pl in enumerate(ua_ta_plotters(qbo_dir, **kw).values()):
        pl.plot_lagcorr_persistence(ref, sources, init_months=init_months, max_lag=max_lag, font_scale=font_scale,
                                    fig_id=2 * c, skip_bad=skip_bad, axes=axes[:, c])
    for ax in axes[:, 1]:
        ax.set_ylabel("")
    fig.tight_layout()
    _legend_right(fig, axes[0, 0], font_scale)
    return fig


def plot_hovmoller_ua_ta(qbo_dir, ref: str, sources, *, period: str = PERIOD, ta_field: Field = "qbo_anom",
                         ta_anom_mode: AnomMode = "ref_clim", ylim=(100, 1), w_per_col=4.2, h_per_row=3.4,
                         font_scale=12, time_every=6, cmap="RdBu_r", skip_bad=True) -> plt.Figure:
    """Time-pressure sections of U (top row) and T (bottom row) for ``ref`` and every source.

    T is drawn as ``ta_field``/``ta_anom_mode`` (default: anomalies from the ``ref`` monthly
    climatology); each source panel gives RMSE/PCC vs ``ref`` of the field shown. A source must
    have both variables on the reference grid, so the two rows share their columns.
    """
    plotters = ua_ta_plotters(qbo_dir, period=period, ta_field=ta_field, ta_anom_mode=ta_anom_mode)
    loaded = {v: pl.load_set(ref, sources, skip_bad=skip_bad) for v, pl in plotters.items()}
    names = [ref, *(n for n in sources if all(any(t.name == n for t in loaded[v][1]) for v in plotters))]
    fs = dict(title=round(font_scale), label=round(font_scale * 0.95), tick=round(font_scale * 0.9),
              anno=round(font_scale * 0.75))
    fig, axes = plt.subplots(2, len(names), figsize=(w_per_col * len(names), 2 * h_per_row),
                             sharex=True, sharey=True, squeeze=False, layout="constrained")
    for r, (v, pl) in enumerate(plotters.items()):
        ref_d, tests = loaded[v]
        by_name = {d.name: d for d in [ref_d, *tests]}
        ref_f = pl.to_field(ref_d.qbo, ref_d.qbo)
        levels, ticks = pl._default_levels(diff=False)
        p = pl.parameter
        cs = None
        for c, n in enumerate(names):
            ax = axes[r, c]
            da = ref_f if n == ref else pl.to_field(by_name[n].qbo, ref_d.qbo)
            cs = ax.contourf(np.arange(da.sizes["time"]), da.plev.values, da.transpose("plev", "time").values,
                             levels=levels, cmap=cmap, extend="both")
            ax.set_title(f"({string.ascii_lowercase[r * len(names) + c]}) {n}", fontsize=fs["title"], loc="left")
            _pressure_axis(ax, ylim, labelsize=fs["tick"])
            if n != ref:
                rmse, pcc = pl.rmse_and_pcc(da, ref_f)
                ax.text(0.97, 0.04, f"RMSE {rmse:.2f} {p.var_units}\nPCC {pcc:.2f}", transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=fs["anno"],
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.7", alpha=0.85))
            if r == 1:
                QboPlotter._time_ticks(ax, da.time.values, time_every, fs["tick"])
        anom = pl._field_tag(None, None) != "raw"
        axes[r, 0].set_ylabel(f"{QBO_VARIABLES[v]['long_name'].capitalize()}\nPressure (hPa)",
                              fontsize=fs["label"])
        cbar = fig.colorbar(cs, ax=axes[r, :], ticks=ticks, shrink=0.9, pad=0.01, aspect=15)
        cbar.ax.tick_params(labelsize=fs["tick"])
        cbar.set_label(f"{p.var_label}{' anomaly' if anom else ''} ({p.var_units})", fontsize=fs["label"])
    fig.suptitle(f"QBO 5°S–5°N: U and {'T anomalies' if ta_field == 'qbo_anom' else 'T'} vs {ref}",
                 fontsize=fs["title"])
    return fig


# ---------------------------------------------------------------------------
# Memory and phase timing
# ---------------------------------------------------------------------------
SCALE_HEIGHT_KM = 7.0


def first_crossing(x: np.ndarray, direction: str = "E2W", persist: int = 2, start: float = 0.0) -> float:
    """First month (linearly interpolated) at or after ``start`` at which ``x`` changes sign in
    ``direction`` and keeps the new sign for ``persist`` months; ``E2W``: westerly onset."""
    x = np.asarray(x, dtype=float)
    sgn = 1.0 if direction == "E2W" else -1.0
    for k in range(max(1, int(np.floor(start)) + 1), len(x) - persist + 1):
        a, b = sgn * x[k - 1], sgn * x[k]
        if a < 0 <= b and np.all(sgn * x[k:k + persist] > 0):
            t = float(k - 1 + a / (a - b))
            if t >= start:
                return t
    return np.nan


def descending_onsets(q: xr.DataArray, levels_hpa, direction: str = "E2W", persist: int = 2) -> dict:
    """Onset month of a descending shear zone at each level: the first onset at the top level, then
    at each lower level the first onset after the one above (flickers left by the previous phase
    near the tropopause are skipped). NaN below a level without an onset."""
    out, start = {}, 0.0
    for p in sorted(levels_hpa):
        t = first_crossing(q.sel(plev=p).values, direction, persist, start) if np.isfinite(start) else np.nan
        out[p], start = t, t
    return out


def memory_months(test: np.ndarray, ref: np.ndarray, threshold: float = 0.5, persist: int = 2) -> float:
    """Lead months before the profile correlation with the reference stays below ``threshold``
    for ``persist`` months (``test``/``ref``: (time, plev)); the record length if it never does."""
    r = np.array([_pcorr(test[t], ref[t]) for t in range(test.shape[0])])
    low = r < threshold
    for t in range(len(r) - persist + 1):
        if low[t:t + persist].all():
            return float(t)
    return float(len(r))


def descent_rate(onsets: Mapping[float, float]) -> float:
    """Descent rate (km month⁻¹, positive downward) of a shear zone from its onset month per level (hPa)."""
    p = np.array([k for k, v in onsets.items() if np.isfinite(v)], dtype=float)
    t = np.array([onsets[k] for k in p])
    if p.size < 2:
        return np.nan
    z = SCALE_HEIGHT_KM * np.log(1000.0 / p)
    return float(-np.polyfit(t, z, 1)[0])


def qbo_timing_table(qbo_dir, refs, sources, *, period: str = PERIOD, levels_hpa=(10.0, 20.0, 30.0, 50.0, 70.0),
                     memory_band_hpa=(10.0, 100.0), descent_band_hpa=(10.0, 50.0), memory_threshold: float = 0.5,
                     persist: int = 2):
    """Memory and phase timing of every source (``ua`` sections) against the mean of ``refs``.

    Columns: ``memory_months`` (profile correlation over ``memory_band_hpa``),
    ``W_onset_<p>`` / ``E_onset_<p>`` (months since the first month, following the
    descending zone from the top level, ``descending_onsets``) and their
    errors ``*_err`` (source − reference), ``descent_km_per_month`` of the
    westerly zone over ``descent_band_hpa`` and its error ``descent_err`` (the model zone stalls near 70 hPa, so a fit
    through the lowest levels would hide its fast descent above), and ``abs_W_onset_err_<p>``.
    """
    plotter = QboPlotter(QboPlotParameter.for_variable("ua", period=period), qbo_dir)
    refs = [refs] if isinstance(refs, str) else list(refs)
    lo, hi = memory_band_hpa
    ref_q = sum(plotter.load(r).qbo for r in refs) / len(refs)
    band = ref_q.sel(plev=slice(lo, hi))

    def timing(q):
        row = {}
        for kind, direction in (("W", "E2W"), ("E", "W2E")):
            for p, t in descending_onsets(q, levels_hpa, direction, persist).items():
                row[f"{kind}_onset_{p:g}"] = t
        row["descent_km_per_month"] = descent_rate({p: row[f"W_onset_{p:g}"] for p in levels_hpa
                                                    if descent_band_hpa[0] <= p <= descent_band_hpa[1]})
        return row

    ref_row = timing(ref_q)
    rows = {}
    for s in [*refs, *sources]:
        q = plotter.load(s).qbo
        row = dict(memory_months=memory_months(q.sel(plev=slice(lo, hi)).values, band.values,
                                               memory_threshold, persist), **timing(q))
        for p in levels_hpa:
            for kind in ("W", "E"):
                row[f"{kind}_onset_{p:g}_err"] = row[f"{kind}_onset_{p:g}"] - ref_row[f"{kind}_onset_{p:g}"]
            row[f"abs_W_onset_err_{p:g}"] = abs(row[f"W_onset_{p:g}_err"])
        row["descent_err"] = row["descent_km_per_month"] - ref_row["descent_km_per_month"]
        rows[s] = row
    tab = pd.DataFrame(rows).T
    tab.index.name = "source"
    tab.attrs = dict(reference=" + ".join(refs), memory_threshold=memory_threshold, persist=persist)
    return tab


def plot_qbo_timing(qbo_dir, ref: str, sources, timing, *, period: str = PERIOD, levels_hpa=(30.0, 50.0),
                    style_map=None, font_scale=12, figsize=(16, 5.2)) -> plt.Figure:
    """(a) westerly-onset month by level; (b, c) U at ``levels_hpa`` vs lead month."""
    plotter = QboPlotter(QboPlotParameter.for_variable("ua", period=period), qbo_dir)
    fig, axes = plt.subplots(1, 1 + len(levels_hpa), figsize=figsize,
                             gridspec_kw=dict(width_ratios=(1, *[1.4] * len(levels_hpa))))
    onset_cols = [c for c in timing.columns if c.startswith("W_onset_") and not c.endswith("_err")]
    p = np.array([float(c.split("_")[-1]) for c in onset_cols])
    for s in [*sources, ref]:
        kw = _line_kwargs(s, style_map, ref)
        axes[0].plot(timing.loc[s, onset_cols].astype(float).values, p, marker="o", label=s, **kw)
        q = plotter.load(s).qbo
        for ax, lev in zip(axes[1:], levels_hpa):
            ax.plot(np.arange(q.sizes["time"]), q.sel(plev=lev).values, **kw)
    axes[0].set_yscale("log")
    axes[0].set_ylim(max(p) * 1.1, min(p) / 1.1)
    axes[0].set_yticks(p, [f"{v:g}" for v in p])
    axes[0].minorticks_off()
    axes[0].set_xlabel("Westerly onset (months since start)", fontsize=font_scale)
    axes[0].set_ylabel("Pressure (hPa)", fontsize=font_scale)
    axes[0].set_title("(a) Westerly onset by level", loc="left", fontsize=font_scale)
    for k, (ax, lev) in enumerate(zip(axes[1:], levels_hpa)):
        ax.axhline(0, color="k", linewidth=0.6)
        ax.set_xlabel("Lead month", fontsize=font_scale)
        ax.set_ylabel(r"U (m s$^{-1}$)", fontsize=font_scale)
        ax.set_title(f"({string.ascii_lowercase[k + 1]}) U {lev:g} hPa (5°S–5°N)", loc="left", fontsize=font_scale)
    for ax in axes:
        ax.grid(True, alpha=0.3)
    axes[0].legend(fontsize=font_scale * 0.7, frameon=False, loc="lower right")
    fig.tight_layout()
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
    ax.set_ylabel("Pressure (hPa)", fontsize=font_scale)
    ax.set_xlim(*xlim)
    _pressure_axis(ax, ylim)
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


# ---------------------------------------------------------------------------
# Thermal-wind consistency of the temperature errors (paper figure)
# ---------------------------------------------------------------------------
def thermal_wind_slope(qbo_dir, refs, *, period: str = PERIOD, layer_hpa=(10.0, 70.0)) -> float:
    """Regression slope dT'/d(dU/dln p) (K s/m) of the reanalyses, pooled over the (month, level)
    points of ``layer_hpa`` and averaged over ``refs``; T' from the first reference's climatology."""
    ref_t = QboPlotter(QboPlotParameter.for_variable("ta", period=period), qbo_dir).load(refs[0]).qbo
    lo, hi = min(layer_hpa), max(layer_hpa)
    slopes = []
    for r in refs:
        T, S = _t_and_shear(qbo_dir, r, period, ref_t)
        slopes.append(np.polyfit(S.sel(plev=slice(lo, hi)).values.ravel(),
                                 T.sel(plev=slice(lo, hi)).values.ravel(), 1)[0])
    return float(np.mean(slopes))


def thermal_wind_errors(qbo_dir, ref: str, source: str, slope: float, *, period: str = PERIOD):
    """Errors vs ``ref`` of T and of dU/dln p, and the unbalanced part dT - slope * dS (time, plev)."""
    load = lambda s, v: QboPlotter(QboPlotParameter.for_variable(v, period=period), qbo_dir).load(s).qbo
    dT = load(source, "ta") - load(ref, "ta")
    dS = shear_dlnp(load(source, "ua")) - shear_dlnp(load(ref, "ua"))
    return dT, dS, (dT - slope * dS).rename("dT_unbalanced")


def plot_thermal_wind_consistency(qbo_dir, ref: str, groups: Mapping[str, Any], *, refs=("ERAI", "ERA5"),
                                  sections=("DARTOP", "ERA5OP"), profile_sources=("DARTOP", "ERA5OP", "CTRL"),
                                  profile_months=3, floor: str | None = "ERA5", period: str = PERIOD,
                                  layer_hpa=(10.0, 70.0), ylim=(100.0, 5.0),
                                  colors=("#2a78d6", "#eb6834", "0.55"), linestyles=("-", "--", ":"),
                                  dt_levels=np.arange(-12, 13, 2), shear_contours=(-8, -4, -2, 2, 4, 8),
                                  figsize=(11.5, 8.6), font_scale=11):
    """Are the temperature errors of each initial state consistent with its wind errors?

    Thermal-wind balance on the equator ties T' to dU/dln p with the slope alpha fitted on the
    reanalyses (``thermal_wind_slope``), so the balanced part of a temperature error is
    alpha * dS and dT - alpha * dS is the part that the wind errors do not support.

    (a, b) dT vs ``ref`` (shading) and alpha * dS (contours) for ``sections``;
    (c) RMS over ``layer_hpa`` of dT - alpha * dS vs lead month: line the group mean, shading the
    range over the group's runs (land initial states), ``floor`` vs ``ref`` in black;
    (d) dT (solid) and alpha * dS (dashed) averaged over the first ``profile_months`` months for
    ``profile_sources``, coloured by group.
    ``groups`` maps a label to its run keys; the first ``len(colors)`` groups are drawn.
    Returns (figure, per-run table: corr and slope of dT on dS, unbalanced RMS and mean dT over
    ``layer_hpa`` for months 1-6 and 7-end).
    """
    groups = dict(list(groups.items())[:len(colors)])
    slope = thermal_wind_slope(qbo_dir, list(refs), period=period, layer_hpa=layer_hpa)
    lo, hi = min(layer_hpa), max(layer_hpa)
    fs = font_scale
    names = [*sections, *profile_sources, *(s for v in groups.values() for s in v), *([floor] if floor else [])]
    err = {s: thermal_wind_errors(qbo_dir, ref, s, slope, period=period) for s in dict.fromkeys(names)}
    lead = np.arange(1, next(iter(err.values()))[0].sizes["time"] + 1)
    color_of = {s: c for (g, members), c in zip(groups.items(), colors) for s in members}

    def pressure_axis(ax):
        _pressure_axis(ax, ylim)
        ax.tick_params(labelsize=fs * 0.9)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.36, wspace=0.26)
    # (a, b) sections
    cf = None
    for i, s in enumerate(sections):
        ax = fig.add_subplot(gs[0, i])
        dT, dS, _ = err[s]
        p = dT.plev.values
        cf = ax.contourf(lead, p, dT.transpose("plev", "time").values, levels=dt_levels, cmap="RdBu_r", extend="both")
        cs = ax.contour(lead, p, (slope * dS).transpose("plev", "time").values, levels=list(shear_contours),
                        colors="k", linewidths=0.8)
        ax.clabel(cs, fmt="%g", fontsize=fs * 0.7)
        pressure_axis(ax)
        ax.set_xticks(lead[::3])
        ax.set_title(f"({string.ascii_lowercase[i]}) {s}", fontsize=fs * 1.05, loc="left")
        ax.set_title(r"$\Delta T$ (shading), $\alpha\,\Delta S$ (contours, K)", fontsize=fs * 0.9, loc="right")
        ax.set_xlabel("Lead month", fontsize=fs)
        if i == 0:
            ax.set_ylabel("Pressure (hPa)", fontsize=fs)
    cax = fig.add_axes([0.92, 0.585, 0.012, 0.29])
    cb = fig.colorbar(cf, cax=cax, ticks=dt_levels[::2])
    cb.set_label(f"$\\Delta T$ vs {ref} (K)", fontsize=fs * 0.9)
    cb.ax.tick_params(labelsize=fs * 0.85)

    # (c) unbalanced error vs lead month
    ax = fig.add_subplot(gs[1, 0])
    rms = lambda x: np.sqrt((x.sel(plev=slice(lo, hi)) ** 2).mean("plev")).values
    if floor:
        ax.plot(lead, rms(err[floor][2]), color="k", linewidth=1.6, label=f"{floor} vs {ref}")
    for (g, members), c, ls in zip(groups.items(), colors, linestyles):
        curves = np.array([rms(err[s][2]) for s in members])
        ax.fill_between(lead, curves.min(0), curves.max(0), color=c, alpha=0.22, linewidth=0)
        ax.plot(lead, curves.mean(0), color=c, linestyle=ls, linewidth=2.2, label=f"{g} (n={len(members)})")
    ax.set_xlim(1, lead[-1])
    ax.set_ylim(bottom=0)
    ax.set_xticks(lead[::3])
    ax.set_xlabel("Lead month", fontsize=fs)
    ax.set_ylabel("RMS (K)", fontsize=fs)
    ax.set_title(f"(c) Unbalanced error $\\Delta T-\\alpha\\,\\Delta S$, {lo:g}–{hi:g} hPa",
                 fontsize=fs * 1.05, loc="left")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=fs * 0.82, frameon=False, loc="upper right")
    ax.tick_params(labelsize=fs * 0.9)

    # (d) early-lead profiles of the error and its balanced part
    ax = fig.add_subplot(gs[1, 1])
    early = slice(0, profile_months)
    for s in profile_sources:
        dT, dS, _ = err[s]
        c = color_of.get(s, "k")
        ax.plot(dT.isel(time=early).mean("time"), dT.plev, color=c, linewidth=2.2, marker="o", markersize=5,
                label=rf"{s} $\Delta T$")
        ax.plot(slope * dS.isel(time=early).mean("time"), dS.plev, color=c, linewidth=1.6, linestyle="--",
                label=rf"{s} $\alpha\,\Delta S$")
    ax.axvline(0, color="k", linewidth=0.4)
    pressure_axis(ax)
    ax.set_xlabel("K", fontsize=fs)
    ax.set_title(f"(d) Months 1–{profile_months} mean: $\\Delta T$ (solid) vs $\\alpha\\,\\Delta S$ (dashed)",
                 fontsize=fs * 1.05, loc="left")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=fs * 0.78, frameon=False, loc="upper left")
    fig.text(0.5, 0.015, f"α = {slope:.3f} K s m$^{{-1}}$: regression of T′ on ∂U/∂ln p in "
             f"{', '.join(refs)} ({lo:g}–{hi:g} hPa); ΔS = Δ(∂U/∂ln p); Δ = run − {ref}.",
             ha="center", fontsize=fs * 0.85)

    rows = []
    for g, members in groups.items():
        for s in members:
            dTl, dSl, ul = (z.sel(plev=slice(lo, hi)) for z in err[s])
            x, y = dSl.values.ravel(), dTl.values.ravel()
            rows.append({"group": g, "source": s, "r_dT_dS": _pcorr(x, y, 2), "slope_dT_dS": np.polyfit(x, y, 1)[0],
                         "unbalanced_rms_m1-6": float(np.sqrt((ul.isel(time=slice(0, 6)) ** 2).mean())),
                         f"unbalanced_rms_m7-{lead[-1]}": float(np.sqrt((ul.isel(time=slice(6, None)) ** 2).mean())),
                         "mean_dT_m1-6": float(dTl.isel(time=slice(0, 6)).mean())})
    table = pd.DataFrame(rows).set_index("source")
    table.attrs["alpha"] = slope
    return fig, table
