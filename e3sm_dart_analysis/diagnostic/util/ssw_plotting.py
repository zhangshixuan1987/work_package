"""Figures for the SSW indices written by ``util.ssw_diagnostics``.

Every figure uses two reanalyses (``refs``, ERA5 and ERAI by default): lines are
the mean over the references and shaded bands / whiskers their half-range,
so differences smaller than the reanalysis disagreement are visible as such.

* ``plot_bias_timeseries``: bias of U(55-65N), Z and T (60-90N) at two levels;
* ``plot_bias_and_running_rmse``: bias and running RMSE of one index;
* ``plot_rmse_attribution``: RMSE change vs CTRL by initial-state component,
  per verification window (SSW season or lead months), with the noise floor of a
  difference of two single runs (``util.hind_factorial.pair_noise``: t * sqrt(2) *
  the interaction residual of the additive atm + lnd fit of the RMSEs); only the
  bars whose two runs are among ``experiments`` are drawn, and the noise floor
  needs the full factorial (NaN otherwise);
* ``winter_table`` / ``plot_pseudo_ensemble``: beyond a few months the runs are
  effectively independent realizations of the forced (AMIP) model, so the
  8 runs are a pseudo-ensemble of the 2012-13 vortex; with the QBO wind at
  ``qbo_level`` in the same winter (Holton-Tan: an easterly lower-stratospheric
  QBO favours a weak vortex), the plausible pathway from the initial state;
* ``plot_evolution``: raw index evolution with both references and each
  experiment's correlation with them;
* ``plot_rmse_profiles``: RMSE profiles over a window.
"""
from __future__ import annotations

import string
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from configs.hind_experiment_config import STYLE_MAP
from util.hind_factorial import pair_noise
from util.ssw_diagnostics import lead_key, ssw_path

INDEX_META = {
    "u_55_65N": dict(label="U", units=r"m s$^{-1}$", region="55–65°N", scale=1.0),
    "u_vortex": dict(label="U", units=r"m s$^{-1}$", region="55–75°N", scale=1.0),
    "z_pcap": dict(label="Z", units="dam", region="60–90°N", scale=0.1),
    "t_pcap": dict(label="T", units="K", region="60–90°N", scale=1.0),
}

# (group, label, minuend, subtrahend) of the attribution bars, in drawing order
ATTRIBUTION = (
    ("Total vs CTRL", "ERA5(atm)", "ERA5BL", "CTRL"),
    ("Total vs CTRL", "ERA5(atm+lnd)", "ERA5OP", "CTRL"),
    ("Total vs CTRL", "DART(atm)", "DARTBL", "CTRL"),
    ("Total vs CTRL", "DART(lnd)", "DARTBA", "CTRL"),
    ("Total vs CTRL", "DART(atm+lnd)", "DARTOP", "CTRL"),
    ("Total vs CTRL", "DART(atm)+ERA5(lnd)", "DARTAE", "CTRL"),
    ("Land effect", "+LND | ERA5(atm)", "ERA5OP", "ERA5BL"),
    ("Land effect", "+LND | DART(atm)", "DARTOP", "DARTBL"),
    ("DA − ERA5", "DART − ERA5", "DARTOP", "ERA5OP"),
)


def _month_dates(time) -> np.ndarray:
    return np.array([np.datetime64(f"{t.year:04d}-{t.month:02d}-15") for t in time.values])


class SswPlotter:
    def __init__(self, ssw_dir, experiments, refs=("ERA5", "ERAI"), *, ssw_window=("2012-11", "2013-03"),
                 style_map=None, font_scale: float = 12):
        self.ssw_dir = Path(ssw_dir)
        self.experiments = list(experiments)
        self.refs = tuple(refs)
        self.attribution = [row for row in ATTRIBUTION if row[2] in self.experiments and row[3] in self.experiments]
        self.ssw_window = ssw_window
        self.style_map = style_map or STYLE_MAP
        self.fs = font_scale

    def _open(self, source: str, ref: str) -> xr.Dataset:
        path = ssw_path(self.ssw_dir, source, ref)
        if not path.exists():
            raise FileNotFoundError(f"Missing SSW diagnostics: {path}")
        return xr.load_dataset(path)

    def _series(self, source: str, name: str, p_hpa: float) -> np.ndarray:
        """``name`` (with ``{ref}``) for each reference at one level -> (nref, time)."""
        out = []
        for ref in self.refs:
            ds = self._open(source, ref)
            out.append(ds[name.format(ref=ref)].sel(plev=p_hpa * 100.0, method="nearest").values)
        return np.array(out, dtype=float)

    def _shade_ssw(self, ax, label_xy=(0.43, 0.02)):
        y0, m0 = map(int, self.ssw_window[0].split("-"))
        y1, m1 = map(int, self.ssw_window[1].split("-"))
        end = np.datetime64(f"{y1 + (m1 == 12):04d}-{m1 % 12 + 1:02d}-01")
        ax.axvspan(np.datetime64(f"{y0:04d}-{m0:02d}-01"), end, color="0.85", alpha=0.6, zorder=0, linewidth=0)
        ax.text(*label_xy, "SSW window", transform=ax.transAxes, fontsize=self.fs * 0.8, color="0.35")

    def _dates(self) -> np.ndarray:
        return _month_dates(self._open(self.experiments[0], self.refs[0]).time)

    def _legend(self, fig, ax, ncol=None):
        h, lab = ax.get_legend_handles_labels()
        fig.legend(h, lab, loc="lower center", ncol=ncol or min(len(lab), 5), fontsize=self.fs * 0.9, frameon=True)

    # ------------------------------------------------------------------
    def plot_bias_timeseries(self, indices=("u_55_65N", "z_pcap", "t_pcap"), plevs_hpa=(10.0, 100.0),
                             band_alpha=0.18, ylims=None, fig_id=0) -> plt.Figure:
        """Bias of each index at each level: mean over references ± half-range."""
        t = self._dates()
        fig, axes = plt.subplots(len(indices), len(plevs_hpa), figsize=(6.0 * len(plevs_hpa), 2.9 * len(indices)),
                                 sharex=True, squeeze=False)
        for i, idx in enumerate(indices):
            meta = INDEX_META[idx]
            for j, p in enumerate(plevs_hpa):
                ax = axes[i, j]
                for exp in self.experiments:
                    y = self._series(exp, f"{idx}_bias_{{ref}}", p) * meta["scale"]
                    c, d = y.mean(0), 0.5 * np.ptp(y, axis=0)
                    (line,) = ax.plot(t, c, label=exp, **self.style_map.get(exp, {}))
                    ax.fill_between(t, c - d, c + d, color=line.get_color(), alpha=band_alpha, linewidth=0)
                if ylims and (idx, p) in ylims:
                    ax.set_ylim(*ylims[(idx, p)])
                self._shade_ssw(ax)
                ax.axhline(0, color="k", linewidth=0.4)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"({string.ascii_lowercase[fig_id + i * len(plevs_hpa) + j]}) "
                             f"{meta['label']} {p:g} hPa bias", loc="left", fontsize=self.fs)
                ax.set_title(meta["region"], loc="right", fontsize=self.fs)
                ax.set_ylabel(f"Bias ({meta['units']})", fontsize=self.fs * 0.95)
                ax.tick_params(axis="x", labelrotation=30)
        fig.text(0.01, 0.01, f"Lines: mean bias vs {' and '.join(self.refs)}; bands: half-range across references",
                 fontsize=self.fs * 0.75)
        self._legend(fig, axes[0, 0], ncol=len(self.experiments))
        fig.tight_layout(rect=(0.02, 0.06, 0.98, 0.98))
        return fig

    def plot_bias_and_running_rmse(self, index="u_55_65N", plevs_hpa=(10.0,), window_months=3,
                                   band_alpha=0.18, fig_id=0) -> plt.Figure:
        """Bias (top) and centred ``window_months`` running RMSE (bottom) per level."""
        meta = INDEX_META[index]
        t = self._dates()
        fig, axes = plt.subplots(2, len(plevs_hpa), figsize=(7.0 * len(plevs_hpa), 7.5), sharex=True, squeeze=False)
        for j, p in enumerate(plevs_hpa):
            for exp in self.experiments:
                b = self._series(exp, f"{index}_bias_{{ref}}", p) * meta["scale"]
                st = self.style_map.get(exp, {})
                (line,) = axes[0, j].plot(t, b.mean(0), label=exp, **st)
                axes[0, j].fill_between(t, b.min(0), b.max(0), color=line.get_color(), alpha=band_alpha, linewidth=0)
                r = np.sqrt(xr.DataArray(b ** 2, dims=("ref", "time")).rolling(time=window_months, center=True)
                            .mean().values)
                axes[1, j].plot(t, r.mean(0), **st)
                axes[1, j].fill_between(t, r.min(0), r.max(0), color=line.get_color(), alpha=band_alpha, linewidth=0)
            for k, (ax, what) in enumerate(((axes[0, j], "bias"), (axes[1, j], f"{window_months}-month running RMSE"))):
                self._shade_ssw(ax)
                ax.grid(True, alpha=0.3)
                ax.set_title(f"({string.ascii_lowercase[fig_id + 2 * j + k]}) {meta['label']} {p:g} hPa {what}",
                             loc="left", fontsize=self.fs)
                ax.set_title(meta["region"], loc="right", fontsize=self.fs)
                ax.set_ylabel(f"{what.split()[-1].upper() if k else 'Bias'} ({meta['units']})", fontsize=self.fs * 0.95)
                ax.tick_params(axis="x", labelrotation=30)
            axes[0, j].axhline(0, color="k", linewidth=0.4)
        self._legend(fig, axes[0, 0])
        fig.tight_layout(rect=(0.02, 0.07, 0.98, 0.98))
        return fig

    def rmse(self, source: str, ref: str, index: str, p_hpa: float, window) -> float:
        """RMSE over ``window``: ``"ssw"``, ``"winter"``, ``"all"`` or a lead window (L0, L1)."""
        ds = self._open(source, ref)
        if isinstance(window, str):
            name = f"{index}_rmse_{window}_{ref}"
            return float(ds[name].sel(plev=p_hpa * 100.0, method="nearest"))
        name = f"{index}_rmse_{lead_key(window)}_{ref}"
        if name in ds:
            return float(ds[name].sel(plev=p_hpa * 100.0, method="nearest"))
        bias = ds[f"{index}_bias_{ref}"].sel(plev=p_hpa * 100.0, method="nearest").isel(
            time=slice(window[0], window[1] + 1))
        return float(np.sqrt((bias ** 2).mean()))

    def rmse_attribution_table(self, index="u_55_65N", p_hpa=10.0, windows=((0, 6), "ssw"), as_percent=True):
        """RMSE differences between the experiment pairs of ``ATTRIBUTION`` that were run, per window.

        ``value``: mean over references; ``ref_halfrange``: half-range across them;
        ``noise95``: 95% noise floor of a difference of two single runs (same units).
        ``as_percent`` divides by CTRL's RMSE.
        """
        import pandas as pd
        rows = []
        for window in windows:
            vals, rm = [], []
            for ref in self.refs:
                r = {e: self.rmse(e, ref, index, p_hpa, window) for e in self.experiments}
                den = r["CTRL"] if as_percent else 1.0
                vals.append([(r[a] - r[b]) / den * (100 if as_percent else 1) for _, _, a, b in self.attribution])
                rm.append({e: v / den * (100 if as_percent else 1) for e, v in r.items()})
            vals = np.array(vals)
            noise = float(pair_noise({e: np.mean([x[e] for x in rm]) for e in self.experiments})[0])
            for (group, lab, a, b), m, h in zip(self.attribution, vals.mean(0), 0.5 * np.ptp(vals, axis=0)):
                rows.append(dict(window=str(window), group=group, bar=lab, minuend=a, subtrahend=b, value=m,
                                 ref_halfrange=h, noise95=noise, exceeds_noise=bool(abs(m) > noise)))
        return pd.DataFrame(rows)

    def plot_rmse_attribution(self, index="u_55_65N", p_hpa=10.0, windows=((0, 6), "ssw"), as_percent=True,
                              ylims=None, figsize=(10, 9), fig_id=0, table=None):
        """Bars of ``rmse_attribution_table`` (or ``table``, computed with the same arguments).

        Bars: mean over references (dark when beyond the noise floor); whiskers: half-range
        across references; shading: the 95% noise floor of a two-run difference.
        Returns (figure, table).
        """
        if table is None:
            table = self.rmse_attribution_table(index, p_hpa, windows, as_percent)
        fig, axes = plt.subplots(len(windows), 1, figsize=figsize, squeeze=False)
        for k, window in enumerate(windows):
            ax = axes[k, 0]
            labels = [lab for _, lab, _, _ in self.attribution]
            sub = table[table.window == str(window)].set_index("bar").loc[labels]
            mean, half = sub.value.values, sub.ref_halfrange.values
            x = np.arange(len(labels))
            noise = float(sub.noise95.iloc[0])
            if np.isfinite(noise):
                ax.axhspan(-noise, noise, color="tab:orange", alpha=0.18, zorder=1, linewidth=0)
            ax.bar(x, mean, color=["0.3" if abs(v) > noise else "0.75" for v in mean], edgecolor="0.3", zorder=2)
            ax.errorbar(x, mean, yerr=half, fmt="none", ecolor="k", capsize=4, zorder=3)
            ax.axhline(0, color="k", linewidth=0.5)
            groups = [g for g, _, _, _ in self.attribution]
            for name in dict.fromkeys(groups):
                idx = [i for i, g in enumerate(groups) if g == name]
                if idx[-1] < len(groups) - 1:
                    ax.axvline(idx[-1] + 0.5, color="tab:red", linestyle=":", linewidth=1)
                ax.text((idx[0] + idx[-1]) / 2, 0.97, name, transform=ax.get_xaxis_transform(), ha="center",
                        va="top", fontsize=self.fs * 0.85)
            wlab = "SSW season" if window == "ssw" else (window.capitalize() if isinstance(window, str)
                                                        else f"Lead months {window[0]}–{window[1]}")
            ax.set_title(f"({string.ascii_lowercase[fig_id + k]}) {INDEX_META[index]['label']} {p_hpa:g} hPa "
                         f"RMSE attribution", loc="left", fontsize=self.fs)
            ax.set_title(wlab, loc="right", fontsize=self.fs)
            ax.set_ylabel(r"$\Delta$RMSE / RMSE(CTRL) (%)" if as_percent else r"$\Delta$RMSE", fontsize=self.fs * 0.95)
            if ylims and k < len(ylims) and ylims[k]:
                ax.set_ylim(*ylims[k])
            else:  # headroom for the group labels
                y0, y1 = ax.get_ylim()
                ax.set_ylim(y0, y1 + 0.15 * (y1 - y0))
            ax.set_xticks(x)
            ax.set_xticklabels(labels if k == len(windows) - 1 else [],
                               rotation=30, ha="right")
            ax.grid(True, axis="y", alpha=0.3)
        fig.text(0.01, 0.005, f"Whiskers: half-range across {' and '.join(self.refs)}. Shading: 95% noise floor of a "
                 "difference of two single runs (interaction residual of the atm + lnd fit); dark bars exceed it.",
                 fontsize=self.fs * 0.75)
        fig.tight_layout()
        return fig, table

    def plot_evolution(self, p_hpa=10.0, indices=("u_55_65N", "z_pcap", "t_pcap"), figsize=(10, 9), fig_id=0):
        """Raw index evolution with both references; table of correlations with each reference."""
        t = self._dates()
        fig, axes = plt.subplots(len(indices), 1, figsize=figsize, sharex=True, squeeze=False)
        ref_styles = {self.refs[0]: dict(color="k", linestyle="-", linewidth=1.6),
                      self.refs[-1]: dict(color="k", linestyle="--", linewidth=1.6)}
        for i, idx in enumerate(indices):
            ax = axes[i, 0]
            meta = INDEX_META[idx]
            refs = {r: self._open(r, r)[idx].sel(plev=p_hpa * 100, method="nearest").values * meta["scale"]
                    for r in self.refs}
            lines = [f"{'corr':<7}" + "".join(f"{r:>6}" for r in self.refs)]
            for exp in self.experiments:
                y = self._open(exp, self.refs[0])[idx].sel(plev=p_hpa * 100, method="nearest").values * meta["scale"]
                ax.plot(t, y, label=exp, **self.style_map.get(exp, {}))
                lines.append(f"{exp:<7}" + "".join(f"{np.corrcoef(y, refs[r])[0, 1]:6.2f}" for r in self.refs))
            for r, y in refs.items():
                ax.plot(t, y, label=r, **ref_styles.get(r, {}))
            ax.text(1.01, 1.0, "\n".join(lines), transform=ax.transAxes, va="top", family="monospace",
                    fontsize=self.fs * 0.7)
            self._shade_ssw(ax, (0.46, 0.88))
            ax.grid(True, alpha=0.3)
            ax.set_ylabel(f"{meta['label']} ({meta['units']})", fontsize=self.fs * 0.95)
            ax.set_title(f"({string.ascii_lowercase[fig_id + i]}) {meta['label']} {p_hpa:g} hPa evolution",
                         loc="left", fontsize=self.fs)
            ax.set_title(meta["region"], loc="right", fontsize=self.fs)
        axes[-1, 0].tick_params(axis="x", labelrotation=30)
        self._legend(fig, axes[0, 0])
        fig.tight_layout(rect=(0.02, 0.07, 0.82, 0.98))
        return fig

    def plot_rmse_profiles(self, window="winter", indices=("u_55_65N", "z_pcap", "t_pcap"), figsize=(13, 5.2)):
        """RMSE(p) over ``window`` (mean over references, half-range shaded)."""
        fig, axes = plt.subplots(1, len(indices), figsize=figsize, sharey=True)
        for ax, idx in zip(np.atleast_1d(axes), indices):
            meta = INDEX_META[idx]
            for exp in self.experiments:
                prof = np.array([self._open(exp, r)[f"{idx}_rmse_{window}_{r}"].values for r in self.refs]) * meta["scale"]
                p = self._open(exp, self.refs[0]).plev.values / 100.0
                (line,) = ax.plot(prof.mean(0), p, label=exp, **self.style_map.get(exp, {}))
                ax.fill_betweenx(p, prof.min(0), prof.max(0), color=line.get_color(), alpha=0.15, linewidth=0)
            ax.set_yscale("log")
            ax.set_ylim(1000, 1)
            ax.grid(True, alpha=0.3)
            ax.set_title(f"{meta['label']} {meta['region']} RMSE ({window})", fontsize=self.fs)
            ax.set_xlabel(f"RMSE ({meta['units']})", fontsize=self.fs * 0.95)
        np.atleast_1d(axes)[0].set_ylabel("hPa", fontsize=self.fs * 0.95)
        self._legend(fig, np.atleast_1d(axes)[0])
        fig.tight_layout(rect=(0, 0.1, 1, 1))
        return fig

    # ------------------------------------------------------------------ pseudo-ensemble
    def winter_table(self, qbo_dir, *, index="u_55_65N", p_hpa=10.0, winter=("2012-12", "2013-02"),
                     qbo_level=50.0, period=None):
        """Per source (references and experiments): winter-mean and SSW-window minimum of
        ``index`` at ``p_hpa`` (with its month), months with easterly mean wind in the SSW window,
        and the winter-mean equatorial QBO wind at ``qbo_level`` (from the ``ua`` sections of
        ``util.qbo_diagnostics`` in ``qbo_dir``; NaN when missing)."""
        import pandas as pd
        from configs.hind_experiment_config import PERIOD
        from util.qbo_diagnostics import qbo_section_path
        key = lambda t: t.year * 100 + t.month
        ym = lambda s: int(s[:4]) * 100 + int(s[5:7])
        w0, w1 = ym(winter[0]), ym(winter[1])
        s0, s1 = ym(self.ssw_window[0]), ym(self.ssw_window[1])
        rows = {}
        for src in [*self.refs, *self.experiments]:
            ds = self._open(src, src if src in self.refs else self.refs[0])
            k = np.array([key(t) for t in ds.time.values])
            u = ds[index].sel(plev=p_hpa * 100.0, method="nearest").values * INDEX_META[index]["scale"]
            ssw = (k >= s0) & (k <= s1)
            imin = int(np.flatnonzero(ssw)[np.argmin(u[ssw])])
            row = dict(winter_mean=float(u[(k >= w0) & (k <= w1)].mean()), ssw_min=float(u[imin]),
                       ssw_min_month=str(ds.time.values[imin])[:7], easterly_months=int((u[ssw] < 0).sum()))
            path = qbo_section_path(qbo_dir, src, "ua", period or PERIOD)
            if path.exists():
                q = xr.load_dataset(path).qbo.sel(plev=qbo_level)
                kq = np.array([key(t) for t in q.time.values])
                row[f"qbo{qbo_level:g}_winter"] = float(q.values[(kq >= w0) & (kq <= w1)].mean())
            else:
                row[f"qbo{qbo_level:g}_winter"] = np.nan
            rows[src] = row
        tab = pd.DataFrame(rows).T.infer_objects()
        tab.index.name = "source"
        tab.attrs = dict(index=index, p_hpa=p_hpa, winter=f"{winter[0]}..{winter[1]}", qbo_level=qbo_level)
        return tab

    def plot_pseudo_ensemble(self, table, *, index="u_55_65N", p_hpa=10.0, qbo_level=50.0, figsize=(16, 5.2),
                             fig_id=0):
        """(a) ``index`` of the 8 runs (members thin, range shaded) and the references; (b) winter-mean
        vortex wind vs winter-mean QBO wind at ``qbo_level`` (Holton-Tan)."""
        meta = INDEX_META[index]
        t = self._dates()
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=figsize, gridspec_kw=dict(width_ratios=(1.8, 1)))
        Y = []
        for exp in self.experiments:
            y = self._open(exp, self.refs[0])[index].sel(plev=p_hpa * 100, method="nearest").values * meta["scale"]
            Y.append(y)
            ax.plot(t, y, color=self.style_map.get(exp, {}).get("color", "0.5"), linewidth=0.9, alpha=0.8)
        Y = np.array(Y)
        ax.fill_between(t, Y.min(0), Y.max(0), color="0.6", alpha=0.25, linewidth=0, label="range of the 8 runs")
        ax.plot(t, Y.mean(0), color="0.3", linewidth=2.2, label="mean of the 8 runs")
        for r, ls in zip(self.refs, ("-", "--")):
            y = self._open(r, r)[index].sel(plev=p_hpa * 100, method="nearest").values * meta["scale"]
            ax.plot(t, y, color="k", linestyle=ls, linewidth=2, label=r)
        ax.axhline(0, color="k", linewidth=0.5)
        self._shade_ssw(ax, (0.46, 0.02))
        ax.set_ylabel(f"{meta['label']} ({meta['units']})", fontsize=self.fs * 0.95)
        ax.set_title(f"({string.ascii_lowercase[fig_id]}) {meta['label']} {p_hpa:g} hPa {meta['region']}: "
                     "pseudo-ensemble", loc="left", fontsize=self.fs)
        ax.legend(fontsize=self.fs * 0.75, frameon=False, loc="upper left")
        ax.tick_params(axis="x", labelrotation=30)
        qcol = f"qbo{qbo_level:g}_winter"
        for src in table.index:
            ref = src in self.refs
            ax2.scatter(table.loc[src, qcol], table.loc[src, "winter_mean"], s=90 if ref else 60,
                        marker="*" if ref else "o", color="k" if ref else self.style_map.get(src, {}).get("color", "0.5"),
                        edgecolor="k", linewidth=0.5, zorder=3)
            ax2.annotate(src, (table.loc[src, qcol], table.loc[src, "winter_mean"]), fontsize=self.fs * 0.65,
                         xytext=(4, 3), textcoords="offset points")
        ax2.axvline(0, color="k", linewidth=0.5)
        ax2.set_xlabel(f"QBO U {qbo_level:g} hPa, 5°S–5°N, {table.attrs.get('winter', 'winter')} (m s$^{{-1}}$)",
                       fontsize=self.fs * 0.9)
        ax2.set_ylabel(f"Winter-mean {meta['label']} {p_hpa:g} hPa ({meta['units']})", fontsize=self.fs * 0.9)
        ax2.set_title(f"({string.ascii_lowercase[fig_id + 1]}) Holton–Tan: vortex vs QBO phase", loc="left",
                      fontsize=self.fs)
        for a in (ax, ax2):
            a.grid(True, alpha=0.3)
        fig.tight_layout()
        return fig
