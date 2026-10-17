"""Figures for the SSW indices written by ``util.ssw_diagnostics``.

Every figure uses two reanalyses (``refs``, ERA5 and ERAI by default): lines are
the mean over the references and shaded bands / whiskers their half-range,
so differences smaller than the reanalysis disagreement are visible as such.

* ``plot_bias_timeseries``: bias of U(55-65N), Z and T (60-90N) at two levels;
* ``plot_bias_and_running_rmse``: bias and running RMSE of one index;
* ``plot_rmse_attribution``: RMSE change vs CTRL by initial-state component,
  per verification window (SSW season or lead months);
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
from util.ssw_diagnostics import lead_key, ssw_path

INDEX_META = {
    "u_55_65N": dict(label="U", units=r"m s$^{-1}$", region="55–65°N", scale=1.0),
    "u_vortex": dict(label="U", units=r"m s$^{-1}$", region="55–75°N", scale=1.0),
    "z_pcap": dict(label="Z", units="dam", region="60–90°N", scale=0.1),
    "t_pcap": dict(label="T", units="K", region="60–90°N", scale=1.0),
}

# (label, minuend, subtrahend) of the attribution bars; groups split by ``|`` markers
ATTRIBUTION = (
    ("ERA5(atm)", "ERA5BL", "CTRL"),
    ("ERA5(atm+lnd)", "ERA5OP", "CTRL"),
    ("DART(atm)", "DARTBL", "CTRL"),
    ("DART(lnd)", "DARTBA", "CTRL"),
    ("DART(atm+lnd)", "DARTOP", "CTRL"),
    ("DART(atm)+ERA5(lnd)", "DARTAE", "CTRL"),
    ("+LND | ERA5(atm)", "ERA5OP", "ERA5BL"),
    ("+LND | DART(atm)", "DARTOP", "DARTBL"),
    ("DART − ERA5", "DARTOP", "ERA5OP"),
)
ATTRIBUTION_GROUPS = (("Total vs CTRL", 0, 5), ("Land effect", 6, 7), ("DA − ERA5", 8, 8))


def _month_dates(time) -> np.ndarray:
    return np.array([np.datetime64(f"{t.year:04d}-{t.month:02d}-15") for t in time.values])


class SswPlotter:
    def __init__(self, ssw_dir, experiments, refs=("ERA5", "ERAI"), *, ssw_window=("2012-11", "2013-03"),
                 style_map=None, font_scale: float = 12):
        self.ssw_dir = Path(ssw_dir)
        self.experiments = list(experiments)
        self.refs = tuple(refs)
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
        """RMSE differences between experiment pairs (``ATTRIBUTION``) per window.

        ``value``: mean over references; ``ref_halfrange``: half-range across them.
        ``as_percent`` divides by CTRL's RMSE.
        """
        import pandas as pd
        rows = []
        for window in windows:
            vals = []
            for ref in self.refs:
                r = {e: self.rmse(e, ref, index, p_hpa, window) for e in {x for _, a, b in ATTRIBUTION for x in (a, b)}}
                den = r["CTRL"] if as_percent else 1.0
                vals.append([(r[a] - r[b]) / den * (100 if as_percent else 1) for _, a, b in ATTRIBUTION])
            vals = np.array(vals)
            for (lab, a, b), m, h in zip(ATTRIBUTION, vals.mean(0), 0.5 * np.ptp(vals, axis=0)):
                rows.append(dict(window=str(window), bar=lab, minuend=a, subtrahend=b, value=m, ref_halfrange=h))
        return pd.DataFrame(rows)

    def plot_rmse_attribution(self, index="u_55_65N", p_hpa=10.0, windows=((0, 6), "ssw"), as_percent=True,
                              ylims=None, figsize=(10, 9), fig_id=0, table=None):
        """Bars of ``rmse_attribution_table`` (or ``table``, computed with the same arguments).

        Bars: mean over references; whiskers: half-range across references.
        Returns (figure, table).
        """
        if table is None:
            table = self.rmse_attribution_table(index, p_hpa, windows, as_percent)
        fig, axes = plt.subplots(len(windows), 1, figsize=figsize, squeeze=False)
        for k, window in enumerate(windows):
            ax = axes[k, 0]
            sub = table[table.window == str(window)].set_index("bar").loc[[lab for lab, _, _ in ATTRIBUTION]]
            mean, half = sub.value.values, sub.ref_halfrange.values
            x = np.arange(len(ATTRIBUTION))
            ax.bar(x, mean, color="0.7", edgecolor="0.3", zorder=2)
            ax.errorbar(x, mean, yerr=half, fmt="none", ecolor="k", capsize=4, zorder=3)
            ax.axhline(0, color="k", linewidth=0.5)
            for _, g0, g1 in ATTRIBUTION_GROUPS[:-1]:
                ax.axvline(g1 + 0.5, color="tab:red", linestyle=":", linewidth=1)
            for name, g0, g1 in ATTRIBUTION_GROUPS:
                ax.text((g0 + g1) / 2, 0.97, name, transform=ax.get_xaxis_transform(), ha="center", va="top",
                        fontsize=self.fs * 0.85)
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
            ax.set_xticklabels([lab for lab, _, _ in ATTRIBUTION] if k == len(windows) - 1 else [],
                               rotation=30, ha="right")
            ax.grid(True, axis="y", alpha=0.3)
        fig.text(0.01, 0.005, f"Whiskers: half-range across {' and '.join(self.refs)}", fontsize=self.fs * 0.75)
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
