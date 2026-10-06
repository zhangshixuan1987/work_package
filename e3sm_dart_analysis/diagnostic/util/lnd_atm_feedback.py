import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import xarray as xr


LagSummary = Dict[int, Dict[str, float]]
FeedbackResults = Dict[str, LagSummary]


class LandAtmosphereFeedbackAnalyzer:
    def __init__(
        self,
        land_mask: Optional[xr.DataArray] = None,
        alpha: float = 0.05,
        n_bootstrap: int = 1000,
        min_samples: int = 30,
        random_seed: Optional[int] = None,
    ):
        if not 0 < alpha < 1:
            raise ValueError("alpha must be between 0 and 1.")
        if n_bootstrap < 1:
            raise ValueError("n_bootstrap must be at least 1.")
        if min_samples < 2:
            raise ValueError("min_samples must be at least 2.")

        self.land_mask = land_mask
        self.alpha = alpha
        self.n_bootstrap = n_bootstrap
        self.min_samples = min_samples
        self.rng = np.random.default_rng(random_seed)

    def _bootstrap_corr(self, x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
        if len(x) < self.min_samples or len(y) < self.min_samples:
            return np.nan, np.nan, np.nan

        idx = np.arange(len(x))
        boot_corrs = []
        for _ in range(self.n_bootstrap):
            sample_idx = self.rng.choice(idx, size=len(idx), replace=True)
            x_sample = x[sample_idx]
            y_sample = y[sample_idx]
            if np.std(x_sample) > 0 and np.std(y_sample) > 0:
                boot_corrs.append(np.corrcoef(x_sample, y_sample)[0, 1])

        if not boot_corrs:
            return np.nan, np.nan, np.nan

        boot_corrs = np.asarray(boot_corrs, dtype=float)
        lower = np.nanpercentile(boot_corrs, 100 * (self.alpha / 2))
        upper = np.nanpercentile(boot_corrs, 100 * (1 - self.alpha / 2))
        return float(np.nanmean(boot_corrs)), float(lower), float(upper)

    @staticmethod
    def _require_vars(dataset: xr.Dataset, variables: Sequence[str]) -> None:
        missing = [var for var in variables if var not in dataset]
        if missing:
            available = ", ".join(sorted(dataset.data_vars))
            raise KeyError(f"Missing variables {missing}. Available variables: {available}")

    @staticmethod
    def _apply_lag(x: xr.DataArray, y: xr.DataArray, lag: int) -> Tuple[xr.DataArray, xr.DataArray]:
        if "time" not in x.dims or "time" not in y.dims:
            raise ValueError("Both variables must include a time dimension.")
        if abs(lag) >= min(x.sizes["time"], y.sizes["time"]):
            return x.isel(time=slice(0, 0)), y.isel(time=slice(0, 0))
        if lag > 0:
            return x.isel(time=slice(0, -lag)), y.isel(time=slice(lag, None))
        if lag < 0:
            return x.isel(time=slice(-lag, None)), y.isel(time=slice(0, lag))
        return x, y

    def _apply_land_mask(self, data: xr.DataArray) -> xr.DataArray:
        return data.where(self.land_mask) if self.land_mask is not None else data

    def _valid_flat_values(self, x: xr.DataArray, y: xr.DataArray) -> Tuple[np.ndarray, np.ndarray]:
        x_values = np.asarray(x.values).ravel()
        y_values = np.asarray(y.values).ravel()
        valid = np.isfinite(x_values) & np.isfinite(y_values)
        return x_values[valid], y_values[valid]

    @staticmethod
    def _summarize_member_correlations(correlations: List[Tuple[float, float, float]]) -> Dict[str, float]:
        if not correlations:
            return {"mean": np.nan, "low": np.nan, "high": np.nan}

        values = np.asarray(correlations, dtype=float)
        return {
            "mean": float(np.nanmean(values[:, 0])),
            "low": float(np.nanmean(values[:, 1])),
            "high": float(np.nanmean(values[:, 2])),
        }

    def compute_lag_correlation_curve(
        self,
        model_data: Dict[str, List[xr.Dataset]],
        var_land: str,
        var_atm: str,
        lags: List[int],
    ) -> FeedbackResults:
        if not lags:
            raise ValueError("lags must contain at least one lag value.")

        results: FeedbackResults = {}
        for exp, members in model_data.items():
            if not members:
                raise ValueError(f"Experiment {exp!r} has no members.")

            exp_result: LagSummary = {}
            for lag in lags:
                member_correlations = []
                for member in members:
                    self._require_vars(member, [var_land, var_atm])
                    x, y = self._apply_lag(member[var_land], member[var_atm], lag)
                    x = self._apply_land_mask(x)
                    y = self._apply_land_mask(y)
                    x_values, y_values = self._valid_flat_values(x, y)
                    corr = self._bootstrap_corr(x_values, y_values)
                    if np.isfinite(corr[0]):
                        member_correlations.append(corr)

                exp_result[lag] = self._summarize_member_correlations(member_correlations)
            results[exp] = exp_result

        return results

    def compute_spatial_map(
        self,
        model_data: List[xr.Dataset],
        var_land: str,
        var_atm: str,
        lag: int,
    ) -> xr.DataArray:
        if not model_data:
            raise ValueError("model_data must contain at least one member dataset.")

        maps = []
        for member in model_data:
            self._require_vars(member, [var_land, var_atm])
            x, y = self._apply_lag(member[var_land], member[var_atm], lag)
            rmap = xr.corr(x, y, dim="time")
            rmap = self._apply_land_mask(rmap)
            maps.append(rmap)

        out = xr.concat(maps, dim="member").mean(dim="member", skipna=True)
        out.name = f"{var_land}_{var_atm}_lag{lag}_feedback"
        out.attrs.update(
            {
                "var_land": var_land,
                "var_atm": var_atm,
                "lag": lag,
                "description": "Ensemble mean land-atmosphere lag correlation map",
            }
        )
        return out

    def save_to_netcdf(self, data: FeedbackResults, filename: str) -> None:
        if not data:
            raise ValueError("No feedback results to save.")

        datasets = []
        for exp, lag_dict in data.items():
            lags = sorted(lag_dict)
            ds = xr.Dataset(
                {
                    f"{exp}_mean": xr.DataArray(
                        [lag_dict[lag]["mean"] for lag in lags],
                        coords={"lag": lags},
                        dims=["lag"],
                    ),
                    f"{exp}_low": xr.DataArray(
                        [lag_dict[lag]["low"] for lag in lags],
                        coords={"lag": lags},
                        dims=["lag"],
                    ),
                    f"{exp}_high": xr.DataArray(
                        [lag_dict[lag]["high"] for lag in lags],
                        coords={"lag": lags},
                        dims=["lag"],
                    ),
                }
            )
            datasets.append(ds)

        final_ds = xr.merge(datasets)
        final_ds.attrs.update(
            {
                "alpha": self.alpha,
                "n_bootstrap": self.n_bootstrap,
                "min_samples": self.min_samples,
            }
        )

        output_dir = os.path.dirname(os.path.abspath(filename))
        os.makedirs(output_dir, exist_ok=True)
        final_ds.to_netcdf(filename)


def plot_lag_with_ci(
    results_dict: FeedbackResults,
    title: str = "Land-Atmosphere Lag Correlation",
):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    for exp, lag_r in results_dict.items():
        lags = sorted(lag_r.keys())
        means = [lag_r[lag]["mean"] for lag in lags]
        lows = [lag_r[lag]["low"] for lag in lags]
        highs = [lag_r[lag]["high"] for lag in lags]
        ax.plot(lags, means, label=exp)
        ax.fill_between(lags, lows, highs, alpha=0.2)

    ax.axvline(0, linestyle="--", color="gray")
    ax.set_xlabel("Lag (hours)")
    ax.set_ylabel("Correlation")
    ax.set_title(title)
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    plt.show()
    return fig, ax


def plot_feedback_map_with_sig(
    rmap: xr.DataArray,
    threshold: float = 0.2,
    title: str = "Spatial Feedback Map",
):
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    sig = np.abs(rmap) >= threshold

    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    rmap.plot(ax=ax, transform=ccrs.PlateCarree(), cmap="RdBu_r", vmin=-1, vmax=1, add_colorbar=True)
    sig.plot.contour(ax=ax, levels=[0.5], colors="black", linewidths=0.5, transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(title)
    plt.show()
    return fig, ax
