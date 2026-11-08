NOTEBOOK = "9_ssw_polar_vortex"   # outputs: rc.OUTPUT_ROOT/9_ssw_polar_vortex/{figures,data}/
# "auto": reuse the zonal-mean caches; "rebuild": re-read the model output; "require": fail if a cache is missing
derivation_mode = "auto"
CLIM_METHOD = "direct"            # Feb/Aug climatology: "direct" (Reanalysis Feb/Aug starts 1985-2010), "lead_drift", "target_month"
WORKERS = 4                       # processes reading model output (cache builds only; login node: 30 GB per user)

# Data locations (machine-specific)
OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")                       # reference-data archive
ERA5_ZONAL_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily_zonal"                    # ERA5 daily zonal means
ERA5_DAILY_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily"                          # ERA5 daily maps (PSL)
SIMULATION_ROOT = Path("/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation")         # raw eam.h1 history

YEAR = 2026
INIT_MONTHS = (2, 5, 8)
MONTH_NAMES = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug",
               9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"}
INIT_COLORS = {2: "#1f77b4", 5: "#2ca02c", 8: "#d62728"}
MODEL_PREFIX = rc.MODEL_PREFIX
POOL = {
    "Reanalysis-CLIM": f"{MODEL_PREFIX}_RTENSO2026_clim",
    "Reanalysis-OBS": f"{MODEL_PREFIX}_RTENSO2026_obs",
    "Reanalysis-NSIDC": f"{MODEL_PREFIX}_RTENSO2026_obs_ice_nsidc",
}

CLIM_YEARS = rc.CLIM_YEARS          # 1981-2010: model drift and ERA5 anomaly climatology
ERA5_FIRST_YEAR = 1979              # ERA5 daily record read (SSW climatology from winter 1979/80)
REF_LEVEL = 10.0                    # hPa, SSW wind
QBO_LEVEL = 50.0                    # hPa, QBO phase (Holton-Tan)
QBO_MONTHS = (10, 11)               # QBO phase months before the winter
DJF = ((0, 12), (1, 1), (1, 2))     # (years after the winter year, month)
SSW = {"sep_days": 20, "recover_days": 10}
NINO_PRODUCT_OBS = "HadISST1.1"     # observed DJF Nino3.4 (1870 on)
ENSO_THRESHOLD = 0.5                # degC: El Nino / La Nina winters (DJF Nino3.4)
LAGS = np.arange(-40, 61)           # days around SSW central dates
PCH_LAGS = (-30, 60)

fontz = 13
fig_type = "png"
fig_dpi = 200
plt.rcParams.update({"font.size": fontz, "axes.titlesize": fontz, "axes.labelsize": fontz,
                     "legend.fontsize": fontz * 0.8, "xtick.labelsize": fontz * 0.85, "ytick.labelsize": fontz * 0.85})


def save(fig, *parts):
    path = rc.save_figure(fig, NOTEBOOK, *parts, dpi=fig_dpi, ext=fig_type)
    print(f"[FIGURE] {path}")


DATA_DIR = rc.output_dir(NOTEBOOK, "data")