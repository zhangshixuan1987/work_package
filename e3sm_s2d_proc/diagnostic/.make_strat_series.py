"""Write the 9a-9i stratosphere notebooks (split of the former 9, 10, 11).

Sections reuse the cell sources saved in .split_work/ (9_cNN, 10_cNN, 11_cNN); loading goes through
util/strat_diag.py. Run from the diagnostic folder: python .make_strat_series.py
"""
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
W = HERE / ".split_work"


def src(name):
    return (W / name).read_text()


def write(name, cells):
    nb = nbf.v4.new_notebook()
    nb["cells"] = [nbf.v4.new_markdown_cell(s.strip("\n")) if kind == "md" else nbf.v4.new_code_cell(s.strip("\n"))
                   for kind, s in cells]
    nb["metadata"] = {"kernelspec": {"display_name": "e3sm_analysis", "language": "python", "name": "e3sm_analysis"},
                      "language_info": {"name": "python"}}
    nbf.write(nb, HERE / f"{name}.ipynb")
    print("wrote", name, len(cells), "cells")


SERIES = r"""
**The 9a-9i series** (polar vortex and SSWs in the RTENSO2026 2026 starts): `9a` builds and checks the
caches; `9b` model climatology; `9c` vortex forecast; `9d` SSW occurrence; `9e` ENSO/QBO background;
`9f` SSW dynamics; `9g` SSW mechanism (featured member vs controls); `9h` SSW-event animation; `9i`
seasonal-evolution animation. Shared loading: `util/strat_diag.py`; readers `util/strat_case.py`; maps
`util/vortex_plot.py`.
"""

IMPORTS = r"""
%matplotlib inline

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from IPython.display import Image, display
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from util import rtenso_case as rc
from util import strat_case as sc
from util import strat_diag as sd

warnings.filterwarnings("ignore", message="invalid value encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice")
"""


def controls(notebook, extra="", mode="require"):
    return rf'''
NOTEBOOK = "{notebook}"   # outputs: rc.OUTPUT_ROOT/{notebook}/{{figures,data}}/
# "require": read the caches built by 9a only; "auto": build a missing cache; "rebuild": recompute
derivation_mode = "{mode}"
CLIM_METHOD = "direct"      # model climatology: "direct" (each start's Reanalysis starts; Feb/Aug 1985-2010), "lead_drift"
WORKERS = 4                 # processes reading model output (cache builds only; login node: 30 GB per user)

# Data locations (machine-specific)
OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")                       # reference-data archive
ERA5_ZONAL_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily_zonal"                    # ERA5 daily zonal means
ERA5_DAILY_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily"                          # ERA5 daily maps (PSL)
SIMULATION_ROOT = Path("/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation")         # raw eam.h1 history

cfg = sd.StratConfig(era5_zonal_dir=ERA5_ZONAL_DIR, era5_daily_dir=ERA5_DAILY_DIR, simulation_root=SIMULATION_ROOT,
                     clim_method=CLIM_METHOD, mode=derivation_mode, workers=WORKERS)
YEAR, INIT_MONTHS, POOL, CLIM_YEARS = cfg.year, cfg.init_months, cfg.pool, cfg.clim_years
ERA5_FIRST_YEAR, REF_LEVEL, SSW = cfg.era5_first_year, cfg.ref_level, cfg.ssw
MONTH_NAMES, INIT_COLORS, MODEL_PREFIX = sd.MONTH_NAMES, sd.INIT_COLORS, rc.MODEL_PREFIX
{extra.strip()}

fontz = 13
fig_type = "png"
fig_dpi = 200
plt.rcParams.update({{"font.size": fontz, "axes.titlesize": fontz, "axes.labelsize": fontz,
                     "legend.fontsize": fontz * 0.8, "xtick.labelsize": fontz * 0.85, "ytick.labelsize": fontz * 0.85}})
rng = np.random.default_rng(0)


def save(fig, *parts):
    path = rc.save_figure(fig, NOTEBOOK, *parts, dpi=fig_dpi, ext=fig_type)
    print(f"[FIGURE] {{path}}")


DATA_DIR = rc.output_dir(NOTEBOOK, "data")
'''


LOAD_ERA5 = r"""
%%time
E = sd.load_era5(cfg)                         # ERA5 daily zonal-mean indices and climatologies (cached)
era5, obs_end = E.fields, E.obs_end
e_u60, e_pcz, e_pct, e_vt, e_qbo, e_pch = E.u60, E.pcz, E.pct, E.vt, E.qbo, E.pch
clim_u60, std_u60, clim_pcz, std_pcz = E.clim_u60, E.std_u60, E.clim_pcz, E.std_pcz
e_month_clim, e_u60_month = E.month_clim, E.u60_month
on_days = sd.on_days
"""
LOAD_REFS = r"""
R = sd.load_references(cfg)                   # Reanalysis monthly zonal means of every start month
ref, ref_direct, REF_ALL, MPLEV = R.ref, R.direct, R.all, R.plev
e_month_clim_m = e_month_clim.sel(plev=MPLEV)   # ERA5 on the model monthly levels
"""
LOAD_MODEL = r"""
%%time
M = sd.load_model(cfg, E, R)                  # pooled 2026 starts: zonal means, climatologies, daily indices
monthly, daily, lead_clim, z_offset = M.monthly, M.daily, M.lead_clim, M.z_offset
m_u60, m_u60_clim, m_u60_adj, m_pch, m_vt = M.u60, M.u60_clim, M.u60_adj, M.pch, M.vt
mm_u60, mm_u60_anom, mm_qbo = M.mm_u60, M.mm_u60_anom, M.mm_qbo
"""
LOAD_CAT = r"""
C = sd.ssw_catalogue(cfg, E, M)               # SSW central dates and event properties (raw and bias-adjusted)
obs_events, obs_winters, obs_counts = C.obs_events, C.obs_winters, C.obs_counts
model_ssw, model_ssw_adj, events_all = C.model_ssw, C.model_ssw_adj, C.events
mid_month, case_of, djf_mean, DJF_DAYS = sd.mid_month, sd.case_of, sd.djf_mean, sd.DJF_DAYS
"""

CAUTION_FIRST = r"""
**Read with 9b.** E3SM's DJF vortex is weaker than ERA5's (U(60°N, 10 hPa) about 24 vs 29 m/s, and
about 24 vs 33 m/s in Dec-Jan) and varies far less between winters (σ about 4 vs 10 m/s). Two forecast
winters from one initial state are conditional forecast statistics, not a model climatology.
"""

# =====================================================================================================
# 9a: build and check the caches
# =====================================================================================================
nb9a = [
    ("md", r"""
# 9a. Stratosphere caches: build and check

Builds every cache the 9b-9i notebooks read, so they can run with `derivation_mode = "require"`, and
checks them: what exists, its size and coverage, and the one correction applied to the model data.

* **ERA5**: `ERA5/1x1/daily_zonal` (daily zonal means of U, V, T, Z and eddy fluxes on 37 levels; built
  by `ERA5/scripts/daily_zonal`) and its reduced indices (cached here); ERA5 00Z 10-hPa maps of the
  observed 2025/26 winter (for `WINTER = 2025` in 9h).
* **Model**: daily zonal means from the native `eam.h1` files, monthly zonal means of the post-processed
  U/T/Z3, daily 10-hPa height maps (native `Z010` remapped), for every `POOL` case and 2026 start.
* **Reanalysis climatologies**: monthly zonal means and the `Z10`, `PSL`, `TS` lead-time climatologies of
  every start month (May/Nov 1981-2010; Feb/Aug the six 1985-2010 starts: `CLIM_METHOD = "direct"`).

First build: about 1 h on a login node (most of it the daily model reads); later runs only check.
""" + SERIES),
    ("code", IMPORTS),
    ("md", "### Run controls"),
    ("code", controls("9a_strat_build_caches", r'''
MAP_FILE = Path("/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc")
ERA5_MIRROR = Path("/global/cfs/projectdirs/m3522/cmip6/ERA5")           # hourly 0.25-degree ERA5 (10-hPa maps)
OBSERVED_WINTER = ("2025-10-01", "2026-06-30")                          # ERA5 10-hPa maps of the observed 2025/26 winter
CLIM_FIELDS = ("Z10", "PSL", "TS")                                      # Reanalysis map climatologies (9c, 9e, 9f, 9g-9i)
''', mode="auto")),
    ("md", r"""
### 1. ERA5

The daily zonal-mean archive (years present, last day) and the reduced indices used by every notebook.
"""),
    ("code", r"""
%%time
files = sorted(ERA5_ZONAL_DIR.glob("ERA5.daily_zonal.en00.U.*.nc"))
years = [int(f.name.split(".")[-2][:4]) for f in files]
missing = sorted(set(range(ERA5_FIRST_YEAR, max(years) + 1)) - set(years))
print(f"ERA5 daily zonal: {len(files)} years {min(years)}-{max(years)}" + (f"; missing {missing}" if missing else ""))
E = sd.load_era5(cfg)
print(f"last day {E.obs_end:%d %b %Y}; indices cached")
era5_maps = sc.era5_daily_maps(ERA5_MIRROR, *OBSERVED_WINTER, mode=derivation_mode, workers=WORKERS)
print(f"ERA5 10-hPa maps {pd.Timestamp(era5_maps['time'].values[0]):%d %b %Y} - {pd.Timestamp(era5_maps['time'].values[-1]):%d %b %Y}")
rc.obs_monthly_climatology("Z10", mode=derivation_mode)
"""),
    ("md", r"""
### 2. Model zonal means and Reanalysis climatologies
"""),
    ("code", r"""
%%time
R = sd.load_references(cfg)
for m, ds in R.all.items():
    ok = np.isfinite(ds["U"].sel(plev=REF_LEVEL).isel(lead=0, lat=-1))
    print(f"Reanalysis init{m:02d}: {int(ok.sum())} member-years in {int(ok.any('member').sum())} years")
M = sd.load_model(cfg, E, R)
for field in CLIM_FIELDS:
    for init in INIT_MONTHS:
        c = rc.matched_climatology(field, init, method=CLIM_METHOD, mode=derivation_mode)
    print(f"{field} climatologies ready ({c.attrs.get('climatology_method', '')} ...)")
"""),
    ("md", r"""
### 3. Model daily 10-hPa maps
"""),
    ("code", r"""
%%time
for init in INIT_MONTHS:
    maps = sc.pooled_daily_maps(SIMULATION_ROOT, MAP_FILE, POOL, YEAR, init, mode=derivation_mode, workers=WORKERS)
    print(f"[init {MONTH_NAMES[init]}] {maps.sizes['member']} members, {maps.sizes['time']} days "
          f"({pd.Timestamp(maps['time'].values[0]):%b %Y} - {pd.Timestamp(maps['time'].values[-1]):%b %Y}, Sep-May)")
"""),
    ("md", r"""
### 4. Cache inventory and the daily-minus-monthly height offset

Every cache with its size. Then the one correction applied to the model data: the native daily `Z`
(eam.h1, computed online on pressure levels) is higher than the monthly `Z3` the climatology is made
from (interpolated afterwards; extrapolated below ground at 1000 hPa); U, T and PSL agree to 0.1. A
processing difference is nearly constant across leads (seasons) and starts, a physical signal is not:
the offset is about 38 m at 10 hPa with a standard deviation of about 2 m over 24 leads, identical in
the three starts, so adding it to the climatology removes a processing artifact, not signal.
"""),
    ("code", r"""
root = rc.diag_dir("")
rows = []
for pattern in ("*/strat/*/*.nc", "observations/strat/*.nc", "observations/strat/*/*.nc",
                "Reanalysis/leadtime_clim/*/Z10/*.nc", "Reanalysis/leadtime_clim/*/PSL/*init0[28]*.nc"):
    for f in sorted(root.glob(pattern)):
        rows.append({"cache": str(f.relative_to(root)), "MB": round(f.stat().st_size / 1e6, 1)})
inventory = pd.DataFrame(rows)
inventory.to_csv(DATA_DIR / "strat_cache_inventory.csv", index=False)
print(f"{len(inventory)} caches, {inventory['MB'].sum() / 1e3:.1f} GB under {root}")
display(inventory)
print("Polar-cap (60-90N) daily-minus-monthly Z offset (ensemble mean) by start: mean and std over 24 leads (m)")
display(sd.offset_table(M).round(1))
"""),
]

# =====================================================================================================
# 9b: model climatology
# =====================================================================================================
c10 = src("9_c10.py").replace(
    'e_series = sc.monthly_means(sc.lat_mean(era5[field].sel(plev=REF_LEVEL), band))',
    'e_series = sc.monthly_means(e_u60 if field == "U" else E.pct)')
nb9b = [
    ("md", r"""
# 9b. E3SM stratospheric climatology against ERA5

How does E3SM's stratosphere compare with ERA5 before any forecast is read? The model climatology is
the Reanalysis (BruteForce) 1981-2010 hindcasts at each lead (May and Nov starts, all years and
members), which is also the baseline the forecast anomalies of 9c-9i are measured from.

* Seasonal cycle of U(60°N, 10 hPa) and 10-hPa polar-cap temperature by calendar month.
* DJF zonal-mean U, model minus ERA5 (latitude-pressure).
* DJF vortex strength and its interannual variability (the spread that sets the model tercile bounds in 9c).

**Findings (as of the 2026 runs).** E3SM's vortex is about 10 m/s too weak in Dec-Jan (24 vs 33 m/s)
with no midwinter peak, and breaks down later (stronger than ERA5 in Mar-Apr); its DJF interannual
standard deviation is about 4 m/s against 10 m/s in ERA5; the 10-hPa polar cap is 4-6 K too cold and
the summer easterlies too weak.
""" + SERIES),
    ("code", IMPORTS),
    ("md", "### Run controls"),
    ("code", controls("9b_strat_model_climatology")),
    ("md", "### Load"),
    ("code", LOAD_ERA5),
    ("code", LOAD_REFS),
    ("md", src("9_c09.md").replace("### 3. ", "### ")),
    ("code", c10),
]

# =====================================================================================================
# 9c: vortex forecast
# =====================================================================================================
plume_controls = r'''
NINO_PRODUCT_OBS = "HadISST1.1"     # observed DJF Nino3.4 (El Nino winters marked)
ENSO_THRESHOLD = 0.5                # degC
DJF = sd.DJF
# Vortex index plume (whole forecast period)
INDEX_SPAN = (f"{YEAR}-01-01", f"{YEAR + 2}-07-31")
SMOOTH_DAYS = 5                     # centred running mean of the daily indices (1: none)
INDEX_RANGE = (0.1, 0.9)            # member quantiles shaded
'''
nb9c = [
    ("md", r"""
# 9c. The vortex forecast of the 2026 starts

What do the 1 Feb, 1 May and 1 Aug 2026 starts say about the polar vortex, and how much of it can be checked?

1. **Feb 2026 start against the observed end of winter 2025/26** (the only verifiable part).
2. **Daily U(60°N, 10 hPa) plumes** for 2026/27 and 2027/28, with the model's own climatology.
3. **Vortex index plume** over the whole forecast period: U(60°N), 10-hPa polar-cap temperature and
   height anomaly, polar-cap sea-level pressure.
4. **DJF vortex strength**: drift-corrected anomalies with uncalibrated tercile probabilities, and their
   sensitivity to the climatology method.

**Findings.** The Feb start caught the late-February 2026 weakening 3-4 weeks ahead (about a week early)
but not the reversal. With the direct climatology the DJF 2026/27 anomalies are +1.5 (Feb), +1.1 (May)
and 0.0 m/s (Aug), tercile probabilities near one third each: no predictable vortex anomaly. The
`lead_drift` climatology gave a spurious weak vortex for Feb/Aug (-5 to -7 m/s).
""" + CAUTION_FIRST + SERIES),
    ("code", IMPORTS),
    ("md", "### Run controls"),
    ("code", controls("9c_vortex_forecast", plume_controls)),
    ("md", "### Load"),
    ("code", LOAD_ERA5),
    ("code", LOAD_REFS),
    ("code", LOAD_MODEL),
    ("code", LOAD_CAT + "\nnino_obs, obs_djf_nino = sd.observed_nino(cfg, NINO_PRODUCT_OBS)\n"),
    ("md", src("9_c13.md").replace("### 5. ", "### 1. ")),
    ("code", src("9_c14.py")),
    ("md", src("9_c15.md").replace("### 6. The 2026/27 and 2027/28 winters", "### 2. Daily U(60°N, 10 hPa) plumes, 2026/27 and 2027/28")),
    ("code", src("9_c16.py")),
    ("md", src("10_c06.md").replace("### Vortex index plumes", "### 3. Vortex index plume")),
    ("code", src("10_c07.py")),
    ("code", src("10_c08.py")),
    ("md", "### 4. DJF vortex strength\n\n" + src("9_c20.md")),
    ("code", src("9_c21.py")),
]

# =====================================================================================================
# 9d: SSW occurrence
# =====================================================================================================
nb9d = [
    ("md", r"""
# 9d. SSW occurrence in the 2026 starts

How many SSWs do the forecasts produce, when, and does it depend on the start or the case?

1. **Frequency**: events per member-winter and the probability of at least one SSW, by start × case
   × winter (raw wind primary; bias-adjusted wind as a threshold-sensitivity test), against ERA5
   1979/80 on (all, El Niño and La Niña winters).
2. **Timing**: month of the central dates.
3. **Final warming**: dates with their coverage status (records ending before 31 Jul are censored or
   unresolved, never read as "no recovery").

**Findings.** 2026/27: 0.45 (Feb), 0.27 (May), 0.37 (Aug) SSWs per member-winter, 0.35 pooled; 2027/28:
0.30. Half the model events fall in March (ERA5: Jan-Feb). Case differences (e.g. May NSIDC 0.5 vs OBS
0.1) are within the 10-member sampling intervals. ERA5: 0.68 per winter. The model and ERA5 rates are not
directly comparable (conditional forecasts vs a climatology; a weaker, less variable model vortex).
""" + CAUTION_FIRST + SERIES),
    ("code", IMPORTS),
    ("md", "### Run controls"),
    ("code", controls("9d_ssw_occurrence", 'NINO_PRODUCT_OBS = "HadISST1.1"\nENSO_THRESHOLD = 0.5\nDJF = sd.DJF')),
    ("md", "### Load"),
    ("code", LOAD_ERA5),
    ("code", LOAD_REFS),
    ("code", LOAD_MODEL),
    ("code", LOAD_CAT),
    ("md", "### 1. SSW frequency and the event catalogue\n\n" + src("9_c17.md")),
    ("code", src("9_c18.py")),
    ("md", "### 2. Frequency, timing and final warming"),
    ("code", src("9_c19.py")),
]

# =====================================================================================================
# 9e: ENSO and QBO background
# =====================================================================================================
nb9e = [
    ("md", r"""
# 9e. ENSO and QBO background state of the vortex forecasts

In what tropical background do the forecasts put the 2026/27 vortex, and do the members' vortices vary with it?

1. **QBO**: equatorial zonal wind of each start against ERA5 (time-height and 50 hPa).
2. **Background state by start and case**: DJF Niño3.4 and Oct-Nov equatorial wind at 30 and 50 hPa.
3. **Member covariation within each start and case** (group means removed), so the differences
   between starts do not enter the fit.

This describes the forecast sample; it is not a test of the ENSO and Holton-Tan pathways, which needs a
long hindcast record. No ENSO × QBO interaction is fitted.

**Findings.** DJF Niño3.4 +2.4 to +4.7 °C (a very strong El Niño) with a westerly autumn-2026 QBO in every
member (U50 about +13 m/s): opposing influences, consistent with the near-normal vortex of 9c. The QBO
turns easterly for 2027/28 in nearly every member, without a weaker vortex. Within groups the member
vortex barely covaries with either (Niño3.4 r = +0.13, QBO r = -0.08).
""" + SERIES),
    ("code", IMPORTS),
    ("md", "### Run controls"),
    ("code", controls("9e_enso_qbo_background", 'QBO_LEVEL = 50.0\nQBO_MONTHS = (10, 11)\nENSO_THRESHOLD = 0.5                # degC: El Nino threshold line\nDJF = sd.DJF')),
    ("md", "### Load"),
    ("code", LOAD_ERA5),
    ("code", LOAD_REFS),
    ("code", LOAD_MODEL + "\nmid_month, case_of, djf_mean = sd.mid_month, sd.case_of, sd.djf_mean\n"
             "m_djf = sd.djf_vortex(cfg, M, R, E)          # {(start, winter): DJF U(60N, 10 hPa) anomaly per member}\n"),
    ("md", src("9_c11.md").replace("### 4. QBO", "### 1. QBO")),
    ("code", src("9_c12.py")),
    ("md", src("9_c22.md").replace("### 7. ENSO and QBO background state of the forecasts", "### 2. Background state and member covariation")),
    ("code", src("9_c23.py")),
]

# =====================================================================================================
# 9f: SSW dynamics
# =====================================================================================================
nb9f = [
    ("md", r"""
# 9f. Dynamics of the simulated SSWs

Are E3SM's SSWs realistic? Lag composites of every model SSW (raw wind; every start, case and
member-winter) against ERA5 1979-2025, season-matched to the model's onset months.

1. **Composites**: U(60°N, 10 hPa), the 100-hPa eddy heat flux, polar-cap height (time-height).
2. **Event properties**: pre-onset wave forcing, minimum wind, easterly days, recovery, persistence.
3. **Vertical evolution**: polar-cap anomalies at 10-500 hPa and polar-cap SLP in pre- and post-onset windows.

**Findings.** Model SSWs are late-season (half in March), shallow and brief (median minimum wind about
-3 vs -9 m/s; recovery about 4 vs 8 days, then a rebound to about +14 m/s), driven by weaker 100-hPa
heat flux (about 16 vs 27 K m/s; partly the daily-mean sampling). Their polar-cap anomaly appears at all
levels at once and lasts longer in the lower stratosphere; the simultaneous tropospheric signal suggests
co-variability (a post-onset tropospheric association), not a demonstrated downward influence.
""" + SERIES),
    ("code", IMPORTS),
    ("md", "### Run controls"),
    ("code", controls("9f_ssw_dynamics")),
    ("md", "### Load"),
    ("code", LOAD_ERA5),
    ("code", LOAD_REFS),
    ("code", LOAD_MODEL),
    ("code", LOAD_CAT),
    ("md", src("9_c24.md").replace("### 8. SSW event dynamics", "### 1. Composites")),
    ("code", src("9_c25.py")),
    ("md", "### 2. Event properties and vertical evolution"),
    ("code", src("9_c26.py")),
]

# =====================================================================================================
# 9g: SSW mechanism (featured member vs controls)
# =====================================================================================================
nb9g = [
    ("md", r"""
# 9g. What distinguishes an SSW member from comparable members that do not reverse?

An event-mechanism study of one simulated SSW: the **featured member** (`EVENT`, chosen by hand) against
two controls from the same start and winter, both without an SSW:

* **matched member**: the closest pre-onset vortex state (mean U(60°N, 10 hPa) over `MATCH_WINDOW`
  days before the featured central date);
* **near-miss member**: the strongest weakening that does not reverse (lowest minimum U(60°N, 10 hPa)
  within `NEAR_WINDOW` days of the central date). It is chosen **from its later outcome**, so it is a
  descriptive control, not an independent prediction test.

Then: aligned timelines (raw U(60°N, 10 hPa), 10-hPa polar-cap temperature, 100-hPa heat flux) over
`EVENT_WINDOW`; where the featured member ranks among all members of its start; **whether the featured
trajectory reflects an ensemble relationship** (preconditioning and pre-onset wave activity against the
outcome, every member); pre-onset 10-hPa height maps (`PRE_WINDOWS`); and a side-by-side animation.

**Interpretation caution.** These comparisons identify precursors and a plausible mechanism; they do not
prove what caused the SSW. Similar U(60°N, 10 hPa) does not mean identical preconditioning (vertical wind
structure, vortex geometry and wave-propagation conditions can differ). The 100-hPa heat flux is a proxy for
upward wave activity, not the momentum forcing at 10 hPa (EP-flux convergence would be closer; it needs the
eddy momentum flux, which the model cache does not hold). The model flux is the eddy covariance of daily-mean
fields, so sub-daily contributions are missing: compare members with one another (same sampling), not with
the magnitude of the differently sampled ERA5 flux. Polar-cap temperature is an area mean and U(60°N) one
latitude ring, so they need not move together; vortex geometry matters.

For the broad seasonal view see 9i (fixed members, no event-based selection) and 9h (selected SSW members).
""" + SERIES),
    ("code", IMPORTS + "\nfrom matplotlib import animation\nfrom esp_lab.utils.filename_utils import figure_filename\nfrom util import vortex_plot as vp\n"),
    ("md", "### Run controls"),
    ("code", controls("9g_ssw_mechanism", r'''
MAP_FILE = Path("/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc")
EVENT = {"init": 8, "member": "Reanalysis-CLIM:9", "winter": 2026}   # featured SSW member (its first SSW of the winter)
MATCH_WINDOW = (-40, -21)          # days before the central date: matched member = closest mean U(60N, 10 hPa)
NEAR_WINDOW = (-15, 15)            # days around the central date: near-miss = lowest minimum U60 without reversal
EVENT_WINDOW = (-30, 30)           # timelines and animation
PRE_WINDOWS = [(-20, -11), (-10, -1)]   # pre-onset map averages and ranks
RANK_WINDOW = (-20, -1)            # pre-onset window of the ensemble ranks (and the wave activity of section 4)
PRECOND_WINDOW = (-30, -21)        # preconditioning: mean U(60N, 10 hPa)
OUTCOME_WINDOW = (0, 10)           # outcome: minimum U(60N, 10 hPa)
ROLES = {"featured": "#d62728", "matched": "#1f77b4", "near-miss": "#ff7f0e"}
LAT_MIN = 30.0
STYLE = {"ci": 100.0, "range": 1200.0, "cmap": "blue2red", "units": "m", "name": "10-hPa height anomaly",
         "line": 400.0, "line_levels": np.arange(26000.0, 33000.0, 400.0)}
FPS, HOLD, ANIM_DPI, ANIM_FORMATS = 5, 5, 100, ("gif", "video")
FIGDIR = rc.figure_dir(NOTEBOOK)
''')),
    ("md", "### Load"),
    ("code", LOAD_ERA5),
    ("code", LOAD_REFS),
    ("code", LOAD_MODEL),
    ("code", LOAD_CAT),
    ("md", r"""
### 1. The featured event and its controls

The featured member's first SSW of the winter sets day 0. Controls are drawn from the members of the
same start with no SSW in that winter (raw wind).
"""),
    ("code", r"""
init, feat, w = EVENT["init"], EVENT["member"], EVENT["winter"]
u, t10, vt = m_u60[init], M.pct[init], m_vt[init]
times = pd.DatetimeIndex(u["time"].values)
ev = model_ssw[init]
ev_feat = ev[(ev["series"] == feat) & (ev["winter"] == w)]
if ev_feat.empty:
    raise ValueError(f"{feat} has no SSW in {sd.winter_label(w)} (start {MONTH_NAMES[init]}); choose another EVENT")
day0 = pd.Timestamp(ev_feat["central_date"].iloc[0])
ssw_members = set(ev[ev["winter"] == w]["series"])
controls_pool = [str(m) for m in u["member"].values if str(m) not in ssw_members]


def window(da, lo, hi):
    return da.sel(time=slice(day0 + pd.Timedelta(days=lo), day0 + pd.Timedelta(days=hi)))


pre_match = window(u, *MATCH_WINDOW).mean("time").to_series()
near_min = window(u, *NEAR_WINDOW).min("time").to_series()
matched = (pre_match[controls_pool] - pre_match[feat]).abs().idxmin()
near = near_min[[m for m in controls_pool if m != matched]].idxmin()
MEMBERS = {"featured": feat, "matched": matched, "near-miss": near}
print(f"Featured {feat}, start {MONTH_NAMES[init]} {YEAR}: SSW central date {day0:%d %b %Y}; "
      f"{len(ssw_members)} of {u.sizes['member']} members have an SSW in {sd.winter_label(w)}, {len(controls_pool)} controls")
summary = pd.DataFrame({role: {
    f"U60 mean, days {MATCH_WINDOW[0]} to {MATCH_WINDOW[1]} (m/s)": pre_match[m],
    f"U60 mean, days {RANK_WINDOW[0]} to {RANK_WINDOW[1]} (m/s)": float(window(u.sel(member=m), *RANK_WINDOW).mean()),
    f"[v*T*] 100 hPa mean, days {RANK_WINDOW[0]} to {RANK_WINDOW[1]} (K m/s)": float(window(vt.sel(member=m), *RANK_WINDOW).mean()),
    f"U60 minimum, days {NEAR_WINDOW[0]} to {NEAR_WINDOW[1]} (m/s)": near_min[m],
    "case": sd.case_of(m)} for role, m in MEMBERS.items()}).T.assign(member=list(MEMBERS.values()))
display(summary.round(1))
"""),
    ("md", r"""
### 2. Aligned timelines

Days from the featured central date. Thick lines: the three members; grey: the median and 10-90% range of
all members of the start. The heat flux is the daily-mean eddy covariance (5-day running mean thin line on top).
"""),
    ("code", r"""
lags = np.arange(EVENT_WINDOW[0], EVENT_WINDOW[1] + 1)
lag_days = (day0 + pd.to_timedelta(lags, unit="D")).values


def on_lags(da):
    return da.sel(time=lag_days).assign_coords(time=lags).rename(time="lag")


fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True, gridspec_kw={"hspace": 0.12})
for ax, (series, title, units) in zip(axes, [(u, "U(60°N, 10 hPa)", "m s$^{-1}$"), (t10, "T(60-90°N, 10 hPa)", "K"),
                                             (vt, "[v*T*] 100 hPa, 45-75°N (daily means)", "K m s$^{-1}$")]):
    s = on_lags(series)
    q = s.quantile([0.1, 0.5, 0.9], "member")
    ax.fill_between(lags, q.isel(quantile=0), q.isel(quantile=2), color="0.85", lw=0, label="all members 10-90%")
    ax.plot(lags, q.isel(quantile=1), color="0.45", lw=1.5, label="all members median")
    for role, m in MEMBERS.items():
        v = s.sel(member=m)
        ax.plot(lags, v, color=ROLES[role], lw=2.4 if role == "featured" else 1.8, label=f"{role} ({m.replace('Reanalysis-', '')})")
        if "v*T*" in title:
            ax.plot(lags, v.rolling(lag=5, center=True, min_periods=1).mean(), color=ROLES[role], lw=0.8, ls=":")
    if "U(60" in title:
        ax.axhline(0, color="k", lw=1)
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_ylabel(units)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(alpha=0.3)
axes[0].legend(ncol=5, fontsize=fontz * 0.65, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.08))
axes[-1].set_xlabel(f"days from the featured central date ({day0:%d %b %Y})")
fig.suptitle(f"Featured SSW member and controls, 1 {MONTH_NAMES[init]} {YEAR} start", x=0.12, ha="left", y=0.94, fontweight="bold")
save(fig, "mechanism", "timelines", MONTH_NAMES[init], feat.replace(":", ""))
plt.show()
"""),
    ("md", r"""
### 3. Where the featured member ranks in its ensemble

Pre-onset (`RANK_WINDOW`) mean U(60°N, 10 hPa) and 100-hPa heat flux of every member of the start: SSW
members (filled) and non-SSW members (open), the three compared members highlighted. A featured member at the
weak-wind, high-flux end supports a wave-driven precondition; one in the middle suggests chance.
"""),
    ("code", r"""
pre_u = window(u, *RANK_WINDOW).mean("time").to_series()
pre_vt = window(vt, *RANK_WINDOW).mean("time").to_series()
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
rank_rows = []
for ax, (vals, label) in zip(axes, [(pre_u, "U(60°N, 10 hPa) (m/s)"), (pre_vt, "[v*T*] 100 hPa (K m/s)")]):
    order = vals.sort_values()
    y = rng.uniform(-0.15, 0.15, order.size)
    is_ssw = np.array([m in ssw_members for m in order.index])
    ax.scatter(order.values[is_ssw], y[is_ssw], color="0.3", s=30, label="SSW members")
    ax.scatter(order.values[~is_ssw], y[~is_ssw], facecolor="none", edgecolor="0.3", s=30, label="non-SSW members")
    for role, m in MEMBERS.items():
        ax.scatter([vals[m]], [0.32], marker="v", s=140, color=ROLES[role], edgecolor="k", zorder=5, label=role)
        rank_rows.append({"quantity": label, "role": role, "member": m, "value": vals[m],
                          "rank (1 = lowest)": int(order.index.get_loc(m)) + 1, "of": order.size,
                          "percentile": 100.0 * (order < vals[m]).mean()})
    ax.set_yticks([])
    ax.set_ylim(-0.4, 0.5)
    ax.set_xlabel(f"{label}, days {RANK_WINDOW[0]} to {RANK_WINDOW[1]}")
    ax.grid(alpha=0.3, axis="x")
axes[0].legend(fontsize=fontz * 0.65, frameon=False, ncol=2, loc="lower left")
fig.suptitle(f"Pre-onset state of all {u.sizes['member']} members, 1 {MONTH_NAMES[init]} {YEAR} start", x=0.12, ha="left", fontweight="bold")
save(fig, "mechanism", "ensemble_ranks", MONTH_NAMES[init], feat.replace(":", ""))
plt.show()
display(pd.DataFrame(rank_rows).round(1))
"""),
    ("md", r"""
### 4. Is the featured member representative? The ensemble relationship

Every member of the start, on the featured member's calendar days: preconditioning (mean U(60°N, 10 hPa),
`PRECOND_WINDOW`) against pre-onset wave activity (mean 100-hPa heat flux, `RANK_WINDOW`), coloured by the
outcome (minimum U(60°N, 10 hPa), `OUTCOME_WINDOW`); filled markers have a qualifying SSW within the
outcome window, marker shapes are the cases. If members with more pre-onset wave activity on a weaker vortex
weaken more, the featured member is one instance of an ensemble relationship; otherwise it is one
illustrative trajectory. Spearman rank correlations quantify it (with 20-30 members, only strong
relationships are distinguishable from chance).
"""),
    ("code", r"""
from scipy import stats

pre_c = window(u, *PRECOND_WINDOW).mean("time").to_series()
pre_w = window(vt, *RANK_WINDOW).mean("time").to_series()
out_min = window(u, *OUTCOME_WINDOW).min("time").to_series()
lo, hi = day0 + pd.Timedelta(days=OUTCOME_WINDOW[0] - 10), day0 + pd.Timedelta(days=OUTCOME_WINDOW[1])
ssw_near = set(ev[(ev["central_date"] >= lo) & (ev["central_date"] <= hi)]["series"])
rel = pd.DataFrame({"precond_u60": pre_c, "preonset_vt": pre_w, "outcome_min_u60": out_min})
rel["ssw"] = [m in ssw_near for m in rel.index]
rel["case"] = [sd.case_of(m) for m in rel.index]
rel.to_csv(DATA_DIR / f"mechanism_ensemble_relationship_{MONTH_NAMES[init]}_{feat.replace(':', '')}.csv")

fig, ax = plt.subplots(figsize=(9.5, 6.5))
markers = dict(zip(sorted(rel["case"].unique()), ("o", "s", "^")))
norm_o = plt.Normalize(-10, 30)
for case, g in rel.groupby("case"):
    for filled in (True, False):
        gg = g[g["ssw"] == filled]
        if gg.empty:
            continue
        kw = {"c": gg["outcome_min_u60"], "cmap": "RdBu", "norm": norm_o} if filled else              {"facecolors": "none", "edgecolors": plt.get_cmap("RdBu")(norm_o(gg["outcome_min_u60"].values))}
        ax.scatter(gg["precond_u60"], gg["preonset_vt"], s=90, marker=markers[case], lw=1.6, zorder=3,
                   label=f"{case.replace('Reanalysis-', '')}{', SSW' if filled else ''}", **kw)
for role, m in MEMBERS.items():
    ax.scatter(rel.loc[m, "precond_u60"], rel.loc[m, "preonset_vt"], s=320, facecolors="none", edgecolors=ROLES[role],
               lw=2.5, zorder=4)
    ax.annotate(role, (rel.loc[m, "precond_u60"], rel.loc[m, "preonset_vt"]), xytext=(8, 8), textcoords="offset points",
                color=ROLES[role], fontweight="bold", fontsize=fontz * 0.8)
sm = plt.cm.ScalarMappable(norm=norm_o, cmap="RdBu")
fig.colorbar(sm, ax=ax, label=f"minimum U(60°N, 10 hPa), days {OUTCOME_WINDOW[0]} to {OUTCOME_WINDOW[1]} (m/s)")
ax.set_xlabel(f"preconditioning: mean U(60°N, 10 hPa), days {PRECOND_WINDOW[0]} to {PRECOND_WINDOW[1]} (m/s)")
ax.set_ylabel(f"pre-onset wave activity: mean [v*T*] 100 hPa, days {RANK_WINDOW[0]} to {RANK_WINDOW[1]} (K m/s)")
ax.set_title(f"All {len(rel)} members of the 1 {MONTH_NAMES[init]} {YEAR} start, around {day0:%d %b %Y}", loc="left")
ax.legend(fontsize=fontz * 0.65, frameon=False, loc="upper right")
ax.grid(alpha=0.3)
save(fig, "mechanism", "ensemble_relationship", MONTH_NAMES[init], feat.replace(":", ""))
plt.show()

rows = []
for x, name in (("preonset_vt", "pre-onset heat flux"), ("precond_u60", "preconditioning U60")):
    r, p = stats.spearmanr(rel[x], rel["outcome_min_u60"])
    rows.append({"predictor": name, "outcome": "minimum U60", "Spearman r": r, "p": p, "n": len(rel)})
r, p = stats.spearmanr(rel["preonset_vt"], rel["precond_u60"])
rows.append({"predictor": "pre-onset heat flux", "outcome": "preconditioning U60", "Spearman r": r, "p": p, "n": len(rel)})
print(f"Members with an SSW between days {OUTCOME_WINDOW[0] - 10} and {OUTCOME_WINDOW[1]}: {len(ssw_near)} of {len(rel)}")
display(pd.DataFrame(rows).round(3))
"""),
    ("md", r"""
### 5. Pre-onset 10-hPa height

Mean 10-hPa height of each member over the `PRE_WINDOWS` and the first days after onset: raw E3SM contours
(fixed levels) and the anomaly from the start's climatology (shading; `CLIM_METHOD`). Map orientation:
North America at the bottom, Europe on the right, Asia at the top, the Pacific on the left.
"""),
    ("code", r"""
%%time
maps = sc.pooled_daily_maps(SIMULATION_ROOT, MAP_FILE, POOL, YEAR, init, mode=derivation_mode)
win_days = pd.date_range(day0 + pd.Timedelta(days=EVENT_WINDOW[0]), day0 + pd.Timedelta(days=EVENT_WINDOW[1]), freq="D")
win_days = win_days[win_days.isin(pd.DatetimeIndex(maps["time"].values))]
raw = maps.sel(member=list(MEMBERS.values()), time=win_days.values).sel(lat=slice(LAT_MIN - 0.01, 90.0)).load()
clim = sc.map_lead_climatology(POOL, YEAR, init, maps, win_days.values, method=CLIM_METHOD, lat_min=LAT_MIN, mode="auto").load()
anom = raw - clim
map_windows = PRE_WINDOWS + [(0, 9)]
fig, axes = plt.subplots(len(MEMBERS), len(map_windows), figsize=(3.6 * len(map_windows), 3.7 * len(MEMBERS)),
                         subplot_kw={"projection": vp.PROJ}, squeeze=False, gridspec_kw={"hspace": 0.15, "wspace": 0.04})
for r, (role, m) in enumerate(MEMBERS.items()):
    for c, (lo, hi) in enumerate(map_windows):
        sel = slice(day0 + pd.Timedelta(days=lo), day0 + pd.Timedelta(days=hi))
        st = vp.setup_map(axes[r, c], f"days {lo:+d} to {hi:+d}" if r == 0 else "", lat_min=LAT_MIN, color="k",
                          title_size=fontz * 0.8, badge_size=fontz * 0.7)
        uu = float(u.sel(member=m, time=sel).mean())
        vp.fill_map(axes[r, c], st, anom.sel(member=m, time=sel).mean("time"), raw.sel(member=m, time=sel).mean("time"),
                    style=STYLE, total_levels=STYLE["line_levels"], **vp.u60_badge(uu, "U60 mean"))
        if c == 0:
            axes[r, c].text(-0.1, 0.5, f"{role}\n{m.replace('Reanalysis-', '')}", transform=axes[r, c].transAxes,
                            rotation=90, ha="center", va="center", color=ROLES[role], fontweight="bold", fontsize=fontz * 0.85)
levels, cmap, norm = vp.shading(STYLE)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, orientation="horizontal", fraction=0.03, pad=0.03,
                  shrink=0.6, extend="both", ticks=np.arange(-1200, 1201, 400))
cb.set_label(f"10-hPa height anomaly (m, {CLIM_METHOD} climatology); contours: raw height every {STYLE['line']:g} m")
fig.suptitle(f"Pre-onset 10-hPa height, featured SSW ({day0:%d %b %Y}) and controls", y=0.93, fontweight="bold")
save(fig, "mechanism", "preonset_maps", MONTH_NAMES[init], feat.replace(":", ""))
plt.show()
"""),
    ("md", r"""
### 6. Side-by-side animation

The three members over `EVENT_WINDOW`, one frame per day, above their U(60°N, 10 hPa) timelines.
"""),
    ("code", r"""
%%time
fig = plt.figure(figsize=(13.0, 7.6))
gs = fig.add_gridspec(3, 3, height_ratios=(1.0, 0.045, 0.5), hspace=0.32, wspace=0.05, left=0.05, right=0.98, top=0.86, bottom=0.08)
panels = {}
for k, (role, m) in enumerate(MEMBERS.items()):
    ax = fig.add_subplot(gs[0, k], projection=vp.PROJ)
    panels[role] = (ax, vp.setup_map(ax, f"{role}: {m.replace('Reanalysis-', '')}", lat_min=LAT_MIN, color=ROLES[role]))
cax = fig.add_subplot(gs[1, :])
pos = cax.get_position()
cax.set_position([pos.x0 + pos.width * 0.25, pos.y0, pos.width * 0.5, pos.height])
fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", extend="both",
             ticks=np.arange(-1200, 1201, 400)).set_label("10-hPa height anomaly (m); contours: raw height every 400 m")
axt = fig.add_subplot(gs[2, :])
s = on_lags(u)
q = s.quantile([0.1, 0.9], "member")
axt.fill_between(lags, q.isel(quantile=0), q.isel(quantile=1), color="0.88", lw=0)
axt.axhline(0, color="k", lw=1)
lines = {}
for role, m in MEMBERS.items():
    axt.plot(lags, s.sel(member=m), color=ROLES[role], lw=1.2, alpha=0.25)
    lines[role] = axt.plot([], [], color=ROLES[role], lw=2.4)[0]
cursor = axt.axvline(lags[0], color="0.2", ls="--")
axt.set_xlim(lags[0], lags[-1])
axt.set_xlabel(f"days from the featured central date ({day0:%d %b %Y})")
axt.set_ylabel("U(60°N, 10 hPa) (m s$^{-1}$)")
axt.grid(alpha=0.3)
title = fig.suptitle("", x=0.05, ha="left", y=0.97, fontsize=15, fontweight="bold")


def draw(lag):
    t = day0 + pd.Timedelta(days=int(lag))
    title.set_text(f"Featured SSW and controls, 1 {MONTH_NAMES[init]} {YEAR} start: day {int(lag):+d} ({t:%d %b %Y})")
    for role, m in MEMBERS.items():
        ax, st = panels[role]
        tn = np.datetime64(t, "ns")
        if tn in raw["time"].values:
            vp.fill_map(ax, st, anom.sel(member=m, time=tn), raw.sel(member=m, time=tn), style=STYLE,
                        total_levels=STYLE["line_levels"], **vp.u60_badge(float(u.sel(member=m, time=tn))))
        lines[role].set_data(lags[lags <= lag], s.sel(member=m).values[lags <= lag])
    cursor.set_xdata([lag, lag])
    return []


draw(0)
fig.savefig(FIGDIR / figure_filename("mechanism", "preview", MONTH_NAMES[init], feat.replace(":", ""), ext="png"), dpi=fig_dpi)
plt.show()
video_ext, video = vp.video_writer(FPS)
writers = {"gif": ("gif", lambda: animation.PillowWriter(fps=FPS)), "video": (video_ext, video)}
anim = animation.FuncAnimation(fig, draw, frames=list(lags) + [lags[-1]] * HOLD, blit=False, cache_frame_data=False)
for fmt in ANIM_FORMATS:
    ext, writer = writers[fmt]
    if writer is None:
        continue
    path = FIGDIR / figure_filename("mechanism", "animation", MONTH_NAMES[init], feat.replace(":", ""), ext=ext)
    anim.save(path, writer=writer(), dpi=ANIM_DPI)
    print(f"{path} ({path.stat().st_size / 1e6:.1f} MB)")
plt.close(fig)
"""),
    ("md", r"""
### How to read these results

* **Featured vs matched.** Members that start from a similar vortex (`MATCH_WINDOW`) can follow very
  different wave-activity histories and wind evolutions: the initial U(60°N, 10 hPa) alone does not set the
  outcome. A featured member with repeated heat-flux pulses before onset, against a matched member with
  little, is consistent with enhanced upward planetary-wave activity contributing to the disruption.
* **Featured vs near-miss.** Comparable flux peaks can differ in timing and duration; forcing after the
  featured onset cannot explain its earlier reversal. Compare accumulated pre-onset wave activity and the
  evolving vortex state, not only peak flux.
* **Wording.** The 100-hPa heat flux is a *prominent diagnostic contrast* and a *proxy* for upward wave
  activity. Pulses that *precede* episodes of weakening are consistent with wave-driven deceleration; the
  lag is not quantified here. Ranks ("largest of the members") come from the individual-member values in
  section 3. Polar-cap warming *accompanies* the disruption; the maps (section 5) are needed to tell a
  displacement from broader warming. Magnitudes are not compared with ERA5 here (different sampling and no
  matched averaging). A polar-cap warming without a U(60°N) reversal is described as such; a formal minor-warming
  classification needs further diagnostics.
* **Representativeness.** Section 4 says whether the featured member is one instance of an ensemble
  relationship or one illustrative trajectory.
"""),
]

# =====================================================================================================
# 9h: SSW-event animation (former 10, without the plume, which moved to 9c)
# =====================================================================================================
c03 = src("10_c03.py").replace('NOTEBOOK = "10_polar_vortex_animation"   # outputs: rc.OUTPUT_ROOT/10_polar_vortex_animation/{figures,data}/',
                               'NOTEBOOK = "9h_vortex_ssw_animation"   # outputs: rc.OUTPUT_ROOT/9h_vortex_ssw_animation/{figures,data}/')
c00 = src("10_c00.md").replace("# RTENSO2026: animated evolution of the stratospheric polar vortex",
                               "# 9h. Animated polar vortex with selected SSW members")
i0 = c00.find("* **Vortex index plumes**")
if i0 >= 0:
    j0 = c00.find("* **Timeline**", i0)
    c00 = c00[:i0] + c00[j0:]
c00 += "\n\nThe vortex index plume formerly here is in 9c; the unselected seasonal evolution is 9i; the mechanism of one event is 9g.\n" + SERIES
nb9h = [("md", c00), ("code", src("10_c01.py") + "\nfrom util import strat_diag as sd\n"), ("md", src("10_c02.md")),
        ("code", c03.replace('derivation_mode = "auto"', 'derivation_mode = "require"') + r'''

cfg = sd.StratConfig(era5_zonal_dir=ERA5_ZONAL_DIR, era5_daily_dir=ERA5_DAILY_DIR, simulation_root=SIMULATION_ROOT,
                     clim_method=CLIM_METHOD, mode=derivation_mode)
'''),
        ("md", src("10_c04.md")), ("code", src("10_c05.py")),
        ("code", r'''
# Each start's own lead-time climatology of U(60N, 10 hPa) (dotted lines of the timeline)
R = sd.load_references(cfg)
E = sd.load_era5(cfg) if CLIM_METHOD != "direct" else None
m_idx = {}
for init in fc:
    zd = sc.pooled_daily_zonal(SIMULATION_ROOT, POOL, YEAR, init, mode=derivation_mode)
    lc = sc.lead_climatology(init, R.all if CLIM_METHOD == "direct" else R.ref, None if CLIM_METHOD == "direct" else
                             E.month_clim.sel(plev=R.plev), method=CLIM_METHOD)
    m_idx[init] = {"clim": {"u60": sc.lat_mean(sc.daily_from_lead(lc["U"].sel(plev=10.0), YEAR, init, zd["time"].values),
                                               sc.U60).assign_coords(time=zd["time"])}}
''')]
for k in range(9, 18):
    name = f"10_c{k:02d}"
    p = W / f"{name}.md"
    nb9h.append(("md", src(f"{name}.md")) if p.exists() else ("code", src(f"{name}.py")))

# =====================================================================================================
# 9i: seasonal evolution animation (former 11)
# =====================================================================================================
nb9i = []
for k in range(14):
    name = f"11_c{k:02d}"
    p = W / f"{name}.md"
    s = src(f"{name}.md") if p.exists() else src(f"{name}.py")
    s = (s.replace('NOTEBOOK = "11_polar_vortex_evolution"   # outputs: rc.OUTPUT_ROOT/11_polar_vortex_evolution/{figures,data}/',
                   'NOTEBOOK = "9i_vortex_evolution_animation"   # outputs: rc.OUTPUT_ROOT/9i_vortex_evolution_animation/{figures,data}/')
          .replace("Notebook 10 explains selected SSW events (members picked for their SSW); this notebook shows the\n"
                   "seasonal evolution without event-based selection. Notebook 9 has the statistics (SSW frequency,\n"
                   "vortex terciles, composites) behind both.",
                   "9h animates selected SSW members and 9g studies one event's mechanism; this notebook shows the\n"
                   "seasonal evolution without event-based selection. 9b-9f have the statistics behind them.")
          .replace("# RTENSO2026: seasonal evolution of the 2026/27 polar vortex", "# 9i. Seasonal evolution of the 2026/27 polar vortex"))
    if k == 0:
        s += "\n" + SERIES
    nb9i.append(("md" if p.exists() else "code", s))

for name, cells in (("9a_strat_build_caches", nb9a), ("9b_strat_model_climatology", nb9b), ("9c_vortex_forecast", nb9c),
                    ("9d_ssw_occurrence", nb9d), ("9e_enso_qbo_background", nb9e), ("9f_ssw_dynamics", nb9f),
                    ("9g_ssw_mechanism", nb9g), ("9h_vortex_ssw_animation", nb9h), ("9i_vortex_evolution_animation", nb9i)):
    write(name, cells)
