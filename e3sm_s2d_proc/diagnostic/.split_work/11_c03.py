NOTEBOOK = "11_polar_vortex_evolution"   # outputs: rc.OUTPUT_ROOT/11_polar_vortex_evolution/{figures,data}/
# "require": use the existing caches only (no raw-data build); "auto" builds a missing cache
derivation_mode = "require"
CLIM_METHOD = "direct"              # shading climatology: "direct" (each start's Reanalysis starts) or "lead_drift"

# Data locations (machine-specific)
OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")                       # reference-data archive
ERA5_ZONAL_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily_zonal"                    # ERA5 daily zonal means
SIMULATION_ROOT = Path("/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation")         # raw eam.h1 (cache keys only)
MAP_FILE = Path("/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc")

YEAR = 2026
EVENT_LABEL = "2026/27 polar vortex"
INIT_MONTHS = (2, 5, 8, 11)
MONTH_NAMES = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug",
               9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"}
INIT_COLORS = {2: "#1f77b4", 5: "#2ca02c", 8: "#d62728", 11: "#9467bd"}   # one colour per start (as notebook 8)
MODEL_PREFIX = rc.MODEL_PREFIX
POOL = {
    "Reanalysis-CLIM": f"{MODEL_PREFIX}_RTENSO2026_clim",
    "Reanalysis-OBS": f"{MODEL_PREFIX}_RTENSO2026_obs",
    "Reanalysis-NSIDC": f"{MODEL_PREFIX}_RTENSO2026_obs_ice_nsidc",
}

# One fixed member per start, chosen before looking at the results (here the first member of the clim case,
# present in every start); change it here, not according to what the member does.
MEMBERS = {2: "Reanalysis-CLIM:0", 5: "Reanalysis-CLIM:0", 8: "Reanalysis-CLIM:0", 11: "Reanalysis-CLIM:0"}

FRAMES = ("2026-10-01", "2027-04-30")      # first and last frame day
FRAME_STEP = 2                             # days between frames
PREVIEW_DAY = "2027-01-15"
STORYBOARD = ["2026-10-15", "2026-11-15", "2026-12-15", "2027-01-15", "2027-02-15", "2027-03-15", "2027-04-15"]

SHOW_ENSMEAN = True                 # ensemble-mean row under the member row
SHADING = "anomaly"                 # "anomaly": direct-climatology anomaly shading; None: raw contours only
SHOW_REFERENCE = True               # reference column: ERA5 and E3SM climatologies of the day (raw contours)
SHOW_SSW = False                    # mark the members' SSW central dates on the timeline
SSW = {"sep_days": 20, "recover_days": 10}
CLIM_YEARS = rc.CLIM_YEARS
LAT_MIN = 30.0
STYLE = {"ci": 100.0, "range": 1200.0, "cmap": "blue2red", "units": "m",
         "name": "10-hPa height anomaly",
         "line": 400.0, "line_levels": np.arange(26000.0, 33000.0, 400.0)}   # fixed absolute contour levels (m)
U60_YLIM = (-20.0, 55.0)

# Animation output (as notebook 8)
FPS = 6                    # frames per second (FRAME_STEP days each)
HOLD = 6                   # extra copies of the last frame, so a looping GIF pauses at the end
ANIM_DPI = 100             # 14 x 10 in at 100 dpi: 1400 x 1000 px
ANIM_FORMATS = ("gif", "video")   # "video": H.264 MP4 when ffmpeg has libx264, else VP9 WebM
SHOW_INLINE = True         # show the GIF in the notebook
fig_dpi = 200              # preview frame and storyboard

OUTDIR = rc.output_dir(NOTEBOOK, "data")
FIGDIR = rc.figure_dir(NOTEBOOK)