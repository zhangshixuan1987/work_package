NOTEBOOK = "10_polar_vortex_animation"   # outputs: rc.OUTPUT_ROOT/10_polar_vortex_animation/{figures,data}/
# "auto": reuse the map and climatology caches; "rebuild": re-read the model output and ERA5; "require": fail if missing
derivation_mode = "auto"
CLIM_METHOD = "direct"              # Feb/Aug climatology: "direct" (Reanalysis Feb/Aug starts 1985-2010), "lead_drift", "target_month"
WORKERS = 4                         # processes reading model output and ERA5 (cache builds only)

# Data locations (machine-specific)
OBSERVATION_ROOT = Path("/global/cfs/cdirs/e3sm/zhan391/data")                       # reference-data archive
ERA5_ZONAL_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily_zonal"                    # ERA5 daily zonal means
ERA5_DAILY_DIR = OBSERVATION_ROOT / "ERA5" / "1x1" / "daily"                          # ERA5 daily maps (PSL)
ERA5_MIRROR = Path("/global/cfs/projectdirs/m3522/cmip6/ERA5")                        # hourly 0.25-degree ERA5
SIMULATION_ROOT = Path("/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation")         # raw eam.h1 history
MAP_FILE = Path("/global/cfs/cdirs/e3sm/diagnostics/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc")

YEAR = 2026
INIT_MONTHS = (2, 5, 8, 11)
MONTH_NAMES = {1: "Jan", 2: "Feb", 3: "Mar", 4: "Apr", 5: "May", 6: "Jun", 7: "Jul", 8: "Aug",
               9: "Sep", 10: "Oct", 11: "Nov", 12: "Dec"}
INIT_COLORS = {2: "#1f77b4", 5: "#2ca02c", 8: "#9467bd", 11: "#d62728"}
MODEL_PREFIX = rc.MODEL_PREFIX
POOL = {
    "Reanalysis-CLIM": f"{MODEL_PREFIX}_RTENSO2026_clim",
    "Reanalysis-OBS": f"{MODEL_PREFIX}_RTENSO2026_obs",
    "Reanalysis-NSIDC": f"{MODEL_PREFIX}_RTENSO2026_obs_ice_nsidc",
}

WINTER = 2026                       # year of the November: 2026 -> 2026/27 (forecast); 2025 -> 2025/26 (observed)
FRAMES = (f"{WINTER}-10-01", f"{WINTER + 1}-04-30")   # first and last frame day
FRAME_STEP = 2                      # days between frames
FEATURED = {}                       # {init month: "Reanalysis-CLIM:3"}; a missing start picks its first SSW member
STORYBOARD = [f"{WINTER}-11-15", f"{WINTER}-12-15", f"{WINTER + 1}-01-15", f"{WINTER + 1}-02-15",
              f"{WINTER + 1}-03-15", f"{WINTER + 1}-04-15"]
CLIM_YEARS = rc.CLIM_YEARS
ERA5_FIRST_YEAR = 1979              # ERA5 daily record read for the climatologies
SSW = {"sep_days": 20, "recover_days": 10}

# Vortex index plumes (static figure over the whole forecast period)
INDEX_SPAN = (f"{YEAR}-01-01", f"{YEAR + 2}-07-31")
SMOOTH_DAYS = 5                     # centred running mean of the daily indices (1: none)
INDEX_RANGE = (0.1, 0.9)            # member quantiles shaded

LAT_MIN = 30.0                      # map edge
STYLE = {"ci": 100.0, "range": 1200.0, "line": 400.0, "units": "m", "cmap": "blue2red",
         "name": "10-hPa geopotential height anomaly"}

# Animation output
FPS = 6                             # frames per second (FRAME_STEP days each)
HOLD = 6                            # extra copies of the last frame, so a looping GIF pauses at the end
ANIM_DPI = 100                      # 14 x 9.5 in at 100 dpi: 1400 x 950 px
ANIM_FORMATS = ("gif", "video")     # "video": H.264 MP4 if ffmpeg has libx264 (imageio-ffmpeg), else VP9 WebM
SHOW_INLINE = True                  # show the GIF in the notebook
fig_dpi = 200                       # static storyboard

OUTDIR = rc.output_dir(NOTEBOOK, "data")
FIGDIR = rc.figure_dir(NOTEBOOK)
WINTER_LABEL = f"{WINTER}/{str(WINTER + 1)[2:]}"