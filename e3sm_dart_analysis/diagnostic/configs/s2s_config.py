"""Shared experiment, variable, region, and lead-window catalogs for S2S diagnostics."""

EXPERIMENT_GROUPS = {
    "Dec2011": dict(models=("CTRL10-S0", "DART20-S0", "DART40-S0"), period="201112-201112", season="Analysis", init_date="2011-12-01", init_month=12, run="da", freq="daily"),
    "Jan2012": dict(models=("CTRL10-S0", "CAPT10-S0", "DART20-S0", "DART40-S0"), period="201201-201203", season="Winter", init_date="2012-01-01", init_month=1, run="fc", freq="daily"),
    "Jun2012": dict(models=("CTRL10-S1", "CAPT10-S1", "DART40-S1"), period="201206-201208", season="Summer", init_date="2012-06-01", init_month=6, run="fc", freq="daily"),
}

OBS_VARIABLES = {
    "ERA5": {
        "UBOT":"U10", "VBOT":"V10", "PRECT":"PRECT", "PSL":"PSL", "T925":"T925", "Q925":"Q925",
        "U850":"U850", "V850":"V850", "T850":"T850", "Q850":"Q850", "U500":"U500", "V500":"V500", "T500":"T500", "Q500":"Q500", "U200":"U200", "V200":"V200", "T200":"T200", "Q200":"Q200", "Z200":"Z200", "Z500":"Z500", "Z850":"Z850",
        "OMEGA850":"OMEGA850", "OMEGA500":"OMEGA500", "OMEGA200":"OMEGA200",
        "FLDS":"FLDS", "FLDSC":"FLDSC", "FLNS":"FLNS", "FLNSC":"FLNSC",
        "FSDS":"FSDS", "FSDSC":"FSDSC", "FSNS":"FSNS", "FSNSC":"FSNSC",
        "LHFLX":"LHFLX", "SHFLX":"SHFLX", "QFLX":"QFLX", "TREFHT":"TREFHT", "TS":"TS",
        "PRECSL":"PRECSL", "PRECC":"PRECC", "PRECL":"PRECL",
    },
    "GPCP": {"PRECT":"PRECT"},
    "GPM": {"PRECT":"PRECT"},
    "NOAA-OLR": {"FLUT":"FLUT"},
}

PCC_REGIONS = (
    dict(name="GLOBAL", lat=(-90.,90.), lon=(0.,360.)),
    dict(name="GLOBAL_LAND", lat=(-90.,90.), lon=(0.,360.), mask="land"),
    dict(name="NINO3.4", lat=(-5.,5.), lon=(190.,240.)),
    dict(name="MAR_CONT", lat=(-10.,10.), lon=(95.,150.)),
    dict(name="CONUS", lat=(25.,50.), lon=(235.,294.)),
    dict(name="NH_POLAR", lat=(60.,90.), lon=(0.,360.)),
)
TCC_REGIONS = (
    dict(name="GLOBAL", lat=(-90.,90.), lon=(0.,360.)),
    dict(name="OCEAN", lat=(-90.,90.), lon=(0.,360.), mask="ocean"),
    dict(name="LAND", lat=(-90.,90.), lon=(0.,360.), mask="land"),
    *PCC_REGIONS[2:],
)
LEAD_WINDOWS = {"week1-2":(1,14), "week3-4":(15,28), "week5-6":(29,42), "week7-8":(43,56)}
WEEKLY_WINDOWS = {"week1":(1,7), "week2":(8,14), "week3":(15,21), "week4":(22,28)}
FIFTEEN_DAY_WINDOWS = {"day01-15":(1,15), "day16-30":(16,30), "day31-45":(31,45), "day46-60":(46,60)}
MONTH_WINDOWS = {"month1":(1,31)}
WINDOW_SCHEMES = {"biweekly":LEAD_WINDOWS, "weekly":WEEKLY_WINDOWS, "15day":FIFTEEN_DAY_WINDOWS, "month":MONTH_WINDOWS}
AVAILABLE_PCC_METRICS = ("ACC", "RMSE", "ACC_member", "RMSE_member", "ACC_member_mean", "RMSE_member_mean", "ACC_member_std", "RMSE_member_std", "ACC_member_quantile", "RMSE_member_quantile")
AVAILABLE_PCC_MAP_METRICS = ("ACC_EMEAN", "N_SAMPLES_EMEAN", "ACC_MEMBER_MEAN", "ACC_MEMBER_STD", "P_ACC_POS", "N_SAMPLES_MEMBER_MEAN", "N_ENS_USED", "ACC_MEMBER", "N_SAMPLES_MEMBER")
AVAILABLE_TCC_METRICS = ("TCC", "RMSE", "BIAS", "TCC_ensmean", "RMSE_ensmean", "BIAS_ensmean", "<METRIC>_regional_member", "<METRIC>_member_ens_mean", "<METRIC>_member_ens_std", "<METRIC>_member_qpostXX", "<METRIC>_member_qpreXX", "<METRIC>_regional_ensmean", "<METRIC>_valid_wsum_member", "<METRIC>_valid_frac_member", "<METRIC>_valid_wsum_ensmean", "<METRIC>_valid_frac_ensmean", "region_area_weight_sum")
