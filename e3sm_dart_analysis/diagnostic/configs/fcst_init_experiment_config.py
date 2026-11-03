"""Experiments of the F20TR initialization-strategy hindcasts (``fcst_init_eval`` workflow).

Five ways of initializing the same E3SMv3 LR F20TR (prescribed SST/sea ice) model,
each run as a 20-member, 2-month ensemble from 8 start dates (4 lagged starts per
season, 5 days apart). One case per configuration and start date, one directory
per member::

    INPUT_DATA_ROOT/F20TR_ne30pg2_r05_IcoswISC30E3r5_<CONFIG>_<YYYYMMDDHH>/<ENxx>/
        archive/{atm,lnd,rof}/hist/       raw EAM/ELM/MOSART history
        post/<comp>/180x360_aave/ts/<freq>/0yr/<VAR>_<YYYYMM>_<YYYYMM>.nc

``post`` holds one file per variable over the whole forecast (atm daily and 6-hourly,
lnd and rof daily), written by ``post_process/regen_f20tr_ts.sh``; the zppy-style
``ncclimo --split`` files they replace stopped after 4 days (atm) or at the end of
the start month (lnd). Each member records the check in
``post/scripts/regen_ts_<stream>.status`` (``OK files=... records=...``).

Initial states (from each case's ``user_nl_eam`` ncdata and ``user_nl_elm`` finidat):

==========  ==========================================  ============================
config      atmosphere                                  land
==========  ==========================================  ============================
CLIM_BASE   FREE_CLIM free-running ensemble             FREE_CLIM
HYB_BASE    HYB_DARTPERT20 (ERA5 + DART perturbations)  DATM-forced land (v3-LR.DATM)
DA_ATM      DART EN40 atmosphere analysis               DATM-forced land
DA_LND      ERA5 (v3.LR.ne30pg2.ERA5 eam.i)             DART EN40 land analysis
DA_COUPLED  DART EN40 atmosphere analysis               DART EN40 land analysis
==========  ==========================================  ============================

The BASE runs use members EN01-EN20; the DA runs use the 20 of the 40 DART members
listed in ``DA_MEMBERS``.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

DEFAULT_INPUT_DATA_ROOT = Path("/global/cfs/cdirs/m4849/zhan391/eam_hindcast")
CASE_PREFIX = "F20TR_ne30pg2_r05_IcoswISC30E3r5"
CASE_TEMPLATE = CASE_PREFIX + "_{config}_{init:%Y%m%d}00"

# key -> (atmosphere initial state, land initial state)
CONFIGS = {
    "CLIM_BASE": ("FREE_CLIM", "FREE_CLIM"),
    "HYB_BASE": ("ERA5+DART pert.", "DATM-forced"),
    "DA_ATM": ("DART", "DATM-forced"),
    "DA_LND": ("ERA5", "DART"),
    "DA_COUPLED": ("DART", "DART"),
}
CONFIG_KEYS = tuple(CONFIGS)

BASE_MEMBERS = tuple(f"EN{k:02d}" for k in range(1, 21))
DA_MEMBERS = tuple(f"EN{k:02d}" for k in (1, 3, 4, 7, 8, 9, 10, 14, 17, 18, 20, 21, 25, 26, 28, 31, 35, 36, 37, 40))
MEMBERS = {c: (DA_MEMBERS if c.startswith("DA_") else BASE_MEMBERS) for c in CONFIGS}

# season -> start dates (lagged by 5 days; forecasts run 2 months)
SEASONS = {
    "DJF": ("2011-12-16", "2011-12-21", "2011-12-26", "2011-12-31"),
    "JJA": ("2012-05-16", "2012-05-21", "2012-05-26", "2012-05-31"),
}
FORECAST_MONTHS = 2

# Line styles shared by every fcst_init_eval figure (references in black).
STYLE_MAP = {
    "ERA5": dict(color="k", linestyle="-", linewidth=2.4, zorder=6),
    "OBS": dict(color="k", linestyle="-", linewidth=2.4, zorder=6),
    "CLIM_BASE": dict(color="0.55", linestyle="-", linewidth=2.0, zorder=3),
    "HYB_BASE": dict(color="tab:blue", linestyle="-", linewidth=2.0, zorder=4),
    "DA_ATM": dict(color="tab:orange", linestyle="-", linewidth=2.0, zorder=4),
    "DA_LND": dict(color="tab:green", linestyle="-", linewidth=2.0, zorder=4),
    "DA_COUPLED": dict(color="tab:red", linestyle="-", linewidth=2.4, zorder=5),
}

# post-processed streams: (component, frequency) -> ts subdirectory
STREAMS = {
    ("atm", "daily"): "atm/180x360_aave/ts/daily/0yr",
    ("atm", "6hourly"): "atm/180x360_aave/ts/6hourly/0yr",
    ("lnd", "daily"): "lnd/180x360_aave/ts/daily/0yr",
    ("rof", "daily"): "rof/180x360_aave/ts/daily/0yr",
}
_STATUS = {("atm", "daily"): "atm_daily", ("atm", "6hourly"): "atm_6hourly",
           ("lnd", "daily"): "lnd_daily", ("rof", "daily"): "rof_daily"}


def check_configs(configs) -> tuple[str, ...]:
    unknown = [c for c in configs if c not in CONFIGS]
    if unknown:
        raise ValueError(f"Unknown configurations {unknown}; choose from {list(CONFIGS)}")
    return tuple(configs)


def init_dates(seasons=tuple(SEASONS)) -> dict[str, tuple[np.datetime64, ...]]:
    """{season: start dates as datetime64[D]}."""
    unknown = [s for s in seasons if s not in SEASONS]
    if unknown:
        raise ValueError(f"Unknown seasons {unknown}; choose from {list(SEASONS)}")
    return {s: tuple(np.datetime64(d, "D") for d in SEASONS[s]) for s in seasons}


def case_name(config: str, init) -> str:
    """``F20TR_ne30pg2_r05_IcoswISC30E3r5_<CONFIG>_<YYYYMMDD>00``."""
    check_configs([config])
    return CASE_TEMPLATE.format(config=config, init=np.datetime64(init, "D").astype(object))


def member_dir(config: str, init, member: str, data_root=DEFAULT_INPUT_DATA_ROOT) -> Path:
    return Path(data_root) / case_name(config, init) / member


def ts_dir(config: str, init, member: str, component: str = "atm", freq: str = "daily",
           data_root=DEFAULT_INPUT_DATA_ROOT) -> Path:
    if (component, freq) not in STREAMS:
        raise KeyError(f"No {component} {freq} stream; choose from {list(STREAMS)}")
    return member_dir(config, init, member, data_root) / "post" / STREAMS[component, freq]


_RANGE = re.compile(r"_(\d{6})_(\d{6})\.nc$")


def ts_file(config: str, init, member: str, variable: str, component: str = "atm", freq: str = "daily",
            data_root=DEFAULT_INPUT_DATA_ROOT) -> Path:
    """The single ``<VAR>_<YYYYMM>_<YYYYMM>.nc`` file of a member's stream."""
    d = ts_dir(config, init, member, component, freq, data_root)
    hits = [f for f in d.glob(f"{variable}_*.nc") if _RANGE.search(f.name)
            and f.name[: -len(_RANGE.search(f.name).group(0))] == variable]
    if len(hits) != 1:
        raise FileNotFoundError(f"{'No' if not hits else 'Several'} {variable} files in {d}: {[h.name for h in hits]}")
    return hits[0]


def stream_ready(config: str, init, member: str, component: str = "atm", freq: str = "daily",
                 data_root=DEFAULT_INPUT_DATA_ROOT) -> bool:
    """True when ``regen_f20tr_ts.sh`` finished this member's stream (status OK)."""
    status = member_dir(config, init, member, data_root) / "post/scripts" / f"regen_ts_{_STATUS[component, freq]}.status"
    return status.is_file() and status.read_text().startswith("OK")


def readiness(configs=CONFIG_KEYS, seasons=tuple(SEASONS), streams=(("atm", "daily"),),
              data_root=DEFAULT_INPUT_DATA_ROOT):
    """Table (config, season, init, stream) -> number of members whose stream is regenerated."""
    import pandas as pd

    rows = []
    for c in configs:
        for s, inits in init_dates(seasons).items():
            for i in inits:
                for comp, freq in streams:
                    n = sum(stream_ready(c, i, m, comp, freq, data_root) for m in MEMBERS[c])
                    rows.append(dict(config=c, season=s, init=str(i), stream=f"{comp}_{freq}",
                                     ready=n, members=len(MEMBERS[c])))
    return pd.DataFrame(rows)
