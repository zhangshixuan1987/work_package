"""Experiments of the F20TR initialization-strategy hindcasts (``hind_lead_eval`` workflow).

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

==========  ==========  ==========================================  ============================
config      in figures  atmosphere                                  land
==========  ==========  ==========================================  ============================
CLIM_BASE   CLIM/CLIM   FREE_CLIM free-running ensemble             FREE_CLIM
HYB_BASE    ERA5/DATM   HYB_DARTPERT20 (ERA5 + DART perturbations)  DATM-forced land (v3-LR.DATM)
DA_ATM      DART/DATM   DART EN40 atmosphere analysis               DATM-forced land
DA_LND      ERA5/DART   ERA5 (v3.LR.ne30pg2.ERA5 eam.i)             DART EN40 land analysis
DA_COUPLED  DART/DART   DART EN40 atmosphere analysis               DART EN40 land analysis
==========  ==========  ==========================================  ============================

The BASE runs use members EN01-EN20; the DA runs use 20 of the 40 DART members, a different
20 in each season (``DA_MEMBERS``; use ``members(config, season)``).
"""
from __future__ import annotations

import json
import os
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

# Names shown in figures: <atmosphere initial state>/<land initial state>. The keys above stay in the
# case, member and saved-file names.
CONFIG_LABELS = {
    "CLIM_BASE": "CLIM/CLIM",
    "HYB_BASE": "ERA5/DATM",
    "DA_ATM": "DART/DATM",
    "DA_LND": "ERA5/DART",
    "DA_COUPLED": "DART/DART",
}


def config_label(key: str) -> str:
    """Figure name of a configuration key (``CONFIG_LABELS``); other names unchanged."""
    return CONFIG_LABELS.get(key, key)


def relabel(text: str) -> str:
    """``text`` with every configuration key replaced by its figure name, e.g. a contrast
    "DA_COUPLED - CLIM_BASE" -> "DART/DART - CLIM/CLIM"."""
    for key in sorted(CONFIG_LABELS, key=len, reverse=True):
        text = re.sub(rf"(?<![A-Za-z0-9_]){re.escape(key)}(?![A-Za-z0-9_])", CONFIG_LABELS[key], text)
    return text

BASE_MEMBERS = tuple(f"EN{k:02d}" for k in range(1, 21))
# Every configuration has 20 members. The DA runs keep the numbers of the 20 DART members they use,
# a different subset of the 40 in each season (DARTEN40 in Dec 2011, DARTEN40S1 in May 2012); the
# BASE perturbations come from the same 20 DART members (EN01 = unperturbed control). The numbers are
# labels only: the members of each ensemble are exchangeable and are not paired across configurations.
DA_MEMBERS = {
    "DJF": tuple(f"EN{k:02d}" for k in (1, 3, 4, 7, 8, 9, 10, 14, 17, 18, 20, 21, 25, 26, 28, 31, 35, 36, 37, 40)),
    "JJA": tuple(f"EN{k:02d}" for k in (2, 3, 5, 7, 11, 12, 13, 14, 15, 16, 20, 22, 24, 26, 27, 28, 31, 32, 34, 38)),
}

# season -> start dates (lagged by 5 days; forecasts run 2 months)
SEASONS = {
    "DJF": ("2011-12-16", "2011-12-21", "2011-12-26", "2011-12-31"),
    "JJA": ("2012-05-16", "2012-05-21", "2012-05-26", "2012-05-31"),
}
FORECAST_MONTHS = 2

# Line styles shared by every hind_lead_eval figure (references in black).
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


def members(config: str, season: str) -> tuple[str, ...]:
    """Member directories of ``config`` in ``season``."""
    if season not in SEASONS:
        raise ValueError(f"Unknown season {season!r}; choose from {list(SEASONS)}")
    return DA_MEMBERS[season] if config.startswith("DA_") else BASE_MEMBERS


def season_of(init) -> str:
    """Season whose start dates include ``init``."""
    d = np.datetime64(init, "D")
    for s, dates in SEASONS.items():
        if any(np.datetime64(x, "D") == d for x in dates):
            return s
    raise ValueError(f"{init} is not a start date of any season")


def init_dates(seasons=tuple(SEASONS)) -> dict[str, tuple[np.datetime64, ...]]:
    """{season: start dates as datetime64[D]}."""
    unknown = [s for s in seasons if s not in SEASONS]
    if unknown:
        raise ValueError(f"Unknown seasons {unknown}; choose from {list(SEASONS)}")
    return {s: tuple(np.datetime64(d, "D") for d in SEASONS[s]) for s in seasons}


# configuration key -> raw case name, without the start-date suffix: the case of a start is
# ``<case name>_<YYYYMMDD>00``. Notebooks set their own with ``set_case_names(CONFIGS)``; the mapping is kept in
# the environment (``CASE_NAMES_ENV``) so that reader processes started afterwards use it too.
CASE_NAMES = {c: f"{CASE_PREFIX}_{c}" for c in CONFIGS}
CASE_NAMES_ENV = "HIND_LEAD_EVAL_CASE_NAMES"


def case_names() -> dict[str, str]:
    """The current key -> raw case name mapping (``CASE_NAMES`` updated by ``set_case_names``)."""
    return {**CASE_NAMES, **json.loads(os.environ.get(CASE_NAMES_ENV, "{}"))}


def set_case_names(names: dict[str, str]) -> dict[str, str]:
    """Use the raw case names ``names`` ({configuration key: case name without ``_<YYYYMMDD>00``}) for the
    configurations they list, the defaults for the others. Returns the full mapping."""
    if not isinstance(names, dict):
        raise TypeError("set_case_names expects a dict {configuration key: raw case name}")
    check_configs(list(names))
    bad = [k for k, v in names.items() if not isinstance(v, str) or not v or "/" in v]
    if bad:
        raise ValueError(f"Invalid case names for {bad}: give the case directory name without the start date")
    os.environ[CASE_NAMES_ENV] = json.dumps({**json.loads(os.environ.get(CASE_NAMES_ENV, "{}")), **names})
    return case_names()


def case_attrs(config: str) -> dict[str, str]:
    """Attributes recording the raw case a saved product of ``config`` was built from."""
    return {"case_name": case_names()[config]}


def case_current(attrs, config: str) -> bool:
    """True when ``attrs`` (of a saved product) name the current case of ``config``. Products saved before
    case names were recorded count as built from the default case."""
    return attrs.get("case_name", CASE_NAMES[config]) == case_names()[config]


def case_stamp(configs) -> str:
    """The non-default case names of ``configs`` (empty with the defaults), for the identity of saved results."""
    changed = {c: case_names()[c] for c in configs if case_names()[c] != CASE_NAMES[c]}
    return f"|cases={json.dumps(changed, sort_keys=True)}" if changed else ""


def cache_current(path, config: str) -> bool:
    """True when the NetCDF file ``path`` exists and was built from the current case of ``config``."""
    path = Path(path)
    if not path.exists():
        return False
    import xarray as xr

    with xr.open_dataset(path) as ds:
        ok = case_current(ds.attrs, config)
    if not ok:
        print(f"[STALE] {path}: built from another case of {config}, recomputing")
    return ok


def case_name(config: str, init) -> str:
    """``<case name of config>_<YYYYMMDD>00``, by default ``F20TR_ne30pg2_r05_IcoswISC30E3r5_<CONFIG>_<YYYYMMDD>00``."""
    check_configs([config])
    return f"{case_names()[config]}_{np.datetime64(init, 'D').astype(object):%Y%m%d}00"


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


def balanced_runs(status, require_all_members: bool = True) -> list[tuple[str, str]]:
    """(config, season) pairs of a ``readiness`` table, the same configurations in every season.

    Both start sets (Dec 2011, May 2012) hold the same experiments, so a configuration missing in one
    season would make the pooled and cross-season comparisons unbalanced; this raises instead.
    """
    summary = status.groupby(["config", "season"])[["ready", "members"]].sum()
    ok = summary.ready >= (summary.members if require_all_members else 1)
    missing = [f"{c} {s} ({r.ready}/{r.members})" for (c, s), r in summary[~ok].iterrows()]
    if missing:
        raise RuntimeError("Not every configuration is ready in every season: " + ", ".join(missing))
    return list(summary.index)


def readiness(configs=CONFIG_KEYS, seasons=tuple(SEASONS), streams=(("atm", "daily"),),
              data_root=DEFAULT_INPUT_DATA_ROOT):
    """Table (config, season, init, stream) -> number of members whose stream is regenerated."""
    import pandas as pd

    rows = []
    for c in configs:
        for s, inits in init_dates(seasons).items():
            for i in inits:
                for comp, freq in streams:
                    n = sum(stream_ready(c, i, m, comp, freq, data_root) for m in members(c, s))
                    rows.append(dict(config=c, season=s, init=str(i), stream=f"{comp}_{freq}",
                                     ready=n, members=len(members(c, s))))
    return pd.DataFrame(rows)
