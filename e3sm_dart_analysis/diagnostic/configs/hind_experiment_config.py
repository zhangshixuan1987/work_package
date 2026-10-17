"""Experiments of the 2012-2013 AMIP hindcasts (``hind_*`` workflows).

Each experiment is one continuous E3SMv3 LR AMIP run started on 2012-01-01
from a different atmosphere/land initial state; ``CTRL`` starts both from
the piControl. Runs live in ``INPUT_DATA_ROOT/<case>`` with zppy output
under ``post/atm/180x360_aave`` (``ts/daily/1yr``, ``ts/6hourly/1yr``,
``ts/monthly/1yr/<YYYY>_<YYYY>`` on 19 pressure levels, and CMIP-named
``cmip_ts/monthly``).
"""
from __future__ import annotations

from pathlib import Path

DEFAULT_INPUT_DATA_ROOT = Path("/global/cfs/cdirs/m4849/zhan391/eam_hindcast")
CASE_TEMPLATE = "20260114.v3-LR-AMIP-{tag}.ne30pg2_r05_IcoswISC30E3r5"
PERIOD = "2012-2013"
START_DATE = "2012-01-01"
END_DATE = "2013-12-31"

# key -> (case tag, initial state)
EXPERIMENTS = {
    "CTRL": ("CTRL1", "piControl atm + piControl lnd"),
    "DARTOP": ("DART1", "DART atm + DART lnd"),
    "ERA5OP": ("ERA51", "ERA5 atm + ERA5-forced lnd"),
    "ERA5DL": ("ERA5DL", "ERA5 atm + DART lnd"),
    "ERA5BL": ("ERA5BL", "ERA5 atm + piControl lnd"),
    "DARTBL": ("DARTBL", "DART atm + piControl lnd"),
    "DARTBA": ("DARTBA", "piControl atm + DART lnd"),
    "DARTAE": ("DARTAE", "DART atm + ERA5-forced lnd"),
}
EXPERIMENT_KEYS = tuple(EXPERIMENTS)

# Line styles shared by every hind_* figure (references in black).
STYLE_MAP = {
    "ERA5": dict(color="k", linestyle="-", linewidth=2.6, zorder=6),
    "ERAI": dict(color="k", linestyle="--", linewidth=2.6, zorder=6),
    "OBS": dict(color="k", linestyle="-", linewidth=2.6, zorder=6),
    "CTRL": dict(color="0.55", linestyle="-", linewidth=2.2, zorder=5),
    "DARTOP": dict(color="tab:red", linestyle="-", linewidth=2.4, zorder=6),
    "ERA5OP": dict(color="tab:red", linestyle="--", linewidth=2.4, zorder=6),
    "ERA5DL": dict(color="tab:blue", linestyle="-", linewidth=1.8, zorder=4),
    "ERA5BL": dict(color="tab:blue", linestyle=":", linewidth=2.0, zorder=4),
    "DARTBL": dict(color="tab:purple", linestyle="-", linewidth=1.8, zorder=4),
    "DARTBA": dict(color="tab:purple", linestyle="--", linewidth=1.8, zorder=4),
    "DARTAE": dict(color="tab:purple", linestyle=":", linewidth=2.0, zorder=4),
}

# Pressure levels (Pa) of the monthly model output; references are subset to them.
MODEL_PLEV_PA = (
    100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000,
    15000, 10000, 7000, 5000, 3000, 2000, 1000, 500, 100,
)

ATM_POST = Path("post/atm/180x360_aave")
SUBDIRS = {
    "daily": ATM_POST / "ts/daily/1yr",
    "6hourly": ATM_POST / "ts/6hourly/1yr",
    "monthly": ATM_POST / "ts/monthly/1yr",
    "cmip_monthly": ATM_POST / "cmip_ts/monthly",
}


def case_name(experiment: str) -> str:
    """Full case name of an experiment key (``CTRL`` -> ``20260114.v3-LR-AMIP-CTRL1...``)."""
    if experiment not in EXPERIMENTS:
        raise KeyError(f"Unknown experiment {experiment!r}; choose from {list(EXPERIMENTS)}")
    return CASE_TEMPLATE.format(tag=EXPERIMENTS[experiment][0])


def experiment_dir(experiment: str, data_root=DEFAULT_INPUT_DATA_ROOT) -> Path:
    return Path(data_root) / case_name(experiment)


def output_dir(experiment: str, kind: str, data_root=DEFAULT_INPUT_DATA_ROOT) -> Path:
    """zppy output directory of ``kind`` (daily, 6hourly, monthly, cmip_monthly)."""
    if kind not in SUBDIRS:
        raise KeyError(f"Unknown output kind {kind!r}; choose from {list(SUBDIRS)}")
    return experiment_dir(experiment, data_root) / SUBDIRS[kind]


def check_experiments(experiments) -> tuple[str, ...]:
    unknown = [e for e in experiments if e not in EXPERIMENTS]
    if unknown:
        raise ValueError(f"Unknown experiments {unknown}; choose from {list(EXPERIMENTS)}")
    return tuple(experiments)
