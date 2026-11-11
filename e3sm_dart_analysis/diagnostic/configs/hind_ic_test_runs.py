"""Runs and comparisons of the ERA5 initial-condition tests (single members started 2011-12-16 00Z).

Every run is one member of an F20TR_ne30pg2_r05_IcoswISC30E3r5 hindcast; its directory is
``<root>/<case>_<YYYYMMDD>00/<member>``, where ``<root>`` is the notebook's ``TEST_DATA_ROOT`` for the
test runs and ``INPUT_DATA_ROOT`` for the production runs they are compared with. Atmosphere initial
files are given relative to the notebook's ``IC_ROOT``. Only names and relative layouts are kept here;
the roots depend on the machine and are notebook inputs.

The production members are the ones sharing the tests' other initial states:
* HYB_BASE EN01 is the unperturbed bare-ERA5 reference of the HYB_BASE ensemble (DATM land EN01), the
  state that ERA5_DARTFILL and ERA5_CLIMFILL complete;
* DA_COUPLED, DA_ATM and DA_LND EN40 share the DART EN40 atmosphere and/or land with CLIMO_ATM and
  CLIMO_LND, which replace one component of DA_COUPLED EN40 with the CLIM_BASE (FREE_CLIM) state.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

CASE_PREFIX = "F20TR_ne30pg2_r05_IcoswISC30E3r5"
INIT = np.datetime64("2011-12-16")
DATE_TAG = "{init:%Y-%m-%d}-00000"


@dataclass(frozen=True)
class Run:
    case: str          # case name without the _<YYYYMMDD>00 start suffix
    member: str
    source: str        # "test" (TEST_DATA_ROOT) or "production" (INPUT_DATA_ROOT)
    atm_ic: str        # atmosphere initial state
    lnd_ic: str        # land initial state
    ic_file: str       # atmosphere initial file, relative to IC_ROOT ({date} = YYYY-MM-DD-00000)
    label: str         # short name in figures
    color: str
    linestyle: str = "-"


_DART40 = "e3sm_dart/DARTEN40_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/EN40/archive/rest/{date}/" \
          "DARTEN40_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy.EN40.eam.i.{date}.nc"

RUNS: dict[str, Run] = {
    # ERA5 atmosphere over the DATM-forced land (EN01, unperturbed)
    "ERA5_BARE": Run(f"{CASE_PREFIX}_HYB_BASE", "EN01", "production",
                     "bare ERA5 (T, U, V, Q, PS, CLDLIQ, CLDICE only)", "DATM EN01",
                     "e3sm_dart/ic_atm_era5_dartpert20/{date}/HYB_DARTPERT20.EN01.eam.i.{date}.nc",
                     "ERA5 bare (HYB_BASE)", "0.35"),
    "ERA5_DARTFILL": Run(f"{CASE_PREFIX}_ERA5COMPLETE_TEST", "EN01", "test",
                         "ERA5 + DART-mean fields", "DATM EN01",
                         "e3sm_dart/ic_atm_era5complete_dartpert20/{date}/ERA5COMPLETE_DARTPERT20.EN01.eam.i.{date}.nc",
                         "ERA5 + DART fill", "tab:blue"),
    "ERA5_CLIMFILL": Run(f"{CASE_PREFIX}_ERA5CLIMFILL_TEST", "EN01", "test",
                         "ERA5 + CLIM_BASE fields", "DATM EN01",
                         "e3sm_dart/ic_atm_era5climfill_dartpert20/{date}/ERA5CLIMFILL_DARTPERT20.EN01.eam.i.{date}.nc",
                         "ERA5 + CLIM fill", "tab:cyan", "--"),
    # DART EN40 joint state and its one-component replacements
    "DA_COUPLED": Run(f"{CASE_PREFIX}_DA_COUPLED", "EN40", "production", "DART EN40", "DART EN40",
                      _DART40, "DA_COUPLED (DART/DART)", "tab:red"),
    "DA_ATM": Run(f"{CASE_PREFIX}_DA_ATM", "EN40", "production", "DART EN40", "DATM EN40",
                  _DART40, "DA_ATM (DART/DATM)", "tab:orange"),
    "DA_LND": Run(f"{CASE_PREFIX}_DA_LND", "EN40", "production",
                  "bare ERA5 (T, U, V, Q, PS, CLDLIQ, CLDICE only)", "DART EN40",
                  "DART_INIT/ATM/v3.LR.ne30pg2.ERA5.eam.i.{date}.nc", "DA_LND (bare ERA5/DART)", "tab:green"),
    "CLIMO_ATM": Run(f"{CASE_PREFIX}_CLIMO_ATM_TEST", "EN40", "test",
                     "CLIM_BASE + (DART EN40 - DART mean)", "DART EN40",
                     "e3sm_dart/ic_climoatm_dartpert20/{date}/CLIMOATM_DARTPERT20.EN40.eam.i.{date}.nc",
                     "CLIMO_ATM (CLIM/DART)", "tab:purple"),
    "CLIMO_LND": Run(f"{CASE_PREFIX}_CLIMO_LND_TEST", "EN40", "test", "DART EN40", "CLIM_BASE (FREE_CLIM)",
                     _DART40, "CLIMO_LND (DART/CLIM)", "tab:brown"),
}

# Extra atmosphere initial states shown in the initial-state check only (no run of their own here)
EXTRA_IC_FILES = {
    "CLIM_BASE state": "e3sm_dart/F20TR_ne30pg2_r05_IcoswISC30E3r5.FREE_CLIM/archive/rest/{date}/"
                       "F20TR_ne30pg2_r05_IcoswISC30E3r5.FREE_CLIM.eam.i.{date}.nc",
}
# ERA5 fields the completion keeps (H2OLNZ is rescaled to the ERA5 Q, so it is not among them)
ERA5_FIELDS = ("T", "U", "V", "Q", "PS", "CLDLIQ", "CLDICE")

# Comparisons: (run A, run B, question). Differences are A - B.
COMPARISONS: dict[str, tuple[str, str, str]] = {
    "tracer_fix": ("ERA5_DARTFILL", "ERA5_BARE",
                   "Does completing the ERA5 state remove the cold start? (ERA5 + DART fill - bare ERA5, DATM land)"),
    "fill_source": ("ERA5_CLIMFILL", "ERA5_DARTFILL",
                    "Does the source of the filled fields matter? (CLIM_BASE fill - DART-mean fill)"),
    "atm_clim": ("CLIMO_ATM", "DA_COUPLED",
                 "Atmosphere initial state, complete states: CLIM_BASE - DART (both over DART land)"),
    "atm_bare_era5": ("DA_LND", "DA_COUPLED",
                      "Atmosphere initial state as run in production DA_LND: bare ERA5 - DART (DART land)"),
    "lnd_clim": ("CLIMO_LND", "DA_COUPLED",
                 "Land initial state: CLIM_BASE - DART (both under the DART atmosphere)"),
    "lnd_datm": ("DA_ATM", "DA_COUPLED",
                 "Land initial state as run in production DA_ATM: DATM-forced - DART (DART atmosphere)"),
}

# Run groups plotted together (one figure per group)
GROUPS = {
    "era5_completion": ("ERA5_BARE", "ERA5_DARTFILL", "ERA5_CLIMFILL"),
    "component_swap": ("DA_COUPLED", "CLIMO_ATM", "DA_LND", "CLIMO_LND", "DA_ATM"),
}


def check_runs(keys) -> tuple[str, ...]:
    unknown = [k for k in keys if k not in RUNS]
    if unknown:
        raise KeyError(f"unknown runs {unknown}; choose from {list(RUNS)}")
    return tuple(keys)


def date_tag(init=INIT) -> str:
    return DATE_TAG.format(init=np.datetime64(init, "D").astype(object))


def case_dir(case: str, init, root) -> Path:
    return Path(root) / f"{case}_{np.datetime64(init, 'D').astype(object):%Y%m%d}00"


def member_dir(key: str, roots: dict, init=INIT) -> Path:
    """``roots`` = {"test": TEST_DATA_ROOT, "production": INPUT_DATA_ROOT}."""
    r = RUNS[key]
    return case_dir(r.case, init, roots[r.source]) / r.member


def ic_path(key: str, ic_root, init=INIT) -> Path:
    rel = RUNS[key].ic_file if key in RUNS else EXTRA_IC_FILES[key]
    return Path(ic_root) / rel.format(date=date_tag(init))


def ensemble_members(case: str, init, root) -> list[str]:
    """Members present in a production case (the DA cases keep their DART member numbers)."""
    return sorted(p.name for p in case_dir(case, init, root).glob("EN[0-9][0-9]") if p.is_dir())
