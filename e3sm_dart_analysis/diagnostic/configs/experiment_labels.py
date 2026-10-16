"""Experiment names shown in figures, following the analysis_atm_da convention.

Every experiment is labelled ``<TYPE><members>-S<cycle>`` (``CTRL10-S0``,
``CAPT10-S1``, ``DART40-S0``, ``CAM80-S1``): S0 is the boreal-winter cycle
(DA Dec 2011, hindcasts Jan-Feb 2012), S1 the boreal-summer cycle (DA May 2012,
hindcasts Jun-Jul 2012). The DA run and the hindcasts started from it share
one name. Workflows that key their inputs on case names (``DARTEN20_15day``,
``CAPTEN10``) keep those keys for paths and caches and use
``experiment_label`` only for display.
"""
from __future__ import annotations

import re

# Case-name stems (DA runs, 15-day hindcasts, short aliases) -> canonical name.
_CASE_ALIASES = {
    "CTRLEN10": "CTRL10-S0", "CTRLEN10_15day": "CTRL10-S0",
    "CTRLEN10s1_15day": "CTRL10-S1", "CTRLEN10S1_15day": "CTRL10-S1",
    "CTRFEN10": "CTRF10-S0", "CTRFEN10S1": "CTRF10-S1",
    "CAPTEN10": "CAPT10-S0", "CAPTEN10_15day": "CAPT10-S0",
    "CAPTEN10S1": "CAPT10-S1", "CAPTEN10S1_15day": "CAPT10-S1",
    "DARTEN10": "DART10-S0", "DARTEN10_15day": "DART10-S0",
    "DARTEN20": "DART20-S0", "DARTEN20_15day": "DART20-S0",
    "DARTEN40": "DART40-S0", "DARTEN40_15day": "DART40-S0",
    "DARTEN40S1": "DART40-S1", "DARTEN40S1_15day": "DART40-S1",
    "CAM6EN80": "CAM80-S0", "CAM6EN80S1": "CAM80-S1",
}
_CANONICAL = re.compile(r"(CTRL|CTRF|CAPT|DART|CAM)\d+-S\d")
_CASE_SUFFIX = re.compile(r"_(F20TR|WCYCL20TR)_.*$")


def experiment_label(name: str) -> str:
    """Canonical display name for an experiment key or case name.

    Accepts canonical names (returned unchanged), case-name stems
    (``DARTEN20_15day``) and full run ids
    (``DARTEN20_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy``). Names that are
    not experiments (``OBS``, ``ERA5``) are returned unchanged.
    """
    name = str(name)
    if _CANONICAL.fullmatch(name):
        return name
    return _CASE_ALIASES.get(_CASE_SUFFIX.sub("", name), name)
