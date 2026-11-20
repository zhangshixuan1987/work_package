"""Experiment names shown in figures, following the analysis_atm_da convention.

Every experiment is labelled ``<TYPE><members>-S<cycle>`` (``CTRL10-S0``,
``RHYB10-S1``, ``DART40-S0``, ``CAM80-S1``): S0 is the boreal-winter cycle
(DA Dec 2011, hindcasts Jan-Feb 2012), S1 the boreal-summer cycle (DA May 2012,
hindcasts Jun-Jul 2012). The DA run and the hindcasts started from it share
one name. Workflows that key their inputs on case names (``DARTEN20_15day``,
``CAPTEN10``) keep those keys for paths and caches and use
``experiment_label`` only for display.

Hindcast figures (``fcst_*`` and ``hind_*`` workflows) write the cycle as
``H<cycle>`` instead (``CTRL10-H0``, ``DART40-H1``) so they are not mistaken
for the DA runs; ``hindcast_label`` gives that name. The ``S<cycle>`` names
stay the keys of every registry, path and cache.

Two hindcast sets start from ERA5 (2026-11-19 naming):

- ``REAN10`` (keys ``REAN10-S0``/``REAN10-S1``, case names ``ERA5EN10*_15day``):
  the ERA5 reanalysis initial state.
- ``RHYB10`` (keys ``CAPT10-S0``/``CAPT10-S1``, case names ``CAPTEN10*``):
  the hybrid ERA5-based state, shown as ``RHYB`` (``RHYB10-H0``) and not as CAPT
  so it is not mistaken for the official CAPT protocol; the ``CAPT`` keys stay
  the keys of registries, paths and caches.
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
    "ERA5EN10": "REAN10-S0", "ERA5EN10_15day": "REAN10-S0",
    "ERA5EN10S1": "REAN10-S1", "ERA5EN10S1_15day": "REAN10-S1",
    "DARTEN10": "DART10-S0", "DARTEN10_15day": "DART10-S0",
    "DARTEN20": "DART20-S0", "DARTEN20_15day": "DART20-S0",
    "DARTEN40": "DART40-S0", "DARTEN40_15day": "DART40-S0",
    "DARTEN40S1": "DART40-S1", "DARTEN40S1_15day": "DART40-S1",
    "CAM6EN80": "CAM80-S0", "CAM6EN80S1": "CAM80-S1",
}
_CANONICAL = re.compile(r"(CTRL|CTRF|CAPT|REAN|DART|CAM)\d+-S\d")
_CASE_SUFFIX = re.compile(r"_(F20TR|WCYCL20TR)_.*$")

# Experiment types renamed for display (key type -> figure type).
_DISPLAY_TYPE = {"CAPT": "RHYB"}
_DISPLAY_TYPE_RE = re.compile(r"\b(" + "|".join(_DISPLAY_TYPE) + r")(?=\d+-[SH]\d\b|[AL]-)")


def display_type(text: str) -> str:
    """Rename key experiment types in ``text`` for display (``CAPT10-S0`` -> ``RHYB10-S0``,
    ``CAPTA-DARTL-H0`` -> ``RHYBA-DARTL-H0``); other text is unchanged."""
    return _DISPLAY_TYPE_RE.sub(lambda m: _DISPLAY_TYPE[m.group(1)], str(text))


def experiment_label(name: str) -> str:
    """Display name for an experiment key or case name.

    Accepts canonical names (``CAPT10-S0`` -> ``RHYB10-S0``, others unchanged), case-name stems
    (``DARTEN20_15day``) and full run ids
    (``DARTEN20_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy``). Names that are
    not experiments (``OBS``, ``ERA5``) are returned unchanged.
    """
    name = str(name)
    if not _CANONICAL.fullmatch(name):
        name = _CASE_ALIASES.get(_CASE_SUFFIX.sub("", name), name)
    return display_type(name)


_CYCLE = re.compile(r"(?<=\d)-S(\d)\b")


def hindcast_label(name: str) -> str:
    """Hindcast display name: ``experiment_label`` with the cycle as ``H<n>``.

    ``CTRL10-S0`` and ``CTRLEN10_15day`` give ``CTRL10-H0``; ``DART40-S1``
    gives ``DART40-H1``, ``CAPT10-S0`` gives ``RHYB10-H0``, ``REAN10-S0`` gives ``REAN10-H0``. Also rewrites names
    inside longer text (``"Jan2012:CTRL10-S0"`` -> ``"Jan2012:CTRL10-H0"``). Names that are not
    experiments are returned unchanged.
    """
    return display_type(_CYCLE.sub(r"-H\1", experiment_label(name)))
