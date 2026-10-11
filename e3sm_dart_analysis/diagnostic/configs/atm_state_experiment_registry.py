"""Canonical DART state-space experiment metadata.

Keep native-summary path information here so profile and spatial diagnostics do
not silently use different aliases, run IDs, or DART directory keys.
"""

from __future__ import annotations


ATM_STATE_EXPERIMENTS = {
    # CTRL runs DART without assimilating anything: its posterior equals its
    # prior, so it is the no-DA baseline on the same files and grid.
    "CTRL10-S0": {
        "run_id": "CTRLEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        "dart_key": "dart_en10",
        "label": "CTRL10-S0",
        "color": "tab:blue",
        "marker": "o",
    },
    "DART10-S0": {
        "run_id": "DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        "dart_key": "dart_en10",
        "label": "DART10-S0",
        "color": "tab:orange",
        "marker": "s",
    },
    "DART20-S0": {
        "run_id": "DARTEN20_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        "dart_key": "dart_en20",
        "label": "DART20-S0",
        "color": "tab:green",
        "marker": "^",
    },
    "DART40-S0": {
        "run_id": "DARTEN40_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        "dart_key": "dart_en40",
        "label": "DART40-S0",
        "color": "tab:red",
        "marker": "D",
    },
    "DART40-S1": {
        "run_id": "DARTEN40S1_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        "dart_key": "dart_en40",
        "label": "DART40-S1",
        "color": "tab:red",
        "marker": "D",
    },
}


def state_experiments(*names: str) -> dict[str, dict[str, str]]:
    """Return independent metadata dictionaries for the requested experiments."""
    missing = [name for name in names if name not in ATM_STATE_EXPERIMENTS]
    if missing:
        raise KeyError(f"Unknown state-space experiment(s): {', '.join(missing)}")
    return {name: dict(ATM_STATE_EXPERIMENTS[name]) for name in names}
