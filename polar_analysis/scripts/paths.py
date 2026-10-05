"""Locations of the v3 LE paper data, diagnostic products, and figures.

The root is read from ``config/paths.json`` in this package and can be
overridden for a session with the ``V3LE_ROOT`` environment variable::

    V3LE_ROOT=/some/other/v3LE_paper jupyter lab

Importing this module has no side effects on the filesystem.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Union

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
CONFIG_FILE = PACKAGE_ROOT / "config" / "paths.json"

with CONFIG_FILE.open() as _f:
    _CONFIG = json.load(_f)

V3LE_ROOT = Path(os.environ.get("V3LE_ROOT") or _CONFIG["v3le_root"]).expanduser()
V3LE_DATA_DIR = V3LE_ROOT / _CONFIG["subdirs"]["data"]
V3LE_DIAG_DIR = V3LE_ROOT / _CONFIG["subdirs"]["diag_data"]
V3LE_FIG_ROOT = V3LE_ROOT / _CONFIG["subdirs"]["figures"]

PathLike = Union[str, os.PathLike]


def fig_dir(*parts: PathLike) -> Path:
    """Return a figure subdirectory under ``V3LE_FIG_ROOT`` (not created)."""
    return V3LE_FIG_ROOT.joinpath(*map(str, parts))


def require(*paths: PathLike) -> None:
    """Raise one FileNotFoundError listing every path that does not exist."""
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required input(s):\n  " + "\n  ".join(missing)
            + f"\n(V3LE_ROOT = {V3LE_ROOT})"
        )


def describe() -> str:
    """Human-readable summary of the resolved locations."""
    source = "V3LE_ROOT env" if os.environ.get("V3LE_ROOT") else str(CONFIG_FILE)
    return (
        f"V3LE_ROOT     = {V3LE_ROOT}  [{source}]\n"
        f"V3LE_DATA_DIR = {V3LE_DATA_DIR}\n"
        f"V3LE_DIAG_DIR = {V3LE_DIAG_DIR}\n"
        f"V3LE_FIG_ROOT = {V3LE_FIG_ROOT}"
    )


__all__ = [
    "V3LE_ROOT", "V3LE_DATA_DIR", "V3LE_DIAG_DIR", "V3LE_FIG_ROOT",
    "fig_dir", "require", "describe",
]
