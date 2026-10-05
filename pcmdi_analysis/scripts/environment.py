"""Runtime environment helpers shared by plotting workflows."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional


def configure_proj_data(*, verbose: bool = False) -> Optional[Path]:
    """Select a readable, compatible PROJ database for the active Python env.

    Conda kernels can inherit ``PROJ_LIB`` from a different login/compute
    environment.  Prefer pyproj's own data directory, validate it with a CRS
    lookup, and update both modern and legacy environment variables.
    """
    try:
        import pyproj
    except ImportError:
        if verbose:
            print("[WARN] pyproj is unavailable; PROJ data were not configured.")
        return None

    candidates = [
        Path(pyproj.datadir.get_data_dir()),
        Path(sys.prefix) / "share" / "proj",
    ]
    for variable in ("CONDA_PREFIX", "PROJ_DATA", "PROJ_LIB"):
        value = os.environ.get(variable)
        if value:
            base = Path(value)
            candidates.append(base / "share" / "proj" if variable == "CONDA_PREFIX" else base)

    seen: set[Path] = set()
    errors = []
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        if not (candidate / "proj.db").is_file():
            continue
        try:
            pyproj.datadir.set_data_dir(str(candidate))
            os.environ["PROJ_DATA"] = str(candidate)
            os.environ["PROJ_LIB"] = str(candidate)
            pyproj.CRS.from_epsg(4326)
        except Exception as exc:  # Try the next candidate when databases differ.
            errors.append(f"{candidate}: {exc}")
            continue
        if verbose:
            print(f"[INFO] Using PROJ data: {candidate}")
        return candidate

    detail = "; ".join(errors) if errors else "no candidate contained proj.db"
    raise RuntimeError(f"Could not configure a compatible PROJ database: {detail}")
