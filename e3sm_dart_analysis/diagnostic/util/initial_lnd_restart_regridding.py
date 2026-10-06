"""Conservative regridding utilities for aggregated ELM restart fields."""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable

import numpy as np
from netCDF4 import Dataset


CORNER_PERMUTATIONS = {
    "keep": None,
    "0132": (0, 1, 3, 2),
    "0321": (0, 3, 2, 1),
}


def _as_path(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()


def _require_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")


def _temporary_output(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent, delete=False
    )
    handle.close()
    temporary = Path(handle.name)
    temporary.unlink()
    return temporary


def write_unstructured_scrip(
    aggregated_file: str | Path,
    scrip_file: str | Path,
    *,
    corner_permutation: str = "keep",
    overwrite: bool = False,
) -> Path:
    """Write an unstructured SCRIP grid from an aggregated restart file."""
    source = _as_path(aggregated_file)
    destination = _as_path(scrip_file)
    _require_file(source, "Aggregated restart")
    if corner_permutation not in CORNER_PERMUTATIONS:
        choices = ", ".join(CORNER_PERMUTATIONS)
        raise ValueError(f"corner_permutation must be one of: {choices}")
    if destination.exists() and not overwrite:
        return destination

    with Dataset(source, "r") as dataset:
        missing = [name for name in ("lat", "lon", "lat_b", "lon_b") if name not in dataset.variables]
        if missing:
            raise KeyError(
                f"{source} cannot be conservatively regridded; missing variables: {', '.join(missing)}"
            )
        lat = np.asarray(dataset.variables["lat"][:], dtype=np.float64)
        lon = np.asarray(dataset.variables["lon"][:], dtype=np.float64)
        lat_b = np.asarray(dataset.variables["lat_b"][:], dtype=np.float64)
        lon_b = np.asarray(dataset.variables["lon_b"][:], dtype=np.float64)
        frac = (
            np.asarray(dataset.variables["frac"][:], dtype=np.float64)
            if "frac" in dataset.variables
            else None
        )

    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape:
        raise ValueError(f"lat/lon must be matching 1-D arrays; got {lat.shape} and {lon.shape}")
    if lat_b.ndim != 2 or lon_b.shape != lat_b.shape or lat_b.shape[0] != lat.size:
        raise ValueError(
            "lat_b/lon_b must have matching (gridcell, corner) shapes and agree with lat/lon; "
            f"got {lat_b.shape}, {lon_b.shape}, and {lat.shape}"
        )

    permutation = CORNER_PERMUTATIONS[corner_permutation]
    if permutation is not None:
        if lat_b.shape[1] != len(permutation):
            raise ValueError(
                f"Corner permutation {corner_permutation} requires four corners; got {lat_b.shape[1]}"
            )
        lat_b = lat_b[:, permutation]
        lon_b = lon_b[:, permutation]

    temporary = _temporary_output(destination)
    try:
        with Dataset(temporary, "w", format="NETCDF3_CLASSIC") as output:
            output.createDimension("grid_rank", 1)
            output.createDimension("grid_corners", lat_b.shape[1])
            output.createDimension("grid_size", lat.size)

            grid_dims = output.createVariable("grid_dims", "i4", ("grid_rank",))
            grid_dims[:] = lat.size
            center_lat = output.createVariable("grid_center_lat", "f8", ("grid_size",))
            center_lon = output.createVariable("grid_center_lon", "f8", ("grid_size",))
            corner_lat = output.createVariable(
                "grid_corner_lat", "f8", ("grid_size", "grid_corners")
            )
            corner_lon = output.createVariable(
                "grid_corner_lon", "f8", ("grid_size", "grid_corners")
            )
            mask = output.createVariable("grid_imask", "i4", ("grid_size",))

            center_lat.units = corner_lat.units = "degrees"
            center_lon.units = corner_lon.units = "degrees"
            center_lat[:] = lat
            center_lon[:] = lon
            corner_lat[:] = lat_b
            corner_lon[:] = lon_b
            if frac is None:
                mask[:] = np.ones(lat.size, dtype=np.int32)
            else:
                if frac.shape != lat.shape:
                    raise ValueError(f"frac must have shape {lat.shape}; got {frac.shape}")
                mask[:] = np.isfinite(frac) & (frac > 0.0)

        temporary.replace(destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


def _run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def conservative_regrid_restart(
    aggregated_file: str | Path,
    destination_scrip: str | Path,
    source_scrip: str | Path,
    map_file: str | Path,
    output_file: str | Path,
    *,
    corner_permutation: str = "keep",
    variables: Iterable[str] | None = None,
    force: bool = False,
    debug_level: int = 1,
) -> Path:
    """Create/reuse conservative weights and regrid an aggregated restart."""
    aggregated = _as_path(aggregated_file)
    destination_grid = _as_path(destination_scrip)
    source_grid = _as_path(source_scrip)
    weights = _as_path(map_file)
    output = _as_path(output_file)
    _require_file(aggregated, "Aggregated restart")
    _require_file(destination_grid, "Destination SCRIP grid")

    ncremap = shutil.which("ncremap")
    if ncremap is None:
        raise RuntimeError("ncremap is not available on PATH. Load an NCO environment first.")

    write_unstructured_scrip(
        aggregated,
        source_grid,
        corner_permutation=corner_permutation,
        overwrite=force,
    )

    if force or not weights.is_file():
        temporary_map = _temporary_output(weights)
        try:
            _run(
                [
                    ncremap,
                    "-a",
                    "aave",
                    "-s",
                    str(source_grid),
                    "-g",
                    str(destination_grid),
                    "-m",
                    str(temporary_map),
                ]
            )
            if not temporary_map.is_file():
                raise RuntimeError(f"ncremap did not create the expected map: {temporary_map}")
            temporary_map.replace(weights)
        finally:
            if temporary_map.exists():
                temporary_map.unlink()

    if output.is_file() and not force:
        print(f"Reusing existing regridded restart: {output}")
        return output

    temporary_output = _temporary_output(output)
    command = [
        ncremap,
        "-D",
        str(debug_level),
        "-R",
        "--rgr col_nm=gridcell,lat_nm_in=lat,lon_nm_in=lon",
    ]
    selected = list(variables or [])
    if selected:
        command.extend(["-v", ",".join(selected)])
    command.extend(
        [
            "-m",
            str(weights),
            "-i",
            str(aggregated),
            "-o",
            str(temporary_output),
        ]
    )
    try:
        _run(command)
        if not temporary_output.is_file():
            raise RuntimeError(f"ncremap did not create the expected output: {temporary_output}")
        temporary_output.replace(output)
    finally:
        if temporary_output.exists():
            temporary_output.unlink()
    return output
