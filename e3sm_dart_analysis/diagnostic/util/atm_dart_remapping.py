"""Consistent native-to-pressure-grid remapping for EAM-DART summaries.

Every experiment and summary product follows the same two-stage NCO path:

1. conservative horizontal remapping from the native ne30 grid to 1 degree;
2. logarithmic interpolation to the canonical 19 pressure levels.

Targets are written atomically and validated before they become visible.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import xarray as xr


PRESSURE_LEVELS_PA = np.asarray(
    [
        100000.0, 92500.0, 85000.0, 70000.0, 60000.0, 50000.0,
        40000.0, 30000.0, 25000.0, 20000.0, 15000.0, 10000.0,
        7000.0, 5000.0, 3000.0, 2000.0, 1000.0, 500.0, 100.0,
    ],
    dtype=np.float64,
)

SUMMARY_PRODUCTS = {
    "forecast_mean": ("e", "forecast_mean"),
    "forecast_sd": ("e", "forecast_sd"),
    "output_mean": ("i", "output_mean"),
    "output_sd": ("i", "output_sd"),
}

REQUIRED_STATE_VARIABLES = {"PS", "U", "V", "T", "Q", "CLDLIQ", "CLDICE"}
_NETCDF_LOCK = threading.Lock()


@dataclass(frozen=True)
class RemapRequest:
    experiment: str
    timestamp: str
    product: str
    source: Path
    target: Path


def native_summary_path(
    input_root: Path,
    experiment: Mapping[str, str],
    timestamp: str,
    product: str,
) -> Path:
    """Return the native DART summary path for one product and cycle."""
    try:
        stream, native_product = SUMMARY_PRODUCTS[product]
    except KeyError as error:
        raise ValueError(f"Unsupported DART summary product: {product}") from error
    run_id = experiment["run_id"]
    return (
        Path(input_root) / run_id / experiment["dart_key"] / "eam" / timestamp
        / f"{run_id}.dart.{stream}.eam_{native_product}.{timestamp}.nc"
    )


def ensemble_mean_source(source: Path) -> Path:
    """Native *_mean summary that pairs with a native *_sd summary."""
    source = Path(source)
    name = source.name.replace("_forecast_sd.", "_forecast_mean.").replace("_output_sd.", "_output_mean.")
    if name == source.name:
        raise ValueError(f"Not a DART *_sd summary: {source}")
    mean_source = source.with_name(name)
    if not mean_source.is_file():
        raise FileNotFoundError(f"Ensemble-mean summary needed for {source.name} is missing: {mean_source}")
    return mean_source


def remapped_summary_path(
    output_root: Path, experiment: str, timestamp: str, product: str
) -> Path:
    """Return the standardized pressure-grid target path."""
    return Path(output_root) / experiment / f"{experiment}.{product}.{timestamp}.nc"


def build_remap_plan(
    input_root: Path,
    output_root: Path,
    experiments: Mapping[str, Mapping[str, str]],
    timestamps: Iterable[str],
    *,
    products: Sequence[str] = tuple(SUMMARY_PRODUCTS),
    include_existing: bool = False,
) -> list[RemapRequest]:
    """Build a deterministic plan and fail early when a native input is absent."""
    unknown = sorted(set(products) - set(SUMMARY_PRODUCTS))
    if unknown:
        raise ValueError(f"Unsupported DART summary product(s): {', '.join(unknown)}")

    requests: list[RemapRequest] = []
    missing_sources: list[Path] = []
    for experiment_name, experiment in experiments.items():
        for timestamp in timestamps:
            for product in products:
                source = native_summary_path(
                    Path(input_root), experiment, timestamp, product
                )
                target = remapped_summary_path(
                    Path(output_root), experiment_name, timestamp, product
                )
                if not source.is_file():
                    missing_sources.append(source)
                elif include_existing or not _is_netcdf_file(target):
                    requests.append(
                        RemapRequest(
                            experiment_name, timestamp, product, source, target
                        )
                    )

    if missing_sources:
        preview = "\n".join(f"  - {path}" for path in missing_sources[:30])
        suffix = (
            "" if len(missing_sources) <= 30
            else f"\n  ... and {len(missing_sources) - 30} more"
        )
        raise FileNotFoundError(
            f"Missing {len(missing_sources)} native DART summary file(s):\n"
            f"{preview}{suffix}"
        )
    return requests


def write_pressure_grid(path: Path) -> Path:
    """Write the canonical NCO vertical-remapping grid."""
    path = Path(path)
    dataset = xr.Dataset(
        coords={
            "plev": xr.DataArray(
                PRESSURE_LEVELS_PA,
                dims=("plev",),
                attrs={
                    "long_name": "pressure",
                    "standard_name": "air_pressure",
                    "units": "Pa",
                    "positive": "down",
                    "axis": "Z",
                },
            )
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with _NETCDF_LOCK:
        dataset.to_netcdf(temporary)
        dataset.close()
        os.replace(temporary, path)
    return path


def validate_pressure_grid(path: Path) -> None:
    """Check that a vertical-grid file contains the canonical pressure levels."""
    with _NETCDF_LOCK:
        with xr.open_dataset(path) as dataset:
            if "plev" not in dataset:
                raise ValueError(f"Pressure coordinate 'plev' is missing from {path}")
            if not np.array_equal(dataset["plev"].values, PRESSURE_LEVELS_PA):
                raise ValueError(f"Unexpected pressure levels in {path}")


def _resolve_executable(name: str, nco_bin: Path | None) -> str:
    candidate = str(Path(nco_bin) / name) if nco_bin else shutil.which(name)
    if not candidate or not Path(candidate).is_file():
        location = f" in {nco_bin}" if nco_bin else " on PATH"
        raise RuntimeError(f"Required NCO executable '{name}' was not found{location}.")
    return candidate


def _run(command: Sequence[str], *, environment: Mapping[str, str] | None = None) -> None:
    subprocess.run(list(map(str, command)), check=True, env=environment)


def validate_remapped_summary(path: Path) -> None:
    """Validate the grid and core variables of a completed remap."""
    path = Path(path)
    with _NETCDF_LOCK:
        with xr.open_dataset(path) as dataset:
            expected_sizes = {"plev": 19, "lat": 180, "lon": 360}
            wrong_sizes = {
                name: (dataset.sizes.get(name), size)
                for name, size in expected_sizes.items()
                if dataset.sizes.get(name) != size
            }
            if wrong_sizes:
                raise ValueError(f"Unexpected dimensions in {path}: {wrong_sizes}")
            if not np.array_equal(dataset["plev"].values, PRESSURE_LEVELS_PA):
                raise ValueError(f"Unexpected pressure levels in {path}")
            missing = sorted(REQUIRED_STATE_VARIABLES - set(dataset.data_vars))
            if missing:
                raise ValueError(f"Missing state variables in {path}: {', '.join(missing)}")
            # Mid-troposphere is above ground almost everywhere; an empty level
            # means the vertical interpolation used a wrong surface pressure.
            valid = float(dataset["T"].sel(plev=50000.0).notnull().mean())
            if valid < 0.95:
                raise ValueError(f"T at 500 hPa is only {valid:.0%} valid in {path}")
            # Mean products must carry a real surface pressure; a spread file
            # mislabelled as a mean has PS of only tens to hundreds of Pa.
            if "PS_ENSMEAN" not in dataset and float(dataset["PS"].min()) < 30000.0:
                raise ValueError(f"PS in {path} is not a surface pressure (min {float(dataset['PS'].min()):.0f} Pa)")


NETCDF_SIGNATURES = (b"CDF\x01", b"CDF\x02", b"CDF\x05", b"\x89HDF\r\n\x1a\n")


def _is_netcdf_file(path: Path) -> bool:
    """True when ``path`` exists and starts with a NetCDF-3 or NetCDF-4 signature.

    Truncated or zero-filled targets (e.g. left by an interrupted copy) are then
    treated as missing and regenerated instead of being trusted as complete.
    """
    try:
        with open(path, "rb") as handle:
            header = handle.read(8)
    except OSError:
        return False
    return any(header.startswith(signature) for signature in NETCDF_SIGNATURES)


def remap_summary(
    request: RemapRequest,
    *,
    map_file: Path,
    nco_bin: Path | None = None,
    force: bool = False,
    pressure_grid_file: Path | None = None,
) -> Path:
    """Create and atomically publish one horizontally/vertically remapped file."""
    map_file = Path(map_file)
    if not map_file.is_file():
        raise FileNotFoundError(f"Horizontal map file not found: {map_file}")
    if _is_netcdf_file(request.target) and not force:
        validate_remapped_summary(request.target)
        return request.target

    ncrename = _resolve_executable("ncrename", nco_bin)
    ncremap = _resolve_executable("ncremap", nco_bin)
    environment = os.environ.copy()
    if nco_bin:
        # ncremap is a driver script that locates ncks/ncap2 and other NCO
        # binaries through PATH. Pin the whole suite to one installation.
        environment["PATH"] = os.pathsep.join(
            [str(Path(nco_bin)), environment.get("PATH", "")]
        )
    request.target.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix=".dart-remap-", dir=request.target.parent
    ) as temporary_directory:
        temporary_directory = Path(temporary_directory)
        native_copy = temporary_directory / "native.nc"
        horizontal = temporary_directory / "horizontal.nc"
        pressure_grid = (
            Path(pressure_grid_file) if pressure_grid_file
            else write_pressure_grid(temporary_directory / "plev19.nc")
        )
        completed = temporary_directory / "completed.nc"

        shutil.copy2(request.source, native_copy)
        spread_product = request.product.endswith("_sd")
        if spread_product:
            # In *_sd files PS is the ensemble spread of surface pressure (tens
            # to hundreds of Pa). Vertical interpolation needs the actual
            # surface pressure, so borrow the ensemble-mean PS of the matching
            # *_mean product and keep the PS spread as PS_SD until the end.
            mean_source = ensemble_mean_source(request.source)
            ncks = _resolve_executable("ncks", nco_bin)
            _run([ncrename, "-v", "PS,PS_SD", native_copy], environment=environment)
            _run([ncks, "-A", "-v", "PS", mean_source, native_copy], environment=environment)
        _run(
            [ncrename, "-d", ".ncol_d,ncol", native_copy],
            environment=environment,
        )
        _run(
            [
                ncremap, "-P", "eam", "-m", map_file,
                "-i", native_copy, "-o", horizontal,
            ],
            environment=environment,
        )
        _run(
            [
                ncremap,
                "--vrt_ntp=log",
                "--vrt_xtr=mss_val",
                f"--vrt_out={pressure_grid}",
                horizontal,
                completed,
            ],
            environment=environment,
        )
        if spread_product:
            # Restore the file's meaning: PS is the PS spread; keep the mean
            # PS used for the interpolation as PS_ENSMEAN.
            _run([ncrename, "-v", "PS,PS_ENSMEAN", "-v", "PS_SD,PS", completed],
                 environment=environment)
        validate_remapped_summary(completed)
        os.replace(completed, request.target)
    return request.target


def ensure_remapped_summaries(
    input_root: Path,
    output_root: Path,
    experiments: Mapping[str, Mapping[str, str]],
    timestamps: Iterable[str],
    *,
    map_file: Path,
    products: Sequence[str] = tuple(SUMMARY_PRODUCTS),
    nco_bin: Path | None = None,
    generate_missing: bool = False,
    workers: int = 1,
) -> list[Path]:
    """Validate availability and optionally generate every missing remap."""
    plan = build_remap_plan(
        input_root, output_root, experiments, timestamps, products=products
    )
    if plan and not generate_missing:
        preview = "\n".join(f"  - {request.target}" for request in plan[:30])
        suffix = "" if len(plan) <= 30 else f"\n  ... and {len(plan) - 30} more"
        raise FileNotFoundError(
            f"Missing {len(plan)} remapped DART summary file(s). "
            "Set generate_missing=True or run the shared remapping command:\n"
            f"{preview}{suffix}"
        )

    if workers < 1:
        raise ValueError("workers must be at least 1")

    pressure_grid = Path(output_root) / "_provenance" / "vrt_remap_plev19.nc"
    if pressure_grid.is_file():
        validate_pressure_grid(pressure_grid)
    else:
        write_pressure_grid(pressure_grid)

    def process(item: tuple[int, RemapRequest]) -> Path:
        index, request = item
        print(
            f"[REMAP {index}/{len(plan)}] {request.experiment} "
            f"{request.timestamp} {request.product}",
            flush=True,
        )
        return remap_summary(
            request,
            map_file=map_file,
            nco_bin=nco_bin,
            pressure_grid_file=pressure_grid,
        )

    indexed_plan = list(enumerate(plan, start=1))
    if workers == 1:
        return [process(item) for item in indexed_plan]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(process, indexed_plan))
