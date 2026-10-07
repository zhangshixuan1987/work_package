#!/usr/bin/env python3
"""Move the legacy ``diag_dart`` tree into workflow-oriented output folders.

The command is a dry run unless ``--execute`` is supplied.  It refuses unknown
files, duplicate destinations, existing destination files, cross-filesystem
moves, and source/destination nesting.  An execution manifest is written before
the first atomic rename so an interrupted migration remains auditable.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_SOURCE = Path("/compyfs/zhan391/v3_dart_cda_scratch/diag_dart")
DEFAULT_DESTINATION = Path("/compyfs/www/zhan391/e3sm_dart/diag_dart_2026/data")

FORECAST_EXPERIMENT_DIRS = {"CAPTEN10", "CTRLEN10", "DARTEN20", "DARTEN40"}
PROVENANCE_FILES = {
    "diagnostic_run_config.txt",
    "run_e3sm_regrid_ic_s0.bash",
    "run_e3sm_regrid_ic_s1.bash",
    "sh0.log",
    "sh1.log",
}


def destination_relative(source_relative: Path) -> Path:
    """Return the new path below ``diag_dart_2026/data`` for one legacy file."""
    parts = source_relative.parts
    top = parts[0]
    name = source_relative.name

    if top in FORECAST_EXPERIMENT_DIRS:
        if len(parts) != 2:
            raise ValueError(f"Unexpected nested experiment product: {source_relative}")
        if name.startswith("seb_") and name.endswith(".nc"):
            return Path(
                "fcst_atm_model_evaluation",
                "surface_energy_budget",
                "legacy",
                top,
                name,
            )
        if (
            "_ensemble_mean_bias." in name
            or "_ensemble_bias_summary." in name
            or name.startswith("run_config_")
        ):
            return Path("fcst_atm_bias", "bias_metrics", "legacy", top, name)
        raise ValueError(f"Unknown experiment-level product: {source_relative}")

    if len(parts) == 1:
        if name in PROVENANCE_FILES:
            return Path(
                "analysis_atm_init",
                "remapped_ensemble_diagnostics",
                "_provenance",
                name,
            )
        if name.endswith("_ensemble_metrics.nc") or name.endswith("_ensemble_summary.csv"):
            return Path("fcst_atm_model_evaluation", "legacy_metrics", name)
        if "_bias_spread_" in name and name.endswith(".nc"):
            return Path("analysis_lnd_init", "surface_variables", name)
        raise ValueError(f"Unknown top-level product: {source_relative}")

    remainder = Path(*parts[1:])
    if top == "analysis_bias":
        return Path("analysis_atm_model_evaluation", "legacy_metrics", remainder)
    if top == "analysis_ic":
        return Path("analysis_atm_init", "remapped_ensemble_diagnostics", remainder)
    if top == "som_error":
        return Path("analysis_lnd_init", "initial_condition_error", remainder)
    if top == "tci":
        return Path("analysis_lac", "dirmeyer_tci", remainder)
    if top == "tcc_rmse_base":
        return Path("fcst_s2s_tcc", "biweekly_metrics", remainder)
    if top == "tcc_rmse_sigvar":
        return Path("fcst_s2s_tcc", "weekly_metrics", remainder)
    if top == "acc_rmse":
        if name.startswith("s2s_skill_"):
            return Path("fcst_s2s_pcc", name)
        if name.startswith("s2s_accmap_"):
            return Path("fcst_s2s_pcc", "acc_maps", name)
        raise ValueError(f"Unknown PCC product: {source_relative}")

    raise ValueError(f"Unknown legacy product: {source_relative}")


def build_plan(source: Path, destination: Path) -> list[dict[str, object]]:
    files = sorted(path for path in source.rglob("*") if path.is_file() or path.is_symlink())
    plan: list[dict[str, object]] = []
    destinations: dict[Path, Path] = {}
    errors: list[str] = []

    for source_path in files:
        relative = source_path.relative_to(source)
        try:
            target = destination / destination_relative(relative)
        except ValueError as error:
            errors.append(str(error))
            continue
        previous = destinations.get(target)
        if previous is not None:
            errors.append(f"Duplicate destination {target}: {previous} and {source_path}")
        destinations[target] = source_path
        if target.exists() or target.is_symlink():
            errors.append(f"Destination already exists: {target}")
        stat = source_path.lstat()
        plan.append(
            {
                "source": str(source_path),
                "destination": str(target),
                "size": stat.st_size,
            }
        )

    if errors:
        preview = "\n".join(f"  - {message}" for message in errors[:50])
        suffix = f"\n  ... and {len(errors) - 50} more" if len(errors) > 50 else ""
        raise RuntimeError(f"Migration plan has {len(errors)} error(s):\n{preview}{suffix}")
    return plan


def remove_empty_directories(root: Path) -> None:
    for directory in sorted(
        (path for path in root.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    ):
        directory.rmdir()
    root.rmdir()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--execute", action="store_true", help="perform atomic moves")
    parser.add_argument("--manifest", type=Path, help="execution manifest path")
    args = parser.parse_args()

    source = args.source.expanduser().resolve()
    destination = args.destination.expanduser().resolve()
    if not source.is_dir():
        raise FileNotFoundError(f"Legacy source directory does not exist: {source}")
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError("Source and destination must not be equal or nested")
    if source.stat().st_dev != destination.stat().st_dev:
        raise OSError("Source and destination are on different filesystems; atomic moves are required")

    plan = build_plan(source, destination)
    total_bytes = sum(int(item["size"]) for item in plan)
    groups = Counter(
        str(Path(str(item["destination"])).relative_to(destination).parts[0]) for item in plan
    )
    print(f"Mode: {'EXECUTE' if args.execute else 'DRY RUN'}")
    print(f"Source: {source}")
    print(f"Destination: {destination}")
    print(f"Files: {len(plan)}")
    print(f"Logical bytes: {total_bytes}")
    for group, count in sorted(groups.items()):
        print(f"  {group}: {count}")

    if not args.execute:
        print("No files moved. Re-run with --execute after reviewing this summary.")
        return

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest = args.manifest or destination.parent / "manifests" / f"diag_dart_{timestamp}.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "created_utc": timestamp,
        "source": str(source),
        "destination": str(destination),
        "file_count": len(plan),
        "logical_bytes": total_bytes,
        "status": "planned",
        "files": plan,
    }
    manifest.write_text(json.dumps(payload, indent=2) + "\n")

    moved = 0
    try:
        for item in plan:
            source_path = Path(str(item["source"]))
            target = Path(str(item["destination"]))
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source_path, target)
            moved += 1
    except Exception:
        payload["status"] = "interrupted"
        payload["moved_count"] = moved
        manifest.write_text(json.dumps(payload, indent=2) + "\n")
        raise

    remove_empty_directories(source)
    payload["status"] = "complete"
    payload["moved_count"] = moved
    payload["completed_utc"] = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Moved {moved} files and removed the empty source tree.")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
