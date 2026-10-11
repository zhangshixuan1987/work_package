#!/usr/bin/env python3
"""Generate standardized pressure-grid EAM-DART summary files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


DIAGNOSTIC_DIR = Path(__file__).resolve().parents[1]
if str(DIAGNOSTIC_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTIC_DIR))

from configs.atm_state_experiment_registry import state_experiments
from util.atm_dart_remapping import build_remap_plan, ensure_remapped_summaries


def timestamps(start: str, end: str) -> tuple[str, ...]:
    return tuple(
        f"{value:%Y-%m-%d}-{value.hour * 3600:05d}"
        for value in pd.date_range(start, end, freq="6h")
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root", type=Path,
        default=Path("/compyfs/zhan391/v3_dart_cda_scratch"),
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--map-file", type=Path, required=True)
    parser.add_argument("--start", default="2011-12-16 00:00:00")
    parser.add_argument("--end", default="2012-01-01 00:00:00")
    parser.add_argument(
        "--experiments", nargs="+",
        default=["DART10-S0", "DART20-S0", "DART40-S0"],
    )
    parser.add_argument("--nco-bin", type=Path)
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Number of independent remaps to run concurrently (default: 1).",
    )
    parser.add_argument(
        "--execute", action="store_true",
        help="Generate missing files. Without this flag, only print the plan.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected_experiments = state_experiments(*args.experiments)
    selected_timestamps = timestamps(args.start, args.end)
    plan = build_remap_plan(
        args.input_root, args.output_root, selected_experiments, selected_timestamps
    )
    print(
        f"{len(plan)} missing remap(s) across {len(selected_experiments)} "
        f"experiment(s) and {len(selected_timestamps)} cycle(s)."
    )
    if not args.execute:
        for request in plan:
            print(f"[PLAN] {request.source} -> {request.target}")
        return
    ensure_remapped_summaries(
        args.input_root,
        args.output_root,
        selected_experiments,
        selected_timestamps,
        map_file=args.map_file,
        nco_bin=args.nco_bin,
        generate_missing=True,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
