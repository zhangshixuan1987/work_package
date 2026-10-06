#!/usr/bin/env python3
"""Command-line entry point for conservative initial-land restart regridding."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


DIAGNOSTIC_DIR = Path(__file__).resolve().parents[2]
if str(DIAGNOSTIC_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTIC_DIR))

from util.initial_lnd_restart_regridding import CORNER_PERMUTATIONS, conservative_regrid_restart


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Conservatively regrid an ELM restart aggregated to gridcells."
    )
    parser.add_argument("--input", required=True, help="Aggregated gridcell NetCDF file.")
    parser.add_argument("--destination-scrip", required=True, help="Destination SCRIP grid.")
    parser.add_argument("--source-scrip", required=True, help="Generated source SCRIP grid.")
    parser.add_argument("--map-file", required=True, help="Conservative weight file to create/reuse.")
    parser.add_argument("--output", required=True, help="Regridded output NetCDF file.")
    parser.add_argument(
        "--corner-permutation",
        choices=tuple(CORNER_PERMUTATIONS),
        default="keep",
        help="Optional ordering adjustment for four-corner source polygons.",
    )
    parser.add_argument(
        "--variables",
        nargs="*",
        help="Optional variable subset. By default ncremap processes every compatible variable.",
    )
    parser.add_argument("--debug-level", type=int, default=1, help="ncremap debug level.")
    parser.add_argument("--force", action="store_true", help="Recreate SCRIP, weights, and output.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = conservative_regrid_restart(
        args.input,
        args.destination_scrip,
        args.source_scrip,
        args.map_file,
        args.output,
        corner_permutation=args.corner_permutation,
        variables=args.variables,
        force=args.force,
        debug_level=args.debug_level,
    )
    print(f"Regridded restart: {output}")


if __name__ == "__main__":
    main()
