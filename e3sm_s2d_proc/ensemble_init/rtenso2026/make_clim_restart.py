#!/usr/bin/env python
"""Select a BruteForce year and copy its land, river, and sea-ice restarts.

Branches: ``clim`` (year closest to the climatology) and ``obs`` (year closest
to the ERA5 anomaly on the target date). Outputs go to
<out-root>/YYYYMMDD/<branch>/.

Examples:
    python make_clim_restart.py --date 2015-02-01 --branch obs --dry-run
    python make_clim_restart.py --date 2026-05-01 --branch clim
    python make_clim_restart.py --date 2026-05-01 --branch obs --mode per_component
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "util"))

import pandas as pd  # noqa: E402

from clim_restart import ALL_COMPONENTS, make_plan, run_selection  # noqa: E402

SRC_ROOT = "/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/BruteForce"
OUT_ROOT = "/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--date", required=True, help="target initialization date, YYYY-MM-DD")
    p.add_argument("--branch", default="clim", choices=["clim", "obs"])
    p.add_argument("--components", nargs="+", default=list(ALL_COMPONENTS), choices=ALL_COMPONENTS)
    p.add_argument("--mode", default="joint", choices=["joint", "per_component"])
    p.add_argument("--box", type=float, nargs=2, default=[30.0, 60.0], metavar=("DLAT", "DLON"),
                   help="compare means over DLAT x DLON degree boxes (default 30 60; 0 0 = model points)")
    p.add_argument("--component-weights", nargs="*", default=[], metavar="COMP=W",
                   help="weights a_c in the joint distance, e.g. mpassi=2 (default 1 each)")
    p.add_argument("--src-root", default=SRC_ROOT, help="initial-condition archive (read-only)")
    p.add_argument("--out-root", default=OUT_ROOT, help="files go to <out-root>/YYYYMMDD/<branch>/")
    p.add_argument("--src-prefix", default="20251024_s2d_spinup")
    p.add_argument("--out-prefix", default="20251024_s2d_spinup")
    p.add_argument("--first-year", type=int)
    p.add_argument("--last-year", type=int)
    p.add_argument("--exclude-years", type=int, nargs="*", default=[])
    p.add_argument("--keep-target-year", action="store_true",
                   help="allow the target year itself as a candidate (excluded by default)")
    p.add_argument("--cache-dir", help="key-field cache (default <out-root>/cache; '' to disable)")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="compute distances and selection only")
    args = p.parse_args()

    weights = {k: float(v) for k, v in (item.split("=") for item in args.component_weights)}
    try:
        plan = make_plan(
            args.date, args.src_root, args.out_root, branch=args.branch,
            components=args.components, src_prefix=args.src_prefix, out_prefix=args.out_prefix,
            first_year=args.first_year, last_year=args.last_year,
            exclude_years=args.exclude_years, exclude_target_year=not args.keep_target_year,
        )
    except FileNotFoundError as err:      # e.g. ERA5 not yet available for an obs date
        print(f"skip {args.date} [{args.branch}]: {err}")
        return
    distances, _ = run_selection(
        plan, box=tuple(args.box) if all(args.box) else None, mode=args.mode, component_weights=weights,
        cache_dir=(Path(args.out_root) / "cache" if args.cache_dir is None
                   else Path(args.cache_dir) if args.cache_dir else None),
        write=not args.dry_run, overwrite=args.overwrite,
        log=lambda m: print(m, flush=True),
    )
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(distances.sort_values("joint").round(3).head(10).to_string())


if __name__ == "__main__":
    main()
