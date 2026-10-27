#!/usr/bin/env python
"""Build the Reanalysis lead-time climatologies and observation anomalies used by the notebooks.

    python build_caches.py                  # default fields, reuse valid caches
    python build_caches.py TS PRECT --mode rebuild

Each cache is written once under S2D_DIAG_ROOT and reused by the notebooks
(``derivation_mode="auto"``), so running this ahead of time keeps the
notebooks fast.
"""

import argparse
import time

import rtenso_case as rc

DEFAULT_FIELDS = ["TS", "PRECT", "TREFHT", "PSL", "U850"]

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("fields", nargs="*", default=DEFAULT_FIELDS)
parser.add_argument("--mode", default="auto", choices=["auto", "rebuild"])
args = parser.parse_args()

for field in args.fields:
    for init_month in rc.REFERENCE["init_months"]:
        t0 = time.time()
        rc.leadtime_climatology(field, init_month, mode=args.mode)
        print(f"{field} init{init_month:02d} climatology ready ({time.time() - t0:.0f} s)", flush=True)
    if rc.FIELDS[field]["obs"] is not None:
        t0 = time.time()
        rc.obs_anomaly(field, mode=args.mode)
        print(f"{field} observation anomalies ready ({time.time() - t0:.0f} s)", flush=True)
