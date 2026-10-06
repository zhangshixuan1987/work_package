#!/usr/bin/env python3
"""Command-line entry point for column-to-gridcell restart aggregation."""

from pathlib import Path
import sys


DIAGNOSTIC_DIR = Path(__file__).resolve().parents[2]
if str(DIAGNOSTIC_DIR) not in sys.path:
    sys.path.insert(0, str(DIAGNOSTIC_DIR))

from util.initial_lnd_restart_aggregation import main


if __name__ == "__main__":
    main()
