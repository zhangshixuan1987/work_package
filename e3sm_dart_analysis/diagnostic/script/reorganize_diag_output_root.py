#!/usr/bin/env python3
"""Tidy the diagnostic output root into data/, figure/, paper/, archive/, manifests/.

Dry run unless ``--execute`` is given. Nothing is deleted: each move is an
atomic same-filesystem rename, and a JSON manifest is written before the first
move so the change can be audited or reversed.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_ROOT = Path("/compyfs/www/zhan391/e3sm_dart/diag_dart_2026")

MOVES = (
    ("migration_manifests", "manifests"),
    ("figure/ncl_paper_figures", "paper/ncl_paper_figures"),
    ("figure.backup", "archive/figure_backup"),
    ("scp_test", "archive/scp_test"),
    (
        "figure/analysis_atm_da/legacy_analysis_da",
        "archive/figure/analysis_atm_da/legacy_analysis_da",
    ),
    (
        "figure/analysis_atm_da/legacy_comparisons",
        "archive/figure/analysis_atm_da/legacy_comparisons",
    ),
)

README = """# E3SM-DART diagnostic outputs

| Path | Contents | Written by |
| --- | --- | --- |
| `data/<workflow>/` | NetCDF/CSV products, one folder per workflow | `diagnostic/jupyter/<workflow>` notebooks |
| `figure/<workflow>/` | Figures, one folder per workflow | `diagnostic/jupyter/<workflow>` notebooks |
| `data/ncl/`, `figure/ncl/` | NCL products (`active/`, `legacy/`) | `diagnostic/script/run_process_ncl.bash` |
| `paper/ncl_paper_figures/` | Historical NCL paper-figure sources | read-only archive |
| `archive/` | Superseded outputs kept for reference | nothing writes here |
| `manifests/` | JSON records of every reorganization/migration | migration scripts |

Workflow names are defined in `diagnostic/configs/output_paths.py`.
"""


def plan(root: Path) -> list[tuple[Path, Path]]:
    planned = []
    for src_rel, dst_rel in MOVES:
        src, dst = root / src_rel, root / dst_rel
        if not src.exists():
            print(f"[SKIP] missing: {src_rel}")
            continue
        if dst.exists():
            raise FileExistsError(f"Destination already exists: {dst}")
        planned.append((src, dst))
    return planned


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--execute", action="store_true", help="perform the moves")
    args = parser.parse_args()
    root = args.root.resolve()

    planned = plan(root)
    for src, dst in planned:
        print(f"[MOVE] {src.relative_to(root)} -> {dst.relative_to(root)}")
    readme = root / "README.md"
    print(f"[README] {'exists, unchanged' if readme.exists() else 'create'}: README.md")

    if not args.execute:
        print("Dry run only; rerun with --execute to apply.")
        return

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_dir = root / "manifests"
    if (root / "migration_manifests").exists():
        manifest_dir = root / "migration_manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest = manifest_dir / f"reorganize_{timestamp}.json"
    manifest.write_text(
        json.dumps(
            {
                "root": str(root),
                "created_utc": timestamp,
                "moves": [
                    {"source": str(s.relative_to(root)), "destination": str(d.relative_to(root))}
                    for s, d in planned
                ],
            },
            indent=2,
        )
    )
    print(f"[MANIFEST] {manifest}")

    for src, dst in planned:
        dst.parent.mkdir(parents=True, exist_ok=True)
        os.rename(src, dst)
        print(f"[DONE] {src.relative_to(root)} -> {dst.relative_to(root)}")

    if not readme.exists():
        readme.write_text(README)
        print(f"[DONE] wrote {readme}")


if __name__ == "__main__":
    main()
