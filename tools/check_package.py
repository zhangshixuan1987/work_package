#!/usr/bin/env python3
"""Consistency and path checks for a work_package analysis package.

Usage::

    python tools/check_package.py polar_analysis [pcmdi_analysis ...]

Errors (exit status 1):
  * a notebook is not valid JSON or a code cell does not compile;
  * the v3 LE root from config/paths.json is written literally anywhere
    outside config/ (it must come from scripts/paths.py or shell_script/paths.sh);
  * a notebook uses V3LE_* / fig_dir / require before importing them from paths;
  * notebooks that use the standard setup cell do not all share the same one;
  * a shell script fails ``bash -n`` or does not source paths.sh.

Warnings (reported, exit status unaffected):
  * an absolute /lcrc path or a V3LE_* sub-path used in code does not exist.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import warnings
from collections import defaultdict
from pathlib import Path

WORK_PACKAGE = Path(__file__).resolve().parents[1]
PATH_NAMES = ("V3LE_ROOT", "V3LE_DATA_DIR", "V3LE_DIAG_DIR", "V3LE_FIG_ROOT", "fig_dir", "require")
SETUP_MARKER = "# Locate this package and load the shared v3 LE locations"
LCRC_RE = re.compile(r'["\'](/lcrc/[^"\'{}*?\[\]\s]+)')
V3LE_SUB_RE = re.compile(
    r'(?:\{(V3LE_\w+)\}/([^"\'{}*?\s]+)|(V3LE_\w+)\s*/\s*"([^"{}*?]+)")')

try:  # Strip IPython magics/shell escapes the way Jupyter does.
    from IPython.core.inputtransformer2 import TransformerManager
    _TM = TransformerManager()

    def to_python(src: str) -> str:
        return _TM.transform_cell(src)
except ImportError:  # pragma: no cover - fallback without IPython
    def to_python(src: str) -> str:
        return "\n".join("pass" if l.lstrip().startswith(("%", "!")) else l
                         for l in src.splitlines())


def code_cells(nb: dict) -> list[str]:
    return ["".join(c["source"]) for c in nb.get("cells", []) if c.get("cell_type") == "code"]


def load_paths(pkg: Path) -> dict:
    """Import the package's scripts/paths.py in isolation."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"_paths_{pkg.name}", pkg / "scripts" / "paths.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return {n: getattr(mod, n) for n in ("V3LE_ROOT", "V3LE_DATA_DIR", "V3LE_DIAG_DIR", "V3LE_FIG_ROOT")}


def check(pkg: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warns: list[str] = []
    rel = lambda p: p.relative_to(WORK_PACKAGE)  # noqa: E731

    config = pkg / "config" / "paths.json"
    if not config.is_file() or not (pkg / "scripts" / "paths.py").is_file():
        return [f"{pkg.name}: missing config/paths.json or scripts/paths.py"], []
    literal_root = json.loads(config.read_text())["v3le_root"]
    resolved = load_paths(pkg)
    if not resolved["V3LE_ROOT"].is_dir():
        errors.append(f"V3LE_ROOT does not exist: {resolved['V3LE_ROOT']}")

    missing_paths: dict[str, set[str]] = defaultdict(set)
    setup_variants: dict[str, list[str]] = defaultdict(list)

    for nb_path in sorted((pkg / "jupyter").glob("*.ipynb")):
        try:
            cells = code_cells(json.loads(nb_path.read_text()))
        except json.JSONDecodeError as exc:
            errors.append(f"{rel(nb_path)}: invalid JSON ({exc})")
            continue
        imported = False
        for idx, src in enumerate(cells):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", SyntaxWarning)
                    compile(to_python(src), f"{nb_path.name}[{idx}]", "exec")
            except SyntaxError as exc:
                errors.append(f"{rel(nb_path)} cell {idx}: {exc.msg} (line {exc.lineno})")
            if literal_root in src:
                errors.append(f"{rel(nb_path)} cell {idx}: hard-coded v3LE root; use paths.py")
            if re.search(r"^\s*from paths import", src, re.M):
                imported = True
            body = re.sub(r"^\s*from paths import.*$", "", src, flags=re.M)
            # Constants count when referenced; helpers only when called, since
            # names like fig_dir are also common parameter/attribute names.
            used = [n for n in PATH_NAMES
                    if re.search(rf"(?<![\w.]){n}" + (r"\s*\(" if n.islower() else r"\b"), body)
                    and not re.search(rf"^\s*{n}\s*=", body, re.M)
                    and not re.search(rf"^\s*def {n}\b", body, re.M)]
            if used and not imported:
                errors.append(f"{rel(nb_path)} cell {idx}: uses {used} before 'from paths import'")
            if SETUP_MARKER in src:
                key = re.sub(r'\n*FIG_DIR_ROOT = fig_dir\("[^"]*"\)\s*$', "", src)
                setup_variants[key].append(nb_path.name)
            for m in LCRC_RE.finditer(src):
                if not Path(m.group(1)).exists():
                    missing_paths[m.group(1)].add(nb_path.name)
            for m in V3LE_SUB_RE.finditer(src):
                var, sub = (m.group(1), m.group(2)) if m.group(1) else (m.group(3), m.group(4))
                if var in resolved and not (resolved[var] / sub).exists():
                    missing_paths[f"{var}/{sub}"].add(nb_path.name)

    if len(setup_variants) > 1:
        sizes = sorted(((len(v), v[:3]) for v in setup_variants.values()), reverse=True)
        errors.append(f"{len(setup_variants)} different setup cells: {sizes}")

    for py in sorted((pkg / "scripts").glob("*.py")):
        if literal_root in py.read_text():
            errors.append(f"{rel(py)}: hard-coded v3LE root; use paths.py")

    shell_dir = pkg / "shell_script"
    for sh in sorted(shell_dir.glob("*.bash")) if shell_dir.is_dir() else []:
        text = sh.read_text()
        if subprocess.run(["bash", "-n", str(sh)], capture_output=True).returncode:
            errors.append(f"{rel(sh)}: bash -n failed")
        if "paths.sh" not in text:
            errors.append(f"{rel(sh)}: does not source paths.sh")
        if literal_root in text:
            errors.append(f"{rel(sh)}: hard-coded v3LE root; use paths.sh")
        for m in LCRC_RE.finditer(text):
            if not Path(m.group(1)).exists():
                missing_paths[m.group(1)].add(sh.name)

    for path, users in sorted(missing_paths.items()):
        names = sorted(users)
        more = f" (+{len(names) - 4} more)" if len(names) > 4 else ""
        warns.append(f"missing: {path}  <- {', '.join(names[:4])}{more}")
    return errors, warns


def main(argv: list[str]) -> int:
    names = argv or ["polar_analysis", "pcmdi_analysis"]
    status = 0
    for name in names:
        pkg = (WORK_PACKAGE / name).resolve()
        errors, warns = check(pkg)
        print(f"== {name}: {len(errors)} error(s), {len(warns)} warning(s)")
        for e in errors:
            print(f"  ERROR {e}")
        for w in warns:
            print(f"  WARN  {w}")
        status |= bool(errors)
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
