#!/usr/bin/env python3
"""Consistency and path checks for a work_package analysis package.

Usage::

    python tools/check_package.py polar_analysis [pcmdi_analysis ...]

Errors (exit status 1):
  * a notebook is not valid JSON or a code cell does not compile;
  * a notebook still uses the removed ``paths`` module;
  * in packages whose notebooks use a ``parameters`` cell (papermill tag):
    - every notebook has exactly one, directly after the shared setup cell,
    - all notebooks share the same setup cell,
    - absolute /lcrc paths appear only in the parameters cell;
  * config/regions.json (if present) does not resolve against exp_info;
  * a shell script fails ``bash -n`` or lacks a ``DATA_DIR="${DATA_DIR:-...}"``
    parameter.

Warnings (reported, exit status unaffected):
  * a parameter is never used outside the parameters cell;
  * an absolute /lcrc input path (parameter default or literal) does not exist
    (output parameters named FIG_*, DIAG_*, OUT* are skipped).
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
SETUP_MARKER = "PROJECT_ROOT = next("
LCRC_RE = re.compile(r'["\'](/lcrc/[^"\'{}*?\[\]\s]+)')
PARAM_RE = re.compile(r"^([A-Za-z_]\w*)\s*=", re.M)
REMOVED = re.compile(r"from paths import|\bproject_paths\b|\bV3LE_(ROOT|DATA_DIR|DIAG_DIR|FIG_ROOT)\b")
SHELL_PARAM = re.compile(r'^DATA_DIR="\$\{DATA_DIR:-[^}]+\}"', re.M)

try:  # Strip IPython magics/shell escapes the way Jupyter does.
    from IPython.core.inputtransformer2 import TransformerManager
    _TM = TransformerManager()

    def to_python(src: str) -> str:
        return _TM.transform_cell(src)
except ImportError:  # pragma: no cover - fallback without IPython
    def to_python(src: str) -> str:
        return "\n".join("pass" if l.lstrip().startswith(("%", "!")) else l
                         for l in src.splitlines())


def is_parameters(cell: dict) -> bool:
    return "parameters" in cell.get("metadata", {}).get("tags", [])


def check_regions(pkg: Path) -> list[str]:
    if not (pkg / "config" / "regions.json").is_file():
        return []
    sys.path.insert(0, str(pkg / "scripts"))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SyntaxWarning)
            import regions  # package-local module
            for region in regions.REGIONS:
                regions.region_catalogs(region)
        return []
    except Exception as exc:
        return [f"config/regions.json: {exc}"]
    finally:
        sys.path.pop(0)
        for mod in ("regions", "exp_info"):
            sys.modules.pop(mod, None)


def check(pkg: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = check_regions(pkg)
    warns: list[str] = []
    rel = lambda p: p.relative_to(WORK_PACKAGE)  # noqa: E731
    notebooks = sorted((pkg / "jupyter").rglob("*.ipynb"))
    notebooks = [p for p in notebooks if ".ipynb_checkpoints" not in p.parts]

    loaded = {}
    for nb_path in notebooks:
        try:
            loaded[nb_path] = json.loads(nb_path.read_text())
        except json.JSONDecodeError as exc:
            errors.append(f"{rel(nb_path)}: invalid JSON ({exc})")
    use_params = any(is_parameters(c) for nb in loaded.values() for c in nb.get("cells", []))

    missing_paths: dict[str, set[str]] = defaultdict(set)
    setup_variants: dict[str, list[str]] = defaultdict(list)

    for nb_path, nb in loaded.items():
        name = str(nb_path.relative_to(pkg / "jupyter"))
        code = [c for c in nb.get("cells", []) if c.get("cell_type") == "code"]
        for idx, cell in enumerate(code):
            src = "".join(cell["source"])
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", SyntaxWarning)
                    compile(to_python(src), f"{name}[{idx}]", "exec")
            except SyntaxError as exc:
                errors.append(f"{rel(nb_path)} cell {idx}: {exc.msg} (line {exc.lineno})")
            if REMOVED.search(src):
                errors.append(f"{rel(nb_path)} cell {idx}: uses the removed paths module")
            for line in src.splitlines():
                # outputs (FIG_*/DIAG_*/OUT*) are created by the notebook itself
                if is_parameters(cell) and re.match(r"\s*(FIG|DIAG|OUT)\w*\s*=", line):
                    continue
                for m in LCRC_RE.finditer(line):
                    if not Path(m.group(1)).exists():
                        missing_paths[m.group(1)].add(name)

        if not use_params:
            continue
        params = [i for i, c in enumerate(code) if is_parameters(c)]
        setups = [i for i, c in enumerate(code) if SETUP_MARKER in "".join(c["source"])]
        if len(params) != 1:
            errors.append(f"{rel(nb_path)}: needs exactly one 'parameters' cell (found {len(params)})")
            continue
        if len(setups) != 1 or params[0] != setups[0] + 1:
            errors.append(f"{rel(nb_path)}: the 'parameters' cell must directly follow the setup cell")
        if setups:
            setup_variants["".join(code[setups[0]]["source"])].append(name)
        psrc = "".join(code[params[0]]["source"])
        rest = "\n".join("".join(c["source"]) for i, c in enumerate(code) if i != params[0])
        for i, c in enumerate(code):
            if i != params[0] and LCRC_RE.search("".join(c["source"])):
                errors.append(f"{rel(nb_path)} cell {i}: absolute /lcrc path outside the parameters cell")
        for pname in PARAM_RE.findall(psrc):
            if not re.search(rf"\b{pname}\b", rest):
                warns.append(f"unused parameter {pname} in {name}")

    if len(setup_variants) > 1:
        sizes = sorted(((len(v), v[:2]) for v in setup_variants.values()), reverse=True)
        errors.append(f"{len(setup_variants)} different setup cells: {sizes}")

    shell_dir = pkg / "shell_script"
    for sh in sorted(shell_dir.glob("*.*sh")) if shell_dir.is_dir() else []:
        text = sh.read_text()
        if subprocess.run(["bash", "-n", str(sh)], capture_output=True).returncode:
            errors.append(f"{rel(sh)}: bash -n failed")
        if not SHELL_PARAM.search(text):
            errors.append(f"{rel(sh)}: missing DATA_DIR=\"${{DATA_DIR:-<default>}}\" parameter")
        for m in LCRC_RE.finditer(text):
            if not Path(m.group(1)).exists():
                missing_paths[m.group(1)].add(sh.name)
        default = re.search(r'DATA_DIR:-([^}]+)\}', text)
        if default and not Path(default.group(1)).exists():
            missing_paths[default.group(1)].add(sh.name)

    for path, users in sorted(missing_paths.items()):
        names = sorted(users)
        more = f" (+{len(names) - 4} more)" if len(names) > 4 else ""
        warns.append(f"missing: {path}  <- {', '.join(names[:4])}{more}")
    return errors, warns


def main(argv: list[str]) -> int:
    names = argv or ["polar_analysis", "pcmdi_analysis"]
    status = 0
    for name in names:
        errors, warns = check((WORK_PACKAGE / name).resolve())
        print(f"== {name}: {len(errors)} error(s), {len(warns)} warning(s)")
        for e in errors:
            print(f"  ERROR {e}")
        for w in warns:
            print(f"  WARN  {w}")
        status |= bool(errors)
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
