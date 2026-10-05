"""Fast checks that need no model data:  python -m pytest tests/   (from coupled_spinup_eval/)"""
import ast
import glob
import importlib
import json
import os
import re

import pytest

from util.experiments import build_exp_subdirs_from_run_dict, make_run_dict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = sorted(os.path.basename(p)[:-3] for p in glob.glob(f"{ROOT}/util/*.py") if not p.endswith("__init__.py"))
NOTEBOOKS = sorted(glob.glob(f"{ROOT}/jupyter/[0-9]*_*.ipynb"))


@pytest.mark.parametrize("mod", MODULES)
def test_util_module_imports(mod):
    importlib.import_module(f"util.{mod}")


def _code_cells(path):
    nb = json.load(open(path))
    return ["".join(c["source"]).replace("%matplotlib inline", "") for c in nb["cells"] if c["cell_type"] == "code"]


@pytest.mark.parametrize("nb", NOTEBOOKS, ids=os.path.basename)
def test_notebook_cells_compile_and_util_imports_exist(nb):
    cells = _code_cells(nb)
    for i, c in enumerate(cells):
        compile(c, f"{os.path.basename(nb)}[{i}]", "exec")
    for node in ast.walk(ast.parse(cells[0])):
        if isinstance(node, ast.ImportFrom) and (node.module or "").startswith("util."):
            mod = importlib.import_module(node.module)
            for a in node.names:
                assert hasattr(mod, a.name), f"{node.module}.{a.name} missing"


def test_make_run_dict_and_subdirs():
    exps = {"A": {"spin-up": {"name": "caseA", "period": "0001-0220"},
                  "piControl": {"name": "caseA", "period": "0221-0350"}}}
    rd = make_run_dict(exps, "/root")
    assert rd["A"]["spin-up"] == {"name": "caseA", "period": "0001-0220", "path": "/root"}
    sub = build_exp_subdirs_from_run_dict(rd, ["A"], exp_type="spin-up", ts_window=(1, 350))
    assert sub == {"caseA": "post/analysis/mpas_analysis/ts_0001-0350_climo_0301-0350/timeseries"}


def test_subclasses_share_their_base():
    from util import alt_global_timeseries as ag, ctrl_plotting as cp, ctrl_processing as cpr, land_bgc as lb
    from util import regional_timeseries as rt
    assert issubclass(cpr.SSTBuilder, cpr.MPASDiagnosticsBuilder)
    assert issubclass(cpr.OceanBudgetBuilder, cpr.MPASDiagnosticsBuilder)
    assert issubclass(cp.OHCDriftPlotter, cp.MPASTimeDiagnosticsPlotter)
    assert issubclass(rt.NinoIndexPlotter, rt.OCNTimeSeriesPlotter)
    assert issubclass(ag.OHUDiagnosticsPlotter, ag.OCNDiagnosticsPlotter)
    assert issubclass(lb.LandCompareStability, lb.SpinupStabilityBGC)


def test_no_hardcoded_paper_paths_in_notebook_code():
    """Notebook analysis code takes roots from the setup cell (only the setup cell may hold paths)."""
    for nb in NOTEBOOKS:
        cells = _code_cells(nb)
        setup = [c for c in cells if c.lstrip().startswith("# ---- paths ----")]
        assert len(setup) == 1, nb
        for c in cells:
            if c is setup[0]:
                continue
            assert not re.search(r'"/lcrc/', c), f"hard-coded /lcrc path outside setup in {os.path.basename(nb)}"


def test_piControl_sigma_and_z_ref():
    """σ of a metric across piControl windows ignores a linear drift; RMSE_z_ref = (X - mean)/σ."""
    import numpy as np
    import pandas as pd
    from util.model_error import SpinupMetricAnalyzer

    periods = [f"{y:04d}-{y + 49:04d}" for y in range(1, 500, 50)]
    noise = np.array([1, -1, 1, -1, 1, -1, 1, -1, 1, -1], float) * 0.1
    df = pd.DataFrame({
        "drift": 10.0 + 0.01 * np.arange(10) * 50 + noise,   # trend + noise
        "flat": np.full(10, 5.0),                           # σ = 0 -> floor
        "short": [1.0, 2.0, 3.0] + [np.nan] * 7,            # too few windows
    }, index=periods)
    sig = SpinupMetricAnalyzer.piControl_sigma_vector(df, floor_frac=0.001)
    resid = df["drift"] - np.polyval(np.polyfit(np.arange(10), df["drift"], 1), np.arange(10))
    assert np.isclose(sig["drift"], resid.std(ddof=2))
    assert sig["drift"] < df["drift"].std()                 # drift removed
    assert np.isclose(sig["flat"], 0.005)
    assert np.isnan(sig["short"])

    an = SpinupMetricAnalyzer(base_dir=".", spinup_name="x")
    an.variables = ["drift", "flat"]
    an.set_external_reference_vectors({"RMSE": {"drift": 12.0, "flat": 5.0}, "RMSE_sigma": sig})
    z = an.z_to_ref(np.array([[12.0 + 2 * sig["drift"], 5.01]]), "RMSE")
    assert np.allclose(z, [[2.0, 2.0]])
