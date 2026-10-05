"""Experiment table helpers.

Notebooks define their experiments as
    EXPERIMENTS = {tag: {"spin-up": {"name": case_dir, "period": "YYYY-YYYY"},
                         "piControl": {"name": case_dir, "period": "YYYY-YYYY"}}}
and turn them into the `run_dict` structure the analysis classes expect.
"""
from typing import Dict, Iterable, Optional, Tuple


def make_run_entry(name_spinup, period_spinup, name_picontrol, period_picontrol, path):
    return {
        'spin-up': {'name': name_spinup, 'period': period_spinup, 'path': f"{path}"},
        'piControl': {'name': name_picontrol, 'period': period_picontrol, 'path': f"{path}"},
    }


def make_run_dict(experiments, path):
    """{tag: {"spin-up"|"piControl": {"name", "period", "path"}}} from a notebook's experiment table.

    path: directory containing the case directories (an entry may override it with "path").
    """
    return {
        tag: make_run_entry(
            name_spinup=segs["spin-up"]["name"],
            period_spinup=segs["spin-up"]["period"],
            name_picontrol=segs["piControl"]["name"],
            period_picontrol=segs["piControl"]["period"],
            path=segs.get("path", path),
        )
        for tag, segs in experiments.items()
    }


def build_exp_subdirs_from_run_dict(
        run_dict: Dict[str, object],
        exp_tags: Iterable[str],
        exp_type: str = "piControl",            # or "spin-up"
        ts_window: Optional[Tuple[int, int]] = None,
        climo_len: int = 50,
        mpas_analysis: bool = True,
        ts_base: str = "post/analysis/mpas_analysis"
    ) -> Dict[str, dict]:
    """
    Construct {exp_name: '<ts_base>/ts_YYYY-YYYY_climo_YYYY-YYYY/timeseries'}.

    Parameters
    ----------
    run_dict : dict
      Output of make_run_dict(...).
    exp_tags : iterable[str]
      Keys into run_dict, e.g. ["Full-CPL","AltBG1"].
    exp_type : str
      "piControl" or "spin-up" (which block of run_dict to use).
    ts_window : (start,end) or None
      If provided, force the TS token to this window (e.g., (1, 350)).
      If None, use the period from run_dict[exp_tag][exp_type]["period"].
    climo_len : int
      Length of the CLIMO window, ending at TS end (e.g., last 50 years).
    ts_base : str
      Base directory prefix used by mpas_analysis outputs.

    Returns
    -------
    dict : {exp_name: subdir}
    """
    mapping = {}
    for tag in exp_tags:
        meta = run_dict[tag][exp_type]
        exp_name = meta["name"]

        if ts_window is None:
            p = meta["period"]  # "YYYY-YYYY"
            y1, y2 = (int(t) for t in p.split("-"))
        else:
            y1, y2 = ts_window

        if mpas_analysis:
            ts_token = f"{y1:04d}-{y2:04d}"
            clim_start = max(y1, y2 - (climo_len - 1))
            climo_token = f"{clim_start:04d}-{y2:04d}"
            subdir = f"{ts_base}/ts_{ts_token}_climo_{climo_token}/timeseries"
        else:
            subdir = f"{ts_base}"
        mapping[exp_name] = subdir
    return mapping
