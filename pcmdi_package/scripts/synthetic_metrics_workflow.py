from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

try:
    from .logger import _setup_child_logger, _setup_root_logger
except ImportError:
    from logger import _setup_child_logger, _setup_root_logger

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path("/tmp") / f"matplotlib-{os.environ.get('USER', 'user')}")
)

_setup_root_logger()
logger = _setup_child_logger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]
PCMDI_METRICS_ROOT = Path("/lcrc/group/e3sm/diagnostics/pcmdi_data/metrics_data")
DEFAULT_DIAG_ROOT = Path(
    "/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/e3sm-pcmdi-le"
)
DEFAULT_METRIC_FILE = REPO_ROOT / "config" / "synthetic_metrics_list.json"

DEFAULT_CLIM_VARS = (
    "pr,prw,psl,rlds,rldscs,rltcre,rstcre,rsus,rsuscs,rlus,rlut,rlutcs,"
    "rsds,rsdscs,rsut,rsutcs,rtmt,sfcWind,tas,tauu,tauv,ts,ta-200,ta-850,"
    "ua-200,ua-850,va-200,va-850,zg-500"
)

CMIP_EXCLUDE_VARS = {
    "E3SM-1-0": ["ta-850"],
    "E3SM-1-1-ECA": ["ta-850"],
    "CIESM": ["pr"],
    "KIOST-ESM": ["zg-500", "ta-850"],
    "GISS-E2-2-G": ["rlutcs", "zg-500"],
}

CMIP_EXCLUDE_MODELS = [
    "E3SM-1-0",
    "E3SM-1-1",
    "E3SM-1-1-ECA",
    "E3SM-2-0",
    "E3SM-2-1",
]


@dataclass(frozen=True)
class MetricsDataset:
    """Input directories and dataset names for one comparison group."""

    clim_dir: Path
    clim_set: str
    movs_dir: Path
    movs_set: str
    enso_dir: Path
    enso_set: str


@dataclass(frozen=True)
class PlotStyle:
    """Figure layout choices for one synthetic metric panel."""

    title: str
    out_dir: Path
    font_size: float
    figure_size: tuple[float, float]
    legend_lw: float = 1.5


@dataclass
class SyntheticPlotsParameters:
    """Container for model comparison inputs and derived paths."""

    test_group: str
    test_prefix: str = "e3sm.historical.v3-LR"
    test_pattern: str = "v3.LR.historical"
    test_case_id: str = "v20260212"
    test_model_only: bool = False
    test_combined: bool = False
    show_mean_columns: bool = True
    test_dataset: Optional[MetricsDataset] = None
    test_highlight_models: Optional[List[str]] = None

    ref_group: str = "CMIP"
    ref_dataset: MetricsDataset = field(
        default_factory=lambda: dataset_for_group("CMIP", "v20260212")
    )

    clim_vars: Iterable[str] | str = DEFAULT_CLIM_VARS
    clim_regions: Iterable[str] | str = "global"
    movs_group: str = "cbf"
    error_norm: str = "reference"
    exclude_vars: Optional[Dict[str, List[str]]] = field(
        default_factory=lambda: dict(CMIP_EXCLUDE_VARS)
    )
    exclude_models: Optional[List[str]] = field(
        default_factory=lambda: list(CMIP_EXCLUDE_MODELS)
    )
    atm_modes: Iterable[str] | str = "NAM,NAO,PNA,NPO"
    atm_obs: str = "NOAA-20C"
    cpl_modes: Iterable[str] | str = "PDO,NPGO"
    cpl_obs: str = "HadISST"

    diag_path: Path = DEFAULT_DIAG_ROOT
    data_dir: Path = field(default_factory=lambda: DEFAULT_DIAG_ROOT / "climo")
    out_dir: Path = field(default_factory=lambda: DEFAULT_DIAG_ROOT / "CLIM_Metrics")
    run_type: str = "model_vs_obs"
    save_all_data: bool = True
    figure_format: str = "pdf"

    clim_viewer: bool = False
    mova_viewer: bool = False
    movc_viewer: bool = False
    enso_viewer: bool = True

    mean_group1_name: Optional[str] = None
    mean_group2_name: Optional[str] = None
    extra_groups_name: Optional[List[str]] = None

    test_table_id: str = field(init=False, default="Amon")
    test_name: List[str] = field(init=False, default_factory=list)
    test_mip: List[str] = field(init=False, default_factory=list)
    test_case_ids: List[str] = field(init=False, default_factory=list)
    results_dir_full: Path = field(init=False)
    base_test_input_path: str = field(init=False)

    def __post_init__(self) -> None:
        self.diag_path = Path(self.diag_path)
        self.data_dir = Path(self.data_dir)
        self.out_dir = Path(self.out_dir)
        self.results_dir_full = self.out_dir / self.run_type
        self.base_test_input_path = str(
            self.data_dir
            / "put_model_here"
            / "pcmdi_diags"
            / self.run_type
            / "metrics_data"
            / "%(group_type)"
        )
        self.clim_vars = as_list(self.clim_vars)
        self.clim_regions = as_list(self.clim_regions)
        self.atm_modes = as_list(self.atm_modes)
        self.cpl_modes = as_list(self.cpl_modes)
        self.mean_group1_name = self.mean_group1_name or f"{self.ref_group} (mean)"
        self.mean_group2_name = self.mean_group2_name or f"{self.test_group} (mean)"
        self._discover_models()

    def _discover_models(self) -> None:
        model_dirs = sorted(self.data_dir.glob(f"{self.test_pattern}_*"))
        self.test_name = [path.name for path in model_dirs if path.is_dir()]
        self.test_mip = [
            f"{self.test_prefix}.{name.split(f'{self.test_pattern}_')[-1]}"
            for name in self.test_name
        ]
        self.test_case_ids = [self.test_case_id] * len(self.test_name)

    @property
    def ref_clim_dir(self) -> str:
        return str(self.ref_dataset.clim_dir)

    @property
    def ref_movs_dir(self) -> str:
        return str(self.ref_dataset.movs_dir)

    @property
    def ref_enso_dir(self) -> str:
        return str(self.ref_dataset.enso_dir)

    @property
    def test_clim_dir(self) -> Optional[str]:
        return str(self.test_dataset.clim_dir) if self.test_dataset else None

    @property
    def test_movs_dir(self) -> Optional[str]:
        return str(self.test_dataset.movs_dir) if self.test_dataset else None

    @property
    def test_enso_dir(self) -> Optional[str]:
        return str(self.test_dataset.enso_dir) if self.test_dataset else None

    def summary(self) -> None:
        print(f"Test group: {self.test_group}")
        print(f"Reference group: {self.ref_group}")
        print(f"Models ({len(self.test_name)}): {self.test_name}")
        print(f"Results: {self.results_dir_full}")


def as_list(value: Optional[Iterable[str] | str]) -> List[str]:
    """Normalize None, comma-separated strings, or iterables into a clean list."""
    if value is None:
        return []
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(item).strip() for item in value if str(item).strip()]


def load_metric_dict(metric_file: Path = DEFAULT_METRIC_FILE) -> Dict[str, Any]:
    with Path(metric_file).open() as handle:
        return json.load(handle)


def dataset_for_group(
    group: str, case_id: str, root: Path = PCMDI_METRICS_ROOT
) -> MetricsDataset:
    """Return metric-data locations for CMIP or an E3SM aggregate group."""
    root = Path(root)
    group_key = group.lower()

    if group_key == "cmip":
        return MetricsDataset(
            clim_dir=root / "mean_climate",
            clim_set="cmip6.historical.v20250707",
            movs_dir=root / "variability_modes",
            movs_set="cmip6.historical.v20220825",
            enso_dir=root / "enso_metric",
            enso_set="cmip6.historical.v20210620",
        )

    if group_key == "e3smv3-future":
        dataset_prefix = "e3sm.future"
    elif group_key == "e3smv3-historical":
        dataset_prefix = "e3sm.historical"
    else:
        raise ValueError(f"Unknown metrics group: {group!r}")

    dataset_name = f"{dataset_prefix}.{case_id}"
    return MetricsDataset(
        clim_dir=root / "mean_climate",
        clim_set=dataset_name,
        movs_dir=root / "variability_modes",
        movs_set=dataset_name,
        enso_dir=root / "enso_metric",
        enso_set=dataset_name,
    )


def dataset_for_merged_group(
    *,
    mip: str,
    group: str,
    case_id: str,
    root: Path = PCMDI_METRICS_ROOT,
) -> MetricsDataset:
    """Return locations for a merged group produced by metrics_group_merger."""
    dataset_name = f"{mip}.{group}.{case_id}"
    root = Path(root)
    return MetricsDataset(
        clim_dir=root / "mean_climate",
        clim_set=dataset_name,
        movs_dir=root / "variability_modes",
        movs_set=dataset_name,
        enso_dir=root / "enso_metric",
        enso_set=dataset_name,
    )


def modes_for_reference(
    ref_group: str, align_with_cmip: bool = True
) -> tuple[List[str], List[str]]:
    if ref_group == "CMIP" or align_with_cmip:
        return as_list("NAM,NAO,PNA,NPO"), as_list("PDO,NPGO")
    return as_list("NAM,NAO,PNA,NPO,SAM,PSA1,PSA2"), as_list("PDO,NPGO,AMO")


def active_figure_sets(parameters: SyntheticPlotsParameters) -> List[str]:
    figure_sets = []
    if parameters.clim_viewer:
        figure_sets.append("mean_climate")
    if parameters.mova_viewer or parameters.movc_viewer:
        figure_sets.append("variability")
    if parameters.enso_viewer:
        figure_sets.append("enso")
    return figure_sets


def figure_styles(
    parameters: SyntheticPlotsParameters,
    title_label: str,
    panel_label: str = " " * 30,
) -> Dict[str, PlotStyle]:
    clim_scope = "mean" if parameters.test_model_only else "all"
    return {
        "mean_climate": PlotStyle(
            title=f"{panel_label} {title_label} (Mean Climate)",
            font_size=40,
            figure_size=(50.0, 20.0),
            out_dir=parameters.out_dir
            / f"clim_vs_{parameters.error_norm}_{clim_scope}_nocmip"
            / parameters.run_type
            / parameters.ref_group,
        ),
        "variability": PlotStyle(
            title=f"{panel_label} {title_label} (Variability Modes)",
            font_size=40,
            figure_size=(80.0, 30.0),
            out_dir=parameters.out_dir
            / f"movs_{parameters.atm_obs}_{parameters.cpl_obs}_{parameters.movs_group}_nocmip"
            / parameters.run_type
            / parameters.ref_group,
        ),
        "enso": PlotStyle(
            title=f"{panel_label} {title_label} (ENSO)",
            font_size=40,
            figure_size=(50.0, 20.0),
            out_dir=parameters.out_dir
            / "enso_with_feedback_nocmip"
            / parameters.run_type
            / parameters.ref_group,
        ),
    }


def make_parameters(
    *,
    case_id: str = "v20260212",
    ref_group: str = "CMIP",
    test_group: str = "E3SMv3-Historical",
    ref_dataset: Optional[MetricsDataset] = None,
    test_dataset: Optional[MetricsDataset] = None,
    test_highlight_models: Optional[List[str]] = None,
    align_with_cmip: bool = True,
    test_combined: bool = False,
    test_model_only: bool = False,
    show_mean_columns: bool = True,
    plot_mean_groups: Optional[bool] = None,
    test_prefix: str = "e3sm.historical.v3-LR",
    test_pattern: str = "v3.LR.historical",
    clim_vars: Iterable[str] | str = DEFAULT_CLIM_VARS,
    clim_regions: Iterable[str] | str = "global",
    movs_group: str = "cbf",
    error_norm: str = "reference",
    exclude_vars: Optional[Dict[str, List[str]]] = None,
    exclude_models: Optional[List[str]] = None,
    atm_modes: Optional[Iterable[str] | str] = None,
    atm_obs: str = "NOAA-20C",
    cpl_modes: Optional[Iterable[str] | str] = None,
    cpl_obs: str = "HadISST",
    clim_viewer: bool = False,
    mova_viewer: bool = False,
    movc_viewer: bool = False,
    enso_viewer: bool = True,
    metrics_root: Path = PCMDI_METRICS_ROOT,
    diag_root: Path = DEFAULT_DIAG_ROOT,
    data_subdir: str = "climo",
    output_subdir: str = "CLIM_Metrics",
    data_dir: Optional[Path] = None,
    out_dir: Optional[Path] = None,
    run_type: str = "model_vs_obs",
    figure_format: str = "pdf",
    mean_group1_name: Optional[str] = None,
    mean_group2_name: Optional[str] = None,
    extra_groups_name: Optional[List[str]] = None,
) -> SyntheticPlotsParameters:
    """Build a parameter bundle from the small set of notebook knobs."""
    default_atm_modes, default_cpl_modes = modes_for_reference(ref_group, align_with_cmip)
    atm_modes = atm_modes if atm_modes is not None else default_atm_modes
    cpl_modes = cpl_modes if cpl_modes is not None else default_cpl_modes
    ref_dataset = ref_dataset or dataset_for_group(ref_group, case_id, root=metrics_root)
    if test_combined and test_dataset is None:
        test_dataset = dataset_for_group(test_group, case_id, root=metrics_root)
    if plot_mean_groups is not None:
        show_mean_columns = plot_mean_groups

    return SyntheticPlotsParameters(
        test_group=test_group,
        test_prefix=test_prefix,
        test_pattern=test_pattern,
        test_case_id=case_id,
        test_model_only=test_model_only,
        test_combined=test_combined,
        show_mean_columns=show_mean_columns,
        test_dataset=test_dataset,
        test_highlight_models=test_highlight_models,
        ref_group=ref_group,
        ref_dataset=ref_dataset,
        clim_vars=clim_vars,
        clim_regions=clim_regions,
        movs_group=movs_group,
        error_norm=error_norm,
        exclude_vars=exclude_vars if exclude_vars is not None else CMIP_EXCLUDE_VARS,
        exclude_models=exclude_models if exclude_models is not None else CMIP_EXCLUDE_MODELS,
        atm_modes=atm_modes,
        atm_obs=atm_obs,
        cpl_modes=cpl_modes,
        cpl_obs=cpl_obs,
        diag_path=diag_root,
        data_dir=Path(data_dir) if data_dir is not None else Path(diag_root) / data_subdir,
        out_dir=Path(out_dir) if out_dir is not None else Path(diag_root) / output_subdir,
        run_type=run_type,
        save_all_data=True,
        figure_format=figure_format,
        clim_viewer=clim_viewer,
        mova_viewer=mova_viewer,
        movc_viewer=movc_viewer,
        enso_viewer=enso_viewer,
        mean_group1_name=mean_group1_name,
        mean_group2_name=mean_group2_name,
        extra_groups_name=extra_groups_name,
    )


def make_comparison_parameters(**kwargs: Any) -> SyntheticPlotsParameters:
    """Generic entry point for any model/data comparison application."""
    return make_parameters(**kwargs)


def make_hist_vs_cmip_parameters(**overrides: Any) -> SyntheticPlotsParameters:
    """Preset for E3SM historical members compared with CMIP."""
    defaults = {
        "case_id": "v20260212",
        "ref_group": "CMIP",
        "test_group": "E3SMv3-Historical",
        "align_with_cmip": True,
        "test_combined": False,
        "test_model_only": False,
        "show_mean_columns": True,
        "data_subdir": "climo",
        "output_subdir": "CLIM_Metrics",
    }
    defaults.update(overrides)
    return make_parameters(**defaults)


def make_hist_vs_future_parameters(**overrides: Any) -> SyntheticPlotsParameters:
    """Preset for E3SM future aggregate compared with historical aggregate."""
    defaults = {
        "case_id": "v20260212",
        "ref_group": "E3SMv3-Historical",
        "test_group": "E3SMv3-Future",
        "align_with_cmip": True,
        "test_combined": True,
        "test_model_only": False,
        "show_mean_columns": True,
        "data_subdir": "future",
        "output_subdir": "Future_Metrics",
    }
    defaults.update(overrides)
    return make_parameters(**defaults)


def make_plotter(
    parameters: SyntheticPlotsParameters,
    metric_dict: Dict[str, Any],
    style: PlotStyle,
):
    try:
        from .synthetic_metrics_plotter import SyntheticMetricsPlotter
    except ImportError:
        from synthetic_metrics_plotter import SyntheticMetricsPlotter

    enso_options = {
        "reduced_set": False,
        "met_order": None,
        "mod_order": None,
        "sort_y_names": False,
        "show_proj_means": False,
        "show_ref_row": False,
        "show_alt_obs_rows": False,
        "highlight_cmip": True,
    }

    return SyntheticMetricsPlotter(
        test_group=parameters.test_group,
        test_mip=parameters.test_mip,
        test_name=parameters.test_name,
        test_highlight_models=parameters.test_highlight_models,
        test_case_id=parameters.test_case_ids,
        test_table_id=parameters.test_table_id,
        test_combined=parameters.test_combined,
        show_mean_columns=parameters.show_mean_columns,
        figure_format=parameters.figure_format,
        metric_dict=metric_dict,
        save_data=parameters.save_all_data,
        base_test_input_path=parameters.base_test_input_path,
        ref_group=parameters.ref_group,
        clim_viewer=parameters.clim_viewer,
        clim_vars=parameters.clim_vars,
        clim_regions=parameters.clim_regions,
        test_clim_dir=parameters.test_clim_dir,
        test_clim_set=parameters.test_dataset.clim_set if parameters.test_dataset else None,
        ref_clim_dir=parameters.ref_clim_dir,
        ref_clim_set=parameters.ref_dataset.clim_set,
        mova_viewer=parameters.mova_viewer,
        mova_modes=parameters.atm_modes,
        mova_obs=parameters.atm_obs,
        movc_viewer=parameters.movc_viewer,
        movc_modes=parameters.cpl_modes,
        movc_obs=parameters.cpl_obs,
        test_movs_dir=parameters.test_movs_dir,
        test_movs_set=parameters.test_dataset.movs_set if parameters.test_dataset else None,
        ref_movs_dir=parameters.ref_movs_dir,
        ref_movs_set=parameters.ref_dataset.movs_set,
        enso_viewer=parameters.enso_viewer,
        ref_enso_dir=parameters.ref_enso_dir,
        ref_enso_set=parameters.ref_dataset.enso_set,
        test_enso_dir=parameters.test_enso_dir,
        test_enso_set=parameters.test_dataset.enso_set if parameters.test_dataset else None,
        test_model_only=parameters.test_model_only,
        movs_group=parameters.movs_group,
        exclude_vars=parameters.exclude_vars,
        exclude_models=parameters.exclude_models,
        error_norm=parameters.error_norm,
        mean_group1_name=parameters.mean_group1_name,
        mean_group2_name=parameters.mean_group2_name,
        extra_groups_name=parameters.extra_groups_name,
        font_size=style.font_size,
        legend_lw=style.legend_lw,
        figure_size=style.figure_size,
        figure_title=style.title,
        out_dir=str(style.out_dir),
        **enso_options,
    )


def run_synthetic_plots(
    parameters: SyntheticPlotsParameters,
    *,
    title_label: str,
    metric_file: Path = DEFAULT_METRIC_FILE,
    debug: bool = True,
) -> None:
    """Generate all enabled synthetic metric figures."""
    metric_dict = load_metric_dict(metric_file)
    figure_sets = active_figure_sets(parameters)
    styles = figure_styles(parameters, title_label)

    parameters.summary()
    logger.info("Generating synthetic metrics plots for case_id=%s", parameters.test_case_id)

    for figure_set in figure_sets:
        style = styles[figure_set]
        logger.info(
            "Generating figure_set=%s for case_id=%s",
            figure_set,
            parameters.test_case_id,
        )
        print(figure_set)
        plotter = make_plotter(parameters, metric_dict, style)
        plotter.generate(figure_sets=[figure_set], debug=debug)
