"""Fixed output locations shared by all diagnostic workflows."""

from pathlib import Path

DIAG_OUTPUT_ROOT = Path("/compyfs/www/zhan391/e3sm_dart/diag_dart_2026")
OUTPUT_DATA_ROOT = DIAG_OUTPUT_ROOT / "data"
OUTPUT_FIGURE_ROOT = DIAG_OUTPUT_ROOT / "figure"

WORKFLOW_NAMES = (
    "analysis_atm_da",
    "analysis_atm_init",
    "analysis_atm_bias",
    "analysis_lnd_init",
    "analysis_lnd_bias",
    "analysis_lac",
    "analysis_atm_model_evaluation",
    "analysis_lnd_model_evaluation",
    "analysis_mjo",
    "fcst_atm_bias",
    "fcst_atm_initial_shock",
    "fcst_atm_timeseries",
    "fcst_lnd_bias",
    "fcst_atm_model_evaluation",
    "fcst_lnd_model_evaluation",
    "fcst_s2s_pcc",
    "fcst_s2s_tcc",
)

def workflow_output_dirs(
    workflow: str,
    *subdirs: str,
    output_root=DIAG_OUTPUT_ROOT,
    create: bool = True,
):
    """Return fixed data and figure directories for one workflow."""
    if workflow not in WORKFLOW_NAMES:
        available = ", ".join(WORKFLOW_NAMES)
        raise ValueError(f"Unknown workflow {workflow!r}; choose one of: {available}")
    output_root = Path(output_root).expanduser().resolve()
    data_dir = output_root / "data" / workflow
    figure_dir = output_root / "figure" / workflow
    for subdir in subdirs:
        data_dir /= subdir
        figure_dir /= subdir
    if create:
        data_dir.mkdir(parents=True, exist_ok=True)
        figure_dir.mkdir(parents=True, exist_ok=True)
    return data_dir, figure_dir
