# PCMDI diagnostics and paper-plot workflows

This repository contains reusable readers and plotting workflows for PCMDI
mean-climate, modes-of-variability (MOV), ENSO, and tropical-cyclone
diagnostics. The notebooks under `jupyter/` provide the user-facing experiment
configuration; reusable loading, merging, analysis, and plotting logic lives
under `scripts/`.

## Notebook organization

Notebook names generally follow:

```text
plot_<diagnostic>_<figure-or-workflow>_<experiment>.ipynb
```

The main notebook groups are:

- Mean climate: `plot_clim_*.ipynb`
- Modes of variability: `plot_mov_*.ipynb`
- ENSO: `plot_enso_*.ipynb`
- Tropical cyclones: `plot_tc_*.ipynb`

### Synthetic metrics: model versus observations

The synthetic model-versus-observations (MVO) workflows are separated by
diagnostic family. Each notebook merges and plots only its associated metrics:

| Experiment | Mean climate | Variability modes | ENSO |
| --- | --- | --- | --- |
| v3 large ensemble | `plot_clim_synthetic_mvo_v3le.ipynb` | `plot_mov_synthetic_mvo_v3le.ipynb` | `plot_enso_synthetic_mvo_v3le.ipynb` |
| v4P standalone runs | `plot_clim_synthetic_mvo_v4p.ipynb` | `plot_mov_synthetic_mvo_v4p.ipynb` | `plot_enso_synthetic_mvo_v4p.ipynb` |

For the v4P notebooks, the standalone `PCMDIRun` definitions are preserved in
each notebook, while the merge switches and plotting selection enable only one
of `mean_climate`, `variability_modes`, or `enso_metric`. This lets each
diagnostic family be configured and executed independently.

### Synthetic metrics: model versus model

The v3 large-ensemble model-versus-model (MVM) workflow is also split into
independent notebooks:

- `plot_clim_synthetic_mvm_v3le.ipynb`
- `plot_mov_synthetic_mvm_v3le.ipynb`
- `plot_enso_synthetic_mvm_v3le.ipynb`

### Feature and pattern plots

- ENSO feature metrics: `plot_enso_feature_metric_*.ipynb`
- ENSO spatial and temporal patterns: `plot_enso_pattern_*.ipynb`
- MOV pattern maps: `plot_mov_pattern_map_rrm.ipynb`
- MOV pathway metrics: `plot_mov_pathway_metric_rrm.ipynb`
- Mean-climate bias maps: `plot_clim_bias_map.ipynb`

### Data preparation utilities

- `fix_cmip_metrics_data.ipynb`: patch selected CMIP mean-climate JSON metrics.
- `gen_merge_metrics_data.ipynb`: merge configured climate, MOV, and ENSO
  metric groups.
- `gen_tc_track.ipynb`: validate inputs, preview commands, and run
  TempestExtremes tropical-cyclone tracking.

These notebooks default to non-writing or dry-run behavior. Review their
resolved paths and explicitly enable the corresponding `RUN_*` switch before
writing outputs or launching production processing.

### Tropical-cyclone workflow

The TC workflow is separated into generation and analysis stages:

1. `gen_tc_track.ipynb` validates the configured simulations and invokes the
   repository-local `scripts/tc_track.py` TempestExtremes driver.
2. `plot_tc_track_diagnostics.ipynb` constructs lead-time, observational, and
   ENSO-regression diagnostics.
3. `plot_tc_track_density_map.ipynb` compares track-density and ENSO-regression
   maps across configured runs.

TC production operations are explicitly controlled by `RUN_TRACKING`,
`WRITE_DIAGNOSTICS`, `RUN_PLOTS`, `RUN_ENSO_REGRESSION`, and `SAVE_FIGURES`.
Their defaults favor command preview, input inspection, and in-memory results
instead of modifying output products.

## Running a notebook

1. Open the notebook for the desired diagnostic family and experiment.
2. Review its configuration cell, especially input/output paths, run or
   ensemble definitions, metric selections, and viewer switches.
3. Keep the plotting switch disabled while checking resolved inputs when the
   notebook provides a dry-run stage.
4. Enable plotting and execute the workflow cells in order.

The configured output directory is printed by the workflow. Many production
examples reference LCRC or NERSC diagnostic locations and therefore require
access to the corresponding filesystem.

## Reusable modules

Important shared modules include:

- `scripts/metrics_group_merger.py`: assemble synthetic metric groups.
- `scripts/synthetic_metrics_workflow.py`: configure and run synthetic plots.
- `scripts/synthetic_metrics_plotter.py`: shared synthetic figure drivers.
- `scripts/clim_metrics_reader.py`: mean-climate metric input.
- `scripts/movs_metrics_reader.py`: variability-mode metric input.
- `scripts/enso_metrics_reader.py`: ENSO metric input.
- `scripts/enso_feature_sources.py`: ensemble and standalone ENSO sources.
- `scripts/pcmdi_mov_workflow.py`: MOV pattern processing and plotting.
- `scripts/mov_pathway_workflow.py`: MOV pathway metric workflow.
- `scripts/tc_track.py`: TempestExtremes tracking command-line driver.
- `scripts/tc_track_density.py`: reusable box-count and radius-based track
  density calculations.

Experiment-specific paths, labels, and plotting choices should remain in the
notebooks; generally reusable data handling and plotting logic belongs in
`scripts/`.
