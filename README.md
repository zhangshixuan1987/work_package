# E3SM research workflows

This repository collects scripts, configuration, reference data, and diagnostic
outputs used for E3SM atmosphere experiments. Most workflows are independent;
run scripts from their own directory unless their documentation says otherwise.

## Repository layout

| Directory | Purpose |
| --- | --- |
| `CMIP_PGW_forcing/` | Generate and diagnose CMIP6 pseudo-global-warming SST and sea-ice forcing. |
| `Gen_Nudge_Data/` | Generate reanalysis-based nudging and initial-condition data. |
| `Tempest_Extreme/` | Detect and analyze tropical cyclones, extratropical cyclones, and atmospheric rivers. |
| `analysis_script_ncl/` | Shared NCL analysis, plotting, animation, and regridding scripts. |
| `betacast/` | BetaCast configuration and regridding resources. |
| `circulation_pause/` | Circulation, ozone, and related diagnostics. |
| `nudging_data_process/` | Post-process and visualize nudging experiments. |
| `rrm_domain/` | Regionally refined model domain resources. |
| `sample_run_script/` | Example E3SM run scripts and case-specific resources. |
| `tcTrack_betadrift/` | Tropical-cyclone beta-drift analysis. |

## Working with the repository

- Treat paths, machine names, case names, and date ranges in scripts as
  experiment-specific configuration; review them before submitting a job.
- Keep source scripts and small, reproducible configuration files in Git.
- Do not commit editor swap files, caches, scheduler logs, temporary output, or
  generated figures. PNG, GIF, PDF, and EPS outputs are ignored repository-wide;
  use `git add -f` only when a figure is an intentional reference artifact.
- Large NetCDF files are managed by Git LFS, as configured in
  `.gitattributes`.
- Put workflow-specific usage notes beside the relevant scripts when adding or
  changing a workflow.

There is no repository-wide test suite. Validate changes with the interpreter
or model workflow used by the directory you modify.
