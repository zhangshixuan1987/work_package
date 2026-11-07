# e3sm_dart_analysis

## Ensemble Perturbation Workflows

The top-level [ensemble_perturbations](ensemble_perturbations/) directory builds the 20-member initial ensembles in four numbered steps: select 20 representative DART members (01), then one notebook per initial-condition set: CLIM_BASE + DART perturbations (02), ERA5 completed with DART-mean fields + DART perturbations (03), and ERA5 completed with CLIM_BASE fields + DART perturbations (04). See [ensemble_perturbations/README.md](ensemble_perturbations/README.md); the shared code lives in `ensemble_perturbations/util/`.
