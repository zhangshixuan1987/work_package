# e3sm_dart_analysis

## Ensemble Perturbation Workflows

The top-level [ensemble_perturbations](ensemble_perturbations/) directory builds the 20-member initial ensembles in four numbered steps: select 20 representative DART members (01), then one notebook per initial-condition set: CLIM_BASE + DART perturbations (02), ERA5 completed with DART-mean fields + DART perturbations (03), and ERA5 completed with CLIM_BASE fields + DART perturbations (04). See [ensemble_perturbations/README.md](ensemble_perturbations/README.md); the shared code lives in `ensemble_perturbations/util/`.

## ERA5 Initial-Condition Tests

The single-member test hindcasts (start 16 Dec 2011) that address the issues found in the ERA5-initialized F20TR hindcasts are evaluated in [diagnostic/jupyter/hind_ic_test](diagnostic/jupyter/hind_ic_test/): (01) completeness of the atmosphere initial files, (02) daily regional drift and bias against ERA5, and (03) lead-window difference maps. See [diagnostic/jupyter/hind_ic_test/METHODS.md](diagnostic/jupyter/hind_ic_test/METHODS.md); the code is `diagnostic/configs/hind_ic_test_runs.py` and `diagnostic/util/hind_ic_test_eval.py`.
