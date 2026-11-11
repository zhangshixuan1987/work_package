# e3sm_dart_analysis

## Ensemble Perturbation Workflows

The top-level [ensemble_perturbations](ensemble_perturbations/) directory builds the 20-member initial ensembles in four numbered steps: select 20 representative DART members (01), then one notebook per initial-condition set: CLIM_BASE + DART perturbations (02), ERA5 completed with DART-mean fields + DART perturbations (03), and ERA5 completed with CLIM_BASE fields + DART perturbations (04). See [ensemble_perturbations/README.md](ensemble_perturbations/README.md); the shared code lives in `ensemble_perturbations/util/`.

## ERA5 Initial-Condition Tests

The top-level [era5_ic_test](era5_ic_test/) directory evaluates the single-member test hindcasts (start 16 Dec 2011) that address the issues found in the ERA5-initialized F20TR hindcasts. It holds three notebooks: (01) completeness of the atmosphere initial files, (02) daily regional drift and bias against ERA5, and (03) lead-window difference maps. The tests are the completed ERA5 states and the one-component CLIM_BASE swaps of DA_COUPLED. See [era5_ic_test/README.md](era5_ic_test/README.md); the code lives in `era5_ic_test/ic_util/`.
