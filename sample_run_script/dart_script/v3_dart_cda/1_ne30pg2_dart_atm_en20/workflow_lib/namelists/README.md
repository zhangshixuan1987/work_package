# Workflow-owned DART namelists

These files are the authoritative runtime templates used by the numbered
workflow stages:

- `eam/perturb.nml`: Step 3 initial EAM ensemble perturbation.
- `eam/filter.nml`: Step 4 EAM filter assimilation.
- `eam/diagnostics.nml`: Step 6 EAM DART diagnostics.
- `elm/filter.nml`: Step 4 ELM filter assimilation.

Workers copy a template into a stage- or cycle-specific run directory as
`input.nml`, then apply validated runtime values such as ensemble size,
localization, inflation, task layout, and optional coupled-DA settings. Edit
the template here, not a generated runtime copy.

The `input.nml` files beneath DART model `work/` directories remain build-time
inputs for utilities such as `preprocess`; they are not the authoritative
runtime filter templates for this workflow.
