# Shared workflow behavior contract

The `atm_en20` workflow uses the same active restart-safe workflow engine as the
coupled experiments. Its `create_and_setup_case.sh` retains AMIP-specific values
while defining the same E3SM, EAM, ELM, environment, and coupled-DA interfaces.

The shared E3SM timeline is authoritative. Before Step 4 runs, the configured
completed-cycle count, matching completion record, restart archive, and every
ensemble case XML timestamp must resolve to the same valid time.

ELM execution mode is derived once from both coupling switches:

```text
strongly_coupled_on=on AND lnd_da_use_sequential_prior_post=.true.
    -> sequential mode (EAM then ELM; each receives all nodes)
otherwise
    -> direct mode (simultaneous due components split nodes)
```

`my_eam_dart_da` and `my_elm_dart_da` independently enable the component
analyses. Turning strongly coupled DA off does not implicitly turn ELM DA off.
When both workflows receive equivalent configuration, they must produce the
same cadence, component status, execution mode, node allocation, preflight,
failure recovery, handoff, counter update, and continuation decisions.

Engine changes should remain synchronized with the maintained coupled workflow.
Experiment configuration, Step 4 Slurm resources, runtime state, archives, and
README files remain specific to this atmospheric experiment.
