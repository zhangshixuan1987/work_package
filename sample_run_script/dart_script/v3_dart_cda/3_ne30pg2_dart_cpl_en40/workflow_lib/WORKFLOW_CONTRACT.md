# Shared workflow behavior contract

The `scp_en10` and `cpl_en40` experiments use the same active workflow engine.
Their `create_and_setup_case.sh` files contain experiment-specific values, but
must define the same E3SM, EAM, ELM, environment, and strongly coupled DA
interfaces.

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

The only active source difference permitted between experiment trees is the
Step 4 Slurm walltime/node request. Runtime state, archives, configurations,
and README files are not synchronized as engine source.
