from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import xarray as xr

from metrics_group_merger import PCMDIRun
from pcmdi_enso_reader import ENSODiagReader


@dataclass(frozen=True)
class EnsembleSource:
    data_dir: Path
    model: str
    group: str
    period: tuple
    nens: int
    members: Optional[Sequence[int]] = None
    case_id: Optional[str] = None


def load_ensemble_variable(
    source: EnsembleSource,
    metric_name: str,
    metric_config: dict,
    variable: str,
    obs_tag: str,
    default_case_id: str,
    *,
    verbose: bool = False,
    pool_members: bool = False,
):
    case_id = source.case_id or default_case_id
    reader = ENSODiagReader(
        data_dir=str(source.data_dir),
        model=source.model,
        groups=[source.group],
        period_list=[source.period],
        nens=[source.nens],
        members=source.members,
        verbose=verbose,
    )
    model_data, observation_data = reader.load_metric_data(
        enso_group=metric_config["group"],
        var_name=metric_name,
        nc_var=variable,
        ref_dict={source.group: obs_tag},
        case_id=case_id,
    )
    reference = reader.validate_constant_observation(
        observation_data,
        ref_group=source.group,
        ref_member="00",
        sample_dim=None,
        use_allclose=True,
        rtol=metric_config.get("rtol", 1e-2),
        atol=0.0,
        pool_ensemble=pool_members,
    )
    if pool_members:
        model = reader.pool_members_to_samples(model_data[source.group], sample_dim=None)
    else:
        model = reader.combine_members_to_array(model_data[source.group], sample_dim=None)
    return model, reference, case_id


def find_pcmdi_metric_file(run: PCMDIRun, metric_config: dict) -> Path:
    metric_dir = run.metrics_data_dir() / "enso_metric" / metric_config["group"]
    token = run.metrics_case_id or ""
    pattern = f"*{token}*_{metric_config['suffix']}.nc"
    matches = sorted(metric_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {metric_dir / pattern}")
    return matches[-1]


def _pick_variable(ds: xr.Dataset, base_name: str, tag: str) -> str:
    exact = f"{base_name}__{tag}"
    if exact in ds.data_vars:
        return exact
    matches = [
        name for name in ds.data_vars
        if name.startswith(f"{base_name}__") and tag in name
    ]
    if len(matches) != 1:
        raise KeyError(
            f"Could not uniquely identify {base_name!r} tagged {tag!r}; "
            f"candidates={matches}"
        )
    return matches[0]


def load_pcmdi_variable(
    run: PCMDIRun,
    metric_config: dict,
    variable: str,
    obs_tag: Optional[str] = None,
):
    path = find_pcmdi_metric_file(run, metric_config)
    with xr.open_dataset(path, decode_times=False) as ds:
        model_name = _pick_variable(ds, variable, run.model_name)
        if obs_tag is None:
            candidates = [
                name for name in ds.data_vars
                if name.startswith(f"{variable}__") and name != model_name
            ]
            if len(candidates) != 1:
                raise KeyError(
                    f"Set obs_tag for {run.output_name or run.model_name}; "
                    f"candidates={candidates}"
                )
            obs_name = candidates[0]
        else:
            obs_name = _pick_variable(ds, variable, obs_tag)
        model = ds[model_name].squeeze().load()
        reference = ds[obs_name].squeeze().load()
    model = model.expand_dims(member=[run.output_name or run.model_name])
    return model, reference, run.metrics_case_id
