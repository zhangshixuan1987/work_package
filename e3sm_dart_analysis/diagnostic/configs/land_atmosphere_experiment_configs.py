from pathlib import Path

# Experiment paths are case roots; member series live under
# <case>/<member>/archive/post/<component>/180x360_aave/ts (see member_ts_root).

DEFAULT_DATA_ROOT = Path("/compyfs/zhan391/v3_dart_cda_scratch")
TS_SUBDIR = "archive/post/{component}/180x360_aave/ts"

exp_dict1 = {
    'CTRLEN10': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/CTRLEN10_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 10
    },
    'CAPTEN10': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/CAPTEN10_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 10
    },
    'DARTEN20': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/DARTEN20_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 20
    },
    'DARTEN40': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/DARTEN40_15day_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 40
    }
}

exp_dict2 = {
    'CTRLEN10': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/CTRLEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 10
    },
    'DARTEN20': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/DARTEN20_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 20
    },
    'DARTEN40': {
        'path': '/compyfs/zhan391/v3_dart_cda_scratch/DARTEN40_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy',
        'template': '%(variable)s.%(ensemble)s.%(year)s.nc',
        'nens': 40
    }
}

EXPERIMENT_GROUPS = {
    "v3_spinup": exp_dict2,
    "v3_hindcast": exp_dict1,
}


def get_experiment_dict(key, *, data_root=DEFAULT_DATA_ROOT):
    try:
        experiments = EXPERIMENT_GROUPS[key]
    except KeyError as exc:
        available = ", ".join(sorted(EXPERIMENT_GROUPS))
        raise ValueError(f"Unknown experiment key: {key}. Available: {available}") from exc

    data_root = Path(data_root)
    rebased = {}
    for name, spec in experiments.items():
        item = dict(spec)
        source_path = Path(item["path"])
        item["path"] = str(data_root / source_path.relative_to(DEFAULT_DATA_ROOT))
        rebased[name] = item
    return rebased



def member_ts_root(case_root, member, component="atm"):
    """``ts`` directory of one member; util.model_daily picks the subfolder."""
    return Path(case_root) / member / TS_SUBDIR.format(component=component)
