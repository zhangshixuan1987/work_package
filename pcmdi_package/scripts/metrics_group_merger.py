from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

from post_process_merge_jsons import CLIMJSONMerger, ENSOJSONMerger, MOVSJSONMerger


@dataclass(frozen=True)
class PCMDIRun:
    """One raw PCMDI diagnostic run."""

    case: str
    model_name: str
    www: Path
    metrics_case_id: Optional[str] = None
    output_name: Optional[str] = None

    def metrics_data_dir(self, run_type: str = "model_vs_obs") -> Path:
        return Path(self.www) / self.case / "pcmdi_diags" / run_type / "metrics_data"


@dataclass(frozen=True)
class MetricsGroup:
    """Definition of one merged metrics group."""

    name: str
    model_pattern: str
    mips: List[str]
    exps: List[str]
    case_id: str
    source_root: Optional[Path] = None
    runs: Optional[List[PCMDIRun]] = None
    member_names: Optional[List[str]] = None
    model_rename: Optional[Tuple[str, str, bool]] = None
    movs_years: Tuple[int, int] = (1900, 2014)
    movs_modes: str = "NAM,NAO,PNA,NPO,SAM,PSA1,PSA2,PDO,NPGO,AMO"
    movs_obses: str = (
        "NOAA-20C,NOAA-20C,NOAA-20C,NOAA-20C,NOAA-20C,NOAA-20C,NOAA-20C,"
        "HadISST,HadISST,HadISST"
    )
    enso_collections: str = "ENSO_perf,ENSO_tel,ENSO_proc"
    enso_obses: str = "ERA-Interim,ERA-Interim,ERA-Interim"


def discover_member_dirs(group: MetricsGroup, run_type: str = "model_vs_obs") -> List[str]:
    """Find per-run PCMDI metrics_data directories for a group."""
    if group.runs is not None:
        return [str(run.metrics_data_dir(run_type=run_type)) for run in group.runs]

    if group.source_root is None:
        raise ValueError("source_root is required when runs is not provided")

    source_root = Path(group.source_root)
    if group.member_names is not None:
        model_names = group.member_names
    else:
        pattern = group.model_pattern
        if not any(token in pattern for token in "*?[]"):
            pattern = f"{pattern}_*"
        model_names = sorted(path.name for path in source_root.glob(pattern) if path.is_dir())
    return [
        str(source_root / model / "pcmdi_diags" / run_type / "metrics_data")
        for model in model_names
    ]


def discover_member_model_names(group: MetricsGroup) -> List[str]:
    """Return configured PCMDI model_name values when available."""
    if group.runs is None:
        return []
    return [run.model_name for run in group.runs]


def discover_member_case_ids(group: MetricsGroup) -> List[str]:
    """Return raw per-run metric JSON version ids."""
    if group.runs is None:
        return []
    return [run.metrics_case_id or group.case_id for run in group.runs]


def discover_member_output_names(group: MetricsGroup) -> List[str]:
    """Return per-run labels used as model keys in merged output."""
    if group.runs is None:
        return []
    return [run.output_name or run.model_name for run in group.runs]


class PCMDJSONMerger:
    """Orchestrate JSON merging for mean climate, variability modes, and ENSO metrics."""

    def __init__(
        self,
        pmprdir: str | Path,
        mips: Optional[List[str]] = None,
        exps: Optional[List[str]] = None,
        member_dirs: Optional[List[str]] = None,
        model_rename: Optional[Tuple[str, str, bool]] = None,
        *,
        case_id_clim: Optional[str] = None,
        out_path_clim: Optional[str] = None,
        case_id_movs: Optional[str] = None,
        out_path_movs: Optional[str] = None,
        case_id_enso: Optional[str] = None,
        out_path_enso: Optional[str] = None,
        movs_years: Tuple[int, int] = (1900, 2014),
        movs_obses: Optional[str] = None,
        movs_modes: Optional[str] = None,
        enso_collections: Optional[str] = None,
        enso_obses: Optional[str] = None,
        enable_clim: bool = True,
        enable_movs: bool = True,
        enable_enso: bool = True,
        strict: bool = False,
        verbose: bool = True,
        dry_run: bool = False,
    ) -> None:
        self.pmprdir = Path(pmprdir)
        self.mips = mips or ["cmip6"]
        self.exps = exps or ["historical"]
        self.member_dirs = member_dirs or []
        self.member_model_names = []
        self.member_case_ids = []
        self.member_output_names = []
        self.model_rename = model_rename

        self.case_id_clim = case_id_clim
        self.case_id_movs = case_id_movs
        self.case_id_enso = case_id_enso

        self.out_path_clim = out_path_clim
        self.out_path_movs = out_path_movs
        self.out_path_enso = out_path_enso

        self.movs_years = movs_years
        self.movs_obses = movs_obses
        self.movs_modes = movs_modes
        self.enso_collections = enso_collections
        self.enso_obses = enso_obses

        self.enable_clim = enable_clim
        self.enable_movs = enable_movs
        self.enable_enso = enable_enso

        self.strict = strict
        self.verbose = verbose
        self.dry_run = dry_run

    def _log(self, *msg) -> None:
        if self.verbose:
            print(*msg)

    def run_all(self) -> None:
        if self.enable_clim:
            self.run_clim()
        if self.enable_movs:
            self.run_movs()
        if self.enable_enso:
            self.run_enso()

    def run_clim(self) -> None:
        if not self.case_id_clim:
            self._log("[CLIM] No case_id provided - skipping.")
            return
        out_path = self.out_path_clim or str(self.pmprdir)
        member_meanclim_dirs = [f"{d}/mean_climate" for d in self.member_dirs]
        CLIMJSONMerger(
            mips=self.mips,
            exps=self.exps,
            case_id=self.case_id_clim,
            model_rename=self.model_rename,
            pmprdir=str(self.pmprdir),
            out_path=out_path,
            member_dirs=member_meanclim_dirs,
            member_model_names=self.member_model_names,
            member_case_ids=self.member_case_ids,
            member_output_names=self.member_output_names,
            strict=self.strict,
            verbose=self.verbose,
            dry_run=self.dry_run,
        ).merge_all()

    def run_movs(self) -> None:
        if not self.case_id_movs:
            self._log("[MOVS] No case_id provided - skipping.")
            return
        if not self.movs_obses or not self.movs_modes:
            raise ValueError("movs_obses and movs_modes are required when enable_movs=True")
        out_path = self.out_path_movs or str(self.pmprdir)
        syear, eyear = self.movs_years
        member_movs_dirs = [f"{d}/variability_modes" for d in self.member_dirs]
        MOVSJSONMerger(
            mips=self.mips,
            exps=self.exps,
            case_id=self.case_id_movs,
            model_rename=self.model_rename,
            pmprdir=str(self.pmprdir),
            out_path=out_path,
            member_dirs=member_movs_dirs,
            member_model_names=self.member_model_names,
            member_case_ids=self.member_case_ids,
            member_output_names=self.member_output_names,
            movs_obses=list(self.movs_obses.split(",")),
            movs_modes=list(self.movs_modes.split(",")),
            syear=syear,
            eyear=eyear,
            strict=self.strict,
            verbose=self.verbose,
            dry_run=self.dry_run,
        ).merge_all()

    def run_enso(self) -> None:
        if not self.case_id_enso:
            self._log("[ENSO] No case_id provided - skipping.")
            return
        if not self.enso_collections or not self.enso_obses:
            raise ValueError("enso_collections and enso_obses are required when enable_enso=True")
        out_path = self.out_path_enso or str(self.pmprdir)
        member_enso_dirs = [f"{d}/enso_metric" for d in self.member_dirs]
        ENSOJSONMerger(
            mips=self.mips,
            exps=self.exps,
            case_id=self.case_id_enso,
            model_rename=self.model_rename,
            pmprdir=str(self.pmprdir),
            out_path=out_path,
            member_dirs=member_enso_dirs,
            member_model_names=self.member_model_names,
            member_case_ids=self.member_case_ids,
            member_output_names=self.member_output_names,
            collections=list(self.enso_collections.split(",")),
            observations=list(self.enso_obses.split(",")),
            strict=self.strict,
            verbose=self.verbose,
            dry_run=self.dry_run,
        ).merge_all()


def merge_metrics_group(
    group: MetricsGroup,
    *,
    metrics_root: str | Path,
    run_type: str = "model_vs_obs",
    enable_clim: bool = True,
    enable_movs: bool = True,
    enable_enso: bool = True,
    strict: bool = True,
    verbose: bool = True,
    dry_run: bool = False,
    clean: bool = False,
) -> None:
    """Merge raw per-run JSON metrics into one named group under metrics_root."""
    metrics_root = Path(metrics_root)
    member_dirs = discover_member_dirs(group, run_type=run_type)
    member_model_names = discover_member_model_names(group)
    member_case_ids = discover_member_case_ids(group)
    member_output_names = discover_member_output_names(group)

    clim_dir = metrics_root / "mean_climate"
    movs_dir = metrics_root / "variability_modes"
    enso_dir = metrics_root / "enso_metric"

    if clean:
        clean_group_outputs(
            metrics_dirs=[clim_dir, movs_dir, enso_dir],
            group=group,
        )

    merger = PCMDJSONMerger(
        pmprdir=group.source_root or (group.runs[0].www if group.runs else "."),
        mips=group.mips,
        exps=group.exps,
        member_dirs=member_dirs,
        model_rename=group.model_rename,
        case_id_clim=group.case_id,
        out_path_clim=str(clim_dir),
        case_id_movs=group.case_id,
        out_path_movs=str(movs_dir),
        case_id_enso=group.case_id,
        out_path_enso=str(enso_dir),
        movs_years=group.movs_years,
        movs_modes=group.movs_modes,
        movs_obses=group.movs_obses,
        enso_collections=group.enso_collections,
        enso_obses=group.enso_obses,
        enable_clim=enable_clim,
        enable_movs=enable_movs,
        enable_enso=enable_enso,
        strict=strict,
        verbose=verbose,
        dry_run=dry_run,
    )
    merger.member_model_names = member_model_names
    merger.member_case_ids = member_case_ids
    merger.member_output_names = member_output_names
    merger.run_all()
    move_exp_outputs_to_group(
        metrics_dirs=[clim_dir, movs_dir, enso_dir],
        group=group,
        dry_run=dry_run,
    )


def clean_group_outputs(metrics_dirs: List[Path], group: MetricsGroup) -> None:
    for mip in group.mips:
        for exp in group.exps:
            for metrics_dir in metrics_dirs:
                for path in [
                    metrics_dir / mip / exp / group.case_id,
                    metrics_dir / mip / group.name / group.case_id,
                ]:
                    if path.is_dir():
                        shutil.rmtree(path)
                    elif path.exists():
                        path.unlink()


def move_exp_outputs_to_group(
    metrics_dirs: List[Path],
    group: MetricsGroup,
    *,
    dry_run: bool = False,
) -> None:
    for mip in group.mips:
        for exp in group.exps:
            for metrics_dir in metrics_dirs:
                source = metrics_dir / mip / exp / group.case_id
                target = metrics_dir / mip / group.name / group.case_id
                if source.resolve() == target.resolve() or not source.exists():
                    continue
                if dry_run:
                    print(f"[dry-run] move {source} -> {target}")
                    continue
                if target.is_dir():
                    shutil.rmtree(target)
                elif target.exists():
                    target.unlink()
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(source), str(target))
                rename_group_files(target, old_token=exp, new_token=group.name)


def rename_group_files(root: Path, *, old_token: str, new_token: str) -> None:
    for dirpath, _, files in os.walk(root):
        for filename in files:
            if old_token not in filename:
                continue
            old_path = Path(dirpath) / filename
            new_path = Path(dirpath) / filename.replace(old_token, new_token)
            if new_path.exists():
                new_path.unlink()
            old_path.rename(new_path)
