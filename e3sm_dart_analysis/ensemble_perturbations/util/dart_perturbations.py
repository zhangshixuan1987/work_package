"""DART-derived atmospheric perturbations of a reference EAM initial state (notebooks 02-04).

For the 20 selected DART members S of a season and a member s_i in S,

    dA_i = A_DART,s_i - mean_S(A_DART)
    A_out,EN(i+1) = A_ref + alpha * dA_i,  i = 1..19;  EN01 = A_ref (exact copy).

The 20th selected member is part of the mean but gives no output member, so the output ensemble
mean is offset from the reference by -alpha * dA_omitted / 20; ``validate_written_files`` checks that
closure. Only the perturbed variables change; every other variable, the netCDF structure and the
attributes of the reference file are kept and checked. The non-atmosphere restart files of the
reference are linked into each date directory.

Paths are inputs: a :class:`ReferenceCase` holds templates with ``{date}`` for the reference
atmosphere file and the component restart directory, and :func:`dart_date_sets` takes the DART root.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

# --------------------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------------------

# DART case of each season and its 20 selected members (in order; the last one is omitted), as
# chosen by 01_select_dart_20_member_subset.ipynb and recorded in the first build's metadata.
DART_SEASONS = {
    "dec2011": dict(
        dates=["2011-12-16-00000", "2011-12-21-00000", "2011-12-26-00000", "2011-12-31-00000"],
        case="DARTEN40_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        members=[1, 3, 4, 7, 8, 9, 10, 14, 17, 18, 20, 21, 25, 26, 28, 31, 35, 36, 37, 40],
    ),
    "may2012": dict(
        dates=["2012-05-16-00000", "2012-05-21-00000", "2012-05-26-00000", "2012-05-31-00000"],
        case="DARTEN40S1_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy",
        members=[2, 3, 5, 7, 11, 12, 13, 14, 15, 16, 20, 22, 24, 26, 27, 28, 31, 32, 34, 38],
    ),
}


@dataclass
class DartSeason:
    """One season: its start dates, DART case, DART root and selected members."""
    name: str
    dates: list
    case: str
    dart_root: Path
    members: list
    selection_csv: Path | None = None

    @property
    def rest_root(self) -> Path:
        return Path(self.dart_root) / self.case / "archive" / "rest"

    @property
    def diag_root(self) -> Path:
        return Path(self.dart_root) / self.case / "archive" / "dart_en40"

    def member_file(self, date: str, member: int) -> Path:
        name = f"{self.case}.EN{member:02d}.eam.i.{date}.nc"
        shared = self.rest_root / date / name
        if shared.is_file():
            return shared
        return Path(self.dart_root) / self.case / f"EN{member:02d}" / "archive" / "rest" / date / name

    def forecast_file(self, date: str, member: int) -> Path:
        return self.diag_root / date / f"{self.case}.eam_{member:04d}.e.forecast.{date}"


def dart_date_sets(dart_root: Path, names=None, selection_root: Path | None = None) -> list[DartSeason]:
    """The seasons in ``DART_SEASONS`` (all, or ``names``) with their DART files under ``dart_root``.

    ``selection_root/<first date>/selected_20_members.csv``, when given and present, must agree
    with the recorded member list.
    """
    names = list(DART_SEASONS) if names is None else ([names] if isinstance(names, str) else list(names))
    out = []
    for n in names:
        s = DART_SEASONS[n]
        csv = Path(selection_root) / s["dates"][0] / "selected_20_members.csv" if selection_root else None
        out.append(DartSeason(n, list(s["dates"]), s["case"], Path(dart_root), list(s["members"]), csv))
    return out


@dataclass
class ReferenceCase:
    """The reference state that receives the perturbations.

    ``atm_file`` and ``rest_dir`` are path templates with ``{date}``; ``output_root/<date>/`` gets
    ``<out_prefix>.EN01-EN20.eam.i.<date>.nc`` and links to the files of ``rest_dir``.
    """
    name: str
    atm_file: str
    rest_dir: str
    output_root: Path
    out_prefix: str
    out_stream: str = "eam.i"

    def atm(self, date: str) -> Path:
        return Path(self.atm_file.format(date=date))

    def rest(self, date: str) -> Path:
        return Path(self.rest_dir.format(date=date))


@dataclass
class Settings:
    variables: tuple = ("T", "U", "V", "Q", "PS", "H2OLNZ")
    clipped_vars: tuple = ("Q", "H2OLNZ")      # clipped at q_min when clip_q
    alpha: float = 1.0
    clip_q: bool = True
    q_min: float = 0.0
    n_selected: int = 20
    n_perturbed: int = 19
    write_output_files: bool = True
    overwrite: bool = False                  # rewrite existing outputs
    resume: bool = True                      # validate and reuse existing outputs
    link_components: bool = True
    link_mode: str = "symlink"               # "symlink" or "copy"
    require_ncol_d: bool = True
    expected_ncol_d: int = 48602
    mean_tolerance: dict = field(default_factory=lambda: {
        "T": 1.0e-5, "U": 1.0e-5, "V": 1.0e-5, "Q": 1.0e-10, "PS": 1.0e-2, "H2OLNZ": 1.0e-10})
    physical_limits: dict = field(default_factory=lambda: {
        "T": (150.0, 400.0), "U": (-300.0, 300.0), "V": (-300.0, 300.0), "Q": (0.0, 0.1),
        "PS": (4.0e4, 1.2e5), "H2OLNZ": (0.0, 0.1)})
    warn_q_clip_fraction: float = 1.0e-2
    enforce_max_q_clip_fraction: bool = False
    max_q_clip_fraction: float = 2.0e-2
    component_patterns: tuple = ("*.elm.r.*.nc", "*.mosart.r.*.nc", "*.mpaso.rst.*.nc", "*.mpassi.rst.*.nc",
                                 "*.cpl.r.*.nc", "rpointer.drv", "rpointer.ice", "rpointer.lnd",
                                 "rpointer.ocn", "rpointer.rof")
    optional_component_patterns: tuple = ("*.mpaso.rst.*.nc",)   # DOCN cases have no MPAS-O restart
    grid_atol_degrees: float = 2.0e-5
    float32_rtol: float = 1.0e-6
    float32_atol: float = 1.0e-7

    @property
    def perturbed_members(self) -> list[int]:
        return list(range(2, 2 + self.n_perturbed))

    @property
    def all_members(self) -> list[int]:
        return [1] + self.perturbed_members


# --------------------------------------------------------------------------------------------------
# Member selection and run configuration
# --------------------------------------------------------------------------------------------------

def read_selection_csv(csv_path: Path, expected: int) -> list[int]:
    table = pd.read_csv(csv_path)
    if "member" in table.columns:
        members = table["member"].astype(int).tolist()
    elif "member_label" in table.columns:
        members = [int(str(v).upper().removeprefix("EN")) for v in table["member_label"]]
    else:
        raise ValueError(f"{csv_path}: needs a 'member' or 'member_label' column")
    if len(members) != expected or len(set(members)) != expected:
        raise ValueError(f"{csv_path}: expected {expected} distinct members, found {members}")
    return members


def member_selection(season: DartSeason, st: Settings) -> dict:
    """Selected members (CSV checked against the recorded list), the 19 sources, the omitted one."""
    selected = [int(m) for m in season.members]
    if season.selection_csv is not None and Path(season.selection_csv).is_file():
        from_csv = read_selection_csv(season.selection_csv, st.n_selected)
        if from_csv != selected:
            raise ValueError(f"{season.selection_csv} {from_csv} differs from the recorded selection {selected}")
    if len(selected) != st.n_selected or len(set(selected)) != st.n_selected:
        raise ValueError(f"{season.name}: need {st.n_selected} distinct members, got {selected}")
    if any(m < 1 or m > 40 for m in selected):
        raise ValueError(f"{season.name}: members outside EN01-EN40: {selected}")
    sources, omitted = selected[:st.n_perturbed], selected[st.n_perturbed:]
    if len(omitted) != 1:
        raise ValueError("Exactly one selected DART member must be omitted.")
    mapping = pd.DataFrame({
        "output_member": [1] + st.perturbed_members,
        "output_member_label": [f"EN{m:02d}" for m in [1] + st.perturbed_members],
        "dart_member": [pd.NA] + sources,
        "dart_member_label": ["CONTROL"] + [f"EN{m:02d}" for m in sources],
    })
    return dict(selected=selected, sources=sources, omitted=omitted[0],
                output_to_dart=dict(zip(st.perturbed_members, sources)), mapping=mapping)


def _has_vars(path: Path, variables) -> bool:
    if not path.is_file():
        return False
    try:
        with xr.open_dataset(path, decode_times=False) as ds:
            return all(v in ds and 0 not in ds[v].shape for v in variables)
    except Exception:
        return False


def dart_member_files(season: DartSeason, date: str, members, variables) -> tuple[str, dict]:
    """The DART member files: the eam.i restarts, else the DART forecast files."""
    rest = {m: season.member_file(date, m) for m in members}
    if all(_has_vars(p, variables) for p in rest.values()):
        return "rest_eam_i", rest
    fcst = {m: season.forecast_file(date, m) for m in members}
    if all(_has_vars(p, variables) for p in fcst.values()):
        return "dart_en40_forecast", fcst
    raise FileNotFoundError(f"No usable DART atmosphere files for {season.case} {date}")


def build_run_configs(ref: ReferenceCase, seasons: list[DartSeason], st: Settings) -> list[dict]:
    """One configuration per start date; checks that every input exists."""
    configs = []
    for season in seasons:
        sel = member_selection(season, st)
        for date in season.dates:
            source, files = dart_member_files(season, date, sel["selected"], st.variables)
            configs.append(dict(ref=ref, season=season, date=date, dart_source_type=source,
                                dart_files=files, ref_atm_file=ref.atm(date), ref_rest_dir=ref.rest(date),
                                output_dir=Path(ref.output_root) / date, variables=list(st.variables), **sel))
    missing = [f"{c['date']}: {p}" for c in configs for p in (c["ref_atm_file"], c["ref_rest_dir"]) if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing reference inputs:\n" + "\n".join(missing))
    return configs


def config_table(configs: list[dict]) -> pd.DataFrame:
    return pd.DataFrame([dict(reference=c["ref"].name, season=c["season"].name, date=c["date"],
                              dart_source=c["dart_source_type"], selected=c["selected"], omitted=c["omitted"],
                              ref_atm_file=str(c["ref_atm_file"]), output_dir=str(c["output_dir"]))
                         for c in configs])


# --------------------------------------------------------------------------------------------------
# Compatibility checks
# --------------------------------------------------------------------------------------------------

def _horizontal(ds):
    for dim, lat, lon in (("ncol_d", "lat_d", "lon_d"), ("ncol", "lat", "lon")):
        if dim in ds.sizes:
            if dim == "ncol_d" and not ("lat_d" in ds and "lon_d" in ds):
                lat, lon = "lat", "lon"
            return dim, lat, lon
    return None, None, None


def check_ncol_d(path: Path, variables, label: str, expected: int):
    """The variables must sit on ``ncol_d``, or on ``ncol`` of the ``ncol_d`` size (ERA5 files)."""
    bad = {}
    with xr.open_dataset(path, decode_times=False) as ds:
        for v in variables:
            if v not in ds:
                continue
            dims = tuple(d for d in ds[v].dims if d in {"ncol", "ncol_d"})
            if "ncol_d" in dims or ("ncol" in dims and ds.sizes["ncol"] in {ds.sizes.get("ncol_d"), expected}):
                continue
            bad[v] = dims
    if bad:
        raise ValueError(f"{label}: variables not on the dynamics columns ({expected}): {bad} in {path}")


def _align(da, ref_da, label, st: Settings):
    """``da`` renamed (ncol_d -> ncol when the reference uses ncol), transposed and coordinate-checked."""
    rename = {d: "ncol" for d in da.dims if d == "ncol_d" and "ncol" in ref_da.dims and "ncol" not in da.dims}
    da = da.rename(rename)
    if set(da.dims) != set(ref_da.dims):
        raise ValueError(f"{label}: DART dims {da.dims} vs reference {ref_da.dims}")
    da = da.transpose(*ref_da.dims)
    if da.shape != ref_da.shape:
        raise ValueError(f"{label}: shape {da.shape} vs reference {ref_da.shape}")
    for dim in ref_da.dims:
        if dim not in da.coords or dim not in ref_da.coords or ref_da.sizes[dim] == 1:
            continue
        a, b = np.asarray(da.coords[dim].values), np.asarray(ref_da.coords[dim].values)
        if np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number):
            f32 = any(np.issubdtype(x.dtype, np.floating) and x.dtype.itemsize <= 4 for x in (a, b))
            ok = np.allclose(a, b, rtol=st.float32_rtol if f32 else 1e-10, atol=st.float32_atol if f32 else 1e-12,
                             equal_nan=True)
        else:
            ok = np.array_equal(a, b)
        if not ok:
            raise ValueError(f"{label}: coordinate mismatch for {dim}")
    return da.assign_coords({d: ref_da.coords[d] for d in ref_da.dims if d in ref_da.coords})


def _check_grid(dart_ds, ref_ds, label, st: Settings):
    dd, dlat, dlon = _horizontal(dart_ds)
    rd, rlat, rlon = _horizontal(ref_ds)
    if dd is None or rd is None:
        return
    if dart_ds.sizes[dd] != ref_ds.sizes[rd]:
        raise ValueError(f"{label}: DART {dd}={dart_ds.sizes[dd]} vs reference {rd}={ref_ds.sizes[rd]}")
    if not np.allclose(dart_ds[dlat].values, ref_ds[rlat].values, rtol=0, atol=st.grid_atol_degrees):
        raise ValueError(f"{label}: latitude grid mismatch")
    dlon_diff = (dart_ds[dlon].values - ref_ds[rlon].values + 180.0) % 360.0 - 180.0
    if not np.allclose(dlon_diff, 0.0, rtol=0, atol=st.grid_atol_degrees):
        raise ValueError(f"{label}: longitude grid mismatch")


def preflight(cfg: dict, st: Settings):
    """Every selected member has every variable on the reference grid; dynamics columns used."""
    if st.require_ncol_d:
        first = next(iter(cfg["dart_files"].values()))
        check_ncol_d(cfg["ref_atm_file"], st.variables, f"{cfg['date']} reference", st.expected_ncol_d)
        check_ncol_d(first, st.variables, f"{cfg['date']} DART", st.expected_ncol_d)
    errors = []
    with xr.open_dataset(cfg["ref_atm_file"], decode_times=False, cache=False) as ref_ds:
        errors += [f"reference lacks {v}" for v in st.variables if v not in ref_ds]
        for m in cfg["selected"]:
            path = cfg["dart_files"][m]
            try:
                with xr.open_dataset(path, decode_times=False, cache=False) as ds:
                    _check_grid(ds, ref_ds, f"EN{m:02d}", st)
                    for v in st.variables:
                        if v not in ds:
                            errors.append(f"EN{m:02d} lacks {v}: {path}")
                        elif v in ref_ds:
                            _align(ds[v], ref_ds[v], f"EN{m:02d} {v}", st)
            except Exception as exc:
                errors.append(f"EN{m:02d} {path}: {exc}")
    if errors:
        raise ValueError(f"{cfg['date']}: DART member check failed:\n" + "\n".join(errors))


# --------------------------------------------------------------------------------------------------
# Perturbation, writing and validation
# --------------------------------------------------------------------------------------------------

def _load(path, var, ref_ds, st):
    with xr.open_dataset(path, decode_times=False, cache=False) as ds:
        da = ds[var].load().astype("float64")
    return _align(da, ref_ds[var], f"{Path(path).name} {var}", st)


def _finite(da, label):
    v = np.asarray(da.values)
    if not np.isfinite(v).all():
        raise ValueError(f"{label}: {int(np.isnan(v).sum())} NaN, {int(np.isinf(v).sum())} inf")


def dart_mean(cfg: dict, st: Settings) -> xr.Dataset:
    """Mean of the 20 selected members, on the reference layout."""
    out = {}
    with xr.open_dataset(cfg["ref_atm_file"], decode_times=False, cache=False) as ref_ds:
        for v in st.variables:
            acc = np.zeros(ref_ds[v].shape)
            for m in cfg["selected"]:
                da = _load(cfg["dart_files"][m], v, ref_ds, st)
                _finite(da, f"EN{m:02d} {v}")
                acc += da.values
            out[v] = xr.DataArray(acc / len(cfg["selected"]), dims=ref_ds[v].dims,
                                  coords={d: ref_ds[v].coords[d] for d in ref_ds[v].dims if d in ref_ds[v].coords},
                                  attrs=dict(ref_ds[v].attrs), name=v)
    return xr.Dataset(out)


def perturbation_diagnostics(cfg: dict, mean: xr.Dataset, st: Settings) -> pd.DataFrame:
    """Per variable: spread of the member-mean perturbation and the member RMS perturbation."""
    rows = []
    with xr.open_dataset(cfg["ref_atm_file"], decode_times=False, cache=False) as ref_ds:
        for v in st.variables:
            means, rms = [], []
            for m in cfg["selected"]:
                d = (_load(cfg["dart_files"][m], v, ref_ds, st) - mean[v]).values
                means.append(float(d.mean()))
                rms.append(float(np.sqrt((d ** 2).mean())))
            rows.append(dict(ref_name=cfg["ref"].name, date=cfg["date"], variable=v,
                             mean_of_member_mean_perturbations=float(np.mean(means)),
                             std_of_member_mean_perturbations=float(np.std(means, ddof=1)),
                             mean_member_rms_perturbation=float(np.mean(rms)),
                             min_member_rms_perturbation=float(np.min(rms)),
                             max_member_rms_perturbation=float(np.max(rms))))
    return pd.DataFrame(rows)


def output_file(cfg: dict, member: int) -> Path:
    r = cfg["ref"]
    return cfg["output_dir"] / f"{r.out_prefix}.EN{member:02d}.{r.out_stream}.{cfg['date']}.nc"


def _same_attr(a, b):
    try:
        return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)
    except TypeError:
        return np.array_equal(np.asarray(a), np.asarray(b))


def check_same_structure(ref_file: Path, out_file: Path):
    """Dimensions, variables, dtypes, chunking, filters and attributes identical to the reference."""
    errors = []
    with Dataset(ref_file) as a, Dataset(out_file) as b:
        if a.data_model != b.data_model:
            errors.append(f"data model {a.data_model} -> {b.data_model}")
        if list(a.dimensions) != list(b.dimensions):
            errors.append("dimension names/order changed")
        for d in a.dimensions:
            if d in b.dimensions and (len(a.dimensions[d]) != len(b.dimensions[d])
                                      or a.dimensions[d].isunlimited() != b.dimensions[d].isunlimited()):
                errors.append(f"dimension {d} changed")
        if list(a.variables) != list(b.variables):
            errors.append("variable names/order changed")
        if list(a.ncattrs()) != list(b.ncattrs()):
            errors.append("global attribute names/order changed")
        errors += [f"global attribute {k!r} changed" for k in a.ncattrs()
                   if k in b.ncattrs() and not _same_attr(a.getncattr(k), b.getncattr(k))]
        for name, va in a.variables.items():
            if name not in b.variables:
                continue
            vb = b.variables[name]
            for prop, f in (("dtype", lambda x: str(x.dtype)), ("dimensions", lambda x: tuple(x.dimensions)),
                            ("shape", lambda x: tuple(x.shape)), ("filters", lambda x: x.filters()),
                            ("chunking", lambda x: x.chunking()), ("endian", lambda x: x.endian()),
                            ("attributes", lambda x: list(x.ncattrs()))):
                if f(va) != f(vb):
                    errors.append(f"{name} {prop} changed")
            errors += [f"{name} attribute {k!r} changed" for k in va.ncattrs()
                       if k in vb.ncattrs() and not _same_attr(va.getncattr(k), vb.getncattr(k))]
    if errors:
        raise RuntimeError(f"Structure check failed for {out_file}:\n" + "\n".join(errors))


def check_unperturbed_identical(ref_file: Path, out_file: Path, perturbed):
    with Dataset(ref_file) as a, Dataset(out_file) as b:
        for name, va in a.variables.items():
            if name in perturbed:
                continue
            vb = b.variables[name]
            va.set_auto_maskandscale(False)
            vb.set_auto_maskandscale(False)
            if not np.array_equal(va[:], vb[:]):
                raise RuntimeError(f"Unperturbed variable {name} changed in {out_file}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _exists_policy(out: Path, st: Settings) -> bool:
    """True when an existing output is to be reused (validated) instead of rewritten."""
    if not out.exists() or st.overwrite:
        return False
    if not st.resume:
        raise FileExistsError(f"{out} exists; set resume=True to reuse or overwrite=True to rewrite")
    return True


def write_control(cfg: dict, st: Settings) -> tuple[Path, str]:
    """EN01: an exact copy of the reference (checksum verified)."""
    ref, out = cfg["ref_atm_file"], output_file(cfg, 1)
    ref_sum = sha256(ref)
    if _exists_policy(out, st):
        if sha256(out) != ref_sum:
            raise RuntimeError(f"Existing EN01 differs from the reference: {out}")
        print(f"[REUSE] {out.name} (checksum verified)")
        return out, ref_sum
    tmp = out.with_suffix(out.suffix + ".tmp")
    tmp.unlink(missing_ok=True)
    try:
        shutil.copy2(ref, tmp)
        if sha256(tmp) != ref_sum:
            raise RuntimeError("EN01 checksum differs from the reference")
        tmp.replace(out)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    return out, ref_sum


def write_member(cfg: dict, out_member: int, dart_member: int, mean: xr.Dataset, st: Settings):
    """A perturbed member: the reference file with ``A_ref + alpha * (A_DART - mean)`` written in."""
    ref, out = cfg["ref_atm_file"], output_file(cfg, out_member)
    if _exists_policy(out, st):
        check_same_structure(ref, out)
        check_unperturbed_identical(ref, out, set(st.variables))
        prior = cfg["output_dir"] / "q_clipping_diagnostics.csv"
        rows = []
        if prior.is_file() and prior.stat().st_size:
            q = pd.read_csv(prior)
            if "output_member" in q:
                rows = q.loc[q["output_member"] == out_member].to_dict(orient="records")
        print(f"[REUSE] {out.name} (structure verified)")
        return out, rows

    clip_rows, values = [], {}
    with xr.open_dataset(ref, decode_times=False, cache=False) as ref_ds:
        for v in st.variables:
            ref_da = ref_ds[v].load().astype("float64")
            pert = _load(cfg["dart_files"][dart_member], v, ref_ds, st) - mean[v]
            _finite(pert, f"DART EN{dart_member:02d} {v} perturbation")
            new = ref_da + st.alpha * pert
            if st.clip_q and v in st.clipped_vars:
                mask = new < st.q_min
                frac = float(mask.mean())
                vmin = float(new.min())
                clip_rows.append(dict(ref_name=cfg["ref"].name, date=cfg["date"], variable=v,
                                      output_member=out_member, output_member_label=f"EN{out_member:02d}",
                                      dart_member=dart_member, dart_member_label=f"EN{dart_member:02d}",
                                      number_q_clipped=int(mask.sum()), fraction_q_clipped=frac,
                                      minimum_q_before_clip=vmin, maximum_q_correction=max(0.0, st.q_min - vmin),
                                      exceeds_q_clip_warning=frac > st.warn_q_clip_fraction))
                if frac > st.warn_q_clip_fraction:
                    print(f"[WARN] EN{out_member:02d} {v} clipped fraction {frac:.3e}")
                if st.enforce_max_q_clip_fraction and frac > st.max_q_clip_fraction:
                    raise ValueError(f"EN{out_member:02d} {v} clipped fraction {frac:.3e} > {st.max_q_clip_fraction}")
                new = new.clip(min=st.q_min)
            if new.dims != ref_da.dims or new.shape != ref_da.shape:
                raise ValueError(f"EN{out_member:02d} {v} changed structure")
            if (np.isfinite(ref_da.values) & ~np.isfinite(new.values)).any():
                raise ValueError(f"EN{out_member:02d} {v} has new non-finite values")
            lo, hi = st.physical_limits[v]
            if not (lo <= float(new.min()) and float(new.max()) <= hi):
                raise ValueError(f"EN{out_member:02d} {v} outside ({lo}, {hi}): "
                                 f"{float(new.min())}, {float(new.max())}")
            values[v] = np.asarray(new.values)

    tmp = out.with_suffix(out.suffix + ".tmp")
    tmp.unlink(missing_ok=True)
    try:
        shutil.copy2(ref, tmp)
        with Dataset(tmp, mode="r+") as nc:
            for v, arr in values.items():
                if nc.variables[v].shape != arr.shape:
                    raise ValueError(f"{v}: file shape {nc.variables[v].shape} vs {arr.shape}")
                nc.variables[v][:] = arr
        check_same_structure(ref, tmp)
        check_unperturbed_identical(ref, tmp, set(st.variables))
        tmp.replace(out)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    return out, clip_rows


def validate_written_files(cfg: dict, mean: xr.Dataset, st: Settings) -> pd.DataFrame:
    """Per variable: member deltas from the reference and the closure of the ensemble mean."""
    rows = []
    with xr.open_dataset(cfg["ref_atm_file"], decode_times=False, cache=False) as ref_ds:
        for v in st.variables:
            ref_v = ref_ds[v].load().astype("float64")
            means, rms, mins, maxs, total = [], [], [], [], None
            for m in st.perturbed_members:
                with xr.open_dataset(output_file(cfg, m), decode_times=False, cache=False) as ds:
                    out_v = ds[v].load().astype("float64")
                if (np.isfinite(ref_v.values) & ~np.isfinite(out_v.values)).any():
                    raise ValueError(f"written EN{m:02d} {v} has new non-finite values")
                d = out_v - ref_v
                means.append(float(d.mean())); rms.append(float(np.sqrt((d ** 2).mean())))
                mins.append(float(d.min())); maxs.append(float(d.max()))
                total = d if total is None else total + d
            pert_mean = (total / len(st.perturbed_members)).values
            full_mean = (total / len(st.all_members)).values
            expected = (-st.alpha * (_load(cfg["dart_files"][cfg["omitted"]], v, ref_ds, st) - mean[v])
                        / len(st.all_members)).values
            err = full_mean - expected
            rms_ = lambda a: float(np.sqrt(np.mean(a ** 2)))
            rows.append(dict(ref_name=cfg["ref"].name, date=cfg["date"], variable=v,
                             mean_member_mean_delta=float(np.mean(means)), mean_member_rms_delta=float(np.mean(rms)),
                             min_delta_across_members=min(mins), max_delta_across_members=max(maxs),
                             rms_perturbed_mean_delta=rms_(pert_mean),
                             max_abs_perturbed_mean_delta=float(np.abs(pert_mean).max()),
                             rms_complete_mean_delta=rms_(full_mean),
                             max_abs_complete_mean_delta=float(np.abs(full_mean).max()),
                             rms_expected_complete_mean_delta=rms_(expected),
                             max_abs_expected_complete_mean_delta=float(np.abs(expected).max()),
                             rms_mean_closure_error=rms_(err), max_abs_mean_closure_error=float(np.abs(err).max())))
    out = pd.DataFrame(rows)
    out["ensemble_mean_tolerance"] = out["variable"].map(st.mean_tolerance)
    # clipping breaks the closure for the clipped variables, so their tolerance is not enforced
    out["ensemble_mean_tolerance_enforced"] = ~(out["variable"].isin(st.clipped_vars) & st.clip_q)
    out["passes_tolerance"] = (~out["ensemble_mean_tolerance_enforced"]
                               | (out["max_abs_mean_closure_error"] <= out["ensemble_mean_tolerance"]))
    return out


def link_components(cfg: dict, st: Settings) -> list[Path]:
    """Link (or copy) the reference's non-atmosphere restarts and rpointers into the date directory."""
    if not st.link_components:
        return []
    srcs = []
    for pat in st.component_patterns:
        found = list(cfg["ref_rest_dir"].glob(pat))
        if not found and pat not in st.optional_component_patterns:
            raise FileNotFoundError(f"No {pat!r} in {cfg['ref_rest_dir']}")
        srcs += found
    linked = []
    for src in sorted(set(srcs)):
        out = cfg["output_dir"] / src.name
        if out.exists() or out.is_symlink():
            if (st.resume and not st.overwrite and st.link_mode == "symlink" and out.is_symlink()
                    and out.resolve() == src.resolve()):
                linked.append(out)
                continue
            if not st.overwrite:
                raise FileExistsError(f"Component output exists: {out}")
            out.unlink()
        if st.link_mode == "symlink":
            out.symlink_to(src)
        elif st.link_mode == "copy":
            shutil.copy2(src, out)
        else:
            raise ValueError("link_mode must be 'symlink' or 'copy'")
        linked.append(out)
    return linked


def write_metadata(cfg: dict, st: Settings, written, linked, control_sha) -> Path:
    r, s = cfg["ref"], cfg["season"]
    meta = dict(
        ref_name=r.name, date_tag=cfg["date"], set_name=s.name, dart_case_prefix=s.case,
        dart_rest_dir=str(s.rest_root / cfg["date"]), dart_source_type=cfg["dart_source_type"],
        ref_atm_file=str(cfg["ref_atm_file"]), ref_rest_dir=str(cfg["ref_rest_dir"]),
        method="selected-subset-centered DART perturbations",
        selected_dart_member_csv=str(s.selection_csv) if s.selection_csv else None,
        control_output_member="EN01", control_sha256=control_sha,
        selected_dart_members=cfg["selected"], centering_members=cfg["selected"],
        perturbation_source_members=cfg["sources"], omitted_source_member=cfg["omitted"],
        output_to_dart_member={f"EN{o:02d}": f"EN{d:02d}" for o, d in cfg["output_to_dart"].items()},
        compatible_vars=list(st.variables), output_dir=str(cfg["output_dir"].resolve()),
        out_case_prefix=r.out_prefix, out_atm_stream=r.out_stream,
        settings={k: (list(v) if isinstance(v, tuple) else v) for k, v in asdict(st).items()},
        created_atmosphere_files=[str(p) for p in written], created_component_links=[str(p) for p in linked],
        formula="A_output = A_REF + alpha*(A_DART_selected - mean(A_DART_selected_20_member_subset))",
        expected_complete_mean_offset="-alpha*(A_DART_omitted - mean(A_DART_selected_20_member_subset))/20",
    )
    path = cfg["output_dir"] / "dart_perturbation_metadata.json"
    path.write_text(json.dumps(meta, indent=2))
    return path


# --------------------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------------------

def build_date(cfg: dict, st: Settings, display=print) -> dict:
    """All steps for one start date; returns its diagnostics, validation and clipping tables."""
    print(f"\n{'=' * 80}\n[DATE] {cfg['ref'].name} {cfg['season'].name} {cfg['date']} ({cfg['dart_source_type']})")
    cfg["output_dir"].mkdir(parents=True, exist_ok=True)
    preflight(cfg, st)
    cfg["mapping"].to_csv(cfg["output_dir"] / "output_to_dart_member_mapping.csv", index=False)
    mean = dart_mean(cfg, st)
    diag = perturbation_diagnostics(cfg, mean, st)
    diag.to_csv(cfg["output_dir"] / "perturbation_diagnostics.csv", index=False)
    written, clip, valid, control_sha = [], pd.DataFrame(), pd.DataFrame(), None
    if st.write_output_files:
        f, control_sha = write_control(cfg, st)
        written.append(f)
        rows = []
        for o, d in cfg["output_to_dart"].items():
            print(f"[WRITE] EN{o:02d} <- DART EN{d:02d}")
            f, r = write_member(cfg, o, d, mean, st)
            written.append(f)
            rows += r
        clip = pd.DataFrame(rows)
        clip.to_csv(cfg["output_dir"] / "q_clipping_diagnostics.csv", index=False)
        valid = validate_written_files(cfg, mean, st)
        valid.to_csv(cfg["output_dir"] / "written_file_validation.csv", index=False)
        failed = valid.loc[~valid["passes_tolerance"]]
        if not failed.empty:
            display(failed)
            raise RuntimeError(f"{cfg['date']}: ensemble-mean closure outside tolerance")
    linked = link_components(cfg, st)
    meta = write_metadata(cfg, st, written, linked, control_sha)
    print(f"[DONE] {len(written)} atmosphere files, {len(linked)} component links, {meta.name}")
    return dict(diagnostics=diag, validation=valid, q_clipping=clip)


def build_library(configs: list[dict], st: Settings, display=print) -> dict:
    """Run ``build_date`` for every configuration; write the ``all_*.csv`` summaries per output root."""
    results = [build_date(c, st, display) for c in configs]
    tables = {k: pd.concat([r[k] for r in results if not r[k].empty], ignore_index=True)
              if any(not r[k].empty for r in results) else pd.DataFrame()
              for k in ("diagnostics", "validation", "q_clipping")}
    names = {"diagnostics": "all_perturbation_diagnostics.csv", "validation": "all_written_file_validation.csv",
             "q_clipping": "all_q_clipping_diagnostics.csv"}
    for root in {Path(c["ref"].output_root) for c in configs}:
        for k, t in tables.items():
            if not t.empty:
                t.to_csv(root / names[k], index=False)
    return tables


__all__ = ["DART_SEASONS", "DartSeason", "dart_date_sets", "ReferenceCase", "Settings", "member_selection",
           "build_run_configs", "config_table", "preflight", "dart_mean", "build_date", "build_library",
           "output_file", "validate_written_files"]
