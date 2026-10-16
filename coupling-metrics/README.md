# coupling-metrics

Fortran 90 routines for land–atmosphere coupling metrics: the Coupling Metrics
Toolkit (CoMeT) by Ahmed B. Tawfik. Use them to diagnose land–atmosphere
coupling in E3SM output and observations, from NCL, Python (f2py), or a Fortran
program.

## Source

This is a copy of an external repository, not code written here.

| | |
| --- | --- |
| Upstream | <https://github.com/abtawfik/coupling-metrics> |
| Branch | `master` |
| Commit | `73f70d2` (2020-10-20, "fix(hcf): incorrectly specifying the edge case where the BCLP is greater than the surface pressure") |
| License | GNU GPL v2 (`LICENSE`); keep it with the code. |

The upstream Fortran sources, stubs, sample data, and `LICENSE` are unchanged.
Local changes when importing:

- `README.md` (this file) was added. The original upstream README is kept as
  `UPSTREAM_README.md`.
- The upstream `.gitignore`, a generic Python template that also ignored
  `*.csv` and so would have dropped `examples/Sample_Obs/Sample_tcp.csv`, was
  replaced with one that ignores only Fortran build products.
- The upstream `.git` directory, an empty `ReadMe.txt`, and
  `.ipynb_checkpoints/` were removed.

The upstream repository also has the branches `new_version_1`,
`python_version`, and `python_mixing_diagram`. They were not imported. To update
this copy, clone upstream, copy the metric directories over this one, and
update the commit above.

## Metrics

| Directory | Metric | Routine (`static/` module) | Reference |
| --- | --- | --- | --- |
| `convect_trigger_pot/` | Convective Triggering Potential (CTP) and low-level humidity index (HI-low): whether a morning profile favors convection over dry or wet soils. | `ctp_hi_low` (`Conv_Trig_Pot_Mod`) | Findell and Eltahir (2003), J. Hydrometeor. |
| `heated_condensation/` | Heated Condensation Framework (HCF): buoyant condensation level, potential-temperature and energy deficits for local versus non-local convective initiation. | `hcfcalc` (`HCF_vars_calc`) | Tawfik and Dirmeyer (2014), GRL |
| `mixing_diagram/` | Mixing diagrams: surface and entrainment heat and moisture fluxes into the PBL from the diurnal covariation of 2-m temperature and humidity. | `mixing_diag` (`Mixing_Diag_Mod`) | Santanello et al. (2009, 2011) |
| `rh_tendency/` | RH tendency at the top of the boundary layer, split into surface-flux, heating, PBL-growth, and dry-air entrainment terms. | `rh_tend_calc` (`RH_Tend_Mod`) | Ek and Mahrt (1994); Ek and Holtslag (2004) |
| `soil_moisture_memory/` | Soil-moisture memory: timescale from the lagged autocorrelation of soil moisture. | `soilm_memory` (`Soil_Memory_Mod`) | none given in source |
| `terrestrial_coupling_index/` | Terrestrial coupling index: standard deviation of soil moisture times the regression slope of a surface flux on soil moisture, in flux units. | `terra_coupling` (`Terrestrial_Coupling_Mod`); `terra_wrapper` for `(nlon,nlat,ntim)` arrays (`ncl_portable/` only) | Dirmeyer (2011) |

The header of each source file documents argument units and dimensions. All
routines take single-precision (`real(4)`) arrays and a `missing` value, which
they use to mark missing inputs and undefined outputs.

`mixing_diagram/static/` has three alternative files that all define
`Mixing_Diag_Mod`: `mixing_diagram.f90` (`mixing_diag`, the default),
`mixing_diagram_with_advect.f90` (`mixing_diag_advect`, adds advection terms),
and `center_difference.f90` (`mixing_diag_step`). Compile only one of them into
a given program.

## Layout

Each metric directory has two versions of the same calculation:

- `static/`: Fortran 90 modules to `use` from a Fortran program.
- `ncl_portable/`: plain subroutines and `.stub` files for building shared
  libraries with NCL `WRAPIT` or with `f2py`.

`examples/Sample_Obs/` runs the routines on observed profiles and fluxes from
the ARM Southern Great Plains Central Facility (IOP, 6 June 2002).
`Program_run_sample.f90` is written for these input files only; it is not a
general driver.

## Usage

### Fortran sample program

```bash
cd examples/Sample_Obs
bash compile_and_run.sh
```

The script compiles the `static/` modules and `Program_run_sample.f90` with
`gfortran` (debug flags) and runs `sample_data_test.exe`, which prints each
metric. Build products are written into `examples/Sample_Obs/` and ignored by
Git. The script starts with `rm *.o *.mod *.exe`, so a first run prints a
harmless "No such file" message.

On Perlmutter, this works with the default `gfortran`
(`/opt/cray/pe/gcc-native/14/bin/gfortran`) and no extra modules (tested
2026-10-16 at upstream commit `73f70d2`).

### NCL

```bash
cd convect_trigger_pot/ncl_portable
WRAPIT ctp_hilow.stub ctp_hilow.f90     # writes ctp_hilow.so
```

```ncl
external CTPHiLow "./ctp_hilow.so"
CTP   = new(1, float)
HILOW = new(1, float)
CTPHiLow::ctp_hi_low(nlev, t, q, p, t2m, q2m, psfc, CTP, HILOW, missing)
```

Allocate the output variables with `new` before the call. The other metrics
work the same way with their own `.stub` and `.f90` pair.

### Python

```bash
cd soil_moisture_memory/ncl_portable
f2py -c -m soil_moisture_memory soil_moisture_memory.f90
```

This builds an importable extension module. Pass `float32` arrays in Fortran
order to match the `real(4)` arguments. The f2py build was not tested here:
it needs a Python environment with NumPy.
