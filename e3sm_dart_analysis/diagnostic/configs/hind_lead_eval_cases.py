"""Case settings of the hind_lead_eval case notebooks: 05 (initialization benefits), 06 (summer 2012
heat wave), 07 (January 2012 Pacific Northwest storm) and 08 (tracer initialization).

Every notebook of a series builds its settings from these defaults plus its parameter cell (paths, run
options and ``CASE_OVERRIDES``), so a definition changed here changes it in every notebook of the series.
Cached products record the settings they depend on (``stamp``) and are rebuilt when those change; paths
are not part of the stamps (moving an archive does not change the data).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class RunPaths:
    """Machine-specific locations (notebook parameters, see ``jupyter/README.md``)."""
    input_data_root: Path
    work_dir: Path
    diagnostic_output_root: Path
    observation_root: Path
    land_support_dir: Path
    amip_root: Path
    workflow_name: str = "hind_lead_eval"


class _Settings:
    def stamp(self, *names: str) -> str:
        """JSON of the named settings: the identity of a product that depends on them."""
        return json.dumps({n: getattr(self, n) for n in names}, sort_keys=True, default=str)


# Notebook series 05: the 2 x 2 designs (roles: first digit atmosphere, second land; 1 = DART), one per
# question of CONFIG_GROUPS (configs/fcst_init_experiment_config.py). "reference": the E00 initial states the
# DART benefits are measured against.
INIT_BENEFIT_DESIGNS = {
    "dart_vs_capt": dict(   # (c) DART vs a CAPT-like initialization
        roles={"E00": "ERA5CLIMFILL", "E10": "DA_ATM", "E01": "DA_LND_CLIMFILL", "E11": "DA_COUPLED"},
        reference="CAPT-like (ERA5 atmosphere + DATM-forced land)"),
    "lnd_vs_atm": dict(     # (b) land vs atmosphere initialization against the free-running climate states
        roles={"E00": "CLIM_BASE", "E10": "CLIMO_LND", "E01": "CLIMO_ATM", "E11": "DA_COUPLED"},
        reference="free-running climate (CLIM_BASE atmosphere + land)"),
}


@dataclass(frozen=True)
class InitBenefitSettings(_Settings):
    """Notebook series 05: DART initialization benefits for subseasonal 2-m temperature.

    A 2 x 2 of the atmosphere and land initial states (first digit atmosphere, second land; 1 = DART), one of
    ``INIT_BENEFIT_DESIGNS`` (``design``; build with ``for_design``):

    * ``dart_vs_capt`` (default, question c): E00 ERA5CLIMFILL (CAPT-like reference: ERA5 atmosphere +
      DATM-forced land), E11 DA_COUPLED (coupled DART analyses), and the two mixed, land-atmosphere-inconsistent
      experiments E10 DA_ATM (DART atmosphere + DATM-forced land) and E01 DA_LND_CLIMFILL (ERA5 atmosphere +
      DART land). CLIM_BASE (no assimilation) is outside it.
    * ``lnd_vs_atm`` (question b): E00 CLIM_BASE (free-running climate states), E10 CLIMO_LND (DART atmosphere
      + CLIM_BASE land), E01 CLIMO_ATM (CLIM_BASE atmosphere + DART perturbations, DART land), E11 DA_COUPLED.

    The scores, maps and figures of a design go to its own subdirectory (``<design>/``); the daily member
    caches are per configuration and shared.
    """
    paths: RunPaths
    design: str = "dart_vs_capt"            # name of the product set (subdirectory); roles set the experiments
    roles: dict = field(default_factory=lambda: dict(INIT_BENEFIT_DESIGNS["dart_vs_capt"]["roles"]))
    # name: (roles added, roles subtracted). Centred on the full (coupled) initialization:
    # atmosphere S11 - S01, land S11 - S10, interaction (land-atmosphere consistency) S11 - S10 - S01 + S00,
    # total S11 - S00 = atmosphere + land - interaction. Without E01 (e.g. roles["E01"] = None) use the
    # sequential set {"atmosphere": (("E10",), ("E00",)), "land": (("E11",), ("E10",)), "total": (("E11",), ("E00",))}.
    contrasts: dict = field(default_factory=lambda: {"atmosphere": (("E11",), ("E01",)),
                                                     "land": (("E11",), ("E10",)),
                                                     "interaction": (("E11", "E00"), ("E10", "E01")),
                                                     "total": (("E11",), ("E00",))})
    max_members: int | None = None          # e.g. 3 for a quick test
    require_all_members: bool = True        # use a season only when every role's members are post-processed
    variable: str = "TREFHT"
    region: dict = field(default_factory=lambda: dict(name="NH_MIDLAT_LAND", lat=(30.0, 60.0), lon=(0.0, 360.0),
                                                      mask="land"))
    seasons: tuple = ("DJF", "JJA")
    max_lead: int = 56                      # lead days (8 weeks)
    clim_years: tuple = (2001, 2014)        # as notebook 02
    amip_cases: tuple = ("v3.LR.amip_0101", "v3.LR.amip_0151")
    biases: tuple = ("amip", "none")        # forecast anomaly vs the AMIP or the observed climatology
    model_truth_role: str = "E11"           # panel (c): held-out members of this role are the truth
    n_boot: int = 1000
    seed: int = 42
    smooth_days: int = 7                    # running mean over lead days for the daily curves
    window_weeks: tuple = (2, 4)            # shaded subseasonal window
    map_weeks: tuple = (3, 4)               # weeks averaged in the gridpoint-contrast maps

    @classmethod
    def for_design(cls, design: str, paths: RunPaths, **overrides) -> "InitBenefitSettings":
        """Settings of a design of ``INIT_BENEFIT_DESIGNS`` (its roles), with ``overrides`` applied on top."""
        if design not in INIT_BENEFIT_DESIGNS:
            raise ValueError(f"No 05-series design {design!r}; choose from {list(INIT_BENEFIT_DESIGNS)} "
                             f"(CONFIG_GROUP of notebooks 05a-05c)")
        return cls(paths=paths, **{"design": design, "roles": dict(INIT_BENEFIT_DESIGNS[design]["roles"]),
                                   **overrides})

    @property
    def reference(self) -> str:
        """The E00 initial states the benefits are measured against (``INIT_BENEFIT_DESIGNS``)."""
        d = INIT_BENEFIT_DESIGNS.get(self.design)
        return d["reference"] if d and d["roles"] == self.roles else f"E00 {self.roles.get('E00')}"


@dataclass(frozen=True)
class HeatwaveSettings(_Settings):
    """Notebook series 06: the summer 2012 central U.S. heat wave (May 2012 starts)."""
    paths: RunPaths
    configs: tuple = ("CLIM_BASE", "ERA5CLIMFILL", "DA_COUPLED", "DA_ATM", "DA_LND_CLIMFILL")  # Q1/Q2, then Q3
    primary: tuple = ("CLIM_BASE", "ERA5CLIMFILL", "DA_COUPLED")                               # Q1 and Q2
    ic_defect: tuple = ()                       # configs started without radiative tracers: QC only (none now)
    season: str = "JJA"
    n_days: int = 61                            # simulated days used (lead 1 = start date)
    max_members: int | None = None              # e.g. 4 for a quick test
    require_all_members: bool = True
    map_box: dict = field(default_factory=lambda: dict(lat=(24.5, 50.5), lon=(234.5, 294.5)))   # CONUS
    region: dict = field(default_factory=lambda: dict(name="CENTRAL_US", lat=(35.0, 45.0), lon=(255.0, 275.0),
                                                      mask="land"))                         # 35-45N, 105-85W
    obs_period: tuple = ("2012-05-16", "2012-07-31")   # covers every start's simulated days
    clim_years: tuple = (2001, 2014)            # climatologies and threshold
    exclude_years: tuple = (2012,)              # left out of clim_years: the event year
    amip_cases: tuple = ("v3.LR.amip_0101", "v3.LR.amip_0151")
    percentile: float = 0.9
    half_window: int = 7                        # days on either side for the calendar-day percentile
    min_event_days: int = 3
    event_window: tuple | None = None           # e.g. ("2012-06-27", "2012-07-07"); None: largest event from June on
    peak_tolerance_days: int = 7
    biases: tuple = ("none", "amip")
    spread: tuple = (0.1, 0.9)                  # member quantiles shaded in the series
    lead_windows: dict = field(default_factory=lambda: {"week 1": (1, 7), "week 2": (8, 14), "weeks 3-4": (15, 28),
                                                        "weeks 5-6": (29, 42), "weeks 7-8": (43, 56)})
    contrasts: tuple = (("DA_COUPLED", "CLIM_BASE"), ("DA_ATM", "CLIM_BASE"), ("DA_COUPLED", "DA_ATM"),
                        ("ERA5CLIMFILL", "CLIM_BASE"), ("DA_COUPLED", "ERA5CLIMFILL"),
                        ("DA_LND_CLIMFILL", "CLIM_BASE"))  # a minus b
    n_boot: int = 2000
    sens_threshold_years: tuple = ((2001, 2014), (1991, 2020))   # threshold periods tested (exclude_years left out)
    sens_percentiles: tuple = (0.85, 0.9, 0.95)
    workers: int = 4                            # processes reading model output (cache builds only)


@dataclass(frozen=True)
class WinterStormSettings(_Settings):
    """Notebook series 07: the January 2012 Pacific Northwest storm (Dec 2011 starts)."""
    paths: RunPaths
    configs: tuple = ("CLIM_BASE", "ERA5CLIMFILL", "DA_COUPLED", "DA_ATM", "DA_LND_CLIMFILL")
    primary: tuple = ("CLIM_BASE", "ERA5CLIMFILL", "DA_COUPLED")
    ic_defect: tuple = ()                       # configs started without radiative tracers (none now)
    season: str = "DJF"                         # the Dec 2011 starts
    n_days: int = 61
    max_members: int | None = None
    require_all_members: bool = True
    event_window: tuple = ("2012-01-16", "2012-01-20")
    series_xlim: tuple = ("2012-01-06", "2012-01-27")
    obs_period: tuple = ("2011-12-16", "2012-02-29")    # covers every start's simulated days
    map_box: dict = field(default_factory=lambda: dict(lat=(39.5, 52.5), lon=(227.5, 250.5)))   # Pacific Northwest
    env_box: dict = field(default_factory=lambda: dict(lat=(20.5, 65.5), lon=(180.5, 255.5)))   # N Pacific to Rockies
    map_extent: tuple = (-130.0, -110.0, 40.0, 52.0)
    env_extent: tuple = (-178.0, -106.0, 22.0, 64.0)
    regions: tuple = (dict(name="WEST_WA_OR", lat=(42.0, 49.0), lon=(235.5, 238.5), mask="land"),
                      dict(name="INLAND_PNW", lat=(42.0, 49.0), lon=(238.5, 249.0), mask="land"))
    peak_region: str = "WEST_WA_OR"             # region whose observed peak sets the timing reference
    env_days: tuple = ("2012-01-14", "2012-01-16", "2012-01-18", "2012-01-20")
    clim_years: tuple = (2001, 2014)
    amip_cases: tuple = ("v3.LR.amip_0101", "v3.LR.amip_0151")
    peak_tolerance_days: int = 3
    biases: tuple = ("none", "amip")
    spread: tuple = (0.1, 0.9)
    workers: int = 4


# Notebook 08: the variables of each step of the chain from the missing tracers to the forecast (variable, region
# shown in the summary figures). Every variable is read whole; the regions only select what the figures show.
TRACER_CHAIN = {
    "tracers": ("TOZ", "SCO", "TCO"),                               # total, stratospheric, tropospheric ozone (DU)
    "radiation": ("FLDSC", "FSDSC", "FLUTC", "FSNTOAC"),            # clear-sky fluxes
    "stratosphere": ("T010", "T050", "T100"),
    "troposphere": ("T850", "Z500", "TMQ", "PRECT"),
    "surface": ("TREFHT", "TS"),
}


@dataclass(frozen=True)
class TracerInitSettings(_Settings):
    """Notebook 08: the radiative/chemical tracer initialization (question a, ``CONFIG_GROUPS["tracer_init"]``).

    ``test`` (HYB_BASE: the ERA5 eam.i without ozone, linoz gases, aerosols and chemistry) against ``control``
    (ERA5CLIMFILL: the same ERA5 state, DART perturbations and DATM-forced land, the missing fields taken from
    CLIM_BASE). ``contrasts`` (name: (a, b), a minus b) put the tracer defect next to the initialization signals it
    would contaminate. Members are exchangeable and not paired across configurations (``METHODS.md``).
    """
    paths: RunPaths
    configs: tuple = ("CLIM_BASE", "HYB_BASE", "ERA5CLIMFILL", "DA_ATM")
    test: str = "HYB_BASE"
    control: str = "ERA5CLIMFILL"
    contrasts: dict = field(default_factory=lambda: {
        "tracer_defect": ("HYB_BASE", "ERA5CLIMFILL"),     # incomplete - completed ERA5 state: the defect
        "dart_atm": ("DA_ATM", "ERA5CLIMFILL"),            # DART - ERA5 atmosphere (question c, same land)
        "era5_vs_clim": ("ERA5CLIMFILL", "CLIM_BASE"),     # ERA5 atmosphere + DATM land - free-running states
    })
    variables: tuple = tuple(v for vs in TRACER_CHAIN.values() for v in vs)
    seasons: tuple = ("DJF", "JJA")
    max_lead: int = 56
    max_members: int | None = None              # e.g. 3 for a quick test
    require_all_members: bool = True
    regions: tuple = ("GLOBAL", "LAND", "TROPICS", "NH_EXTRATROPICS", "SH_EXTRATROPICS", "NH_MIDLAT_LAND")
    custom_regions: dict = field(default_factory=lambda: {
        "TROPICS": dict(lat=(-30.0, 30.0), lon=(0.0, 360.0)),
        "NH_EXTRATROPICS": dict(lat=(30.0, 90.0), lon=(0.0, 360.0)),
        "SH_EXTRATROPICS": dict(lat=(-90.0, -30.0), lon=(0.0, 360.0)),
        "NH_MIDLAT_LAND": dict(lat=(30.0, 60.0), lon=(0.0, 360.0), mask="land"),   # as notebooks 05
    })
    lead_windows: dict = field(default_factory=lambda: {"days 1-3": (1, 3), "days 4-10": (4, 10),
                                                        "days 11-28": (11, 28), "days 29-56": (29, 56)})
    recovery_k: float = 1.0                     # recovered once |difference| stays within k x the member spread
    z_crit: float = 2.0                         # |difference| / standard error counted as significant
    # (variable, region) of the summary (normalized recovery) figure: one per step of TRACER_CHAIN
    chain: tuple = (("TOZ", "GLOBAL"), ("FLDSC", "GLOBAL"), ("T050", "GLOBAL"), ("TREFHT", "LAND"),
                    ("T850", "NH_EXTRATROPICS"), ("Z500", "NH_EXTRATROPICS"))
    workers: int = 4                            # processes reading model output (summary builds only)


__all__ = ["RunPaths", "INIT_BENEFIT_DESIGNS", "InitBenefitSettings", "HeatwaveSettings", "WinterStormSettings",
           "TRACER_CHAIN", "TracerInitSettings"]
