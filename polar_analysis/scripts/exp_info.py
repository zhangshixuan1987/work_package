from dataclasses import dataclass
from typing import Dict, Tuple

# ---------------- Dataclasses ----------------

@dataclass(frozen=True)
class RunInfo:
    """Metadata for a simulation or observation run."""
    name: str
    period: str


@dataclass(frozen=True)
class Region:
    """Geographical region defined by latitude and longitude bounds."""
    lat_range: Tuple[float, float]
    lon_range: Tuple[float, float]


@dataclass(frozen=True)
class VariableInfo:
    """Metadata for a variable."""
    unit: str
    fac: float
    min: float
    max: float
    nlev: int
    ref: str  # key from RUN_CATALOG

@dataclass(frozen=True)
class BudgetInfo:
    """Metadata for a variable."""
    unit: str
    fac: float
    min: float
    max: float
    nlev: int
    alias: str # alias name
    ref: str  # key from RUN_CATALOG


# ---------------- Run Catalog ----------------

RUN_CATALOG: Dict[str, RunInfo] = {
    "HadISST.analysis":       RunInfo("HadISST",          "2001–2018"),
    "MODIS.observation":      RunInfo("MODIS",            "2001–2018"),
    "CERES_EBAF.observation": RunInfo("CERES_EBAF",       "2001–2018"),
    "MERRA2.analysis":        RunInfo("MERRA2 Analysis",  "2001–2018"),
    "ERA5.analysis":          RunInfo("ERA5 Analysis",    "2001–2018"),
    "ERA5.ensmn":             RunInfo("ERA5 Ensemble",    "2001–2018"),
    "v3.LR.amip_ensmn":       RunInfo("E3SMv3 AMIP",      "2001–2018"),
    "v3.LR.historical_ensmn": RunInfo("E3SMv3 Historical","2001–2018"),
    "v3.LR.piControl":        RunInfo("E3SMv3 piControl", "0101–0200"),
    "v3.LR.piControl.hexagonal": RunInfo("E3SMv3 Hexagonal", "0001–0050"),
    "20250327.v3.SORRME3r3.piControl.alfred2.chrysalis":
        RunInfo("E3SMv3 SORRM", "0041–0050"),
}


# ---------------- Region Catalog ----------------

REGION_CATALOG: Dict[str, Region] = {
    "global":    Region((-90,  90), (-180, 180)),
    "Atlantic":  Region((10,   75), (-75, 0)),
    "CONUS":     Region((25,   50), (-125, -65)),
    "Arctic":    Region((50, 90), (-180, 180)),
    "Antarctic": Region((-90, -50), (-180, 180)),
    "PolarN":    Region((50,   90), (-180, 180)),
    "Greenland": Region((55,   90), (-75,  -10)),
}


# ---------------- Variable Catalog ----------------

VARIABLE_CATALOG: Dict[str, VariableInfo] = {
    "TS":         VariableInfo("K",          1.0, -8, 8,    17, "MODIS.observation"),
    "TURFLX":     VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "MERRA2.analysis"),
    "LHFLX":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "MERRA2.analysis"),
    "SHFLX":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "MERRA2.analysis"),
    "RADNET":     VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FSNS":       VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FSDS":       VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FSUS":       VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FSNSC":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FSDSC":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FSUSC":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FLNS":       VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FLDS":       VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FLUS":       VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FLNSC":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FLDSC":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "FLUSC":      VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "SWCF_SRF":   VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "LWCF_SRF":   VariableInfo("W m$^{-2}$", 1.0, -20, 20,  21, "CERES_EBAF.observation"),
    "ALBEDOC_SRF":VariableInfo("%",        100.0, -8, 8,     9, "CERES_EBAF.observation"),
    "ALBEDO_SRF": VariableInfo("%",        100.0, -8, 8,     9, "CERES_EBAF.observation"),
}

# ---------------- Budget Catalog ----------------

Budget_CATALOG: Dict[str, BudgetInfo] = {
    "TS_ALL":     BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (Model - Ref)",  "best_obs"),
    "TS_Sum":     BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (Sum)",        "best_obs"),
    "TS_nSAF":    BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (FSW_clr)",    "best_obs"),
    "TS_SAF":     BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (Albedo_clr)", "best_obs"),
    "TS_DLW":     BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (FLW_clr)",    "best_obs"),
    "TS_CRE":     BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (CRE_sfc)",    "best_obs"),
    "TS_TURB":    BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (SH+LE)",      "best_obs"),
    "TS_Qn":      BudgetInfo("K", 1.0, -8,  8,  17,  r"$\Delta$TS (Qn)",         "best_obs"),
}