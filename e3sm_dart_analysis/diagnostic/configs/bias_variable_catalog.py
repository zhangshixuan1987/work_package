"""Legacy broad atmosphere/land variable metadata preserved from analysis_bias."""

BIAS_VARIABLE_CATALOG = {
    "lnd": {
        "SOILM": {
            "units": "kg m-2",
            "description": "Total soil moisture content (integrated over depth)",
            "mask": "land"
        },
        "SMOIS": {
            "units": "kg m-2",
            "description": "Volumetric soil moisture content (top 5cm)",
            "mask": "land"
        },
        "SOILLIQ": {
            "units": "kg m-2",
            "description": "Liquid water content in soil layers",
            "mask": "land"
        },
        "SOILICE": {
            "units": "kg m-2",
            "description": "Liquid water content in soil layers",
            "mask": "land"
        },
        "FGEV" : {
            "units": "W m-2",
            "description": "Ground evaporation",
            "mask": "land"
        },
        "QSOIL": {
            "units": "mm s-1",
            "description": "Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)",
            "mask": "land"
        },
        "H2OSFC": {
            "units": "mm",
            "description": "Surface water depth",
            "mask": "land"
        },
        "H2OSNO": {
            "units": "mm",
            "description": "Snow depth (liquid water)",
            "mask": "land"
        },
        "HC": {
            "units": "MJ m-2",
            "description": "heat content of soil/snow/lake",
            "mask": "land"
        },
        "HCSOI": {
            "units": "MJ m-2",
            "description": "soil heat content",
            "mask": "land"
        },
        "H2OSOI": {
            "units": "mm3/mm3",
            "description": "Volumetric soil water (vegetated landunits only)",
            "mask": "land"
        },
        "TSOI": {
            "units": "K",
            "description": "Soil temperature (vegetated landunits only)",
            "mask": "land"
        },
        "TSOI_ICE": {
            "units": "K",
            "description": "Soil temperature (ice landunits only)",
            "mask": "land"
        },
        "TSOI_10CM": {
            "units": "K",
            "description": "Soil temperature in top 10cm of soil",
            "mask": "land"
        },
        "SOILWATER_10CM": {
            "units": "kg m-2",
            "description": "Soil liquid water + ice in top 10cm of soil",
            "mask": "land"
        },
        "RAIN": {
            "units": "mm s-1",
            "description": "Atmospheric rain",
            "mask": "land"
        },
        "SNOW": {
            "units": "mm s-1",
            "description": "Atmospheric snow",
            "mask": "land"
        },
        "SNOWDP": {
            "units": "m",
            "description": "Gridcell mean snow height",
            "mask": "land"
        },
        "TLAI": {
            "units": "1",
            "description": "Total projected leaf area index",
            "mask": "land"
        },
        "TSAI": {
            "units": "1",
            "description": "Total projected stem area index",
            "mask": "land"
        },
        "TAUX": {
            "units": "N m-2",
            "description": "Zonal surface stress",
            "mask": "land"
        },
        "TAUY": {
            "units": "N m-2",
            "description": "Meridional surface stress",
            "mask": "land"
        },
        "TREFMNAV": {
            "units": "K",
            "description": "Daily minimum of average 2-m temperature",
            "mask": "land"
        },
        "TREFMXAV": {
            "units": "K",
            "description": "Daily maximum of average 2-m temperature",
            "mask": "land"
        },
        "U10": {
            "units": "m s-1",
            "description": "10-m wind",
            "mask": "land"
        },
        "TS": {
            "units": "K",
            "description": "Surface temperature (skin temperature)",
            "mask": "land"
        },
        "TSA": {
            "units": "K",
            "description": "Air temperature at 2 meters height",
            "mask": "land"
        },
        "EFLX_LH_TOT": {
            "units": "W m-2",
            "description": "Latent heat flux (evaporation + transpiration)",
            "mask": "land"
        },
        "FSH": {
            "units": "W m-2",
            "description": "Sensible heat flux from surface",
            "mask": "land"
        },
        "FSA": {
            "units": "W m-2",
            "description": "Absorbed solar radiation",
            "mask": "land"
        },
        "FSDS": {
            "units": "W m-2",
            "description": "Atmospheric incident solar radiation",
            "mask": "land"
        },
        "FLDS": {
            "units": "W m-2",
            "description": "Atmospheric longwave radiation",
            "mask": "land"
        },
        "FIRE": {
            "units": "W m-2",
            "description": "Emitted infrared (longwave) radiation",
            "mask": "land" 
        }, 
        "FIRA": {
            "units": "W m-2",
            "description": "Net shortwave radiation at surface",
            "mask": "land"
        },
        "FSR": {
            "units": "W m-2",
            "description": "Reflected solar radiation",
            "mask": "land"
        },
        "TWS": {
            "units": "mm",
            "description": "Total water storage",
            "mask": "land"
        },
        "QVEGT": {
            "units": "mm s-1",
            "description": "Canopy transpiration",
            "mask": "land"
        },
        "QVEGT": {
            "units": "mm s-1",
            "description": "Canopy evaporation",
            "mask": "land"
        }
    },
    "atm": {
        "PRECT": {
            "units": "mm day-1",
            "description": "Total precipitation (convective + large-scale)"
        },
        "TS": {
            "units": "K",
            "description": "Surface temperature (radiative)",
        },
        "TREFHT": {
            "units": "K",
            "description": "Reference height temperature",
        },
        "QREFHT": {
            "units": "kg kg-1",
            "description": "Reference height humidity",
        },
        "RHREFHT": {
            "units": "1",
            "description": "Reference height relative humidity",
        },
        "QFLX": {
            "units": "kg m-2 s-1",
            "description": "Surface water flux",
        },
        "TAUX": {
            "units": "N m-2",
            "description": "Zonal surface stress",
        },
        "TAUY": {
            "units": "N m-2",
            "description": "Meridional surface stress",
        },
        "TMQ": {
            "units": "kg m-2",
            "description": "Total (vertically integrated) precipitable water",
        },
        "TUQ": {
            "units": "kg m-2 s-1",
            "description": "Total (vertically integrated) zonal water flux",
        },
        "TVQ": {
            "units": "kg m-2 s-1",
            "description": "Total (vertically integrated) meridional water flux",
        },
        "TGCLDLWP": {
            "units": "kg m-2",
            "description": "Total grid-box cloud liquid water path",
        },
        "TGCLDLWP": {
            "units": "kg m-2",
            "description": "Total grid-box cloud liquid water path",
        },
        "TGCLDIWP": {
            "units": "kg m-2",
            "description": "Total grid-box cloud ice water path",
        },
        "TGCLDCWP": {
            "units": "kg m-2",
            "description": "Total grid-box cloud water path (ice+liquid)",
        },
        "Z850": {
            "units": "gpm",
            "description": "Geopotential height at 850 hPa"
        },
        "Z700": {
            "units": "gpm",
            "description": "Geopotential height at 700 hPa"
        },
        "Z500": {
            "units": "gpm",
            "description": "Geopotential height at 500 hPa"
        },
        "Z200": {
            "units": "gpm",
            "description": "Geopotential height at 200 hPa"
        },
        "Q850": {
            "units": "kg kg-1",
            "description": "Specific humidity at 850 hPa"
        },
        "Q700": {
            "units": "kg kg-1",
            "description": "Specific humidity at 700 hPa"
        },
        "Q500": {
            "units": "kg kg-1",
            "description": "Specific humidity at 500 hPa"
        },
        "T850": {
            "units": "K",
            "description": "Specific humidity at 850 hPa"
        },
        "T700": {
            "units": "K",
            "description": "Specific humidity at 700 hPa"
        },
        "T500": {
            "units": "K",
            "description": "Specific humidity at 500 hPa"
        },
        "T200": {
            "units": "K",
            "description": "Specific humidity at 200 hPa"
        },
        "U850": {
            "units": "m s-1",
            "description": "Zonal wind at 850 hPa"
        },
        "V850": {
            "units": "m s-1",
            "description": "Meridional wind at 850 hPa"
        },
        "U500": {
            "units": "m s-1",
            "description": "Zonal wind at 500 hPa"
        },
        "V500": {
            "units": "m s-1",
            "description": "Meridional wind at 500 hPa"
        },
        "U200": {
            "units": "m s-1",
            "description": "Zonal wind at 200 hPa"
        },
        "V200": {
            "units": "m s-1",
            "description": "Meridional wind at 200 hPa"
        },
        "OMEGA500": {
            "units": "m s-1",
            "description": "Vertical velocity at 500 mbar pressure surface"
        },
        "U10": {
            "units": "m s-1",
            "description": "10m wind speed"
        },
        "U90M": {
            "units": "m s-1",
            "description": "Zonal wind at turbine hub height (90m above surface)"
        },
        "V90M": {
            "units": "m s-1",
            "description": "Meridional wind at turbine hub height (90m above surface)"
        },
        "PS": {
            "units": "Pa",
            "description": "Surface pressure"
        },
        "PSL": {
            "units": "Pa",
            "description": "Sea level pressure"
        },
        "CLDLOW": {
            "units": "1",
            "description": "Low-level cloud fraction (if available)"
        },
        "CLDTOT": {
            "units": "1",
            "description": "Total cloud cover (if available)"
        },
        "SWCF": {
            "units": "W m-2",
            "description": "Shortwave cloud forcing at TOA"
        },
        "LWCF": {
            "units": "W m-2",
            "description": "Longwave cloud forcing at TOA"
        },
        "FLUT": {
            "units": "W m-2",
            "description": "Upwelling longwave flux at top of model"
        },
        "FSNT": {
            "units": "W m-2",
            "description": "Net shortwave radiation at the top of the atmosphere"
        },
        "FLNT": {
            "units": "W m-2",
            "description": "Net longwave radiation at the top of the atmosphere"
        }
   }
}
