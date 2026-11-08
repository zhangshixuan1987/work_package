%%time
LEVS = [10.0]


def vortex_indices(ds):
    return {"u60": sc.lat_mean(ds["U"].sel(plev=10.0), sc.U60),
            "t10": sc.lat_mean(ds["T"].sel(plev=10.0), sc.PCAP),
            "z10": sc.lat_mean(ds["Z"].sel(plev=10.0), sc.PCAP)}


def on_doy(clim, times):
    d = xr.DataArray(sc.noleap_doy(times), dims="time", coords={"time": times})
    return clim.sel(doy=d).drop_vars("doy")


# ERA5 indices and their smoothed daily climatologies (mean, std)
era5_z = sc.era5_daily_zonal(ERA5_ZONAL_DIR, fields=("U", "T", "Z"), plev=LEVS)
e_idx = vortex_indices(era5_z)
e_idx["slp"] = sc.era5_polar_cap_slp(ERA5_DAILY_DIR, years=(ERA5_FIRST_YEAR, YEAR))
e_iclim = {k: sc.daily_climatology(v, CLIM_YEARS) for k, v in e_idx.items()}

# model indices: raw, the lead-time climatology on the same days, and the anomaly
ref = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", range(CLIM_YEARS[0], CLIM_YEARS[1] + 1), m,
                           mode=derivation_mode, workers=WORKERS).sel(plev=LEVS)
       for m in rc.REFERENCE["init_months"]}
ref_direct = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", range(CLIM_YEARS[0], CLIM_YEARS[1] + 1), m,
                                  mode=derivation_mode, workers=WORKERS).sel(plev=LEVS)
              for m in rc.DIRECT_INIT_MONTHS}      # Feb/Aug starts of every fifth year (1985-2010)
REF_ALL = {**ref, **ref_direct}
e_month = sc.obs_month_climatology(era5_z[["U", "T", "Z"]], CLIM_YEARS)
m_idx = {}
for init in INIT_MONTHS:
    zd = sc.pooled_daily_zonal(SIMULATION_ROOT, POOL, YEAR, init, mode=derivation_mode, workers=WORKERS)
    if zd is None:
        continue
    # daily (native) Z differs from the monthly product by a near-constant offset: add it to the climatology
    mo = sc.pooled_monthly_zonal(POOL, YEAR, init, mode=derivation_mode, workers=WORKERS)
    z_off = sc.daily_monthly_offset(zd, mo, YEAR, init).sel(plev=LEVS)
    zd = zd.sel(plev=LEVS)
    lc = sc.lead_climatology(init, REF_ALL if CLIM_METHOD == "direct" else ref, e_month, method=CLIM_METHOD)
    lc["Z"] = lc["Z"] + z_off
    cd = sc.daily_from_lead(lc, YEAR, init, zd["time"].values).assign_coords(time=zd["time"])
    raw, clim = vortex_indices(zd), vortex_indices(cd)
    # surface: polar-cap sea-level pressure against the drift-corrected PSL climatology (rtenso_case)
    psl_clim = sc.lat_mean(rc.matched_climatology("PSL", init, method=CLIM_METHOD, mode=derivation_mode).mean("lon"), sc.PCAP)
    raw["slp"] = sc.lat_mean(zd["PSL"], sc.PCAP)
    clim["slp"] = sc.daily_from_lead(psl_clim, YEAR, init, zd["time"].values).assign_coords(time=zd["time"])
    alt = None
    ALT = "lead_drift" if CLIM_METHOD == "direct" else "target_month"
    if init not in ref:                       # Feb/Aug: the climatology method matters; show the alternative
        lc_alt = sc.lead_climatology(init, ref, e_month, method=ALT)
        lc_alt["Z"] = lc_alt["Z"] + z_off
        clim_alt = vortex_indices(sc.daily_from_lead(lc_alt, YEAR, init, zd["time"].values).assign_coords(time=zd["time"]))
        pc_alt = sc.lat_mean(rc.matched_climatology("PSL", init, method=ALT, mode=derivation_mode).mean("lon"), sc.PCAP)
        clim_alt["slp"] = sc.daily_from_lead(pc_alt, YEAR, init, zd["time"].values).assign_coords(time=zd["time"])
        alt = {k: raw[k] - clim_alt[k] for k in ("z10", "slp")}
    m_idx[init] = {"raw": raw, "clim": clim, "anom": {k: raw[k] - clim[k] for k in raw}, "anom_alt": alt,
                   "ssw": [d for mem in raw["u60"]["member"].values
                           for d in sc.ssw_central_dates(raw["u60"].sel(member=mem).values, raw["u60"]["time"].values, **SSW)]}
obs_ssw = sc.ssw_central_dates(e_idx["u60"].values, e_idx["u60"]["time"].values, **SSW)