%%time
days = pd.date_range(FRAMES[0], FRAMES[1], freq="D")
frame_days = pd.date_range(FRAMES[0], FRAMES[1], freq=f"{FRAME_STEP}D")


def nh(da):
    return da.sel(lat=slice(LAT_MIN - 0.01, 90.0))


def on_days(da, d):
    return da.sel(time=d.values)


# ERA5: climatology of the 10-hPa height (maps, e3sm_diags zg) and of U(60N, 10 hPa) (daily zonal means)
z_clim_era5 = sc.daily_from_month_clim(nh(rc.obs_monthly_climatology("Z10", mode="auto")).load(), days.values)
e_u60 = sc.lat_mean(sc.era5_daily_zonal(ERA5_ZONAL_DIR, fields=("U",), plev=[10.0])["U"].sel(plev=10.0), sc.U60)
u_clim, u_std = sc.daily_climatology(e_u60, CLIM_YEARS)
doy = xr.DataArray(sc.noleap_doy(days.values), dims="time", coords={"time": days.values})
u60_era5_clim, u60_era5_std = u_clim.sel(doy=doy).drop_vars("doy"), u_std.sel(doy=doy).drop_vars("doy")

# Reanalysis zonal-mean climatologies of every start month (U(60N, 10 hPa) climatology of each start)
REF_ALL = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", range(CLIM_YEARS[0], CLIM_YEARS[1] + 1), m,
                               mode=derivation_mode, verbose=False).sel(plev=[10.0])
           for m in rc.REFERENCE["init_months"] + rc.DIRECT_INIT_MONTHS}

fc = {}
for init in INIT_MONTHS:
    maps = sc.pooled_daily_maps(SIMULATION_ROOT, MAP_FILE, POOL, YEAR, init, mode=derivation_mode)
    zd = sc.pooled_daily_zonal(SIMULATION_ROOT, POOL, YEAR, init, mode=derivation_mode)
    if maps is None or zd is None:
        print(f"[init {MONTH_NAMES[init]}] no cached fields")
        continue
    mem = MEMBERS[init]
    if mem not in maps["member"].values:
        raise ValueError(f"{mem} is not a member of the {MONTH_NAMES[init]} start: {list(maps['member'].values)[:6]} ...")
    d = days[days.isin(pd.DatetimeIndex(maps["time"].values))]
    raw_mem = nh(on_days(maps.sel(member=mem), d)).load()
    raw_mean = nh(on_days(maps, d).mean("member")).load()
    clim = sc.map_lead_climatology(POOL, YEAR, init, maps, d.values, method=CLIM_METHOD, lat_min=LAT_MIN,
                                   mode="auto").load() if (SHADING or SHOW_REFERENCE) else None
    u60 = sc.lat_mean(zd["U"].sel(plev=10.0), sc.U60)
    fc[init] = {"days": d, "n": int(maps.sizes["member"]), "member": mem,
                "raw_mem": raw_mem, "raw_mean": raw_mean, "clim": clim,
                "anom_mem": raw_mem - clim if clim is not None else None,
                "anom_mean": raw_mean - clim if clim is not None else None,
                "u60_mem": u60.sel(member=mem), "u60_mean": u60.mean("member"),
                "u60_lo": u60.quantile(0.1, "member"), "u60_hi": u60.quantile(0.9, "member"),
                "u60_clim": sc.lat_mean(sc.daily_from_lead(
                    sc.lead_climatology(init, REF_ALL, None, method="direct")["U"].sel(plev=10.0), YEAR, init,
                    u60["time"].values), sc.U60).assign_coords(time=u60["time"]),
                "ssw": [t for t in sc.ssw_central_dates(u60.sel(member=mem).values, u60["time"].values, **SSW)
                        if pd.Timestamp(FRAMES[0]) <= t <= pd.Timestamp(FRAMES[1])]}
    print(f"[init {MONTH_NAMES[init]} {YEAR}] member {mem} of {fc[init]['n']}; {len(d)} days; "
          f"climatology: {clim.attrs.get('climatology_method', '') if clim is not None else 'none'}")

# E3SM climatology of the day for the reference column: the May starts (30 years), raw contours
ref_init = 5 if 5 in fc else next(iter(fc))
z_clim_e3sm = fc[ref_init]["clim"]