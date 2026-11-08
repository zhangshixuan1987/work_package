%%time
frame_days = pd.date_range(FRAMES[0], FRAMES[1], freq=f"{FRAME_STEP}D")
all_days = pd.date_range(FRAMES[0], FRAMES[1], freq="D")


def nh(da):
    return da.sel(lat=slice(LAT_MIN - 0.01, 90.0))


# ERA5: U(60N, 10 hPa) and its smoothed daily climatology
era5_u = sc.era5_daily_zonal(ERA5_ZONAL_DIR, fields=("U",), plev=[10.0])["U"].sel(plev=10.0)
e_u60 = sc.lat_mean(era5_u, sc.U60)
clim_u60, std_u60 = sc.daily_climatology(e_u60, CLIM_YEARS)
obs_end = pd.Timestamp(e_u60["time"].values[-1])
doy = xr.DataArray(sc.noleap_doy(all_days.values), dims="time", coords={"time": all_days.values})
u60_clim_days, u60_std_days = clim_u60.sel(doy=doy).drop_vars("doy"), std_u60.sel(doy=doy).drop_vars("doy")

# ERA5 10-hPa height: 1981-2010 monthly climatology (e3sm_diags ERA5 zg) to days, and the observed days
z_clim_month = nh(rc.obs_monthly_climatology("Z10", mode=derivation_mode)).load()
z_clim = sc.daily_from_month_clim(z_clim_month, all_days.values)                     # (time, lat, lon)
obs_days = all_days[all_days <= obs_end]
if len(obs_days):
    z_obs = nh(sc.era5_daily_maps(ERA5_MIRROR, obs_days[0], obs_days[-1], mode=derivation_mode, workers=WORKERS))
    z_obs = z_obs.interp(lat=z_clim["lat"], lon=z_clim["lon"]) if not np.array_equal(z_obs["lon"], z_clim["lon"]) else z_obs
    obs_anom = z_obs - z_clim.sel(time=z_obs["time"])
    print(f"ERA5 10-hPa maps {obs_days[0]:%d %b %Y} - {obs_days[-1]:%d %b %Y}")
else:
    z_obs = obs_anom = None
    print(f"ERA5 ends {obs_end:%d %b %Y}: winter {WINTER_LABEL} not observed yet")

# Model: maps, anomalies and totals of the ensemble mean and the featured member; U(60N, 10 hPa)
fc = {}          # init -> dict of DataArrays on the model days within the frames
for init in INIT_MONTHS:
    zd = sc.pooled_daily_zonal(SIMULATION_ROOT, POOL, YEAR, init, mode=derivation_mode, workers=WORKERS)
    maps = sc.pooled_daily_maps(SIMULATION_ROOT, MAP_FILE, POOL, YEAR, init, mode=derivation_mode, workers=WORKERS)
    if zd is None or maps is None:
        print(f"[init {MONTH_NAMES[init]}] no model output")
        continue
    u60 = sc.lat_mean(zd["U"].sel(plev=10.0), sc.U60)                                   # (member, time)
    days = all_days[all_days.isin(pd.DatetimeIndex(maps["time"].values))]
    if not len(days):
        print(f"[init {MONTH_NAMES[init]}] no days in the frames")
        continue
    z = nh(maps.sel(time=days.values)).load()
    lead_clim = nh(rc.matched_climatology("Z10", init, method=CLIM_METHOD, mode=derivation_mode))
    # native daily Z010 vs the monthly Z3 product the climatology is made from (~40 m): add the
    # ensemble-mean offset of the same runs, month by month (Sep-May, the months the map cache keeps)
    mon = [rc.open_hindcast(prefix, YEAR, init, "Z10") for prefix in POOL.values()]
    mon = nh(xr.concat([m.drop_vars("member") for m in mon if m is not None], dim="member").mean("member"))
    dmean = sc.monthly_means(nh(maps).mean("member")).load()
    targets = sc.to_datetime64(mon["time"].values).astype("datetime64[M]").astype("datetime64[ns]") + np.timedelta64(14, "D")
    z_off = dmean.reindex(time=targets).rename(time="lead").assign_coords(lead=mon["lead"].values) - mon.drop_vars("time")
    clim_days = sc.daily_from_lead(lead_clim + z_off, YEAR, init, days.values)
    anom = z - clim_days
    # featured member: FEATURED, else the first SSW of the winter, else the weakest DJF vortex
    feat = FEATURED.get(init)
    if feat is None:
        first = []
        for mem in u60["member"].values:
            ev = [d for d in sc.ssw_central_dates(u60.sel(member=mem).values, u60["time"].values, **SSW)
                  if sc.winter_of(d) == WINTER]
            if ev:
                first.append((ev[0], str(mem)))
        if first:
            date, feat = min(first)
            why = f"first SSW of the winter ({date:%d %b %Y}); {len(first)} of {u60.sizes['member']} members have an SSW"
        else:
            djf = u60.sel(time=slice(f"{WINTER}-12-01", f"{WINTER + 1}-02-28")).mean("time")
            feat = str(djf["member"].values[int(djf.argmin())])
            why = "weakest DJF vortex (no member has an SSW this winter)"
    else:
        why = "set in FEATURED"
    feat_ssw = [d for d in sc.ssw_central_dates(u60.sel(member=feat).values, u60["time"].values, **SSW) if sc.winter_of(d) == WINTER]
    n_ssw = sum(any(sc.winter_of(d) == WINTER for d in sc.ssw_central_dates(u60.sel(member=m).values, u60["time"].values, **SSW))
                for m in u60["member"].values)
    fc[init] = {
        "days": days, "n": int(z.sizes["member"]), "featured": feat, "why": why,
        "mean_anom": anom.mean("member"), "mean_total": anom.mean("member") + z_clim.sel(time=days.values),
        "feat_anom": anom.sel(member=feat), "feat_total": anom.sel(member=feat) + z_clim.sel(time=days.values),
        "u60_mean": u60.mean("member"), "u60_feat": u60.sel(member=feat), "feat_ssw": feat_ssw, "n_ssw": n_ssw,
        "surrogate": CLIM_METHOD != "direct" and init not in rc.REFERENCE["init_months"],
    }
    print(f"[init {MONTH_NAMES[init]} {YEAR}] {fc[init]['n']} members, {len(days)} days in the frames; "
          f"featured {feat}: {why}")

# All starts pooled (member-weighted), on the days every start covers: the first panel when the winter is not observed
pool_all = None
if fc:
    common = pd.DatetimeIndex(sorted(set.intersection(*(set(f["days"]) for f in fc.values()))))
    n_all = sum(f["n"] for f in fc.values())
    pa = sum(f["n"] * f["mean_anom"].sel(time=common.values) for f in fc.values()) / n_all
    pu = sum(f["n"] * f["u60_mean"].sel(time=common.values) for f in fc.values()) / n_all
    pool_all = {"n": n_all, "anom": pa, "total": pa + z_clim.sel(time=common.values), "u60": pu}
    print(f"[all starts] {n_all} members pooled on {len(common)} common days")
