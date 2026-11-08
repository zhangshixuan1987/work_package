%%time
ref = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", range(CLIM_YEARS[0], CLIM_YEARS[1] + 1), m,
                           mode=derivation_mode, workers=WORKERS)
       for m in rc.REFERENCE["init_months"]}
# Feb/Aug Reanalysis starts exist every fifth year (1985-2010 within CLIM_YEARS): the "direct" climatology
ref_direct = {m: sc.monthly_zonal(rc.REFERENCE["case_prefix"], "Reanalysis", range(CLIM_YEARS[0], CLIM_YEARS[1] + 1), m,
                                  mode=derivation_mode, workers=WORKERS)
              for m in rc.DIRECT_INIT_MONTHS}
REF_ALL = {**ref, **ref_direct}


def ref_for(method):
    return REF_ALL if method == "direct" else ref


MPLEV = ref[rc.REFERENCE["init_months"][0]]["plev"].values
e_month_clim_m = e_month_clim.sel(plev=MPLEV)            # ERA5 on the model monthly levels (all are ERA5 levels)
lead_clim = {init: sc.lead_climatology(init, ref_for(CLIM_METHOD), e_month_clim_m, method=CLIM_METHOD) for init in INIT_MONTHS}


def source_of(prefix):
    return prefix.split(f"{MODEL_PREFIX}_")[-1]


def target_time(init, nlead=rc.NLEAD):
    out = []
    for lead in range(nlead):
        m = init - 1 + lead
        out.append(np.datetime64(f"{YEAR + m // 12:04d}-{m % 12 + 1:02d}-15"))
    return np.asarray(out, dtype="datetime64[ns]")


def pooled_monthly(init):
    # pooled monthly zonal means (member, lead, plev, lat) of one start; repeats counted once
    parts = []
    for name, prefix in POOL.items():
        ds = sc.monthly_zonal(prefix, source_of(prefix), [YEAR], init, mode=derivation_mode, workers=WORKERS)
        if ds is None:
            continue
        ds = ds.isel(year=0, drop=True)
        ds = ds.isel(member=np.flatnonzero(np.isfinite(ds["U"].sel(plev=REF_LEVEL).isel(lead=0, lat=-1)).values))
        u = ds["U"].sel(plev=REF_LEVEL).values
        if any(p.sizes["member"] == ds.sizes["member"] and np.array_equal(u, p["U"].sel(plev=REF_LEVEL).values, equal_nan=True)
               for p in parts):
            print(f"[monthly {name} init{init:02d}] repeats an already pooled case; not pooled again")
            continue
        parts.append(ds.assign_coords(member=[f"{name}:{m}" for m in ds["member"].values]))
    out = xr.concat(parts, dim="member")
    return out.assign_coords(time=("lead", target_time(init, out.sizes["lead"])))


monthly = {init: pooled_monthly(init) for init in INIT_MONTHS}
daily = {init: sc.pooled_daily_zonal(SIMULATION_ROOT, POOL, YEAR, init, mode=derivation_mode, workers=WORKERS)
         for init in INIT_MONTHS}
for init in INIT_MONTHS:
    d = daily[init]
    print(f"[init {MONTH_NAMES[init]} {YEAR}] {d.sizes['member']} members, daily "
          f"{pd.Timestamp(d['time'].values[0]):%d %b %Y} - {pd.Timestamp(d['time'].values[-1]):%d %b %Y}; "
          f"monthly {monthly[init].sizes['member']} members; climatology: {lead_clim[init].attrs['climatology_method']}")