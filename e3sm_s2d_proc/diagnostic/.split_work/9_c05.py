%%time
era5 = sc.era5_daily_zonal(ERA5_ZONAL_DIR, years=(ERA5_FIRST_YEAR, YEAR), load=False)   # lazy: only reductions are loaded
obs_end = pd.Timestamp(era5["time"].values[-1])
print(f"ERA5 daily zonal means {pd.Timestamp(era5['time'].values[0]):%d %b %Y} - {obs_end:%d %b %Y}")


def doy_of(da):
    return xr.DataArray(sc.noleap_doy(da["time"].values), dims="time", coords={"time": da["time"]})


def on_days(clim, da):
    # a (doy, ...) climatology on the days of da
    return clim.sel(doy=doy_of(da)).drop_vars("doy")


# indices
e_u60 = sc.lat_mean(era5["U"].sel(plev=REF_LEVEL), sc.U60).load()       # (time)
e_pcz = sc.lat_mean(era5["Z"], sc.PCAP).load()                           # (time, plev)
e_pct = sc.lat_mean(era5["T"].sel(plev=REF_LEVEL), sc.PCAP).load()
e_vt = sc.lat_mean(era5["VT"].sel(plev=100.0), sc.VT_BAND).load()
e_qbo = sc.monthly_means(sc.lat_mean(era5["U"], sc.EQ)).load()           # (time, plev), monthly

clim_u60, std_u60 = sc.daily_climatology(e_u60, CLIM_YEARS)
clim_pcz, std_pcz = sc.daily_climatology(e_pcz, CLIM_YEARS)
clim_vt, _ = sc.daily_climatology(e_vt, CLIM_YEARS)
clim_pct, _ = sc.daily_climatology(e_pct, CLIM_YEARS)
e_pch = (e_pcz - on_days(clim_pcz, e_pcz)) / on_days(std_pcz, e_pcz)    # standardized polar-cap height

# calendar-month means and their 1981-2010 climatology (U, T, Z) for the model drift correction
e_month_clim = sc.obs_month_climatology(era5[["U", "T", "Z"]], CLIM_YEARS).load()   # (month, plev, lat)
e_u60_month = sc.monthly_means(e_u60)