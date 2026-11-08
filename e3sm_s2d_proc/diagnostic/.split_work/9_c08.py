# Daily model indices, the model climatology on the same days, and the ERA5 climatology on those days
m_u60, m_u60_clim, m_u60_adj, m_pch, m_vt, z_offset = {}, {}, {}, {}, {}, {}
std_pcz_m = sc.interp_logp(std_pcz, sc.DAILY_LEVELS)
for init in INIT_MONTHS:
    d = daily[init]
    lc = sc.interp_logp(lead_clim[init][["U", "Z"]], sc.DAILY_LEVELS)
    # native daily Z is ~40 m (10 hPa) to ~200 m (1000 hPa) above the monthly Z3 the climatology is made
    # from (online vs post-hoc vertical interpolation, extrapolation below ground); U and T agree: add
    # the ensemble-mean daily-minus-monthly offset of the same runs
    z_offset[init] = sc.daily_monthly_offset(d, monthly[init], YEAR, init)
    lc["Z"] = lc["Z"] + z_offset[init]
    cd = sc.daily_from_lead(lc, YEAR, init, d["time"].values)
    cd = cd.assign_coords(time=d["time"])
    m_u60[init] = sc.lat_mean(d["U"].sel(plev=REF_LEVEL), sc.U60)                      # (member, time)
    m_u60_clim[init] = sc.lat_mean(cd["U"].sel(plev=REF_LEVEL), sc.U60)                # (time)
    # bias-adjusted: remove the model-minus-ERA5 climatological difference of each day
    m_u60_adj[init] = m_u60[init] - (m_u60_clim[init] - on_days(clim_u60, m_u60_clim[init]))
    pcz = sc.lat_mean(d["Z"], sc.PCAP) - sc.lat_mean(cd["Z"], sc.PCAP)                  # (member, time, plev)
    m_pch[init] = (pcz / on_days(std_pcz_m, pcz)).dropna("plev", how="all")
    m_vt[init] = sc.lat_mean(d["VT"].sel(plev=100.0), sc.VT_BAND)

# Monthly model indices: U60 at REF_LEVEL (raw and anomaly) and equatorial U (QBO)
mm_u60, mm_u60_anom, mm_qbo = {}, {}, {}
for init in INIT_MONTHS:
    mo = monthly[init]
    mm_u60[init] = sc.lat_mean(mo["U"].sel(plev=REF_LEVEL), sc.U60)
    mm_u60_anom[init] = mm_u60[init] - sc.lat_mean(lead_clim[init]["U"].sel(plev=REF_LEVEL), sc.U60)
    mm_qbo[init] = sc.lat_mean(mo["U"], sc.EQ)                                          # (member, lead, plev)

# Check of the daily-minus-monthly Z offset: a processing difference should be nearly constant across
# leads (seasons) and starts; a physical signal would vary with the season and the initial state
rows = []
for init, off in z_offset.items():
    for p in (10.0, 50.0, 100.0, 500.0, 1000.0):
        v = sc.lat_mean(off.sel(plev=p), sc.PCAP)
        rows.append({"start": MONTH_NAMES[init], "plev (hPa)": p, "mean (m)": float(v.mean()),
                     "min over leads": float(v.min()), "max over leads": float(v.max()), "std over leads": float(v.std())})
print("Polar-cap (60-90N) daily-minus-monthly Z offset by lead (ensemble mean), per start:")
display(pd.DataFrame(rows).pivot_table(index="plev (hPa)", columns="start", values=["mean (m)", "std over leads"]).round(1))
