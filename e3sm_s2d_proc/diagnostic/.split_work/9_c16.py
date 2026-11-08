plume_span = (pd.Timestamp(YEAR, 8, 1), pd.Timestamp(YEAR + 2, 5, 31))
fig, axes = plt.subplots(len(INIT_MONTHS), 1, figsize=(16, 3.6 * len(INIT_MONTHS)), sharex=True, sharey=True,
                         gridspec_kw={"hspace": 0.25})
days = pd.date_range(*plume_span, freq="D")
dummy = xr.DataArray(np.zeros(days.size), dims="time", coords={"time": days.values})
cl, sd = on_days(clim_u60, dummy), on_days(std_u60, dummy)
model_ssw = {}                        # init -> DataFrame (series, winter, central_date) on the raw wind
model_ssw_adj = {}                    # bias-adjusted
for ax, init in zip(np.atleast_1d(axes), INIT_MONTHS):
    u = m_u60[init]
    ax.fill_between(days, cl - sd, cl + sd, color="0.82", label="ERA5 ±1σ")
    ax.plot(days, cl, color="0.35", lw=1.5, label="ERA5 climatology")
    for mem in u["member"].values:
        ax.plot(u["time"], u.sel(member=mem), color=INIT_COLORS[init], lw=0.5, alpha=0.35)
    ax.plot(u["time"], u.mean("member"), color=INIT_COLORS[init], lw=2.5, label="ensemble mean")
    ax.plot(u["time"], m_u60_clim[init], "--", color="k", lw=1.5, label="E3SM lead-time climatology")
    e_part = e_u60.sel(time=slice(*plume_span))
    if e_part.sizes["time"]:
        ax.plot(e_part["time"], e_part, "k-", lw=2.5, label="ERA5")
    tabs, tabs_adj = [], []
    for mem in u["member"].values:
        tabs.append(sc.ssw_table(u.sel(member=mem).values, u["time"].values, str(mem), **SSW))
        tabs_adj.append(sc.ssw_table(m_u60_adj[init].sel(member=mem).values, u["time"].values, str(mem), **SSW))
    model_ssw[init] = pd.concat(tabs, ignore_index=True)
    model_ssw_adj[init] = pd.concat(tabs_adj, ignore_index=True)
    ev = model_ssw[init]
    if len(ev):
        y = [float(u.sel(member=s, time=np.datetime64(d)).values) for s, d in zip(ev["series"], ev["central_date"])]
        ax.plot(ev["central_date"], np.asarray(y), "o", ms=6, mfc="w", mec="k", label="SSW (raw wind)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylim(-35, 70)
    ax.set_ylabel("m s$^{-1}$")
    ax.grid(alpha=0.3)
    ax.set_title(f"E3SM 1 {MONTH_NAMES[init]} {YEAR} start ({u.sizes['member']} members): U(60°N, {REF_LEVEL:g} hPa)", loc="left")
np.atleast_1d(axes)[0].legend(loc="upper left", ncol=3, fontsize=fontz * 0.65, frameon=False)
np.atleast_1d(axes)[-1].set_xlim(*plume_span)
np.atleast_1d(axes)[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
save(fig, "plume", "u60n10", YEAR)
plt.show()