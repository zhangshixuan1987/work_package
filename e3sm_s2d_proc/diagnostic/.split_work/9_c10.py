fig = plt.figure(figsize=(17, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=0.28)
ref_colors = {5: "tab:orange", 11: "tab:purple"}
for k, (field, name, units) in enumerate([("U", "U(60°N, 10 hPa)", "m s$^{-1}$"), ("T", "T(60-90°N, 10 hPa)", "K")]):
    ax = fig.add_subplot(gs[k])
    band = sc.U60 if field == "U" else sc.PCAP
    e_series = sc.monthly_means(sc.lat_mean(era5[field].sel(plev=REF_LEVEL), band))
    e_series = e_series.sel(time=slice(f"{CLIM_YEARS[0]}-01", f"{CLIM_YEARS[1]}-12"))
    g = e_series.groupby("time.month")
    order = [7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6]
    x = np.arange(12)
    mean, std = g.mean().sel(month=order).values, g.std().sel(month=order).values
    ax.fill_between(x, mean - std, mean + std, color="0.8", label="ERA5 ±1σ")
    ax.plot(x, mean, "k-", lw=2.5, label=f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]}")
    for m, ds in ref.items():
        idx = sc.lat_mean(ds[field].sel(plev=REF_LEVEL), band)                       # (year, member, lead)
        months = sc.target_months(m)
        for half in (0, 1):                                                          # leads 0-11, 12-23
            leads = np.arange(12 * half, 12 * half + 12)
            xs = np.array([order.index(months[l]) for l in leads])
            srt = np.argsort(xs)
            mu = idx.isel(lead=leads).mean(("year", "member")).values[srt]
            sd = idx.isel(lead=leads).std(("year", "member")).values[srt]
            ls = "-" if half == 0 else ":"
            ax.plot(xs[srt], mu, ls, color=ref_colors[m], lw=2,
                    label=f"E3SM {MONTH_NAMES[m]} starts, leads {leads[0]}-{leads[-1]}")
            if half == 0:
                ax.fill_between(xs[srt], mu - sd, mu + sd, color=ref_colors[m], alpha=0.15, lw=0)
    if field == "U":
        ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x, [MONTH_NAMES[m][0] for m in order])
    ax.set_title(f"({'ab'[k]}) {name}", loc="left")
    ax.set_ylabel(units)
    ax.grid(alpha=0.3)
    if k == 0:
        ax.legend(loc="upper left", fontsize=fontz * 0.65, frameon=False)

ax = fig.add_subplot(gs[2])
nov = ref[11]["U"].isel(lead=[1, 2, 3]).mean(("year", "member", "lead"))
obs_djf = e_month_clim_m["U"].sel(month=[12, 1, 2]).mean("month")
diff = (nov - obs_djf).sel(lat=slice(0, 90))
levels = np.arange(-12, 13, 2)
cf = ax.contourf(diff["lat"], diff["plev"], diff, levels=levels, cmap="RdBu_r", extend="both")
cs = ax.contour(obs_djf["lat"].sel(lat=slice(0, 90)), obs_djf["plev"], obs_djf.sel(lat=slice(0, 90)),
                levels=np.arange(-40, 81, 10), colors="k", linewidths=0.8)
ax.clabel(cs, fmt="%d", fontsize=fontz * 0.6)
ax.set_yscale("log")
ax.set_ylim(1000, 1)
ax.set_yticks([1000, 300, 100, 30, 10, 3, 1], ["1000", "300", "100", "30", "10", "3", "1"])
ax.set_xlabel("Latitude (°N)")
ax.set_ylabel("hPa")
ax.set_title("(c) DJF zonal-mean U: E3SM (Nov starts, leads 1-3) − ERA5", loc="left", fontsize=fontz * 0.9)
fig.colorbar(cf, ax=ax, label="m s$^{-1}$", shrink=0.9)
save(fig, "climatology", "u60_t60_10hPa", "djf_u_bias")
plt.show()

clim_table = pd.DataFrame({
    "ERA5": e_u60_month.sel(time=slice(f"{CLIM_YEARS[0]}-01", f"{CLIM_YEARS[1]}-12")).groupby("time.month").mean().to_series(),
    **{f"E3SM {MONTH_NAMES[m]} starts (leads 0-11)": pd.Series(
        sc.lat_mean(ds["U"].sel(plev=REF_LEVEL), sc.U60).isel(lead=slice(0, 12)).mean(("year", "member")).values,
        index=sc.target_months(m)[:12]) for m, ds in ref.items()},
}).loc[[10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]].round(1)
clim_table.index = [MONTH_NAMES[m] for m in clim_table.index]
print(f"U(60N, {REF_LEVEL:g} hPa) monthly climatology, m/s")
display(clim_table.T)

# DJF vortex strength and its interannual variability: the model spread sets the width of its tercile bounds (section 6)
W_DJF = xr.DataArray([31.0, 31.0, 28.0], dims="lead")
var_rows = []
for m, ds in ref.items():
    leads = [(12 - m) % 12, (12 - m) % 12 + 1, (12 - m) % 12 + 2]
    u = sc.lat_mean(ds["U"].sel(plev=REF_LEVEL), sc.U60).isel(lead=leads)
    djf = (u * W_DJF.assign_coords(lead=u["lead"])).sum("lead") / 90.0
    var_rows.append({"set": f"E3SM {MONTH_NAMES[m]} starts (leads {leads[0]}-{leads[-1]})", "DJF mean (m/s)": float(djf.mean()),
                     "total sd": float(djf.std()), "sd of yearly ensemble means": float(djf.mean("member").std("year")),
                     "within-year member sd": float(djf.std("member").mean())})
e_d = [float((e_u60_month.sel(time=[np.datetime64(f"{y}-12-15"), np.datetime64(f"{y + 1}-01-15"), np.datetime64(f"{y + 1}-02-15")])
              * np.array([31.0, 31.0, 28.0])).sum() / 90.0) for y in range(CLIM_YEARS[0], CLIM_YEARS[1])]
var_rows.append({"set": f"ERA5 {CLIM_YEARS[0]}/{str(CLIM_YEARS[0] + 1)[2:]}-{CLIM_YEARS[1] - 1}/{str(CLIM_YEARS[1])[2:]}",
                 "DJF mean (m/s)": np.mean(e_d), "total sd": np.std(e_d, ddof=1)})
print(f"DJF U(60N, {REF_LEVEL:g} hPa): mean and interannual standard deviation")
display(pd.DataFrame(var_rows).round(2))
