qbo_span = (np.datetime64(f"{YEAR - 1}-01-01"), np.datetime64(f"{YEAR + 2}-08-01"))
fig, axes = plt.subplots(len(INIT_MONTHS) + 1, 2, figsize=(17, 3.0 * (len(INIT_MONTHS) + 1)),
                         gridspec_kw={"width_ratios": [1.6, 1], "hspace": 0.45, "wspace": 0.18}, sharex="col")
qlev = np.arange(-35, 36, 5)
eq = e_qbo.sel(time=slice(*qbo_span), plev=slice(1, 100))
ax = axes[0, 0]
cf = ax.contourf(eq["time"], eq["plev"], eq.T, levels=qlev, cmap="RdBu_r", extend="both")
ax.set_title(f"(a) ERA5 zonal-mean U, 5°S-5°N (through {obs_end:%b %Y})", loc="left")
axes[0, 1].axis("off")
for k, init in enumerate(INIT_MONTHS, start=1):
    q = mm_qbo[init].sel(plev=slice(1, 100))
    qt = q.assign_coords(time=monthly[init]["time"]).swap_dims(lead="time")
    ax = axes[k, 0]
    ax.contourf(qt["time"], qt["plev"], qt.mean("member").T, levels=qlev, cmap="RdBu_r", extend="both")
    ov = eq.sel(time=slice(qt["time"].values[0], qt["time"].values[-1]))
    if ov.sizes["time"] > 1:
        cs = ax.contour(ov["time"], ov["plev"], ov.T, levels=[-20, -10, 0, 10], colors="k", linewidths=0.8)
        ax.clabel(cs, fmt="%d", fontsize=fontz * 0.55)
    ax.axvline(np.datetime64(f"{YEAR}-{init:02d}-01"), color="k", lw=1.5)
    ax.set_title(f"({'abcdefgh'[2 * k - 1]}) E3SM {MONTH_NAMES[init]} {YEAR} start, ensemble mean "
                 f"({monthly[init].sizes['member']} members; contours: ERA5)", loc="left", fontsize=fontz * 0.85)
    ax2 = axes[k, 1]
    q50 = qt.sel(plev=QBO_LEVEL)
    for mem in q50["member"].values:
        ax2.plot(q50["time"], q50.sel(member=mem), color=INIT_COLORS[init], lw=0.6, alpha=0.5)
    ax2.plot(q50["time"], q50.mean("member"), color=INIT_COLORS[init], lw=2.5)
    e50 = e_qbo.sel(plev=QBO_LEVEL, time=slice(*qbo_span))
    ax2.plot(e50["time"], e50, "k-", lw=2)
    ax2.axhline(0, color="0.4", lw=0.8)
    ax2.set_ylim(-35, 25)
    ax2.set_title(f"({'abcdefgh'[2 * k]}) U {QBO_LEVEL:g} hPa, 5°S-5°N", loc="left", fontsize=fontz * 0.85)
    ax2.grid(alpha=0.3)
for ax in axes[:, 0]:
    ax.set_yscale("log")
    ax.set_ylim(100, 1)
    ax.set_yticks([100, 50, 30, 10, 3, 1], ["100", "50", "30", "10", "3", "1"])
    ax.set_ylabel("hPa")
    ax.set_xlim(*qbo_span)
for ax in axes[-1]:
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
fig.colorbar(cf, ax=axes[:, 0], label="m s$^{-1}$", shrink=0.6, pad=0.01)
save(fig, "qbo", YEAR)
plt.show()

# QBO phase before each forecast winter: member QBO_LEVEL wind over QBO_MONTHS
qbo_rows = []
for init in INIT_MONTHS:
    qt = mm_qbo[init].sel(plev=QBO_LEVEL).assign_coords(time=monthly[init]["time"]).swap_dims(lead="time")
    for w in (YEAR, YEAR + 1):
        sel = [t for t in pd.DatetimeIndex(qt["time"].values) if t.year == w and t.month in QBO_MONTHS]
        if len(sel) == len(QBO_MONTHS):
            v = qt.sel(time=sel).mean("time")
            qbo_rows.append({"start": MONTH_NAMES[init], "winter": f"{w}/{str(w + 1)[2:]}",
                             f"U{QBO_LEVEL:g} mean (m/s)": float(v.mean()), "members easterly": int((v < 0).sum()),
                             "members": int(v.size)})
e_qbo_on = {w: float(e_qbo.sel(plev=QBO_LEVEL).sel(time=[np.datetime64(f"{w}-{m:02d}-15") for m in QBO_MONTHS]).mean())
            for w in range(ERA5_FIRST_YEAR, YEAR + 1)
            if pd.Timestamp(f"{w}-{max(QBO_MONTHS):02d}-15") <= obs_end}
print(f"QBO {QBO_LEVEL:g} hPa over months {QBO_MONTHS}; ERA5 last observed autumn {max(e_qbo_on)}: {e_qbo_on[max(e_qbo_on)]:+.1f} m/s")
display(pd.DataFrame(qbo_rows).round(1))