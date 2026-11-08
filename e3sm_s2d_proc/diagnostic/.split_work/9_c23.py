w = YEAR
groups = {}                                        # (init, case) -> DataFrame of member values
for init in INIT_MONTHS:
    if (init, w) not in m_djf:
        continue
    vortex = m_djf[(init, w)].to_series()
    vortex.index = [str(i) for i in vortex.index]
    nino = {}
    for name, prefix in POOL.items():
        anom = rc.hindcast_anomaly(prefix, init, "TS", year=YEAR, method=CLIM_METHOD, mode=derivation_mode)
        if anom is None:
            continue
        n34 = rc.box_mean(anom, rc.NINO_REGIONS["Nino3.4"]).load()
        n34 = n34.assign_coords(time=("lead", mid_month(n34["time"].values))).swap_dims(lead="time")
        for mem in n34["member"].values:
            v = djf_mean(n34.sel(member=mem), w)
            nino[f"{name}:{mem}"] = np.nan if v is None else float(v)
    q = mm_qbo[init].assign_coords(time=monthly[init]["time"]).swap_dims(lead="time")
    qt = [np.datetime64(f"{w}-{m:02d}-15") for m in QBO_MONTHS]
    qon = q.sel(time=qt).mean("time") if all(t in q["time"].values for t in qt) else None
    df = pd.DataFrame({"vortex": vortex, "nino34": pd.Series(nino).reindex(vortex.index)})
    for p in (30.0, 50.0):
        s = qon.sel(plev=p).to_series() if qon is not None else pd.Series(np.nan, index=vortex.index)
        s.index = [str(i) for i in s.index]
        df[f"qbo{p:g}"] = s.reindex(vortex.index)
    df["case"] = [case_of(m) for m in df.index]
    for case, sub in df.groupby("case"):
        groups[(init, case)] = sub

summary = pd.DataFrame([{"start": MONTH_NAMES[i], "case": c, "members": len(g),
                         "DJF Niño3.4 (°C)": g["nino34"].mean(), "ON U30 (m/s)": g["qbo30"].mean(),
                         "ON U50 (m/s)": g["qbo50"].mean(), "members U50 < 0": int((g["qbo50"] < 0).sum()),
                         "DJF vortex anomaly (m/s)": g["vortex"].mean()} for (i, c), g in groups.items()])
display(summary.round(2))
e_qbo_last = {p: float(e_qbo.sel(plev=p).sel(time=[np.datetime64(f"{max(e_qbo_on)}-{m:02d}-15") for m in QBO_MONTHS]).mean())
              for p in (30.0, 50.0)}

fig, axes = plt.subplots(2, 2, figsize=(16, 10), gridspec_kw={"hspace": 0.45, "wspace": 0.2})
keys = list(groups)
glabel = [f"{MONTH_NAMES[i]}\n{c.replace('Reanalysis-', '')}" for i, c in keys]
gcol = [INIT_COLORS[i] for i, _ in keys]
x = np.arange(len(keys))
ax = axes[0, 0]
for xi, k, c in zip(x, keys, gcol):
    v = groups[k]["nino34"].dropna().values
    ax.scatter(np.full(v.size, xi) + rng.uniform(-0.12, 0.12, v.size), v, color=c, s=16, alpha=0.7)
    ax.plot([xi - 0.3, xi + 0.3], [v.mean()] * 2, color=c, lw=3)
ax.axhline(ENSO_THRESHOLD, color="0.5", ls=":", lw=1)
ax.set_xticks(x, glabel, fontsize=fontz * 0.7)
ax.set_ylabel("DJF Niño3.4 anomaly (°C)")
ax.set_title(f"(a) ENSO state, DJF {w}/{str(w + 1)[2:]}", loc="left")
ax.grid(alpha=0.3, axis="y")
ax = axes[0, 1]
for p, dx, mk in ((30.0, -0.15, "s"), (50.0, 0.15, "o")):
    for xi, k, c in zip(x, keys, gcol):
        v = groups[k][f"qbo{p:g}"].dropna().values
        ax.scatter(np.full(v.size, xi + dx) + rng.uniform(-0.05, 0.05, v.size), v, color=c, marker=mk, s=16, alpha=0.7)
        ax.plot([xi + dx - 0.12, xi + dx + 0.12], [v.mean()] * 2, color=c, lw=3)
    ax.axhline(e_qbo_last[p], color="k", ls="--" if p == 30 else "-", lw=1.2,
               label=f"ERA5 Oct-Nov {max(e_qbo_on)} U{p:g}: {e_qbo_last[p]:+.0f} m/s")
ax.axhline(0, color="0.4", lw=0.8)
ax.set_xticks(x, glabel, fontsize=fontz * 0.7)
ax.set_ylabel("Oct-Nov U, 5°S-5°N (m s$^{-1}$)")
ax.set_title(f"(b) QBO state, Oct-Nov {w} (squares 30 hPa, circles 50 hPa)", loc="left")
ax.legend(frameon=False, fontsize=fontz * 0.65)
ax.grid(alpha=0.3, axis="y")

reg_rows = []
for ax, (xkey, xl) in zip(axes[1], [("nino34", "Niño3.4 − group mean (°C)"), ("qbo50", f"U{QBO_LEVEL:g} − group mean (m s$^{{-1}}$)")]):
    allx, ally = [], []
    for k, c in zip(keys, gcol):
        g = groups[k].dropna(subset=["vortex", xkey])
        if len(g) < 3:
            continue
        dx, dy = g[xkey] - g[xkey].mean(), g["vortex"] - g["vortex"].mean()
        ax.scatter(dx, dy, color=c, s=30, edgecolor="k", lw=0.3)
        allx += list(dx)
        ally += list(dy)
    if len(allx) > 3:
        allx, ally = np.asarray(allx), np.asarray(ally)
        r = np.corrcoef(allx, ally)[0, 1]
        b = np.polyfit(allx, ally, 1)
        xx = np.linspace(allx.min(), allx.max(), 10)
        ax.plot(xx, np.polyval(b, xx), "k-", lw=2, label=f"within-group r={r:+.2f}, slope {b[0]:+.2f} (n={allx.size})")
        reg_rows.append({"predictor": xkey, "set": f"E3SM {w}/{str(w + 1)[2:]}, within start x case", "n": allx.size,
                         "r": r, "slope": b[0]})
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel(xl)
    ax.set_ylabel("DJF U(60°N, 10 hPa) anomaly − group mean (m s$^{-1}$)")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=fontz * 0.7, loc="lower left")
axes[1, 0].set_title("(c) Member covariation with ENSO", loc="left")
axes[1, 1].set_title(f"(d) Member covariation with the QBO ({QBO_LEVEL:g} hPa)", loc="left")
handles = [Line2D([], [], color=INIT_COLORS[i], marker="o", ls="", label=f"{MONTH_NAMES[i]} start") for i in dict.fromkeys(i for i, _ in keys)]
fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, bbox_to_anchor=(0.5, 0.0))
save(fig, "background", "enso", "qbo", YEAR)
plt.show()
display(pd.DataFrame(reg_rows).round(2))
