# Event-level distributions and the vertical evolution table
def event_values(comp, full, lag_range, plev=None):
    c = comp.sel(lag=slice(*lag_range))
    if plev is not None:
        c = c.sel(plev=plev)
    return c.mean("lag").values[full]


cat_raw = events_all[events_all["wind"].isin(["ERA5", "raw"])]
cat_raw = cat_raw.assign(group=np.where(cat_raw["source"] == "ERA5", "ERA5", "E3SM"))
metrics = {
    f"pre-onset [v*T*] 100 hPa, days {PRE_WINDOW[0]} to {PRE_WINDOW[1]} (K m/s)":
        (event_values(e_comp["vt"], e_full, PRE_WINDOW), event_values(m_comp["vt"], m_full, PRE_WINDOW)),
    "minimum U(60°N, 10 hPa) (m/s)": (cat_raw.query("group == 'ERA5'")["min_u"].values, cat_raw.query("group == 'E3SM'")["min_u"].values),
    "easterly days": (cat_raw.query("group == 'ERA5'")["easterly_days"].values, cat_raw.query("group == 'E3SM'")["easterly_days"].values),
    "days to recovery": (cat_raw.query("group == 'ERA5'")["days_to_recovery"].dropna().values,
                         cat_raw.query("group == 'E3SM'")["days_to_recovery"].dropna().values),
    "100-hPa polar-cap height, days 15-44 (σ)":
        (event_values(e_comp["pch"], e_full, (15, 44), 100.0), event_values(m_comp["pch"], m_full, (15, 44), 100.0)),
}
fig, axes = plt.subplots(1, len(metrics), figsize=(18, 4.2), gridspec_kw={"wspace": 0.45})
for ax, (name, (eo, mo)) in zip(axes, metrics.items()):
    bp = ax.boxplot([eo[np.isfinite(eo)], mo[np.isfinite(mo)]], patch_artist=True, widths=0.6, whis=(10, 90))
    for patch, c in zip(bp["boxes"], ("0.6", "tab:red")):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    ax.set_xticks([1, 2], [f"ERA5\n(n={np.isfinite(eo).sum()})", f"E3SM\n(n={np.isfinite(mo).sum()})"])
    ax.set_title(name, fontsize=fontz * 0.7, loc="left", wrap=True)
    ax.grid(alpha=0.3, axis="y")
fig.suptitle("SSW event properties (unweighted ERA5; whiskers 10-90%)", x=0.12, ha="left", y=1.04)
save(fig, "ssw", "event_properties", YEAR)
plt.show()

rows = []
for label, comp, full, w in (("ERA5, season-matched", e_comp, e_full, e_w), ("E3SM", m_comp, m_full, ones_m)):
    for win in PCH_WINDOWS:
        row = {"set": label, "days": f"{win[0]} to {win[1]}"}
        for p in PCH_LEVELS:
            v = comp["pch"].sel(plev=p, lag=slice(*win)).mean("lag").values[full]
            row[f"Z {p:g} hPa (σ)"] = np.nansum(v * w[full]) / w[full].sum()
        v = comp["slp"].sel(lag=slice(*win)).mean("lag").values[full]
        row["SLP (σ)"] = np.nansum(v * w[full]) / w[full].sum()
        rows.append(row)
vertical = pd.DataFrame(rows)
vertical.to_csv(DATA_DIR / f"ssw_vertical_evolution_{YEAR}.csv", index=False)
print("Polar-cap (60-90°N) anomalies around the SSW central dates (composite means; positive SLP: negative NAM):")
display(vertical.round(2))
