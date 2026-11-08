# SSW frequency, central-date months and final-warming dates
fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), gridspec_kw={"width_ratios": [1.3, 1, 1.1], "wspace": 0.3})
ax = axes[0]
sel = freq[freq["set"].str.contains("ERA5") | freq["case"].astype(str).str.startswith("pooled")]
labels = [f"{s}\n({k})" if k != "ERA5" else s for s, k in zip(sel["set"], sel["wind"])]
colors = ["k" if k == "ERA5" else ("tab:red" if k == "raw" else "tab:blue") for k in sel["wind"]]
y = np.arange(len(sel))[::-1]
for yy, (_, r), c in zip(y, sel.iterrows(), colors):   # one call per point: errorbar takes one colour
    ax.errorbar(r["events per winter"], yy, xerr=[[r["events per winter"] - r["5%"]], [r["95%"] - r["events per winter"]]],
                fmt="o", color=c, ms=8, capsize=4, zorder=3)
ax.set_yticks(y, labels, fontsize=fontz * 0.6)
ax.axvline(float(freq.iloc[0]["events per winter"]), color="k", ls=":", lw=1)
ax.set_xlabel("SSWs per winter (90% bootstrap)")
ax.set_title("(a) SSW frequency", loc="left")
ax.grid(alpha=0.3, axis="x")

ax = axes[1]
order = [11, 12, 1, 2, 3]
x = np.arange(len(order))
obs_m = obs_events["central_date"].dt.month.value_counts().reindex(order, fill_value=0) / len(obs_winters)
mod = pd.concat([model_ssw[i] for i in INIT_MONTHS])
n_mw = sum(len(c) for (i, w, k), c in model_counts.items() if k == "raw")
mod_m = (mod["central_date"].dt.month.value_counts().reindex(order, fill_value=0) / max(n_mw, 1)) if len(mod) else 0 * obs_m
ax.bar(x - 0.2, obs_m.values, 0.4, color="k", label=f"ERA5 ({len(obs_winters)} winters)")
ax.bar(x + 0.2, np.asarray(mod_m), 0.4, color="tab:red", label=f"E3SM raw ({n_mw} member-winters)")
ax.set_xticks(x, [MONTH_NAMES[m] for m in order])
ax.set_ylabel("SSWs per winter")
ax.set_title("(b) Month of the central date", loc="left")
ax.legend(frameon=False)
ax.grid(alpha=0.3, axis="y")

ax = axes[2]


def doy_spring(d):
    return (pd.Timestamp(d) - pd.Timestamp(pd.Timestamp(d).year, 1, 1)).days + 1


# only complete records enter the boxes; censored and unresolved springs (record ends before 31 Jul,
# e.g. the May start's 2028) are tabulated but not plotted
obs_fw_status = {y: sc.final_warming(e_u60.values, e_u60["time"].values, y, recover_days=SSW["recover_days"])
                 for y in range(ERA5_FIRST_YEAR, obs_end.year + 1)}
obs_fw_doy = [doy_spring(d) for d, st in obs_fw_status.values() if st == "complete"]
data, labels, cols = [obs_fw_doy], [f"ERA5\n{ERA5_FIRST_YEAR}-{obs_end.year - 1}"], ["0.6"]
fw_rows = [{"source": "ERA5", "member": "ERA5", "spring": y, "final_warming": d, "status": st}
           for y, (d, st) in obs_fw_status.items()]
for init in INIT_MONTHS:
    u = m_u60[init]
    for spring in (YEAR, YEAR + 1, YEAR + 2):
        res = [sc.final_warming(u.sel(member=m).values, u["time"].values, spring, recover_days=SSW["recover_days"])
               for m in u["member"].values]
        if all(st == "no data" for _, st in res) or (spring == YEAR and init != 2):
            continue
        fw_rows += [{"source": f"E3SM init{init:02d}", "member": str(m), "spring": spring, "final_warming": d, "status": st}
                    for m, (d, st) in zip(u["member"].values, res)]
        ok = [d for d, st in res if st == "complete"]
        n_cens = sum(st in ("censored", "unresolved") for _, st in res)
        if n_cens:
            print(f"[final warming] {MONTH_NAMES[init]} start, spring {spring}: {n_cens} of {len(res)} members censored or "
                  f"unresolved (record ends {pd.Timestamp(u['time'].values[-1]):%d %b %Y}); not plotted")
        if len(ok) >= 3:
            data.append([doy_spring(d) for d in ok])
            labels.append(f"{MONTH_NAMES[init]} start\n{spring} (n={len(ok)})")
            cols.append(INIT_COLORS[init])
bp = ax.boxplot(data, vert=False, patch_artist=True, widths=0.6, whis=(5, 95))
for patch, c in zip(bp["boxes"], cols):
    patch.set_facecolor(c)
    patch.set_alpha(0.5)
d_obs, st_obs = obs_fw_status.get(YEAR, (pd.NaT, "no data"))
if pd.notna(d_obs):
    note = "" if st_obs == "complete" else f" ({st_obs}: ERA5 ends {obs_end:%d %b})"
    ax.axvline(doy_spring(d_obs), color="k", ls="--", lw=1.5, label=f"ERA5 {YEAR}: {d_obs:%d %b}{note}")
    ax.legend(frameon=False, loc="lower right", fontsize=fontz * 0.7)
ax.set_yticks(np.arange(1, len(labels) + 1), labels, fontsize=fontz * 0.7)
ticks = [pd.Timestamp(2001, m, 1) for m in (2, 3, 4, 5, 6)]
ax.set_xticks([doy_spring(t) for t in ticks], [f"1 {t:%b}" for t in ticks])
ax.set_title("(c) Final-warming date", loc="left")
ax.grid(alpha=0.3, axis="x")
save(fig, "ssw", "frequency", "final_warming", YEAR)
plt.show()

fw_table = pd.DataFrame(fw_rows)
fw_table.to_csv(DATA_DIR / f"final_warming_dates_{YEAR}.csv", index=False)
print(f"[WRITE] {DATA_DIR / f'final_warming_dates_{YEAR}.csv'}")