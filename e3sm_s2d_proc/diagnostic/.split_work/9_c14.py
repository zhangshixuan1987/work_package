w0 = YEAR - 1
span = (pd.Timestamp(w0, 10, 1), pd.Timestamp(YEAR, 6, 30))
e_obs = e_u60.sel(time=slice(*span))
cl = on_days(clim_u60, e_obs)
sd = on_days(std_u60, e_obs)
obs_ssw_w0 = [d for d in sc.ssw_central_dates(e_u60.values, e_u60["time"].values, **SSW) if sc.winter_of(d) == w0]
obs_fw = {y: sc.final_warming_date(e_u60.values, e_u60["time"].values, y, recover_days=SSW["recover_days"])
          for y in range(ERA5_FIRST_YEAR, obs_end.year + 1)}
fw_date, fw_status = sc.final_warming(e_u60.values, e_u60["time"].values, YEAR, recover_days=SSW["recover_days"])
fw_text = f"{fw_date:%d %b} ({fw_status}: ERA5 ends {obs_end:%d %b %Y})" if pd.notna(fw_date) else f"none yet ({fw_status})"
print(f"ERA5 winter {w0}/{str(YEAR)[2:]}: SSW central dates {[f'{d:%d %b %Y}' for d in obs_ssw_w0] or 'none'}; "
      f"final warming {fw_text}")

init = 2
has_feb = init in INIT_MONTHS
fig, axes = plt.subplots(3 if has_feb else 2, 1, figsize=(15, 12 if has_feb else 8.5),
                         gridspec_kw={"height_ratios": [1.3, 1, 1][: 3 if has_feb else 2], "hspace": 0.32}, sharex=True)
ax = axes[0]
t = e_obs["time"].values
ax.fill_between(t, cl - 2 * sd, cl + 2 * sd, color="0.9", label="ERA5 ±2σ")
ax.fill_between(t, cl - sd, cl + sd, color="0.78", label="ERA5 ±1σ")
ax.plot(t, cl, color="0.4", lw=1.5, label=f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]} mean")
if has_feb:
    u = m_u60[init].sel(time=slice(*span))
    for mem in u["member"].values:
        ax.plot(u["time"], u.sel(member=mem), color=INIT_COLORS[init], lw=0.6, alpha=0.45)
    ax.plot(u["time"], u.mean("member"), color=INIT_COLORS[init], lw=2.5, label=f"E3SM 1 Feb {YEAR} ({u.sizes['member']} members)")
    ax.plot(u["time"], m_u60_clim[init].sel(time=slice(*span)), "--", color=INIT_COLORS[init], lw=1.5,
            label="E3SM lead-time climatology")
ax.plot(t, e_obs, "k-", lw=2.5, label="ERA5")
for d in obs_ssw_w0:
    ax.axvline(d, color="k", ls=":", lw=1.2)
    ax.text(d, 55, f" SSW {d:%d %b}", fontsize=fontz * 0.7, va="top")
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel("m s$^{-1}$")
ax.set_ylim(-30, 60)
ax.set_title(f"(a) U(60°N, {REF_LEVEL:g} hPa), winter {w0}/{str(YEAR)[2:]}", loc="left")
# final-warming rug
if has_feb:
    fw_m = [sc.final_warming_date(m_u60[init].sel(member=mem).values, m_u60[init]["time"].values, YEAR,
                                  recover_days=SSW["recover_days"]) for mem in m_u60[init]["member"].values]
    fw_m = [d for d in fw_m if pd.notna(d)]
    ax.plot(fw_m, np.full(len(fw_m), -27), "|", color=INIT_COLORS[init], ms=14, mew=2, label="final warming (members)")
if pd.notna(obs_fw.get(YEAR, pd.NaT)):
    ax.plot([obs_fw[YEAR]], [-27], "kv", ms=11, label="final warming (ERA5)")
ax.legend(loc="upper right", ncol=2, fontsize=fontz * 0.65, frameon=False)
ax.grid(alpha=0.3)

pch_lev = np.array([-2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5])
ep = e_pch.sel(time=slice(*span)).sel(plev=slice(1, 1000))
ax = axes[1]
cf = ax.contourf(ep["time"], ep["plev"], ep.T, levels=pch_lev, cmap="RdBu_r", extend="both")
ax.set_title("(b) ERA5 polar-cap (60-90°N) height anomaly / σ", loc="left")
if has_feb:
    mp = m_pch[init].sel(time=slice(*span)).mean("member")
    ax = axes[2]
    ax.contourf(mp["time"], mp["plev"], mp.T, levels=pch_lev, cmap="RdBu_r", extend="both")
    ax.set_title(f"(c) E3SM 1 Feb {YEAR} ensemble mean, same quantity (ERA5 σ)", loc="left")
for ax in axes[1:]:
    ax.set_yscale("log")
    ax.set_ylim(1000, 1)
    ax.set_yticks([1000, 300, 100, 30, 10, 1], ["1000", "300", "100", "30", "10", "1"])
    ax.set_ylabel("hPa")
axes[-1].set_xlim(*span)
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%Y"))
fig.colorbar(cf, ax=axes[1:], label="σ", shrink=0.8, pad=0.01)
save(fig, "winter", f"{w0}{YEAR}", "u60n10_pch")
plt.show()