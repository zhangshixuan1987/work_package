PANELS = [("u60", "raw", "U(60°N, 10 hPa)", "m s$^{-1}$"),
          ("t10", "raw", "T(60-90°N, 10 hPa)", "K"),
          ("z10", "anom", "Polar-cap (60-90°N) height anomaly, 10 hPa", "m"),
          ("slp", "std", "Polar-cap (60-90°N) sea-level pressure anomaly (surface)", "σ")]
span = slice(*INDEX_SPAN)
span_days = pd.date_range(*INDEX_SPAN, freq="D").values


def smooth(da):
    return da.rolling(time=SMOOTH_DAYS, center=True, min_periods=1).mean() if SMOOTH_DAYS > 1 else da


fig, axes = plt.subplots(len(PANELS), 1, figsize=(15, 3.3 * len(PANELS)), sharex=True, gridspec_kw={"hspace": 0.18})
for ax, (key, kind, title, units) in zip(axes, PANELS):
    mean, std = e_iclim[key]
    # winters (DJF) shaded lightly
    for y in range(pd.Timestamp(INDEX_SPAN[0]).year - 1, pd.Timestamp(INDEX_SPAN[1]).year + 1):
        ax.axvspan(np.datetime64(f"{y}-12-01"), np.datetime64(f"{y + 1}-03-01"), color="tab:blue", alpha=0.04, lw=0)
    if kind == "raw":
        c, sd = on_doy(mean, span_days), on_doy(std, span_days)
        ax.fill_between(span_days, c - sd, c + sd, color="0.85", lw=0)
        ax.plot(span_days, c, color="0.45", lw=1.2)
    elif kind == "anom":
        sd = on_doy(std, span_days)
        ax.fill_between(span_days, -sd, sd, color="0.85", lw=0)
    else:
        ax.axhspan(-1, 1, color="0.88", lw=0)
    if key == "u60" or kind != "raw":
        ax.axhline(0, color="k", lw=0.8)
    for init, m in m_idx.items():
        color = INIT_COLORS[init]
        s = m["raw" if kind == "raw" else "anom"][key]
        if kind == "std":
            s = s / on_doy(std, s["time"].values)
        s = smooth(s).sel(time=span)
        if not s.sizes["time"]:
            continue
        t = s["time"].values
        q = s.quantile(list(INDEX_RANGE), "member")
        ax.fill_between(t, q.isel(quantile=0), q.isel(quantile=1), color=color, alpha=0.16, lw=0)
        ens = s.mean("member")
        ax.plot(t, ens, color=color, lw=2.3)
        ax.plot(t[0], float(ens[0]), "o", color=color, ms=7, mec="k", mew=0.6, zorder=5)
        if kind == "raw":
            ax.plot(t, m["clim"][key].sel(time=span), color=color, lw=1.1, ls="--", alpha=0.7)
        elif m["anom_alt"] is not None and key in m["anom_alt"]:
            a = m["anom_alt"][key]
            if kind == "std":
                a = a / on_doy(std, a["time"].values)
            ax.plot(t, smooth(a).sel(time=span).mean("member"), color=color, lw=1.6, ls=":")
    o = e_idx[key].sel(time=span)
    if o.sizes["time"]:
        if kind != "raw":
            o = o - on_doy(mean, o["time"].values)
        if kind == "std":
            o = o / on_doy(std, o["time"].values)
        ax.plot(o["time"].values, smooth(o), color="k", lw=2.2, zorder=6)
    ax.set_ylabel(units)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.grid(alpha=0.25)

# SSW central dates: one row of ticks per start, ERA5 in black
ax = axes[0]
rows = [("ERA5", "k", [d for d in obs_ssw if pd.Timestamp(INDEX_SPAN[0]) <= d <= pd.Timestamp(INDEX_SPAN[1])])]
rows += [(f"{MONTH_NAMES[i]} start", INIT_COLORS[i], m["ssw"]) for i, m in m_idx.items()]
for k, (label, color, dates) in enumerate(rows):
    y = 0.03 + 0.05 * k
    ax.plot(dates, np.full(len(dates), y), "|", color=color, ms=9, mew=1.6, transform=ax.get_xaxis_transform(), zorder=7)
n_ssw = ", ".join(f"{label}: {len(d)}" for label, _, d in rows[1:])
ax.text(0.995, 0.97, f"SSWs (ticks): {n_ssw}", transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})

axes[-1].set_xlim(np.datetime64(INDEX_SPAN[0]), np.datetime64(INDEX_SPAN[1]))
axes[-1].xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 4, 7, 10)))
axes[-1].xaxis.set_minor_locator(mdates.MonthLocator())
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
handles = [Line2D([], [], color="k", lw=2.2, label="ERA5"),
           Patch(color="0.85", label=f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]} ±1σ")]
handles += [Line2D([], [], color=INIT_COLORS[i], lw=2.3, marker="o", mec="k", ms=6,
                   label=f"1 {MONTH_NAMES[i]} {YEAR} start ({m_idx[i]['raw']['u60'].sizes['member']})") for i in m_idx]
handles += [Line2D([], [], color="0.3", lw=1.6, ls=":", label=f"Feb/Aug mean, {'lead_drift' if CLIM_METHOD == 'direct' else 'target_month'} climatology"),
            Patch(color="0.5", alpha=0.3, label=f"member {INDEX_RANGE[0]:.0%}-{INDEX_RANGE[1]:.0%}"),
            Line2D([], [], color="0.4", lw=1.1, ls="--", label="E3SM lead-time climatology")]
axes[-1].legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=4, frameon=False, fontsize=10)
fig.suptitle(f"Stratospheric polar vortex, E3SM RTENSO2026 {YEAR} starts and ERA5 (through {obs_end:%d %b %Y}); "
             f"{SMOOTH_DAYS}-day running means", x=0.125, ha="left", y=0.995, fontsize=15, fontweight="bold")
print(rc.save_figure(fig, NOTEBOOK, "RTENSO2026", f"y{YEAR}", "vortex_indices", "plume", dpi=fig_dpi))
plt.show()