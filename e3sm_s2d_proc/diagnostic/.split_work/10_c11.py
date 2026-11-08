def build_figure():
    inits = [i for i in INIT_MONTHS if i in fc]
    ncol = 1 + len(inits)
    fig = plt.figure(figsize=FIG["size"])
    gs = fig.add_gridspec(4, ncol, height_ratios=FIG["height_ratios"], hspace=FIG["hspace"], wspace=FIG["wspace"],
                          left=FIG["left"], right=FIG["right"], top=FIG["top"], bottom=FIG["bottom"])
    maps = {}
    if obs_anom is None and pool_all is not None:        # winter not observed: E3SM all-start mean instead of ERA5
        maps["pool"] = (fig.add_subplot(gs[0, 0], projection=PROJ), f"E3SM all starts: mean of {pool_all['n']}", "k")
    else:
        maps["obs"] = (fig.add_subplot(gs[0, 0], projection=PROJ), "ERA5", "k")
    maps["clim"] = (fig.add_subplot(gs[1, 0], projection=PROJ), f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]} climatology", "0.3")
    for k, init in enumerate(inits, start=1):
        maps[("mean", init)] = (fig.add_subplot(gs[0, k], projection=PROJ),
                                f"1 {MONTH_NAMES[init]} {YEAR} start: mean of {fc[init]['n']}{' †' if fc[init]['surrogate'] else ''}",
                                INIT_COLORS[init])
        maps[("feat", init)] = (fig.add_subplot(gs[1, k], projection=PROJ),
                                f"member {fc[init]['featured']}{' †' if fc[init]['surrogate'] else ''}", INIT_COLORS[init])
    state = {"maps": {key: (ax, setup_map(ax, title, color)) for key, (ax, title, color) in maps.items()}}

    levels, cmap, norm = levels_cmap()
    cax = fig.add_subplot(gs[2, :])
    pos = cax.get_position()
    cax.set_position([pos.x0 + pos.width * 0.25, pos.y0, pos.width * 0.5, pos.height])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", extend="both",
                        ticks=np.arange(-STYLE["range"], STYLE["range"] + 1e-6, 4 * STYLE["ci"]))
    cbar.set_label(f"{STYLE['name']} ({STYLE['units']}); contours: total height every {STYLE['line']:g} m",
                   fontsize=FIG["label"])
    cbar.ax.tick_params(labelsize=FIG["tick"])

    # timeline: static faded series, dynamic solid parts up to the frame day
    ax = fig.add_subplot(gs[3, :])
    t_all = all_days.values
    ax.fill_between(t_all, u60_clim_days - u60_std_days, u60_clim_days + u60_std_days, color="0.85", lw=0)
    ax.plot(t_all, u60_clim_days, color="0.45", lw=1.2)
    ax.axhline(0, color="k", lw=0.8)
    lines = {}
    for init in inits:
        color = INIT_COLORS[init]
        series = {}
        for key, lw, ls in (("u60_mean", FIG["ens_lw"], "-"), ("u60_feat", FIG["member_lw"], "--")):
            s = fc[init][key].sel(time=slice(all_days[0], all_days[-1]))
            ax.plot(s["time"].values, s.values, color=color, lw=lw * 0.7, ls=ls, alpha=FIG["faded_alpha"])
            series[key] = (s, ax.plot([], [], color=color, lw=lw, ls=ls, zorder=4)[0])
        lines[init] = series
        if init in m_idx:                      # the start's own lead-time climatology: the model's weak vortex
            c = m_idx[init]["clim"]["u60"].sel(time=slice(all_days[0], all_days[-1]))
            ax.plot(c["time"].values, c.values, color=color, lw=1.2, ls=":", alpha=0.8)
        for d in fc[init]["feat_ssw"]:         # SSWs of the featured member
            ax.plot([d], [0], marker="v", color=color, ms=9, mec="k", mew=0.6, zorder=6)
    e_part = e_u60.sel(time=slice(all_days[0], all_days[-1]))
    obs_line = ax.plot([], [], color="k", lw=FIG["obs_lw"], zorder=5)[0]
    cursor = ax.axvline(all_days[0], color="0.2", lw=1.2, ls="--", zorder=3)
    ax.set_xlim(all_days[0], all_days[-1])
    lo = min([-15.0] + [float(v[0].min()) for s in lines.values() for v in s.values()] + ([float(e_part.min())] if e_part.size else []))
    hi = max([45.0] + [float(v[0].max()) for s in lines.values() for v in s.values()] + ([float(e_part.max())] if e_part.size else []))
    ax.set_ylim(np.floor(lo / 5) * 5, np.ceil(hi / 5) * 5)
    ax.set_ylabel("U(60°N, 10 hPa) (m s$^{-1}$)", fontsize=FIG["label"])
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.tick_params(labelsize=FIG["tick"])
    ax.grid(alpha=0.25)
    handles = [Line2D([], [], color="k", lw=FIG["obs_lw"], label="ERA5"),
               Patch(color="0.85", label=f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]} ±1σ")]
    handles += [Line2D([], [], color=INIT_COLORS[i], lw=FIG["ens_lw"], label=f"1 {MONTH_NAMES[i]} start, mean") for i in inits]
    handles += [Line2D([], [], color="0.3", lw=FIG["member_lw"], ls="--", label="featured member"),
                Line2D([], [], color="0.3", lw=1.2, ls=":", label="E3SM lead-time climatology"),
                Line2D([], [], color="0.3", marker="v", ls="", mec="k", label="featured-member SSW")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4, fontsize=FIG["legend"],
              frameon=False)
    counts = "; ".join(f"{MONTH_NAMES[i]} start {fc[i]['n_ssw']} of {fc[i]['n']}" for i in inits)
    ax.text(0.005, 0.97, f"members with an SSW in {WINTER_LABEL}: {counts}", transform=ax.transAxes, va="top",
            fontsize=FIG["legend"], bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})
    if CLIM_METHOD == "direct":
        fig.text(FIG["left"], 0.935, "Anomalies from each start's Reanalysis climatology: May 1981-2010 (30 years), "
                 "Feb and Aug 1985-2010 (the six starts of every fifth year); totals add them to the ERA5 climatology",
                 fontsize=FIG["legend"] - 0.5, color="0.3")
    elif any(fc[i]["surrogate"] for i in inits):
        fig.text(FIG["left"], 0.935, "† anomalies relative to the lead_drift climatology (no Reanalysis start in that month): "
                 "their sign is not robust to the drift correction (notebook 9, section 6)", fontsize=FIG["legend"] - 0.5, color="0.3")
    fig.suptitle(f"Stratospheric polar vortex, winter {WINTER_LABEL}: 10-hPa geopotential height", x=FIG["left"],
                 ha="left", y=0.975, fontsize=FIG["suptitle"], fontweight="bold")
    state.update(lines=lines, e_part=e_part, obs_line=obs_line, cursor=cursor,
                 day_text=fig.text(FIG["right"], 0.975, "", ha="right", va="top", fontsize=FIG["day"], fontweight="bold"))
    return fig, state


def draw_frame(state, t):
    t = pd.Timestamp(t)
    state["day_text"].set_text(f"{t:%d %B %Y}")
    for key, (ax, mstate) in state["maps"].items():
        if key == "obs":
            if obs_anom is None or value_at(obs_anom, t) is None:
                fill_map(ax, mstate, empty=f"Not yet observed\n(ERA5 ends {obs_end:%d %b %Y})")
            else:
                fill_map(ax, mstate, value_at(obs_anom, t), value_at(z_obs, t),
                         **u60_badge(float(e_u60.sel(time=t)) if t <= obs_end else None))
        elif key == "pool":
            anom = value_at(pool_all["anom"], t)
            u = value_at(pool_all["u60"], t)
            fill_map(ax, mstate, anom, value_at(pool_all["total"], t), **u60_badge(None if u is None else float(u)))
        elif key == "clim":
            fill_map(ax, mstate, None, value_at(z_clim, t), **u60_badge(float(u60_clim_days.sel(time=t)), "U60 clim"))
        else:
            kind, init = key
            f = fc[init]
            anom = value_at(f[f"{kind}_anom"], t)
            if anom is None:
                fill_map(ax, mstate, empty="Before start" if t < pd.Timestamp(YEAR, init, 1) else "No output")
                continue
            u = value_at(f["u60_mean" if kind == "mean" else "u60_feat"], t)
            fill_map(ax, mstate, anom, value_at(f[f"{kind}_total"], t), **u60_badge(None if u is None else float(u)))
    tn = np.datetime64(t, "ns")
    for init, series in state["lines"].items():
        for s, line in series.values():
            part = s.sel(time=slice(None, tn))
            line.set_data(part["time"].values, part.values)
    part = state["e_part"].sel(time=slice(None, tn))
    state["obs_line"].set_data(part["time"].values, part.values)
    state["cursor"].set_xdata([t, t])
    return []