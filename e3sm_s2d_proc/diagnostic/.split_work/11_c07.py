FIG = {
    "size": (14.0, 10.0 if SHOW_ENSMEAN else 7.4),
    "hspace": 0.32, "wspace": 0.04,
    "left": 0.05, "right": 0.98, "top": 0.89, "bottom": 0.07,
    "suptitle": 16, "day": 20, "title": 10.5, "label": 10.5, "tick": 9.5, "legend": 9.5, "badge": 9.5,
    "member_lw": 2.4, "mean_lw": 1.8, "faded_alpha": 0.22, "range_alpha": 0.15,
}
INITS = [i for i in INIT_MONTHS if i in fc]
ROWS = [("mem", "member")] + ([("mean", "ensemble mean")] if SHOW_ENSMEAN else [])


def map_kwargs():
    return {"style": STYLE if SHADING else None, "total_levels": STYLE["line_levels"]}


def lead_label(init, t):
    m = (t.year - YEAR) * 12 + t.month - init
    return f"lead {m}"


def panel_fields(key, t):
    # (anom, total, u60) of a map panel on day t; None where the panel has no data
    t = np.datetime64(pd.Timestamp(t), "ns")
    if key == "era5":
        return None, z_clim_era5.sel(time=t), float(u60_era5_clim.sel(time=t))
    if key == "e3sm":
        return None, z_clim_e3sm.sel(time=t) if t in z_clim_e3sm["time"].values else None, \
            float(fc[ref_init]["u60_clim"].sel(time=t)) if t in fc[ref_init]["u60_clim"]["time"].values else None
    row, init = key
    f = fc[init]
    if t not in f["raw_mem"]["time"].values:
        return None, None, None
    raw, anom = (f["raw_mem"], f["anom_mem"]) if row == "mem" else (f["raw_mean"], f["anom_mean"])
    u = f["u60_mem"] if row == "mem" else f["u60_mean"]
    return (anom.sel(time=t) if (SHADING and anom is not None) else None), raw.sel(time=t), float(u.sel(time=t))


def panel_title(key):
    if key == "era5":
        return f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]} climatology", "0.25"
    if key == "e3sm":
        return f"E3SM climatology ({MONTH_NAMES[ref_init]} starts)", "0.25"
    row, init = key
    if row == "mem":
        return f"1 {MONTH_NAMES[init]} start: {fc[init]['member'].replace('Reanalysis-', '')}", INIT_COLORS[init]
    return f"1 {MONTH_NAMES[init]} start: mean of {fc[init]['n']}", INIT_COLORS[init]


def draw_map(ax, state, key, t):
    anom, total, u = panel_fields(key, t)
    if total is None:
        fill_map_empty = "Before start" if key not in ("era5", "e3sm") and pd.Timestamp(t) < pd.Timestamp(YEAR, key[1], 1) else "No data"
        vp.fill_map(ax, state, empty=fill_map_empty)
        return
    lead = lead_label(key[1], pd.Timestamp(t)) if key not in ("era5", "e3sm") else ""
    prefix = "U60 clim" if key in ("era5", "e3sm") else ("U60 mean" if key[0] == "mean" else "U60")
    vp.fill_map(ax, state, anom, total, lead=lead, **map_kwargs(), **vp.u60_badge(u, prefix))


def build_figure():
    ncol = len(INITS) + (1 if SHOW_REFERENCE else 0)
    nrow = len(ROWS)
    fig = plt.figure(figsize=FIG["size"])
    ratios = [1.0] * nrow + [0.045, 0.62]
    gs = fig.add_gridspec(nrow + 2, ncol, height_ratios=ratios, hspace=FIG["hspace"], wspace=FIG["wspace"],
                          left=FIG["left"], right=FIG["right"], top=FIG["top"], bottom=FIG["bottom"])
    maps = {}
    for r, (row, row_label) in enumerate(ROWS):
        c0 = 0
        if SHOW_REFERENCE:
            key = "era5" if r == 0 else "e3sm"
            maps[key] = (fig.add_subplot(gs[r, 0], projection=vp.PROJ), False)
            c0 = 1
        for k, init in enumerate(INITS):
            maps[(row, init)] = (fig.add_subplot(gs[r, c0 + k], projection=vp.PROJ), row == "mean")
    state = {"maps": {}}
    for key, (ax, dashed) in maps.items():
        title, color = panel_title(key)
        state["maps"][key] = (ax, vp.setup_map(ax, title, lat_min=LAT_MIN, color=color, title_size=FIG["title"],
                                               badge_size=FIG["badge"], label_size=FIG["label"], dashed_frame=dashed))
    for r, (row, row_label) in enumerate(ROWS):
        y = gs[r, 0].get_position(fig)
        fig.text(0.012, (y.y0 + y.y1) / 2, "one member\n(a realization)" if row == "mem" else "ensemble mean\n(not a realization)",
                 rotation=90, ha="center", va="center", fontsize=FIG["label"], fontweight="bold",
                 color="k" if row == "mem" else "0.4")

    # colour bar (fixed levels) or a note on the contours
    cax = fig.add_subplot(gs[nrow, :])
    pos = cax.get_position()
    cax.set_position([pos.x0 + pos.width * 0.25, pos.y0, pos.width * 0.5, pos.height])
    if SHADING:
        levels, cmap, norm = vp.shading(STYLE)
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", extend="both",
                            ticks=np.arange(-STYLE["range"], STYLE["range"] + 1e-6, 4 * STYLE["ci"]))
        cbar.set_label(f"{STYLE['name']} from each start's {CLIM_METHOD} climatology ({STYLE['units']}); "
                       f"contours: raw E3SM 10-hPa height every {STYLE['line']:g} m", fontsize=FIG["label"])
        cbar.ax.tick_params(labelsize=FIG["tick"])
    else:
        cax.axis("off")
        cax.text(0.5, 0.5, f"contours: raw E3SM 10-hPa height every {STYLE['line']:g} m", ha="center", va="center",
                 fontsize=FIG["label"], transform=cax.transAxes)

    # U(60N, 10 hPa): faded full series (static), solid up to the frame day (dynamic)
    ax = fig.add_subplot(gs[nrow + 1, :])
    t_all = days.values
    ax.fill_between(t_all, u60_era5_clim - u60_era5_std, u60_era5_clim + u60_era5_std, color="0.88", lw=0, zorder=0)
    ax.plot(t_all, u60_era5_clim, color="0.5", lw=1.2, zorder=1)
    ax.axhline(0.0, color="k", lw=1.0, zorder=2)
    ax.text(0.995, 0.0, "easterly below 0 ", transform=ax.get_yaxis_transform(), ha="right", va="top",
            fontsize=FIG["tick"] - 1, color="0.3")
    lines = {}
    for init in INITS:
        f, color = fc[init], INIT_COLORS[init]
        sl = slice(days[0], days[-1])
        lo, hi = f["u60_lo"].sel(time=sl), f["u60_hi"].sel(time=sl)
        ax.fill_between(lo["time"].values, lo, hi, color=color, alpha=FIG["range_alpha"] * 0.5, lw=0)
        ax.plot(f["u60_clim"].sel(time=sl)["time"].values, f["u60_clim"].sel(time=sl), color=color, lw=1.0, ls=":", alpha=0.7)
        series = {}
        for key, lw, ls in (("u60_mem", FIG["member_lw"], "-"), ("u60_mean", FIG["mean_lw"], "--")):
            s = f[key].sel(time=sl)
            ax.plot(s["time"].values, s.values, color=color, lw=lw * 0.7, ls=ls, alpha=FIG["faded_alpha"])
            series[key] = (s, ax.plot([], [], color=color, lw=lw, ls=ls, zorder=4)[0])
        series["dot"] = ax.plot([], [], "o", color=color, ms=7, mec="k", mew=0.6, zorder=6)[0]
        if SHOW_SSW:
            for tt in f["ssw"]:
                ax.plot([tt], [0.0], marker="v", color=color, ms=9, mec="k", mew=0.6, zorder=6)
        lines[init] = series
    cursor = ax.axvline(days[0], color="0.2", lw=1.2, ls="--", zorder=3)
    ax.set_xlim(days[0], days[-1])
    ax.set_ylim(*U60_YLIM)
    ax.set_ylabel("U(60°N, 10 hPa) (m s$^{-1}$)", fontsize=FIG["label"])
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
    ax.tick_params(labelsize=FIG["tick"])
    ax.grid(alpha=0.25)
    handles = [Line2D([], [], color="0.3", lw=FIG["member_lw"], label="member (map row 1)"),
               Line2D([], [], color="0.3", lw=FIG["mean_lw"], ls="--", label="ensemble mean"),
               Patch(color="0.6", alpha=0.3, label="member 10-90%"),
               Line2D([], [], color="0.3", lw=1.0, ls=":", label="E3SM climatology"),
               Patch(color="0.88", label=f"ERA5 {CLIM_YEARS[0]}-{CLIM_YEARS[1]} ±1σ")]
    handles += [Line2D([], [], color=INIT_COLORS[i], lw=FIG["member_lw"], label=f"1 {MONTH_NAMES[i]} {YEAR} start") for i in INITS]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=4, fontsize=FIG["legend"], frameon=False)

    fig.suptitle(f"E3SM RTENSO2026: {EVENT_LABEL} at 10 hPa", x=FIG["left"], ha="left", y=0.975,
                 fontsize=FIG["suptitle"], fontweight="bold")
    fig.text(FIG["left"], 0.94, f"Fixed members (not selected for events): "
             + ", ".join(f"{MONTH_NAMES[i]} {fc[i]['member'].replace('Reanalysis-', '')}" for i in INITS)
             + (f"; shading: anomaly from each start's {CLIM_METHOD} Reanalysis climatology" if SHADING else ""),
             fontsize=FIG["legend"], color="0.3")
    day_text = fig.text(FIG["right"], 0.975, "", ha="right", va="top", fontsize=FIG["day"], fontweight="bold")
    state.update(lines=lines, cursor=cursor, day_text=day_text)
    return fig, state


def draw_frame(state, t):
    t = pd.Timestamp(t)
    state["day_text"].set_text(f"{t:%d %B %Y}")
    for key, (ax, mstate) in state["maps"].items():
        draw_map(ax, mstate, key, t)
    tn = np.datetime64(t, "ns")
    for init, series in state["lines"].items():
        for key in ("u60_mem", "u60_mean"):
            s, line = series[key]
            part = s.sel(time=slice(None, tn))
            line.set_data(part["time"].values, part.values)
        s = series["u60_mem"][0]
        if tn in s["time"].values:
            series["dot"].set_data([t], [float(s.sel(time=tn))])
    state["cursor"].set_xdata([t, t])
    return []