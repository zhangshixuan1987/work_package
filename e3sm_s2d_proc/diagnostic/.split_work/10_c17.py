def storyboard():
    inits = [i for i in INIT_MONTHS if i in fc]
    columns = [("pool", f"All starts mean ({pool_all['n']})", "k")] if obs_anom is None and pool_all is not None else [("obs", "ERA5", "k")]
    for init in inits:
        columns += [(("mean", init), f"{MONTH_NAMES[init]} start mean", INIT_COLORS[init]),
                    (("feat", init), f"{MONTH_NAMES[init]} start, member {fc[init]['featured'].split(':')[-1]}",
                     INIT_COLORS[init])]
    rows = [pd.Timestamp(d) for d in STORYBOARD]
    fig, axes = plt.subplots(len(rows), len(columns), figsize=(2.3 * len(columns), 2.45 * len(rows) + 1.0),
                             subplot_kw={"projection": PROJ}, squeeze=False, gridspec_kw={"hspace": 0.22, "wspace": 0.04})
    for r, t in enumerate(rows):
        for c, (key, title, color) in enumerate(columns):
            ax = axes[r, c]
            mstate = setup_map(ax, title if r == 0 else "", color)
            mstate["badge"].set_fontsize(FIG["badge"] - 2)
            if key == "pool":
                u = value_at(pool_all["u60"], t)
                fill_map(ax, mstate, value_at(pool_all["anom"], t), value_at(pool_all["total"], t),
                         **u60_badge(None if u is None else float(u)))
            elif key == "obs":
                if obs_anom is None or value_at(obs_anom, t) is None:
                    fill_map(ax, mstate, empty="Not yet\nobserved")
                else:
                    fill_map(ax, mstate, value_at(obs_anom, t), value_at(z_obs, t), **u60_badge(float(e_u60.sel(time=t))))
            else:
                kind, init = key
                anom = value_at(fc[init][f"{kind}_anom"], t)
                u = value_at(fc[init]["u60_mean" if kind == "mean" else "u60_feat"], t)
                fill_map(ax, mstate, anom, value_at(fc[init][f"{kind}_total"], t),
                         **u60_badge(None if u is None else float(u)), **({} if anom is not None else {"empty": "Before start"}))
            if c == 0:
                ax.text(-0.12, 0.5, f"{t:%d %b %Y}", transform=ax.transAxes, rotation=90, ha="center", va="center",
                        fontsize=FIG["label"], fontweight="bold")
    levels, cmap, norm = levels_cmap()
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, orientation="horizontal", extend="both",
                        fraction=0.02, pad=0.03, shrink=0.5, aspect=40,
                        ticks=np.arange(-STYLE["range"], STYLE["range"] + 1e-6, 4 * STYLE["ci"]))
    cbar.set_label(f"{STYLE['name']} ({STYLE['units']}); contours: total every {STYLE['line']:g} m", fontsize=FIG["label"])
    fig.suptitle(f"Polar vortex at 10 hPa, winter {WINTER_LABEL}: E3SM ensemble means and featured members",
                 fontsize=FIG["suptitle"] - 2, y=0.92)
    return fig


fig = storyboard()
print(rc.save_figure(fig, NOTEBOOK, "RTENSO2026", f"winter{WINTER}{WINTER + 1}", "Z10", "storyboard", dpi=fig_dpi))
plt.show()