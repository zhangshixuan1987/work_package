def storyboard():
    columns = (["era5"] if SHOW_REFERENCE else []) + [("mem", i) for i in INITS] + ([("mean", i) for i in INITS] if SHOW_ENSMEAN else [])
    rows = [pd.Timestamp(d) for d in STORYBOARD]
    fig, axes = plt.subplots(len(rows), len(columns), figsize=(2.25 * len(columns), 2.4 * len(rows) + 1.2),
                             subplot_kw={"projection": vp.PROJ}, squeeze=False, gridspec_kw={"hspace": 0.2, "wspace": 0.04})
    for r, t in enumerate(rows):
        for c, key in enumerate(columns):
            ax = axes[r, c]
            title, color = panel_title(key)
            mstate = vp.setup_map(ax, title.replace("1 ", "").replace(" start:", "") if r == 0 else "", lat_min=LAT_MIN,
                                  color=color, title_size=FIG["title"] - 2, badge_size=FIG["badge"] - 2,
                                  dashed_frame=key != "era5" and key[0] == "mean")
            draw_map(ax, mstate, key, t)
            if c == 0:
                ax.text(-0.13, 0.5, f"{t:%d %b %Y}", transform=ax.transAxes, rotation=90, ha="center", va="center",
                        fontsize=FIG["label"], fontweight="bold")
    if SHADING:
        levels, cmap, norm = vp.shading(STYLE)
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, orientation="horizontal", extend="both",
                            fraction=0.02, pad=0.03, shrink=0.5, aspect=40,
                            ticks=np.arange(-STYLE["range"], STYLE["range"] + 1e-6, 4 * STYLE["ci"]))
        cbar.set_label(f"{STYLE['name']} ({STYLE['units']}); contours: raw E3SM 10-hPa height every {STYLE['line']:g} m",
                       fontsize=FIG["label"])
    fig.suptitle(f"E3SM RTENSO2026: {EVENT_LABEL} at 10 hPa (fixed members; dashed circles: ensemble means)",
                 fontsize=FIG["suptitle"] - 2, y=0.92)
    return fig


fig = storyboard()
print(rc.save_figure(fig, NOTEBOOK, "RTENSO2026", f"winter{YEAR}{YEAR + 1}", "Z10", "evolution", "storyboard", dpi=fig_dpi))
plt.show()