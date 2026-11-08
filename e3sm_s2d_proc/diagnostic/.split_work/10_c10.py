FIG = {
    "size": (14.0, 9.6),
    "height_ratios": (1.0, 1.0, 0.045, 0.62),     # maps, maps, colorbar, timeline
    "hspace": 0.30, "wspace": 0.04,
    "left": 0.04, "right": 0.98, "top": 0.89, "bottom": 0.08,
    "suptitle": 16, "day": 20, "title": 10.5, "label": 10.5, "tick": 9.5, "legend": 9.5, "badge": 9.5,
    "coast_lw": 0.5, "empty": "0.93", "ens_lw": 2.6, "member_lw": 1.2, "obs_lw": 2.4, "faded_alpha": 0.22,
}
PROJ = ccrs.NorthPolarStereo(central_longitude=-90)    # North America at the bottom-left, Eurasia top-right
_theta = np.linspace(0, 2 * np.pi, 120)
CIRCLE = mpath.Path(np.vstack([np.sin(_theta), np.cos(_theta)]).T * 0.5 + [0.5, 0.5])


def levels_cmap():
    levels = np.arange(-STYLE["range"], STYLE["range"] + STYLE["ci"] / 2, STYLE["ci"])
    cmap = mycolors.blue2red_cmap(levels) if STYLE["cmap"] == "blue2red" else plt.get_cmap(STYLE["cmap"], len(levels) - 1)
    return levels, cmap, BoundaryNorm(levels, ncolors=cmap.N)


def setup_map(ax, title, color="k"):
    # extent as a square in projected coordinates reaching LAT_MIN: set_extent([-180, 180, ...]) collapses
    # the y range of a polar-stereographic map in cartopy 0.26
    x, y = PROJ.transform_point(0.0, LAT_MIN, ccrs.PlateCarree())
    r = float(np.hypot(x, y))
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_boundary(CIRCLE, transform=ax.transAxes)
    ax.add_feature(cfeature.LAND, facecolor="0.88", zorder=0.5)
    ax.coastlines(lw=FIG["coast_lw"], color="0.35", zorder=3)
    ax.gridlines(crs=ccrs.PlateCarree(), lw=0.3, color="0.5", alpha=0.5, xlocs=range(-180, 181, 45), ylocs=[40, 60, 80])
    ax.set_title(title, loc="center", fontsize=FIG["title"], color=color, fontweight="bold")
    return {"artists": [],
            "empty": ax.text(0.5, 0.5, "", transform=ax.transAxes, ha="center", va="center", fontsize=FIG["label"],
                             color="0.35", zorder=6, bbox={"facecolor": "white", "edgecolor": "0.6", "boxstyle": "round,pad=0.4"}),
            "badge": ax.text(0.5, -0.04, "", transform=ax.transAxes, ha="center", va="top", fontsize=FIG["badge"])}


def fill_map(ax, state, anom=None, total=None, empty="", badge="", badge_color="k"):
    # shade ``anom`` and contour ``total`` (lat, lon); a grey panel with ``empty`` when both are None
    for artist in state["artists"]:
        artist.remove()
    state["artists"] = []
    state["badge"].set_text(badge)
    state["badge"].set_color(badge_color)
    if anom is None and total is None:
        ax.set_facecolor(FIG["empty"])
        state["empty"].set_text(empty)
        state["empty"].set_visible(True)
        return
    ax.set_facecolor("white")
    state["empty"].set_visible(False)
    if anom is not None:
        levels, cmap, norm = levels_cmap()
        values, lon = add_cyclic_point(anom.values, coord=anom["lon"].values)
        state["artists"].append(ax.contourf(lon, anom["lat"].values, values, levels=levels, cmap=cmap, norm=norm,
                                            extend="both", transform=ccrs.PlateCarree(), zorder=1))
    if total is not None:
        values, lon = add_cyclic_point(total.values, coord=total["lon"].values)
        lines = np.arange(np.floor(np.nanmin(values) / STYLE["line"]) * STYLE["line"], np.nanmax(values) + STYLE["line"], STYLE["line"])
        state["artists"].append(ax.contour(lon, total["lat"].values, values, levels=lines, colors="k", linewidths=0.8,
                                           transform=ccrs.PlateCarree(), zorder=2))


def value_at(da, t):
    t = np.datetime64(pd.Timestamp(t), "ns")
    return da.sel(time=t) if t in da["time"].values else None


def u60_badge(value, prefix="U60"):
    # keyword arguments (badge, badge_color) of fill_map
    if value is None or not np.isfinite(value):
        return {"badge": "", "badge_color": "k"}
    return {"badge": f"{prefix}(10 hPa) {value:+.0f} m/s" + ("  easterly" if value < 0 else ""),
            "badge_color": "tab:red" if value < 0 else "k"}