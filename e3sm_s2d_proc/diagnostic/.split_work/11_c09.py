fig, state = build_figure()
draw_frame(state, PREVIEW_DAY)
path = FIGDIR / figure_filename("RTENSO2026", f"winter{YEAR}{YEAR + 1}", "Z10", "evolution", "preview",
                                pd.Timestamp(PREVIEW_DAY).strftime("%Y%m%d"), ext="png")
fig.savefig(path, dpi=fig_dpi)
print(path)
plt.show()