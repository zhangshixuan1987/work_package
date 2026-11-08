preview_day = frame_days[len(frame_days) // 2]
for init in fc:                       # a featured member's SSW day, if any, is more telling
    if "first SSW" in fc[init]["why"]:
        preview_day = pd.Timestamp(fc[init]["why"].split("(")[1].split(")")[0])
        break
fig, state = build_figure()
draw_frame(state, preview_day)
plt.show()