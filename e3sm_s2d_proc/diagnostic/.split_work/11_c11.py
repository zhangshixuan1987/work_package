%%time
video_ext, video = vp.video_writer(FPS)
writers = {"gif": ("gif", lambda: animation.PillowWriter(fps=FPS)), "video": (video_ext, video)}
print(f"video: {video_ext or 'no ffmpeg video encoder; skipped'}")
sequence = list(frame_days) + [frame_days[-1]] * HOLD
fig, state = build_figure()
anim = animation.FuncAnimation(fig, lambda t: draw_frame(state, t), frames=sequence, blit=False, cache_frame_data=False)
anim_files = {}
for fmt in ANIM_FORMATS:
    ext, writer = writers[fmt]
    if writer is None:
        continue
    path = FIGDIR / figure_filename("RTENSO2026", f"winter{YEAR}{YEAR + 1}", "Z10", "evolution", ext=ext)
    anim.save(path, writer=writer(), dpi=ANIM_DPI)
    anim_files[ext] = path
    print(f"{path} ({path.stat().st_size / 1e6:.1f} MB)")
plt.close(fig)
if SHOW_INLINE and "gif" in anim_files:
    display(Image(filename=str(anim_files["gif"])))