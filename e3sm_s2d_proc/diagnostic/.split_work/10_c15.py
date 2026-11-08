%%time
def video_writer():
    # (extension, writer factory) for the best video codec the available ffmpeg has, or (None, None)
    try:
        import imageio_ffmpeg   # bundled ffmpeg with H.264
        plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        pass
    if not animation.writers.is_available("ffmpeg"):
        return None, None
    encoders = subprocess.run([plt.rcParams["animation.ffmpeg_path"], "-hide_banner", "-encoders"],
                              capture_output=True, text=True).stdout
    even = ["-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p"]
    if "libx264" in encoders:
        return "mp4", lambda: animation.FFMpegWriter(fps=FPS, codec="libx264", extra_args=even + ["-crf", "20"])
    if "libvpx-vp9" in encoders:
        return "webm", lambda: animation.FFMpegWriter(fps=FPS, codec="libvpx-vp9", extra_args=even + ["-crf", "30", "-b:v", "0"])
    return None, None


video_ext, video = video_writer()
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
    path = FIGDIR / figure_filename("RTENSO2026", f"winter{WINTER}{WINTER + 1}", "Z10", "polar_vortex", ext=ext)
    anim.save(path, writer=writer(), dpi=ANIM_DPI)
    anim_files[ext] = path
    print(f"{path} ({path.stat().st_size / 1e6:.1f} MB)")
plt.close(fig)
if SHOW_INLINE and "gif" in anim_files:
    display(Image(filename=str(anim_files["gif"])))