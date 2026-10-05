"""Small helpers for saved figures: publishing to public_html and pixel comparison."""
import glob
import os
import shutil
import subprocess
import tempfile


def publish(*paths):
    """Make files/directories world-readable (o+r files, o+rx dirs) for the public_html web server.

    The default umask on Chrysalis (007) leaves new files unreadable by the web server.
    Directories are processed recursively. Only paths owned by the current user are changed.
    """
    def add(p, bits):
        try:
            st = os.stat(p)
            if os.getuid() == st.st_uid:
                os.chmod(p, st.st_mode | bits)
        except OSError:
            pass

    for p in paths:
        if os.path.isdir(p):
            add(p, 0o555)
            for root, dirs, files in os.walk(p):
                for d in dirs:
                    add(os.path.join(root, d), 0o555)
                for f in files:
                    add(os.path.join(root, f), 0o444)
        else:
            add(p, 0o444)


def compare_pdfs(new_pdf, ref_pdf, density=30):
    """Render both PDFs with ImageMagick and count differing pixels over all pages.

    Returns 0 when identical, a positive pixel count when they differ, or None if
    either file is missing or cannot be rendered (including when ImageMagick is not
    installed, e.g. on compute nodes).
    """
    if not (os.path.exists(new_pdf) and os.path.exists(ref_pdf)):
        return None
    if shutil.which("convert") is None or shutil.which("compare") is None:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        pages = {}
        for tag, pdf in (("a", ref_pdf), ("b", new_pdf)):
            r = subprocess.run(["convert", "-density", str(density), pdf,
                                os.path.join(tmp, f"{tag}-%03d.png")], capture_output=True)
            if r.returncode != 0:
                return None
            pages[tag] = sorted(glob.glob(os.path.join(tmp, f"{tag}-*.png")))
        if len(pages["a"]) != len(pages["b"]):
            return -1
        total = 0
        for a, b in zip(pages["a"], pages["b"]):
            r = subprocess.run(["compare", "-metric", "AE", a, b, "null:"], capture_output=True, text=True)
            # ImageMagick prints e.g. "5.37e+07 (820)"; the number in parentheses is the pixel count
            out = r.stderr.strip()
            count = out[out.find("(") + 1:out.find(")")] if "(" in out else out.split()[0]
            total += int(float(count))
        return total
