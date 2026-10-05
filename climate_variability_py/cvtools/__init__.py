"""Helper modules used by the climate-variability notebooks.

harmonic         Fourier / harmonic analysis helpers (was HA_helpers.py)
ssa              Singular Spectrum Analysis class ``mySSA`` (was mySSA.py)
spi              Standardized Precipitation Index ``dim_spi_n`` (was dim_spi_n.py)
wavelet          Torrence & Compo continuous wavelet transform
wavelet_inverse  inverse wavelet reconstruction
wave_signif      wavelet significance testing
wave_bases       wavelet mother functions (used by ``wavelet``)

Modules are not imported here so that optional dependencies (e.g. lmoments
for ``spi``) are only needed by the notebooks that use them.
"""
