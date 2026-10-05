# Climate variability with Python (snapshot)

Snapshot of [royalosyin/Python-Practical-Application-on-Climate-Variability-Studies](https://github.com/royalosyin/Python-Practical-Application-on-Climate-Variability-Studies),
reorganized for easier maintenance; see `UPSTREAM.txt` for the source commit.
The original upstream introduction follows the index below.

## Layout

| Path | Contents |
|------|----------|
| `notebooks/` | Tutorial notebooks `exNN_<topic>.ipynb`; run them from this folder |
| `cvtools/` | Helper modules (wavelet, SSA, SPI, harmonic analysis) imported by the notebooks |
| `data/` | Input data and intermediate `.npz` files (paths in notebooks: `../data/...`) |
| `figures/` | Plot output (`../figures/...`); image files are git-ignored |

Changes from upstream: files renamed to snake_case and grouped into the folders
above; Windows paths (`data\x`, `image\x`) replaced with portable `../data/x`
and `../figures/x`; helper scripts moved into the `cvtools` package
(`HA_helpers` -> `cvtools.harmonic`, `mySSA` -> `cvtools.ssa`,
`dim_spi_n` -> `cvtools.spi`); Python 2 `print` statements and removed NumPy
aliases (`np.complex`) fixed in `cvtools`. Notebook code is otherwise unchanged
and some notebooks still use Python 2 or Basemap.

## Notebook index

| Topic | Notebooks |
|-------|-----------|
| Basics and data I/O | ex00_intro_python, ex01_sst_netcdf_subsample, ex02_nino3_index, ex03_sst_land_mask_global_mean, ex04_nino3_ssta_plot |
| Fields, interpolation and maps | ex05_uwind_mean_std, ex06_uwind_zonal_mean_interp, ex07_interp_regular_irregular_grids, ex08_sst_monthly_climatology, ex09_sst_map_projections, ex22_flood_hazard_map |
| Precipitation and drought | ex11_gpcc_precip, ex12_gpcp_precip, ex13_africa_rainfall_hovmoller, ex14_spi_drought_index |
| Trends and time series | ex10_acw_hovmoller, ex15_sst_trend_anomaly, ex23_co2_vs_temperature_anomaly, ex24_co2_time_series, ex25_temperature_anomaly_heatmap, ex26_marine_heatwaves, ex27_wind_rose, ex35_sea_ice_extent |
| Correlation and coupled modes | ex16_nino3_slp_correlation, ex28_mca_slp_sst, ex29_cca_slp_sst, ex34_soi_correlations |
| EOF and decomposition | ex17_eof_hgt500, ex18_eof_global_sst, ex19_eof_central_pacific_sst, ex33_eemd_ne_pacific_sst, ex36_ssa_nh_land_temperature |
| Spectral and wavelet analysis | ex20_power_spectral_density, ex21_wavelet_nino3, ex31_harmonic_analysis_temperature |
| Weather regimes / clustering | ex30_weather_regimes_kmeans, ex32_weather_regimes_som |

ex35 (was `SeaIce.ipynb`) and ex36 (was `Singular Spectrum Analysis for NH Monthly
Land Temperature.ipynb`) were unnumbered upstream.

## Data not included

Upstream ships only small inputs. Download these into `data/` before running
the listed notebooks (mostly NOAA PSL gridded products):

| File | Used by |
|------|---------|
| `skt.mon.mean.nc` | ex01, ex02, ex03, ex07, ex09 |
| `sst.mnmean.nc`, `lsmask.nc` | ex08, ex15, ex18 |
| `sst.mnmean.v5.nc` | ex33 |
| `uwnd3.mon.mean.nc` | ex05, ex06 |
| `precip.mon.total.v7.nc` | ex11 |
| `V22_GPCP.1979-2010.nc` | ex12, ex14 |
| `chirps-v2.0.2016.days_p25.nc` | ex13 |
| `slp.mon.mean.1970.1999.nc` | ex16 |
| `hgt500.mon.mean.nc` | ex17 |
| `eof_data/sst.sw.AVHRR.l4.1982.2000.nc` | ex19 |
| `Hazard_AUS__1000.grd` | ex22 |
| `slp.mnmean.hadslp2.nc`, `sst.mon.anom.kaplan.nc` | ex28, ex29 |
| `z500.DJF.anom.1979.2010.nc` | ex30, ex32 |
| `precip.mon.mean.nc`, `prmsl.mon.mean.nc`, `air.sig995.mon.mean.nc` | ex34 |
| `NH.Ts.csv` | ex36 |

`skt.so.mon.mean.npz` and `annualmeanpr.npz` / `monthlylmeanpr.npz` are written
by ex01 and ex11. `ssta.nino3.30y.npz` (ex02 output) is included.

---

## Upstream introduction

### Python & Practical Application on Climate Variability Studies
Main objective of this tutorial is the transference of know-how in practical applications and management of statistical tools commonly used to explore meteorological time series, focusing on applications to study issues related with the climate variability and climate change. This tutorial starts with some basic statistic for time series analysis as estimation of means, anomalies, standard deviation, correlations, arriving the estimation of particular climate indexes (Niño 3), detrending single time series and decomposition of time series, filtering, interpolation of climate variables on regular or irregular grids, leading modes of climate variability (EOF or HHT), signal processing in the climate system (spectral and wavelet analysis). In addition, this tutorial also deals with different data formats such as CSV, NetCDF, Binary, and matlab'mat, etc. It is assumed that you have basic knowledge and understanding of statistics and Python.

## Generic libraries for scientific analysis
The default Python library for dealing with large arrays of numeric data (e.g. four dimensional latitude/longitude/altitude/time data arrays) is numpy, while the netCDF data format is commonly used in Atmospheric science and Oceanography as it is convinient to store various variables of many dimensions. The default for reading and writing netCDF files is netCDF4 (the capability to read and write text files, including .csv, is built into numpy).
### What is NetCDF?
NetCDF is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data.
NetCDF was developed and is maintained at Unidata. Unidata provides data and software tools for use in geoscience education and research. The NetCDF homepage may be found at http://www.unidata.ucar.edu/software/netcdf/. The NetCDF source-code is hosted at GitHub, and may be found directly at http://github.com/Unidata/netcdf-c.
### How to deal with NetCDF and other data with Python?
we mainly use netCDF4-python, NumPy and SciPy to process NetCDF and other data formats.
* netCDF4-python
> netcdf4-python is a Python interface to the netCDF C library. netCDF version 4 has many features not found in earlier versions of the library and is implemented on top of HDF5. This module can read and write files in both the new netCDF 4 and the old netCDF 3 format, and can create files that are readable by HDF5 clients. The API modelled after Scientific.IO.NetCDF, and should be familiar to users of that module (see more http://unidata.github.io/netcdf4-python/).
* NumPy
> NumPy is the fundamental package for scientific computing with Python (see more http://www.numpy.org/). 
It contains among other things:
a powerful N-dimensional array object
sophisticated (broadcasting) functions
tools for integrating C/C++ and Fortran code
useful linear algebra, Fourier transform, and random number capabilities
Besides its obvious scientific uses, NumPy can also be used as an efficient multi-dimensional container of generic data. Arbitrary data-types can be defined. This allows NumPy to seamlessly and speedily integrate with a wide variety of databases.
* SciPy
> SciPy is a collection of mathematical algorithms and convenience functions built on the Numpy extension of Python (https://www.scipy.org/index.html). It adds significant power to the interactive Python session by providing the user with high-level commands and classes for manipulating and visualizing data. With SciPy an interactive Python session becomes a data-processing and system-prototyping environment rivaling systems such as MATLAB, IDL, Octave, R-Lab, and SciLab.
>The additional benefit of basing SciPy on Python is that this also makes a powerful programming language available for use in developing sophisticated programs and specialized applications. Scientific applications using SciPy benefit from the development of additional modules in numerous niches of the software landscape by developers across the world. Everything from parallel programming to web and data-base subroutines and classes have been made available to the Python programmer. All of this power is available in addition to the mathematical libraries in SciPy.
## Data Sources
We will use the data publicly available as possible.
The data are mainly downloaded from https://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.derived.surfaceflux.html
## Visualization
### Matplotlib
Matplotlib is a Python 2D plotting library which produces publication quality figures in a variety of hardcopy formats and interactive environments across platforms. Matplotlib can be used in Python scripts, the Python and IPython shell, the jupyter notebook, web application servers, and four graphical user interface toolkits.
Matplotlib tries to make easy things easy and hard things possible. You can generate plots, histograms, power spectra, bar charts, errorcharts, scatterplots, etc., with just a few lines of code. For a sampling, see the screenshots, thumbnail gallery, and examples directory
For simple plotting the pyplot module provides a MATLAB-like interface, particularly when combined with IPython. For the power user, you have full control of line styles, font properties, axes properties, etc, via an object oriented interface or via a set of functions familiar to MATLAB users.
See more from https://matplotlib.org/
### Basemap
Basemap is a great tool for creating maps using python in a simple way. It’s a matplotlib extension, so it has got all its features to create data visualizations, and adds the geographical projections and some datasets to be able to plot coast lines, countries, and so on directly from the library.
Basemap has got some documentation, but some things are a bit more difficult to find. I started this documentation to extend a little the original documentation and examples, but it grew a little, and now covers many of the basemap possibilities.
See more from https://basemaptutorial.readthedocs.io/en/latest/
## More
We will mainly apply these mostly generic libraries to carry out data analysis step by step. The procedures or steps are common to the atmospheric and ocean sciences. Although other advanced libraries will simpilfy the procedures, the underlying ideas should be the same in essense.
We will also introduce more other libraries such as xarray and iris in the following parts.
