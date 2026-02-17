#!/bin/csh

#TC track data from the International Best Track Archive for Climate Stewardship (IBTrACS) (Knapp et al., 2010)
#Data are obtained for the 40-year satellite period of 1979-2018 online (https://www.ncdc.noaa.gov/ibtracs/)

set workdir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/OBS/tc_track"

cd $workdir 

rm -rvf *.nc 

wget --recursive --no-parent https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/netcdf/


