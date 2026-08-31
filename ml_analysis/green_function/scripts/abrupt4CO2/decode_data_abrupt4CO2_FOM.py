'''
Description: Python script using PyNGL Python module

 - contour plot on map (rectilinear data)

shixuan.zhang@pnnl.gov
'''

import argparse
from pathlib import Path

import numpy as np
import xarray as xr
import pickle
import pandas as pd

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

target  = "Anomalies_FOM"
DEFAULT_OUTPUT_ROOT = Path(
    "/global/cfs/cdirs/m1867/zhan391/paper_work/green_function"
)
#There are 13 variables in total 
#The following variables have monthly mean data (150 year x 12 months) 
varList1 = ['TS','toa','sfc','LHFLX','SHFLX','PRECT']
#The following variables have annual mean data (150 year) 
varList2 = ['EIS','R_plk','R_alb','R_lr','R_q','R_SWcld','R_LWcld']

def map_plot(lats,lons,var):
    #map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
    #map =  Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
    map = Basemap(projection='robin',lon_0=0,resolution='c')

    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.fillcontinents(color='coral',lake_color='aqua')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='aqua')
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30))
    map.drawparallels(np.arange(-90,90,30))

    # compute native map projection coordinates of lat/lon grid.
    #x, y = map(lons*180./np.pi, lats*180./np.pi)
    x, y = map(*np.meshgrid(lons,lats))

    # contour data over the map.
    #cs = map.contour(x,y,var,15,linewidths=1.5)
    cs = map.contourf(x,y,var,30,cmap=plt.cm.jet)
    map.drawcoastlines() # draw coastlines
    map.drawmapboundary() # draw a line around the map region
    map.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
    map.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
    map.colorbar()

    plt.title('contour plot for selected varibale')
    plt.show()

#-----------------------------------------------------------------------
#-- Function: main
#-----------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="External directory containing the input pickle and CAM grid file.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="External root for generated outputs, figures, and logs.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    data_dir = args.data_dir.expanduser().resolve()
    run_root = args.output_root.expanduser() / "abrupt4CO2"
    for directory in (run_root / "outputs", run_root / "figures", run_root / "logs"):
        directory.mkdir(parents=True, exist_ok=True)

    input_file = data_dir / "Anomalies_of_abrupt4CO2_in_FOM.pkl"
    with input_file.open("rb") as stream:
        mydict = pickle.load(stream)
    #print(mydict[target])
    #for var in varList1[:]:
    #    print(var, mydict['Anomalies_FOM'][var].shape)
    #for var in varList2[:]:
    #    print(var, mydict['Anomalies_FOM'][var].shape)

    #construct the time and grid info for CAM finite volume (FV) grid (96 x 144).
    grid  = xr.open_dataset(data_dir / "cam_grid_info.nc", decode_times=False)
    lat   = grid['lat'].data
    lon   = grid['lon'].data
    nyr   = 150 # note: total number of years
    nmo   = 12  # note: total number of month 
    year  = np.arange(1, nyr+1, 1, dtype=int)
    month = np.arange(1, nmo+1, 1, dtype=int)
    #construct xarray for the extracted variable 
    #note: 
    ds = xr.Dataset(
            data_vars=dict(
                TS=(["year", "month", "lat","lon"], mydict['Anomalies_FOM']['TS']),
                toa=(["year", "month", "lat","lon"], mydict['Anomalies_FOM']['toa']),
                sfc=(["year", "month", "lat","lon"], mydict['Anomalies_FOM']['sfc']), 
                LHFLX=(["year", "month", "lat","lon"], mydict['Anomalies_FOM']['LHFLX']),
                SHFLX=(["year", "month", "lat","lon"], mydict['Anomalies_FOM']['SHFLX']),
                PRECT=(["year", "month", "lat","lon"], mydict['Anomalies_FOM']['PRECT']),
                EIS=(["year", "lat", "lon"], mydict['Anomalies_FOM']['EIS']),
                R_plk=(["year", "lat", "lon"], mydict['Anomalies_FOM']['R_plk']),
                R_alb=(["year", "lat", "lon"], mydict['Anomalies_FOM']['R_alb']),
                R_lr=(["year", "lat", "lon"], mydict['Anomalies_FOM']['R_lr']),
                R_q=(["year", "lat", "lon"], mydict['Anomalies_FOM']['R_q']),
                R_SWcld=(["year", "lat", "lon"], mydict['Anomalies_FOM']['R_SWcld']),
                R_LWcld=(["year", "lat", "lon"], mydict['Anomalies_FOM']['R_LWcld']),
                ),
            coords=dict(
                lon=(["lon"], lon),
                lat=(["lat"], lat),
                year=(["year"], year),
                month=(["month"], month),
                ),
            attrs=dict(description="CESM1 fully-coupled model Abrupt4×CO2 simulation"),)

    # check the map distribution for randomly select year and month 
    map_plot(lat,lon,ds.TS[0,0,:,:])

if __name__ == '__main__':
    main()
