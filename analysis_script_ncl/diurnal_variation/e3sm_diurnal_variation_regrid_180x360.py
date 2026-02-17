'''
Description: Python script using PyNGL Python module

 - contour plot on map (rectilinear data)

shixuan.zhang@pnnl.gov
'''
from datetime import datetime
from scipy import stats
from sklearn.metrics import mean_squared_error
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import numpy as np
import xarray as xr
import Ngl, Nio, os

#-----------------------------------------------------------------------
#-- Function: main
#-----------------------------------------------------------------------
def main():

    diri = '/global/cscratch1/sd/zhan391/DARPA_project/post_process/trainning_data/FV_180x360' #-- data directory
    varName = 'T'
    varDesc = 'DTR (diurnal temperature range, ~S~o~N~C)'
    varUnit = '~S~o~N~C'
    expndg  = 'NDGUV' # nudged experiments for plotting
    expndg  = 'NDGUVT' # nudged experiments for plotting
    expndg  = 'NDGUVTQ' # nudged experiments for plotting

    # Create time slice from dates
    season     = 'July'
    ik         = 71 # bottom model layer
    ntpd       = 8  # 3-hourly data sample per day     
    start_time = '2010-07-01'
    end_time   = '2010-07-31'

    # Create slice variables subset domain
    time_slice = slice(start_time, end_time)
    region     = "CONUS"           #"NH_MidLat"
    regstr     = "CONUS (20-55N)"  #"NH_MidLat (40 - 60N)"
    minlat     = 20                # 40
    maxlat     = 55                # 60
    minlon     = 215               # 0
    maxlon     = 315               # 360
    cenlon     = 180
    lat_slice  = slice(minlat, maxlat) #slice(40, 60)
    lon_slice  = slice(minlon, maxlon) #slice(0, 360)

    # minimum and maximum range for plot
    minval  =  0
    maxval  =  10
    incval  =  0.5

    # data for the plot
    expname = []
    filname = []
    expname.append('ERA5')
    filname.append(diri+'/reference_data/ML_REF_ERA5_3hourly_200911302100-201011301800.nc')
    expname.append('CLIM')
    filname.append(diri+'/clim/E3SMv2_CLIM_3hourly_200911302100-201011301800.nc')

    if (expndg == 'NDGUV'):
        expname.append('NDGUV')
        filname.append(diri+'/after_nudging/ML_IN_E3SMv2_NDGUV_tau6_3hourly_200911302100-201011301800.nc')
    elif (expndg == 'NDGUVT'):
        expname.append('NDGUVT')
        filname.append(diri+'/after_nudging/ML_IN_E3SMv2_NDGUVT_tau6_3hourly_200911302100-201011301800.nc')
    else:
        expname.append('NDGUVTQ')
        filname.append(diri+'/after_nudging/ML_IN_E3SMv2_NDGUVTQ_tau6_3hourly_200911302100-201011301800.nc')

    #read the reference data
    print ('working on ' + expname[0])
    ds  = xr.open_dataset(filname[0])
    lev = ds.variables["lev"][:]

    # Get data, selecting time, level, lat/lon slicea
    data = ds[varName].sel(time=time_slice,
            lev=lev[ik],
            lat=lat_slice,
            lon=lon_slice)
    nmon = int(data.shape[0]/ntpd)
    nlat = data.shape[1]
    nlon = data.shape[2]

    tmpd = np.reshape(np.ravel(data),(nmon,ntpd,nlat,nlon))
    #print(tmpd.shape)
    var0 = np.mean(tmpd.max(axis=1) - tmpd.min(axis=1),axis=0)
    #print(var0)

    # Specify longitude values for chosen domain
    lats = data.lat.values
    lons = data.lon.values

    # Read the CLIM and Nudged simulation
    print ('working on ' + expname[1])
    ds1  = xr.open_dataset(filname[1])
    # Get data, selecting time, level, lat/lon slicea
    dat1 = ds1[varName].sel(time=time_slice,
              lev=lev[ik],
              lat=lat_slice,
              lon=lon_slice)

    tmpd = np.reshape(np.ravel(dat1),(nmon,ntpd,nlat,nlon))
    #print(tmpd.shape)
    var1 = np.mean(tmpd.max(axis=1) - tmpd.min(axis=1),axis=0)
    #print(var1)

    print ('working on ' + expname[2])
    ds2  = xr.open_dataset(filname[2])
    # Get data, selecting time, level, lat/lon slicea
    dat2 = ds2[varName].sel(time=time_slice,
              lev=lev[ik],
              lat=lat_slice,
              lon=lon_slice)

    tmpd = np.reshape(np.ravel(dat2),(nmon,ntpd,nlat,nlon))
    #print(tmpd.shape)
    var2 = np.mean(tmpd.max(axis=1) - tmpd.min(axis=1),axis=0)
    #print(var2)

    wkres = Ngl.Resources()     #-- generate an res object for workstation
    wkres.wkWidth  = 2500       #-- plot resolution 2500 pixel width
    wkres.wkHeight = 2500       #-- plot resolution 2500 pixel height
    wks_type       = 'png'      #-- graphics output type
    
    figname = "fig_diurnal_variation_"+region+"_"+expndg+"_"+varName+"_"+season+".png"
    wks  = Ngl.open_wks(wks_type,figname,wkres)  #-- open workstation
    plot = []
        
    #-- set resources
    res = Ngl.Resources()       #-- generate an resource object for plot
    res.nglDraw                   = False               #-- don't draw individual plots
    res.nglFrame                  = False               #-- don't advance frame

    #-- viewport resources
    #res.vpXF                      = 0.1       #-- start x-position of viewport
    #res.vpYF                      = 0.9       #-- start y-position of viewport
    #res.vpWidthF                  = 0.7       #-- width of viewport
    #res.vpHeightF                 = 0.7       #-- height of viewport

    fontsize = 0.014
    res.tmXBLabelFontHeightF      = fontsize
    res.tmYLLabelFontHeightF      = fontsize
    res.tmYRLabelFontHeightF      = fontsize
    res.tiMainFontHeightF         = fontsize 
    res.tiXAxisFontHeightF        = fontsize
    res.tiYAxisFontHeightF        = fontsize
    
    #-- labelbar resources
    res.lbLabelBarOn              = False
    res.lbOrientation             = "Horizontal";"vertical" 
    #res.lbBoxMinorExtentF         = 0.15
    res.pmLabelBarDisplayMode     = "Never" #-- turn on the label bar
    res.pmLabelBarHeightF         = 0.12 #0.65
    res.pmLabelBarWidthF          = 0.62 #0.1
    res.pmLabelBarParallelPosF    = 0.55
    res.pmLabelBarOrthogonalPosF  = -0.1
    #res.lbRightMarginF            = -0.28
    #res.lbBottomMarginF           = -0.28
    res.lbLabelFontHeightF        = fontsize * 0.8
    res.lbLabelPosition           = "bottom"       
    res.lbTitleString             = varDesc+" ("+varUnit+")" 
    res.lbTitleFontHeightF        = fontsize * 0.8
    #res.lbTitlePosition           = "Right"                           
    #res.lbTitleDirection          = "Across"                          
    #res.lbTitleAngleF             = 90.                               
    #res.lbTitleFontHeightF        = fontsize
    
    BorderThick = 5.0
    res.tmBorderThicknessF    = BorderThick
    res.tmXBMajorThicknessF   = BorderThick
    res.tmXBMinorThicknessF   = BorderThick*0.5
    res.tmYLMajorThicknessF   = BorderThick
    res.tmYLMinorThicknessF   = BorderThick*0.5
    res.tmXTMajorThicknessF   = BorderThick
    res.tmXTMinorThicknessF   = BorderThick*0.5
    res.tmYRMajorThicknessF   = BorderThick
    res.tmYRMinorThicknessF   = BorderThick*0.5

    #-- contour resources
    res.cnLinesOn                      = False
    res.cnLineLabelsOn                 = False  #-- turn off line labels
    res.cnInfoLabelOn                  = False   #-- turn off info label
    res.cnLineThicknessF               = 0.1
    res.cnLineColor                    = "Black"
    res.cnLineLabelDensityF            = 1.0
    res.cnLineLabelFontHeightF         = fontsize
    res.cnLineLabelBackgroundColor     = -1

    res.cnFillOn                       = True         #-- turn on contour fill
    res.cnFillPalette                  = "WhBlGrYeRe" #"BlGrYeOrReVi200"
    res.cnMissingValFillColor          = 'gray50'       #-- missing value fill color
    res.cnLevelSelectionMode           = "ManualLevels" #-- select manual level selection mode
    res.cnMinLevelValF                 = minval #-- minimum contour value
    res.cnMaxLevelValF                 = maxval #-- maximum contour value
    res.cnLevelSpacingF                = incval   #-- contour increment
    
    res.tmXTOn                         = False
    res.tmYROn                         = True
    res.tmYLLabelJust                  = "CenterRight"
    res.tmYRLabelJust                  = "CenterRight"
    res.tmYUseLeft                     = True
    res.tmYRLabelsOn                   = False
    
    res.tiXAxisString                  = ""
    res.tiYAxisString                  = ""
    res.tiXAxisOffsetYF                = 0.0
    res.tiYAxisOffsetXF                = 0.0
    res.tiMainString                   = "" 
    res.tiMainOffsetYF                 = 0.0
    
    res.mpFillOn                       = False
    res.mpDataBaseVersion              = "MediumRes" #-- map database
    res.mpOceanFillColor               = "Transparent"
    res.mpLandFillColor                = "Gray90"
    res.mpInlandWaterFillColor         = "Gray90"
    
    res.mpPerimOn                      = False
    res.mpPerimLineColor               = -1 #'transparent' 
    res.mpOutlineOn                    = True
    res.mpGridAndLimbOn                = True #False 
    res.mpGridLineColor                = -1 #'transparent' 
    res.mpGeophysicalLineThicknessF    = 2.0

    res.sfXArray                       = lons
    res.sfYArray                       = lats
    res.mpLimitMode                    = 'LatLon'#-- must be set using minLatF/maxLatF/minLonF/maxLonF
    res.mpMinLatF                      = minlat  #-- sub-region minimum latitude
    res.mpMaxLatF                      = maxlat  #-- sub-region maximum latitude
    res.mpMinLonF                      = minlon  #-- sub-region minimum longitude
    res.mpMaxLonF                      = maxlon  #-- sub-region maximum longitude
    
    res.tiMainString    = expname[0]
    res.tiMainOffsetYF  = 0.0
    p = Ngl.contour_map(wks,var0,res)
    plot.append(p)

    res.tiMainString    = expname[1]
    res.tiMainOffsetYF  = 0.0
    p = Ngl.contour_map(wks,var1,res)
    plot.append(p)

    res.tiMainString    = expname[2]
    res.tiMainOffsetYF  = 0.0
    p = Ngl.contour_map(wks,var2,res)
    plot.append(p)

    #-- panel resources
    pnlres = Ngl.Resources()
    pnlres.nglDraw                          = True
    pnlres.nglFrame                         = True
    pnlres.nglPanelLabelBar                 = True     # Turn on panel labelbar
    pnlres.nglPanelLabelBarLabelFontHeightF = fontsize*0.8    # Labelbar font height
    pnlres.nglPanelLabelBarHeightF          = 0.06   # Height of labelbar
    pnlres.nglPanelLabelBarWidthF           = 0.60    # Width of labelbar
    pnlres.lbTitleString                    = varDesc
    pnlres.lbTitlePosition                  = "Bottom"
    pnlres.lbTitleFontHeightF               = fontsize*0.8 
    pnlres.nglPanelTop                      =  0.95          #-- top position of panel
    pnlres.nglPanelYWhiteSpacePercent       =  2   #-- reduce space between the panel plots
    pnlres.nglPanelXWhiteSpacePercent       =  5   #-- reduce space between the panel plots
    Ngl.panel(wks,plot[:],[len(plot),1],pnlres)

if __name__ == '__main__':
    main()
