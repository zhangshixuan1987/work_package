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
import Ngl,Nio,os

#-----------------------------------------------------------------------
#-- Function main
#-----------------------------------------------------------------------
def main():

    diri = '/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/FV_180x360'
           #-- data directory
    varName = 'U'
    varDesc = 'Zonal wind'
    varUnit = 'm s~S~-1~N~'
    vtdUnit = 'm s~S~-2~N~'
    expndg  = 'NDGUV' # nudged experiments for plotting
    #expndg = 'NDGUVT' # nudged experiments for plotting
    #expndg = 'NDGUVTQ' # nudged experiments for plotting

    # Create time slice from dates
    season     = 'July'
    p0mb       = 1000.
    plev       = 850 # pressure level at 850 
    ntpd       = 8  # 3-hourly data sample per day     
    start_time = '2010-07-01'
    end_time   = '2010-12-31'

    # Create slice variables subset domain
    time_slice = slice(start_time, end_time)
    region     = "Global"
    regstr     = "Global"  
    minlat     = -90      
    maxlat     =  90
    minlon     =   0  
    maxlon     = 360
    cenlon     = 180
    lat_slice  = slice(minlat, maxlat) 
    lon_slice  = slice(minlon, maxlon) 

    # minimum and maximum range for plot
    minval  =  -10
    maxval  =   10
    incval  =  1.0

    # data for the plot
    expname = []
    filname = []
    expname.append('ERA5')
    filname.append(diri+'/reference_data/ML_REF_ERA5_3hourly_200911302100-201011301800.nc')
    expname.append('CLIM')
    filname.append(diri+'/clim/E3SMv2_CLIM_3hourly_200911302100-201011301800.nc')

    if (expndg == 'NDGUV'):
        expname.append('NDGUV')
        filname.append(diri+'/nudging_tendency/Nudging_Tendency_E3SMv2_NDGUV_tau6_3hourly_200912010000-201011302100.nc')
    elif (expndg == 'NDGUVT'):
        expname.append('NDGUVT')
        filname.append(diri+'/nudging_tendency/Nudging_Tendency_E3SMv2_NDGUVT_tau6_3hourly_200912010000-201011302100.nc')
    else:
        expname.append('NDGUVTQ')
        filname.append(diri+'/nudging_tendency/Nudging_Tendency_E3SMv2_NDGUVTQ_tau6_3hourly_200912010000-201011302100.nc')
        
    #read the reference data
    print ('working on ' + expname[0])
    ds0  = xr.open_dataset(filname[0])

    # Read the CLIM simulation
    print ('working on ' + expname[1])
    ds1  = xr.open_dataset(filname[1])

    # Read the nudged simulation (nudging tendency)
    print ('working on ' + expname[2])
    ds2  = xr.open_dataset(filname[2])

    # Get data0 selecting time0 level0 lat/lon slicea
    dat0 = ds0[varName].sel(time=time_slice,
            lat=lat_slice,
            lon=lon_slice)

    dat1 = ds1[varName].sel(time=time_slice,
            lat=lat_slice,
            lon=lon_slice)

    vart = 'Nudge_'+varName
    dat2 = ds2[vart].sel(time=time_slice,
              lat=lat_slice,
              lon=lon_slice)

    nlat = dat0.shape[1]
    nlon = dat0.shape[2]

    # Specify longitude values for chosen domain
    lats = dat0.lat.values
    lons = dat0.lon.values

    #
    #  Extract the desired variables.
    #
    hyam = ds0["hyam"][0]
    hybm = ds0["hybm"][0]
    psrf = (ds0["PS"][:,:,:])
    #
    #  Do the interpolation.
    #
    tmp  = Ngl.vinth2p(dat1,hyam,hybm,[plev],psrf,2,p0mb,1,False)
    var0 = np.mean((np.mean(tmp,axis=0,keepdims=False)),axis=0,keepdims=False)
    
    tmp  = Ngl.vinth2p(dat1,hyam,hybm,[plev],psrf,2,p0mb,1,False)
    var1 = np.mean((np.mean(tmp,axis=0,keepdims=False)),axis=0,keepdims=False)

    tmp  = Ngl.vinth2p(dat2,hyam,hybm,[plev],psrf,2,p0mb,1,False)
    var2 = np.mean((np.mean(tmp,axis=0,keepdims=False)),axis=0,keepdims=False)

    var1 = var1 - var0 # model bias in CLIM
    var2 = var2 * 1e5  # amplify the tendency 
    var1 = Ngl.add_cyclic(var1)
    var2 = Ngl.add_cyclic(var2)

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
    res.lbLabelBarOn              = True
    res.lbOrientation             = "Horizontal";"vertical" 
    res.pmLabelBarHeightF         = 0.12 #0.65
    res.pmLabelBarWidthF          = 0.62 #0.1
    #res.pmLabelBarParallelPosF    = 0.55
    #res.pmLabelBarOrthogonalPosF  = -0.1
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
    res.cnFillPalette                  = "BlueDarkRed18" #"BlGrYeOrReVi200"
    res.cnMissingValFillColor          = 'gray50'       #-- missing value fill color
    
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
    
    res.mpProjection                   = "Robinson"        # choose projection
    res.mpFillOn                       = False             # turn on map fill
    res.mpGridAndLimbOn                = True              # turn on lat/lon/limb lines
    res.mpGridLineColor                = "transparent"     # we don't want lat/lon lines
    res.mpPerimOn                      = False             # turn off map perimeter
    res.mpPerimLineColor               = -1 #'transparent' 
    res.mpOutlineOn                    = True
    res.mpGridAndLimbOn                = True #False 
    res.mpGridLineColor                = -1 #'transparent' 
    res.pmTickMarkDisplayMode          = 'Never'#-- don't draw tickmark border (box) around plot

    res.mpGeophysicalLineThicknessF    = 2.0
    res.mpCenterLonF                   = 300.

    res.sfXArray                       = lons
    res.sfYArray                       = lats
    res.mpLimitMode                    = 'LatLon'#-- must be set using minLatF/maxLatF/minLonF/maxLonF
    res.mpMinLatF                      = minlat  #-- sub-region minimum latitude
    res.mpMaxLatF                      = maxlat  #-- sub-region maximum latitude
    res.mpMinLonF                      = minlon  #-- sub-region minimum longitude
    res.mpMaxLonF                      = maxlon  #-- sub-region maximum longitude
    
    res.cnLevelSelectionMode           = "ManualLevels" #-- select manual level selection mode
    res.cnMinLevelValF                 = minval #-- minimum contour value
    res.cnMaxLevelValF                 = maxval #-- maximum contour value
    res.cnLevelSpacingF                = incval   #-- contour increment
    res.tiMainString                   = "Mean biases in lower-resolution EAMv2 (CLIM - ERA5)" 
    res.tiMainOffsetYF                 = 0.0
    res.lbTitleString                  = varDesc+" biase ("+varUnit+")"
    p = Ngl.contour_map(wks,var1,res)
    plot.append(p)

    res.cnLevelSelectionMode           = "ManualLevels" #-- select manual level selection mode
    res.cnMinLevelValF                 = -5     #-- minimum contour value
    res.cnMaxLevelValF                 =  5     #-- maximum contour value
    res.cnLevelSpacingF                =  0.5   #-- contour increment
    res.tiMainString                   = "Nuding tendency in NDG_UVT_tau6" +expname[2]
    res.lbTitleString                  = varDesc+" tendency ("+vtdUnit+")"
    res.tiMainOffsetYF                 = 0.0
    p = Ngl.contour_map(wks,var2,res)
    plot.append(p)

    #-- panel resources
    pnlres = Ngl.Resources()
    pnlres.nglDraw                          = True
    pnlres.nglFrame                         = True
    pnlres.nglPanelLabelBar                 = False     # Turn on panel labelbar
    pnlres.nglPanelLabelBarLabelFontHeightF = fontsize*0.8    # Labelbar font height
    pnlres.nglPanelLabelBarHeightF          = 0.06   # Height of labelbar
    pnlres.nglPanelLabelBarWidthF           = 0.60    # Width of labelbar
    pnlres.lbTitleString                    = varDesc
    pnlres.lbTitlePosition                  = "Bottom"
    pnlres.lbTitleFontHeightF               = fontsize*0.8 
    pnlres.nglPanelTop                      =  0.95          #-- top position of panel
    pnlres.nglPanelYWhiteSpacePercent       =  2   #-- reduce space between the panel plots
    pnlres.nglPanelXWhiteSpacePercent       =  5   #-- reduce space between the panel plots
    Ngl.panel(wks,plot[:],[1,len(plot)],pnlres)

if __name__ == '__main__':
    main()
