'''
Description: Python script using PyNGL Python module

 - contour plot on map (rectilinear data)

shixuan.zhang@pnnl.gov
'''

import numpy as np
import xarray as xr
import Ngl

#-------------------------------------------------------
# Function to attach lat/lon labels to a Robinson plot
#-------------------------------------------------------
def add_labels(wks,map,dlat,dlon):

    #-- generate lat/lon values
    dla = int(180.0/dlat)+1                  #-- number of lat labels
    dlo = int(360.0/dlon)+1                  #-- number of lon lables

    lat_values = np.linspace( -90.0,  90.0, dla, endpoint=True)
    lon_values = np.linspace(-180.0, 180.0, dlo, endpoint=True)
    nlat       = len(lat_values)
    nlon       = len(lon_values)

    #-- assign arrays to hold the labels
    lft, rgt                     = [],[]
    lat_lft_label, lat_rgt_label = [],[]
    lon_bot_label                = []

    #-- text resources
    txres               = Ngl.Resources()
    txres.txFontHeightF = 0.024

    #-- add degree sign and S/N to the latitude labels
    #-- don't write 90S label which would be too close to the lon labels
    for l in lat_values:
        if l == -90.0:
           lat_lft_label.append("".format(l))
           lat_rgt_label.append("".format(l))
        elif l < 0:
           lat_lft_label.append("{}~S~o~N~S    ".format(np.fabs(l)))
           lat_rgt_label.append("    {}~S~o~N~S".format(np.fabs(l)))
        elif l > 0:
           lat_lft_label.append("{}~S~o~N~N    ".format(l))
           lat_rgt_label.append("    {}~S~o~N~N".format(l))
        else:
           lat_lft_label.append("0  ")
           lat_rgt_label.append("   0")

    #-- add degree sign and W/E to the longitude labels
    for l in lon_values:
        if l < 0:
           lon_bot_label.append("{}~S~o~N~W".format(np.fabs(l)))
        elif l > 0:
           lon_bot_label.append("{}~S~o~N~E".format(l))
        else:
           lon_bot_label.append("0")

    #-- add the latitude labels left and right to the plot
    for n in range(0,nlat):
        txres.txJust = "CenterRight"
        lft.append(Ngl.add_text(wks,map,lat_lft_label[n],-180.0,\
                                lat_values[n],txres))
        #txres.txJust = "CenterLeft"
        #rgt.append(Ngl.add_text(wks,map,lat_rgt_label[n],180.0,\
        #                               lat_values[n],txres))
    #-- add the longitude labels at the bottom of the plot
    bot = []
    for n in range(0,nlon):
        txres.txJust = "TopCenter"
        bot.append(Ngl.add_text(wks,map,lon_bot_label[n],lon_values[n],\
                                -90.0,txres))
    return

#-----------------------------------------------------------------------
#-- Function: main
#-----------------------------------------------------------------------
def main():
    diri = '/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/FV_PG2' #-- data directory
    varName = 'U'
    varDesc = 'Zonal wind'
    varUnit = 'm/s'
    itime   = 1
    ilev    = 70

    minval  = -12.0
    maxval  =  12.0
    incval  =   4.0

    expname = []
    filname = []

    expname.append('CLIM')
    filname.append(diri+'/E3SMv2_CLIM.2013010100-2013010118.nc')

    for i in range(len(expname)):
        print ('working on ' + expname[i])
        ds = xr.open_dataset(filname[i])
        lat = ds.lat[:]
        lon = ds.lon[:]
        lev = ds.lev[:]
        
        tmp = ds.variables[varName][itime,ilev,:,:]
        var = Ngl.add_cyclic(tmp[:,:])  #-- add cyclic points
        
        if i == 0:
            wkres = Ngl.Resources()     #-- generate an res object for workstation
            wkres.wkWidth  = 2500       #-- plot resolution 2500 pixel width
            wkres.wkHeight = 2500       #-- plot resolution 2500 pixel height
            wks_type       = 'png'      #-- graphics output type
            figname        = 'fig_e3sm_map_180x360_'+varName
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

        fontsize = 0.024
        res.tmXBLabelFontHeightF      = fontsize
        res.tmYLLabelFontHeightF      = fontsize
        res.tmYRLabelFontHeightF      = fontsize
        res.tiMainFontHeightF         = fontsize * 0.9
        res.tiXAxisFontHeightF        = fontsize
        res.tiYAxisFontHeightF        = fontsize
        
        #-- labelbar resources
        res.lbLabelBarOn              = True      
        res.lbOrientation             = "Horizontal";"vertical" 
        #res.lbBoxMinorExtentF         = 0.15
        res.pmLabelBarDisplayMode     = "Always" #-- turn on the label bar
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
        
        BorderThick = 2.0
        res.tmBorderThicknessF    = BorderThick
        res.tmXBMajorThicknessF   = BorderThick
        res.tmXBMinorThicknessF   = BorderThick*0.5
        res.tmYLMajorThicknessF   = BorderThick
        res.tmYLMinorThicknessF   = BorderThick*0.5
        res.tmXTMajorThicknessF   = BorderThick
        res.tmXTMinorThicknessF   = BorderThick*0.5
        res.tmYRMajorThicknessF   = BorderThick
        res.tmYRMinorThicknessF   = BorderThick*0.5

        #-- grid data information resources
        res.sfXCStartV = float(min(lon)) #-- x-axis location of 1st element lon
        res.sfXCEndV   = float(max(lon)) #-- x-axis location of last element lon
        res.sfYCStartV = float(min(lat)) #-- y-axis location of 1st element lat
        res.sfYCEndV   = float(max(lat)) #-- y-axis location of last element lat
        
        #-- contour resources
        res.cnLinesOn                      = False
        res.cnLineLabelsOn                 = False  #-- turn off line labels
        res.cnInfoLabelOn                  = False   #-- turn off info label
        res.cnLineThicknessF               = 0.1
        res.cnLineColor                    = "Black"
        res.cnLineLabelDensityF            = 1.0
        res.cnLineLabelFontHeightF         = fontsize
        res.cnLineLabelBackgroundColor     = -1
        #res.gsnContourLineThicknessesScale = 1.0
        #res.gsnContourZeroLineThicknessF   = 1.
        #res.gsnContourNegLineDashPattern   = 2

        res.cnFillOn                       = True         #-- turn on contour fill
        res.cnFillPalette                  = "cmp_b2r" #"BlGrYeOrReVi200"
       #res.cnFillMode                     = 'CellFill'     #-- change contour fill mode
       #res.cnCellFillEdgeColor            = 'black'        #-- edges color
       #res.cnCellFillMissingValEdgeColor  = 'gray50'       #-- missing value edges color
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
        res.tiMainString                   = expname[i]+ " (" +str(int(lev[ilev]))+" hPa)"
        res.tiMainOffsetYF                 = 0.0
        
        res.mpFillOn                       = False
        #res.mpGridLatSpacingF              = 30. #-- grid lat spacing
        #res.mpGridLonSpacingF              = 60. #-- grid lon spacing
        res.mpDataBaseVersion              = "MediumRes" #-- map database
        res.mpOceanFillColor               = "Transparent"
        res.mpLandFillColor                = "Gray90"
        res.mpInlandWaterFillColor         = "Gray90"
        
        #res.mpProjection                   = "Orthographic"
        #res.mpProjection                   = "Robinson"               #-- set projection
        res.mpProjection                   = "Mollweide"
        res.mpPerimOn                      = False
        res.mpPerimLineColor               = -1 #'transparent' 
        res.mpOutlineOn                    = True
        res.mpGridAndLimbOn                = True #False 
        res.mpGridLineColor                = -1 #'transparent' 
        res.mpGeophysicalLineThicknessF    = 2.0
        res.pmTickMarkDisplayMode          = 'Never'

        #res.mpLimitMode = 'LatLon' #-- must be set using minLatF/maxLatF/minLonF/maxLonF
        #res.mpMinLatF   = -10.     #-- sub-region minimum latitude
        #res.mpMaxLatF   = 80.      #-- sub-region maximum latitude
        #res.mpMinLonF   = -120.    #-- sub-region minimum longitude
        #res.mpMaxLonF   = 60.      #-- sub-region maximum longitude
        
        res.mpCenterLatF =   0.     #-- center latitude
        res.mpCenterLonF =  180.    #-- center longitude

        p = Ngl.contour_map(wks,var,res)
        
        #-- add labels to the plot
        #tx = add_labels(wks,p,30.,120.)
        
        #Ngl.maximize_plot(wks,p)

        #-- Retrieve some resources from map for adding labels
        #vpx  = Ngl.get_float(p,'vpXF')
        #vpy  = Ngl.get_float(p,'vpYF')
        #vpw  = Ngl.get_float(p,'vpWidthF')

        #-- add title string,long_name and units string to panel
        #txres = Ngl.Resources()
        #txres.txFontHeightF = fontsize
        #txres.txJust  = "CenterLeft"
        #Ngl.text_ndc(wks,expname[i],vpx,vpy+0.02,txres)
 
        #txres.txJust  = "CenterRight"
        #Ngl.text_ndc(wks,varName+" ("+varUnit+")",vpx,vpy+0.02,txres)

        plot.append(p)
    
    #-- panel resources
    pnlres = Ngl.Resources()
    pnlres.nglDraw          = True
    pnlres.nglFrame         = True
    pnlres.nglPanelLabelBar = False          #-- common labelbar
    pnlres.nglPanelTop      =  0.95          #-- top position of panel
    pnlres.txString         = ""
    pnlres.txFontHeightF    = fontsize
    pnlres.nglPanelYWhiteSpacePercent =  5   #-- reduce space between the panel plots
    pnlres.nglPanelXWhiteSpacePercent =  5   #-- reduce space between the panel plots
    Ngl.panel(wks,plot[:],[len(plot),1],pnlres)

if __name__ == '__main__':
    main()
