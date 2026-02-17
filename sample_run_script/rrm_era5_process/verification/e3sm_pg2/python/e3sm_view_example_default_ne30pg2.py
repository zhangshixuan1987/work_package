'''
Description: Python script using PyNGL Python module

 - contour plot on map (rectilinear data)

shixuan.zhang@pnnl.gov
'''

import numpy as np
import math, time, sys, os
import Nio, Ngl

#-- rearrange the longitude values to -180.-180.
def rearrange(vlon):
    less_than    = vlon < -180.
    greater_than = vlon >  180.
    vlon[less_than]    = vlon[less_than] + 360.
    vlon[greater_than] = vlon[greater_than] - 360.
    return vlon

if __name__ == '__rearrange__':
    rearrange()

#-----------------------------------------------------------------------
#-- Function: main
#-----------------------------------------------------------------------
def main():

    t1 = time.time()                                    #-- retrieve start time

    diri = '/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis'
    filname = diri+'/era5_to_e3sm_rrm/ERA5_ne30pg2_L72.2009-04-01-00000.nc'
    filmap  = diri+'/era5_to_e3sm_rrm_map/map_640x1280_northamericax4v1pg2.nc'
    print ('working on ' + filname)

    dm = Nio.open_file(filmap,'r')
    vlon, vlat =  dm.variables['xv_b'][:], dm.variables['yv_b'][:]
    ncells, nv =  vlon.shape #-- ncells: number of cells; nv: number of edges
    
    #-- print information to stdout
    print ('')
    print ('cell points:      ', nv)
    print ('cells:            ', str(ncells))
    print ('')
   
    vlon = rearrange(vlon)      #-- set longitude values to -180.-180. degrees
    print('min/max vlon:     {} {}'.format(np.min(vlon), np.max(vlon)))
    print('min/max vlat:     {} {}'.format(np.min(vlat), np.max(vlat)))
    print('')

    varName = 'U'
    varDesc = 'Zonal wind'
    varUnit = 'm/s'
    itime   = 0
    ilev    = 70

    ds = Nio.open_file(filname,'r')
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    lev = ds.variables['lev'][:]
    vor = ds.variables[varName]
    var = vor[itime,ilev,:]
    print(var)
    lon = rearrange(lon)      #-- set longitude values to -180.-180. degrees
    print('min/max lon:     {} {}'.format(np.min(lon), np.max(lon)))
    print('min/max lat:     {} {}'.format(np.min(lat), np.max(lat)))
    print('')

    #-- define _FillValue and missing_value if not existing
    missing = -1e20
    varM = np.ma.array(var, mask=np.equal(var,missing)) #-- mask array with missing values
    nummissing = np.count_nonzero(varM.mask)
    print(varM)

    #-- set data intervals, levels, labels, color indices
    varMin, varMax, varInt = -50, 50, 10                 #-- set data minimum, maximum, interval
    levels   = range(varMin,varMax,varInt)             #-- set levels array
    nlevs    = len(levels)                             #-- number of levels
    labels   = ['{:.2f}'.format(x) for x in levels]     #-- convert list of floats to list of strings
    colors   = range(2,nlevs+6)                        #-- create color indices
    print(colors)

    wkres = Ngl.Resources()     #-- generate an res object for workstation
    wkres.wkColorMap =  'BlAqGrYeOrReVi200' #'WhiteBlueGreenYellowRed'       #-- choose colormap    
    wkres.wkWidth, wkres.wkHeight  =  1024, 1024
    wks_type       = 'png'      #-- graphics output type
    figname        = 'fig_e3sm_map_NE30pg2_'+varName
    wks  = Ngl.open_wks(wks_type,figname,wkres)  #-- open workstation
    plot = []

    #-- define colormap
    cmap     =  Ngl.retrieve_colormap(wks)             #-- RGB ! [256,3]
    ncmap    =  cmap.shape[0]                           #-- number of colors
    colormap =  cmap[:ncmap:12,:]                       #-- select every 13th color
    ncol     =  colormap.shape[0]
    #colormap[20] = ([1.,1.,1.])                         #-- white for missing values
    
    #-- define colormap
    #colormap =  Ngl.read_colormap_file('BlAqGrYeOrReVi200')[22::12,:]    #-- RGB ! [256,4] -> [20,4] #-- select every 12th color
    #colormap[19,:] = [1.,1.,1.,0.]                    #-- white for missing values

    print ('colors index:     ', colors)
    print ('')
    print('levels:           {}'.format(levels))
    print('labels:           {}'.format(labels))
    print ('')
    print('nlevs:            {:3d}'.format(nlevs))
    print('ncols:            {:3d}'.format(ncol))
    print ('')

    #-- overwrite resources of wks
    setlist                    =  Ngl.Resources()
    setlist.wkColorMap         =  colormap              #-- set color map to new colormap array
    setlist.wkBackgroundColor  = 'white'                #-- has to be set when wkColorMap is set to colormap array
    setlist.wkForegroundColor  = 'black'                #-- has to be set when wkColorMap is set to colormap array
    Ngl.set_values(wks,setlist)

    #BorderThick = 2.0
    #res.tmBorderThicknessF    = BorderThick
    #res.tmXBMajorThicknessF   = BorderThick
    #res.tmXBMinorThicknessF   = BorderThick*0.5
    #res.tmYLMajorThicknessF   = BorderThick
    #res.tmYLMinorThicknessF   = BorderThick*0.5
    #res.tmXTMajorThicknessF   = BorderThick
    #res.tmXTMinorThicknessF   = BorderThick*0.5
    #res.tmYRMajorThicknessF   = BorderThick
    #res.tmYRMinorThicknessF   = BorderThick*0.5
    #-- labelbar resources
    #res.lbLabelBarOn              = True
    #res.lbOrientation             = "Horizontal";"vertical"
    #res.lbBoxMinorExtentF         = 0.15
    #res.pmLabelBarDisplayMode     = "Always" #-- turn on the label bar
    #res.pmLabelBarHeightF         = 0.12 #0.65
    #res.pmLabelBarWidthF          = 0.62 #0.1
    #res.pmLabelBarParallelPosF    = 0.55
    #res.pmLabelBarOrthogonalPosF  = -0.1
    #res.lbRightMarginF            = -0.28
    #res.lbBottomMarginF           = -0.28
    #res.lbLabelFontHeightF        = fontsize * 0.8
    #res.lbLabelPosition           = "bottom"
    #res.lbTitleString             = varDesc+" ("+varUnit+")"
    #res.lbTitleFontHeightF        = fontsize * 0.8
    #res.lbTitlePosition           = "Right"                           
    #res.lbTitleDirection          = "Across"                          
    #res.lbTitleAngleF             = 90.                               
    #res.lbTitleFontHeightF        = fontsize
 
    #res.tmXTOn                         = False
    #res.tmYROn                         = True
    #res.tmYLLabelJust                  = "CenterRight"
    #res.tmYRLabelJust                  = "CenterRight"
    #res.tmYUseLeft                     = True
    #res.tmYRLabelsOn                   = False

    #res.tiXAxisString                  = ""
    #res.tiYAxisString                  = ""
    #res.tiXAxisOffsetYF                = 0.0
    #res.tiYAxisOffsetXF                = 0.0
    #res.tiMainString                   = " (" +str(int(lev[ilev]))+" hPa)"
    #res.tiMainOffsetYF                 = 0.0

    fontsize = 0.018
    #-- set map resources
    mpres                             =  Ngl.Resources()
    mpres.tmXBLabelFontHeightF        = fontsize
    mpres.tmYLLabelFontHeightF        = fontsize
    mpres.tmYRLabelFontHeightF        = fontsize
    mpres.tiMainFontHeightF           = fontsize * 0.9
    mpres.tiXAxisFontHeightF          = fontsize
    mpres.tiYAxisFontHeightF          = fontsize
    mpres.nglDraw                     =  False          #-- turn off plot draw and frame advance. We will
    mpres.nglFrame                    =  False          #-- do it later after adding subtitles.
    mpres.mpGridAndLimbOn             =  False
    mpres.mpGeophysicalLineThicknessF =  2.
    mpres.pmTitleDisplayMode          = 'Always'
    mpres.mpFillOn                       = False
    #mpres.mpGridLatSpacingF              = 30. #-- grid lat spacing
    #mpres.mpGridLonSpacingF              = 60. #-- grid lon spacing
    #mpres.mpDataBaseVersion              = "MediumRes" #-- map database
    #mpres.mpOceanFillColor               = "Transparent"
    #mpres.mpLandFillColor                = "Gray90"
    #mpres.mpInlandWaterFillColor         = "Gray90"

    #mpres.mpProjection                   = "Orthographic"
    #mpres.mpProjection                   = "Robinson"               #-- set projection
    #mpres.mpProjection                   = "Mollweide"
    mpres.mpProjection                   = "CylindricalEquidistant"
    #mpres.mpPerimOn                      = True
    #mpres.mpPerimLineColor               = 'black' #'transparent' 
    #mpres.mpOutlineOn                    = True
    #mpres.mpGridAndLimbOn                = False
    #mpres.mpGridLineColor                = 'black' #'transparent' 
    #mpres.mpGeophysicalLineThicknessF    = 2.0
    #mpres.pmTickMarkDisplayMode          = 'Always'

    #mpres.mpLimitMode = 'LatLon' #-- must be set using minLatF/maxLatF/minLonF/maxLonF
    #mpres.mpMinLatF   = -10.     #-- sub-region minimum latitude
    #mpres.mpMaxLatF   = 80.      #-- sub-region maximum latitude
    #mpres.mpMinLonF   = -120.    #-- sub-region minimum longitude
    #mpres.mpMaxLonF   = 60.      #-- sub-region maximum longitude

    mpres.mpCenterLatF =   0.     #-- center latitude
    #mpres.mpCenterLonF =  180.    #-- center longitude
    #mpres.mpCenterLatF =  40
    #mpres.mpCenterLonF = -130


    #-- create only a map
    map = Ngl.map(wks,mpres)
    Ngl.draw(map)

    #-- assign and initialize array which will hold the color indices of the cells
    gscolors = -1*(np.ones((ncells,),dtype=np.int64))     #-- assign array containing zeros; init to transparent: -1

    #-- set color index of all cells in between levels
    for m in range(0,nlevs-1):
        vind = []                                       #-- empty list for color indices
        for i in range(0,ncells-1):
            if (varM[i] >= levels[m] and varM[i] < levels[m+1]):
               gscolors[i] = colors[m+1]
               vind.append(i)
        print ('finished level %3d' % m , ' -- %5d ' % len(vind) , ' polygons considered - gscolors %3d' % colors[m])        
        del vind

    gscolors[varM < varMin]         =  colors[0]        #-- set color index for cells less than level[0]
    gscolors[varM >= varMax]        =  colors[(nlevs-1)+2] #-- set color index for cells greater than levels[nlevs-1]
    gscolors[np.nonzero(varM.mask)] =  20               #-- set color index for missing values

    #-- set polygon resources
    pgres                   =  Ngl.Resources()
    pgres.gsEdgesOn         =  True                     #-- draw the edges
    pgres.gsFillIndex       =  0                        #-- solid fill
    pgres.gsLineColor       = 'black'                   #-- edge line color
    pgres.gsLineThicknessF  =  0.1                      #-- line thickness
    pgres.gsColors          =  -1 #gscolors                 #-- use color array
    pgres.gsSegments        =  list(range(0,len(vlon[:,0])*3,6))
    x1d, y1d = np.ravel(vlon), np.ravel(vlat)           #-- convert to 1D-arrays

    #-- add polygons to map
    polyg  = Ngl.add_polygon(wks,map,x1d,y1d,pgres)

    #-- add a labelbar
    lbres                   =  Ngl.Resources()
    lbres.vpWidthF          =  0.85
    lbres.vpHeightF         =  0.15
    lbres.lbOrientation     = 'Horizontal'
    lbres.lbFillPattern     = 'SolidFill'
    lbres.lbMonoFillPattern =  21                       #-- must be 21 for color solid fill
    lbres.lbMonoFillColor   =  False                    #-- use multiple colors
    lbres.lbFillColors      =  colormap                   #-- indices from loaded colormap
    lbres.lbBoxCount        =  len(colormap[colors,:])
    lbres.lbLabelFontHeightF=  0.014
    lbres.lbLabelAlignment  = 'InteriorEdges'
    lbres.lbLabelStrings    =  labels

    lb = Ngl.labelbar_ndc(wks,nlevs+1,labels,0.1,0.24,lbres)

    #-- maximize and draw the plot and advance the frame
    Ngl.maximize_plot(wks, map)
    Ngl.draw(map)
    Ngl.frame(wks)

    #-- get wallclock time
    t2 = time.time()
    print ('')
    print ('Wallclock time:  %0.3f seconds' % (t2-t1))
    print ('')

    #-- done
    Ngl.end()

if __name__ == '__main__':
    main()
