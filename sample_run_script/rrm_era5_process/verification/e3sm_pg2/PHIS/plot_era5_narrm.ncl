;This example regrids data from a global rectilinear grid (360 x 180) to a 1.9x2.0 degree lat/lon grid (144 x 96), using the default "bilinear" interpolation method.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file

    grid_tag = "northamericax4v1pg2"
    time_tag = "2009-04-01-00000"
    DatDir   = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis"
    FileName = DatDir+"/era5_to_e3sm_rrm/ERA5_ne30pg2_L72_"+grid_tag+"."+time_tag+".nc"
    mapfile  = DatDir+"/era5_to_e3sm_rrm_map/map_640x1280_northamericax4v1pg2.nc"
    MainStr  = "ERA5 reanalysis on " + grid_tag + " grid"
    figname  = "ERA5_"+grid_tag

    f        = addfile(FileName,"r")
    f0       = addfile(mapfile,"r")

    temp     = f->PHIS(0,:)
    lat      = f->lat
    lon      = f->lon
    ;clatbnds = f0->lat_bnds ;* r2d
    ;clonbnds = f0->lon_bnds ;* r2d

    print(num((temp.lt.0.0)))
    temp     = temp / 9.80665 
    temp@units = "m" 

;----------------------------------------------------------------------
; Plot the original and regridded data on a panel plot
;----------------------------------------------------------------------
    wks = gsn_open_wks("png",figname)             ; ESMF.ps

    res                     = True              ; Plot mods desired.
    res@gsnDraw             = False             ; We will panel later.
    res@gsnFrame            = False

    res@gsnMaximize         = True              ; Maximize plot
    
    res@mpMaxLatF                   = 90                 ; choose map range
    res@mpMinLatF                   = -90 
    res@mpFillOn                    = False          ; turn off gray continents
    res@mpGeophysicalLineColor      = "Grey20";"white"        ; color of outlines
    res@mpGeophysicalLineThicknessF = 2              ; line thickness
    res@mpGridAndLimbOn             = False               ; turn on lat/lon grid lines

    res@cnFillOn             = True              ; color plot desired
    res@cnLinesOn            = False             ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour lines
    res@cnInfoLabelOn        = False
    res@cnFillMode           = "RasterFill" ;"CellFill"
    res@cnGridBoundFillColor = "grey"
    res@cnRasterMinCellSizeF = 0.001
    res@cnFillPalette        = "WhBlGrYeRe";"cmocean_dense"
    
    res@cnLevelSelectionMode = "ManualLevels"  ; manual levels
    res@cnMinLevelValF       = -100             ; min level
    res@cnMaxLevelValF       = 3000            ; max level
    res@cnLevelSpacingF      = 100               ; interval

    res@lbLabelBarOn        = False              ; Labelbar will be in panel
    res@lbBoxLinesOn        = False

;---Plot data on original grid
    res@gsnAddCyclic = False

    res@sfXArray      = lon
    res@sfYArray      = lat
    ;res@sfXCellBounds = clonbnds
    ;res@sfYCellBounds = clatbnds

    dims = tostring(dimsizes(temp))
    res@tiMainString = MainStr 
    plot_orig = gsn_csm_contour_map(wks,temp,res)     

;---Resources for paneling
    pres                  = True
    pres@gsnMaximize      = True
    pres@gsnPanelLabelBar = True
    pres@lbBoxLinesOn     = False 
    gsn_panel(wks,(/plot_orig/),(/1,1/),pres)
end

