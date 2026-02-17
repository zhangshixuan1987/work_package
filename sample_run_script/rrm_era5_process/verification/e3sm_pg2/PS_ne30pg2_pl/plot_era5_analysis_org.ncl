;This example regrids data from a global rectilinear grid (360 x 180) to a 1.9x2.0 degree lat/lon grid (144 x 96), using the default "bilinear" interpolation method.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file

    grid_era5 = "721x1440"
    yymm      = "200904"
    yymmd1    = "20090401"
    yymmd2    = "20090430" 
    time_unt  = "hours since 2009-04-01 00:00:00"
    DatDir    = "/global/cfs/projectdirs/m3522/cmip6/ERA5/e5.oper.an.sfc/"+yymm
    FileName  = DatDir+"/e5.oper.an.sfc.128_134_sp.ll025sc."+yymmd1+"00_"+yymmd2+"23.nc"
    MainStr   = "ERA5 reanalysis on " + grid_era5 + " grid"
    figname   = "ERA5_"+grid_era5

    time_plt  = 2009040100

    f         = addfile(FileName,"r")
    time      = cd_calendar(cd_convert(f->time,time_unt), -3)
    indx      = ind(time.eq.time_plt)
    temp      = f->SP(indx,:,:)
    lat       = f->latitude
    lon       = f->longitude
    print(num(ismissing(temp)))
    temp     = temp / 100.0 
    temp@units = "hPa" 

;----------------------------------------------------------------------
; Plot the original and regridded data on a panel plot
;----------------------------------------------------------------------
    wks = gsn_open_wks("png",figname)             ; ESMF.ps

    res                     = True              ; Plot mods desired.
    res@gsnDraw             = False             ; We will panel later.
    res@gsnFrame            = False
    res@gsnMaximize         = False              ; Maximize plot
 
    res@mpFillOn                    = False          ; turn off gray continents
    res@mpGeophysicalLineColor      = "Grey30";"white"        ; color of outlines
    res@mpGeophysicalLineThicknessF = 2              ; line thickness
    res@mpProjection                = "Orthographic"   ;-- change projection
    res@mpCenterLonF                =  260              ;-- center at lon=10
    res@mpCenterLatF                =  40              ;-- center at lat=40
    res@mpDataBaseVersion           = "MediumRes"      ;-- map resolution
    res@mpGridAndLimbOn             =  False            ;-- plot grid lines
    res@mpGridLineColor             = "grey30"         ;-- set grid line color
    res@mpPerimOn                   =  False           ;-- don't draw the box around the plot

    res@cnFillOn                 = True              ; color plot desired
    res@cnLinesOn                = False              ; turn off contour lines
    res@cnLineThicknessF         = 0.01
    res@cnLineColor              = "Grey"
    res@cnLineLabelsOn           = False             ; turn off contour lines
    res@cnInfoLabelOn            = False
    res@cnFillMode               = "RasterFill" ;"CellFill"
    res@cnGridBoundFillColor     = "black"
    res@cnRasterMinCellSizeF     = 0.001
    res@cnFillPalette            = "WhBlGrYeRe";"cmocean_dense"
 
    res@cnLevelSelectionMode     = "ManualLevels"  ; manual levels
    res@cnMinLevelValF           = 600             ; min level
    res@cnMaxLevelValF           = 1100            ; max level
    res@cnLevelSpacingF          = 2               ; interval

    res@lbLabelBarOn             = False              ; Labelbar will be in panel
    res@lbBoxLinesOn             = False
    res@lbLabelStride            =  1               ;-- every other label
    res@lbBoxMinorExtentF        =  0.10            ;-- decrease the height of the labelbar
    res@pmLabelBarOrthogonalPosF = -0.07            ;-- move the labelbar upward

;---Plot data on original grid
    res@gsnAddCyclic             = True
    res@tiMainString             = MainStr 
    plot_orig = gsn_csm_contour_map(wks,temp,res)     

;---Resources for paneling
    pres                          = True
    pres@gsnMaximize              = False
    pres@gsnPanelLabelBar         = True
    pres@lbBoxLinesOn             = False 
    pres@lbLabelStride            = 1               ;-- every other label
    pres@lbBoxMinorExtentF        = 0.15            ;-- decrease the height of the labelbar
    pres@pmLabelBarOrthogonalPosF = 0.03            ;-- move the labelbar upward
    pres@pmLabelBarHeightF        = 0.15
    pres@pmLabelBarWidthF         = 0.95
    gsn_panel(wks,(/plot_orig/),(/1,1/),pres)
end

