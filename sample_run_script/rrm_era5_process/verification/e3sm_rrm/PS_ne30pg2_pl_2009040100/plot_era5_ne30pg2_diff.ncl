;This example regrids data from a global rectilinear grid (360 x 180) to a 1.9x2.0 degree lat/lon grid (144 x 96), using the default "bilinear" interpolation method.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file

    grid_e3sm = "northamericax4v1pg2"
    time_tag = "2009-04-01-00000"
    DatDir   = "../test_data/"
    FileNam0 = DatDir+"/ERA5_640x1280_to_E3SM_ne30pg2_L72_"+grid_e3sm+"."+time_tag+".nc"

    f0       = addfile(FileNam0,"r")
    tem0     = f0->PS(0,:)
    lat0     = f0->lat
    lon0     = f0->lon

    grid_e3sm = "northamericax4v1pg2"
    time_tag  = "2009-04-01-00000"
    DatDir    = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis/era5_pl_to_e3sm_rrm"
    FileName  = DatDir+"/ERA5_ne30pg2_L72_"+grid_e3sm+"."+time_tag+".nc"

    f         = addfile(FileName,"r")
    temp      = f->PS(0,:)
    lat       = f->lat
    lon       = f->lon

    temp      = temp - tem0 
    temp      = temp / tem0 
    temp@units= "%"

    MainStr  = "Difference between ERA5 ML and PL data after remapping"
    figname  = "ERA5_ML_vs_PL_interpolation_diff"

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
    res@cnFillPalette            = "cmocean_curl"
 
    res@cnLevelSelectionMode     = "ManualLevels"  ; manual levels
    res@cnMinLevelValF           = -10             ; min level
    res@cnMaxLevelValF           =  10            ; max level
    res@cnLevelSpacingF          =  0.5           ; interval

    res@lbLabelBarOn             = False              ; Labelbar will be in panel
    res@lbBoxLinesOn             = False
    res@lbLabelStride            =  1               ;-- every other label
    res@lbBoxMinorExtentF        =  0.10            ;-- decrease the height of the labelbar
    res@pmLabelBarOrthogonalPosF = -0.07            ;-- move the labelbar upward

;---Plot data on original grid
   ;res@gsnAddCyclic = False
    res@sfXArray      = lon
    res@sfYArray      = lat
    res@tiMainString = MainStr
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

