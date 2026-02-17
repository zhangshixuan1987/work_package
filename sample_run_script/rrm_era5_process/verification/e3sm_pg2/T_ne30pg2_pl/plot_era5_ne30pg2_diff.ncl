;This example regrids data from a global rectilinear grid (360 x 180) to a 1.9x2.0 degree lat/lon grid (144 x 96), using the default "bilinear" interpolation method.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file

    grid_e3sm = "northamericax4v1pg2"
    time_tag = "2009-04-01-75600"
    DatDir   = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis/era5_ml_to_e3sm_rrm"
    FileNam0 = DatDir+"/ERA5_ne30pg2_L72_"+grid_e3sm+"."+time_tag+".nc"

    grid_e3sm = "northamericax4v1pg2"
    time_tag  = "2009-04-01-75600"
    DatDir    = "/global/cscratch1/sd/zhan391/SciDAC_m3089/acme_init/era5_reanalysis/era5_pl_to_e3sm_rrm"
    FileName  = DatDir+"/ERA5_ne30pg2_L72_"+grid_e3sm+"."+time_tag+".nc"

    MainStr  = "Difference between ERA5 ML and PL data after remapping"

    ilev = ispan(71,0,1) ;(/71,50,30,10,5,1/)
    nlev = dimsizes(ilev)

    varNam   = "T"
    figname  = "ERA5_ml_vs_pl_"+grid_e3sm+"_"+varNam

   do i = 0,nlev-1,1

    f0       = addfile(FileNam0,"r")
    tem0     = f0->$varNam$(0,ilev(i),:)

    f        = addfile(FileName,"r")
    temp     = f->$varNam$(0,ilev(i),:)
    lat      = f->lat
    lon      = f->lon
    lev      = f->lev

    temp      = temp - tem0
    temp      = temp *100.0 / tem0

    temp@long_name = "~F33~D~N~~F~"+temp@long_name 
    temp@units     = "%"

    substr   = temp@long_name+" ("+temp@units+", "+sprintf("%.2f",lev(ilev(i)))+" hPa"+")"

    ;----------------------------------------------------------------------
    ; Plot the original and regridded data on a panel plot
    ;----------------------------------------------------------------------
    ;if(i.eq.0) then
     fignam1 = str_squeeze(figname+"_"+sprintf("%.2f",lev(ilev(i)))+"hPa")
     wks = gsn_open_wks("png",fignam1)             ; ESMF.ps
     plot = new(nlev,graphic)
    ;end if

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
 
    xmin = -decimalPlaces(max(abs(temp)),2,True)
    xmax =  decimalPlaces(max(abs(temp)),2,True)
    xint = (xmax - xmin)/100.
    res@cnLevelSelectionMode     = "ManualLevels"  ; manual levels
    res@cnMinLevelValF           = xmin            ; min level
    res@cnMaxLevelValF           = xmax            ; max level
    res@cnLevelSpacingF          = xint            ; interval

    res@lbLabelBarOn             = True              ; Labelbar will be in panel
    res@lbBoxLinesOn             = False
    res@lbLabelStride            =  1               ;-- every other label
    res@lbBoxMinorExtentF        =  0.10            ;-- decrease the height of the labelbar
    res@pmLabelBarOrthogonalPosF = -0.07            ;-- move the labelbar upward

;---Plot data on original grid
   ;res@gsnAddCyclic   = False
    res@sfXArray       = lon
    res@sfYArray       = lat
    res@tiMainString   = substr
    res@gsnLeftString  = ""
    res@gsnRightString = ""
    plot(i) = gsn_csm_contour_map(wks,temp,res)

;---Resources for paneling
    pres                          = True
    pres@gsnMaximize              = False
    pres@gsnPanelLabelBar         = False
    pres@lbBoxLinesOn             = False 
    pres@lbLabelStride            = 1               ;-- every other label
    pres@lbBoxMinorExtentF        = 0.15            ;-- decrease the height of the labelbar
    pres@pmLabelBarOrthogonalPosF = 0.03            ;-- move the labelbar upward
    pres@pmLabelBarHeightF        = 0.15
    pres@pmLabelBarWidthF         = 0.95
    gsn_panel(wks,(/plot(i)/),(/1,1/),pres)

   end do

end

