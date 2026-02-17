;This example regrids data from a global rectilinear grid (360 x 180) to a 1.9x2.0 degree lat/lon grid (144 x 96), using the default "bilinear" interpolation method.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file
    srcFileName = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/Blocking/ppe_era5/ERA5_BLOCK_blocking_climo.nc"
    dstFileName = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/Blocking/ppe_clim/CLIM_rgd_objectid_2007-2017.nc"

;---Output (and input) files
    srcGridName = "src_SCRIP.nc"
    dstGridName = "dst_SCRIP.nc"
    wgtFileName = "era5_720x361_to_360x180_Rect_bilinear.nc"

;----------------------------------------------------------------------
; Convert original source NCEP grid to a SCRIP convention file.
;----------------------------------------------------------------------
    src_file = addfile(srcFileName,"r")
    temp     = src_file->seasonalmean_object_id(0,:,:)

    Opt                = True
    Opt@GridMask       = where(.not.ismissing(temp),1,0)   ; use Mask2D in NCL V6.2.1 and earlier
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "ERA5 1440x720 Grid"

    printVarSummary(temp)
    rectilinear_to_SCRIP(srcGridName,temp&latitude,temp&longitude,Opt)

;---Clean up
    delete(Opt)
;----------------------------------------------------------------------
; Convert destination grid to a SCRIP convention file.
;----------------------------------------------------------------------
    dst_file = addfile(dstFileName,"r")
    temp1    = dst_file->seasonalmean_object_id(0,:,:)
    print(num(ismissing(temp1)))

    Opt                = True
    Opt@GridMask       = where(.not.ismissing(temp1),1,0)   ; use Mask2D in NCL V6.2.1 and earlier
    Opt@ForceOverwrite = True
    Opt@PrintTimings   = True
    Opt@Title          = "EAM 360x180 Grid"

    printVarSummary(temp1)
    rectilinear_to_SCRIP(dstGridName,temp1&lat,temp1&lon,Opt)

;---Clean up
    delete(Opt)

;----------------------------------------------------------------------
; Generate the weights that take you from the NCEP grid to a
; 5x5 degree grid.
;----------------------------------------------------------------------
    Opt                      = True
    Opt@InterpMethod         = "bilinear"     ; default
    Opt@ForceOverwrite       = True
    Opt@PrintTimings         = True

    ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileName,Opt)

    delete(Opt)


;----------------------------------------------------------------------
; Apply the weights to a given variable on the NCEP file.
;----------------------------------------------------------------------
    Opt                = True
    Opt@PrintTimings   = True

;---In V6.1.0, coordinates and attributes are copied automatically
    temp_regrid = ESMF_regrid_with_weights(temp,wgtFileName,Opt)
    printVarSummary(temp_regrid)

;----------------------------------------------------------------------
; Plot the original and regridded data on a panel plot
;----------------------------------------------------------------------
    wks = gsn_open_wks("png","ESMF")             ; ESMF.ps
    gsn_define_colormap(wks,"gui_default")      ; choose colormap

    res                     = True              ; Plot mods desired.
    res@gsnDraw             = False             ; We will panel later.
    res@gsnFrame            = False

    res@gsnMaximize         = True              ; Maximize plot
    
    res@mpMaxLatF           = 60                 ; choose map range
    res@mpMinLatF           = -60 

    res@cnFillOn            = True              ; color plot desired
    res@cnLinesOn           = False             ; turn off contour lines
    res@cnLineLabelsOn      = False             ; turn off contour lines

    ;res@cnLevelSelectionMode =  "ManualLevels"   ; manual levels
    ;res@cnMinLevelValF       = 4                 ; min level
    ;res@cnMaxLevelValF       = 32                ; max level
    ;res@cnLevelSpacingF      =  2                ; interval

    res@lbLabelBarOn        = False              ; Labelbar will be in panel

;---Plot data on original grid
    res@gsnAddCyclic = False

    dims = tostring(dimsizes(temp))
    res@tiMainString = "NCEP monthly means temp: original data (" + \
                       str_join(dims," x ") + ")"
    plot_orig = gsn_csm_contour_map(wks,temp,res)     

;---Plot data interpolated onto 5x5 degree grid
    res@gsnAddCyclic = True

    dims = tostring(dimsizes(temp_regrid))
    res@tiMainString = "NCEP monthly means temp: regridded to 5x5 grid (" +\
                       str_join(dims," x ") + ")"
    plot_regrid = gsn_csm_contour_map(wks,temp_regrid,res)     

;---Resources for paneling
    pres                  = True
    pres@gsnMaximize      = True
    pres@gsnPanelLabelBar = True

    gsn_panel(wks,(/plot_orig,plot_regrid/),(/2,1/),pres)
end

