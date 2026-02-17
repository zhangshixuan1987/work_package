;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;This script is used to plot the wrf simulations vs NEXRAD data ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  LOAD FUNCTIONS AND PROCEDURES
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

 ;=======================================================================================
 ; experiments for plot
 ;=======================================================================================
 datdir   = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/Blocking"
 Casdir   = (/"ppe_clim"/)
 Groups   = (/"CLIM"/)
 Expnam   = (/"CLIM"/)
 ngrps    = dimsizes(Groups)
 year     = 2010
 month    = (/1, 7/)
 day      = (/1, 1/)
 timstr   = (/"Jan. 1", "Jul. 1"/)
 nt       = dimsizes(month)
 lstring  = "H(500hPa) Blocking Thresholds" 

 ;********************************
 ; select sub-regions
 ;********************************
 minLat   = -90.             ; max lat to read in (deg)
 maxLat   =  90.             ; min lat to read in (deg)
 minLon   =   0.             ; max lon to read in (deg)
 maxLon   = 360.             ; min lon to read in (deg)
 cenlon   = 180.0 ;

 do i = 0,ngrps-1,1

   fl       = datdir+"/"+Casdir(i)+"/"+Groups(i)+"_BLOCK_rgd_threshold_H_filtered.nc"
   f        = addfile(fl,"r")

   lat      = f->lat
   lon      = f->lon

   do j = 0,nt-1,1

    kd  = day_of_year(year,month(j),day(j))  ; doy = (/61,61,61/)
    var = f->threshold_H(kd,:,:)
   ;var = var/10.0 

    if(j.eq.0) then 
     ;;;generate the zoom in information;;;;
     fig_type  = "png" ;;figure types: pdf or png
     fig_dir   = "./figure/"
     fignam    = "fig_2d_BLOCK_threshold_H_filtered"
     fig_type@wkWidth  = 2500
     fig_type@wkHeight = 2500 
     gridthick         = 2.0
     lthick            = 2.0
     wks               = gsn_open_wks(fig_type, fig_dir + fignam)
     plot      = new (nt, graphic)
     countour  = new (nt, graphic)
     cmap      = read_colormap_file("BlAqGrYeOrReVi200")
     cnlevs    = ispan(5100,5900,100) 
    end if 
   
;---Set common resources for all plots 
  res                  = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@gsnMaximize      = False 
  res@tmXTOn           = False
  res@tmYROn           = False
  res@gsnLeftString    = ""
  res@gsnRightString   = ""
  res@tiMainString         = ""

  BorderThick = 2.0
  res@tmBorderThicknessF    = BorderThick
  res@tmXBMajorThicknessF   = BorderThick
  res@tmXBMinorThicknessF   = BorderThick*0.5
  res@tmYLMajorThicknessF   = BorderThick
  res@tmYLMinorThicknessF   = BorderThick*0.5
  res@tmXTMajorThicknessF   = BorderThick
  res@tmXTMinorThicknessF   = BorderThick*0.5
  res@tmYRMajorThicknessF   = BorderThick
  res@tmYRMinorThicknessF   = BorderThick*0.5

  FontHeightF = 0.024
  res@tiMainFontHeightF    = FontHeightF
  res@tmYLLabelFontHeightF = FontHeightF
  res@tmXBLabelFontHeightF = FontHeightF
  res@gsnStringFontHeightF = FontHeightF*1.2
  res@tiXAxisFontHeightF   = FontHeightF
  res@tiYAxisFontHeightF   = FontHeightF

;---Control appearance of map.
  mresP = res
  mresP@gsnMajorLatSpacing          = 10
  mresP@gsnMajorLonSpacing          = 10
  mresP@gsnMinorLatSpacing          = 5
  mresP@gsnMinorLonSpacing          = 5
 ;mresP@mpProjection                = "Robinson"       ; choose projection
  mresP@mpLabelsOn                  = False
  mresP@mpPerimOn                   = False
  mresP@mpDataBaseVersion           = "MediumRes"     ; Alias 'MediumRes'
  mresP@mpDataSetName               = "Earth..1"
  mresP@mpFillOn                    = True    ; False to turn off gray continents
  mresP@mpFillDrawOrder             = "Draw"
  mresP@mpOutlineOn                 = True    ; turn on continental outline
  mresP@mpOutlineBoundarySets       = "AllBoundaries" ;"Geophysical" ;"AllBoundaries"
  mresP@mpLandFillColor             = "Grey90" ;"Yellow2enRod1"
  mresP@mpInlandWaterFillColor      = "PaleTurquoise3"
 ;mresP@mpOceanFillColor            = "PaleTurquoise3"
  mresP@mpGeophysicalLineColor      = "Grey27"
  mresP@mpGeophysicalLineThicknessF = 0.5
  mresP@mpUSStateLineColor          = "Grey27"
  mresP@mpUSStateLineThicknessF     = 0.5
  mresP@mpNationalLineColor         = "Grey27"
  mresP@mpNationalLineThicknessF    = 0.5
  mresP@mpGridAndLimbOn             = "False"
  mresP@mpGridAndLimbDrawOrder      = "PostDraw"
  mresP@mpGridMaskMode              = "MaskLand"
; mresP@mpGridSpacingF              = 5.0
  mresP@mpGridLineColor             = "Grey27"
  mresP@pmTickMarkDisplayMode       = "Always"
; mresP@mpCenterLonF                = 180.0
; mresP@mpCenterLatF                = -90.     ; This is necessary to get the correct map

 ;---Zoom in on plot
  mresP@mpMinLatF              = minLat
  mresP@mpMaxLatF              = maxLat
  mresP@mpMinLonF              = minLon
  mresP@mpMaxLonF              = maxLon
  mresP@mpCenterLonF           = 0.5*(minLon+maxLon) 
  plot(j)                      = gsn_csm_map(wks,mresP)

;---Temperature filled contour plot
  ref_res                       = res
  ref_res@cnFillOn              = True  
  ref_res@cnLevelSelectionMode  = "ExplicitLevels"
 ;ref_res@cnFillMode            = "RasterFill"
  ref_res@cnInfoLabelOn         = False
  ref_res@cnLineLabelsOn        = False
  ref_res@cnLinesOn             = False
  ref_res@cnLineLabelBackgroundColor = -1    ; transparent
  ref_res@cnFillDrawOrder       = "Draw"

  ref_res@lbLabelBarOn             = False   ; Will turn on in panel later
  ref_res@lbLabelFontHeightF       = FontHeightF*0.6
  ref_res@lbTitleFontHeightF       = FontHeightF*0.8
  ref_res@pmLabelBarParallelPosF   = 0.50 ;61
  ref_res@pmLabelBarOrthogonalPosF = 0.20
 ;ref_res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
  ref_res@pmLabelBarHeightF        = 0.12 ;0.60
  ref_res@pmLabelBarWidthF         = 0.65 ;0.06
 ;ref_res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
 ;ref_res@lbTitlePosition          = "Right"         ; title location
 ;ref_res@lbTitleDirection         = "Across"        ; letter angle
 ;ref_res@lbTitleAngleF            = 90.          ; title angle

  ref_res@cnFillPalette       := cmap 
  ref_res@cnLevels            := cnlevs 
  ref_res@cnLevelSelectionMode:= "ExplicitLevels" ;"ManualLevels"
  ref_res@cnLineThicknessF    := lthick
  ref_res@gsnRightString      := Expnam(i) 
  ref_res@gsnLeftString       := lstring+ " ("+timstr(j)+")" 

  countour(j) = gsn_csm_contour(wks,var,ref_res)
  overlay(plot(j),countour(j))
  delete([/var,res,ref_res,res/])

  ;gres                 = True
  ;gres@gsMarkerIndex   = 16                  ; filled dot
  ;gres@gsMarkerColor   = "Navy";"Magenta"
  ;gres@gsMarkerSizeF   = 0.010        ; Increase marker sizes.
  ;gres@gsMarkerThicknessF = 5.0
  ;indxx =  ind(trktime.le.timcur)
  ;if(.not.all(ismissing(indxx)))then
  ;  dum(i) = gsn_add_polymarker(wks,plot(i),psminlon(indxx),psminlat(indxx),gres)
  ;end if
  ;delete([/indxx,trktime,psminlon,psminlat/])

 end do 

 resP                     = True
 resP@gsnMaximize              = True
 resP@gsnPaperOrientation      = "portrait"
 resP@gsnPanelLabelBar         = True
 resP@pmLabelBarParallelPosF   =  0.01
 resP@pmLabelBarOrthogonalPosF = -0.02
 resP@lbLabelFontHeightF       = FontHeightF*0.4
 resP@lbTitleFontHeightF       = FontHeightF*0.45
 resP@lbTitleString            = "Z500 (gpm)" 
 resP@lbOrientation            = "horizontal"  ; vertical label bar
 resP@pmLabelBarHeightF        = 0.05
 resP@pmLabelBarWidthF         = 0.90
 resP@gsnPanelMainString       = ""
 resP@gsnPanelMainFontColor    = "Black"
 resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6
 resP@gsnPanelMainPosXF        = 0.85
;resP@txPosYF                  = 0.7
 resP@gsnPanelBottom           = 0.3
 resP@gsnPanelYWhiteSpacePercent = 0.
 resP@gsnPanelXWhiteSpacePercent = 3.
 gsn_panel(wks,plot,(/1,nt/),resP)

end do 

end
