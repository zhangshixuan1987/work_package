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
 datdir   = "./data"
 Groups   = (/"ERA5","CLIM","NDGUV","NDGUVTQ"/) ;,"NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
 Expnam   = (/"ERA5","CLIM","NDGUV3","NDGUVTQ3"/) ;,"NDGUVTQ3_SRF1","NDGUVTQ3_SRF2"/) 
 Expnam(2:) = Expnam(2:) + "_tau06"
 ngrps    = dimsizes(Groups)

 tim_tag  = "2007-2017"
 timstr   = (/"DJF","MAM","JJA","SON"/)
 nt       = dimsizes(timstr)
 lstring  = "Blocking Occurance (percent)" 

 ;********************************
 ; select sub-regions
 ;********************************
 minLat   = -90.             ; max lat to read in (deg)
 maxLat   =  90.             ; min lat to read in (deg)
 minLon   =   0.             ; max lon to read in (deg)
 maxLon   = 360.             ; min lon to read in (deg)
 cenlon   = 180.0 ;

 do j = 0,nt-1,1

  do i = 0,ngrps-1,1

    fl       = datdir+"/"+Expnam(i)+"_rgd_objectid_"+tim_tag+".nc"
    f        = addfile(fl,"r")
 
    lat = f->lat
    lon = f->lon

    tim = f->time
    var = f->seasonalmean_object_id(j,:,:)
    var = var*100.0

    var!0   = "lat"
    var!1   = "lon"
    var&lat = lat 
    var&lon = lon 

    if(i.eq.0) then 
      var0     = var 
    else 
      var      = var - var0
    end if 
 
    if(i.eq.0) then 
     ;;;generate the zoom in information;;;;
     fig_type  = "png" ;;figure types: pdf or png
     fig_dir   = "./figure/"
     fignam    = "fig_2d_BLOCK_climo_H_filtered_" + timstr(j)
     fig_type@wkWidth  = 2500
     fig_type@wkHeight = 2500 
     gridthick         = 2.0
     lthick            = 2.0

     wks       = gsn_open_wks(fig_type, fig_dir + fignam)
     plot      = new (ngrps, graphic)
     countour  = new (ngrps, graphic)
     cmap1     = read_colormap_file("wh-bl-gr-ye-re")
     cmap2     = read_colormap_file("BlueDarkRed18")
     cmap2(8:9,:)= 0.0 
     cnlevs1   = (/0.1,1,2,3,4,5,6,8,10,12/)
     cnlevs2   = (/-6,-4,-2,-1,-0.1,0.1,1,2,4,6/)
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
  plot(i)                      = gsn_csm_map(wks,mresP)

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

  ref_res@lbLabelBarOn             = True   ; Will turn on in panel later
  ref_res@lbLabelFontHeightF       = FontHeightF*0.6
  ref_res@lbTitleFontHeightF       = FontHeightF*0.8
  ref_res@pmLabelBarParallelPosF   = 0.50 ;61
  ref_res@pmLabelBarOrthogonalPosF = 0.10
 ;ref_res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
  ref_res@pmLabelBarHeightF        = 0.12 ;0.60
  ref_res@pmLabelBarWidthF         = 0.65 ;0.06
 ;ref_res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
 ;ref_res@lbTitlePosition          = "Right"         ; title location
 ;ref_res@lbTitleDirection         = "Across"        ; letter angle
 ;ref_res@lbTitleAngleF            = 90.          ; title angle

  ref_res@cnLevelSelectionMode := "ExplicitLevels" ;"ManualLevels"
  ref_res@cnLineThicknessF     := lthick
  ref_res@gsnRightString       := lstring+ " ("+timstr(j)+")" 

  if (i.eq.0) then 
    ref_res@cnFillPalette       := cmap1
    ref_res@cnLevels            := cnlevs1
    ref_res@gsnLeftString     := Groups(i)
  else
    ref_res@cnFillPalette       := cmap2
    ref_res@cnLevels            := cnlevs2
    ref_res@gsnLeftString      := Groups(i) + " - " + Groups(0)
  end if 

  countour(i) = gsn_csm_contour(wks,var,ref_res)
  overlay(plot(i),countour(i))
  delete([/var,lat,lon,tim,res,ref_res,res/])

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
 resP@gsnPanelLabelBar         = False
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
 gsn_panel(wks,plot,(/2,2/),resP)

end do 

end
