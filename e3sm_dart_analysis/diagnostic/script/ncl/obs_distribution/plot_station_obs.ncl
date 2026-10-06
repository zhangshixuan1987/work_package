;----------------------------------------------------------------------
; station_1.ncl
;----------------------------------------------------------------------
; Concepts illustrated:
;   - Contouring one-dimensional X, Y, Z data
;   - Reading an ASCII file with several columns of data
;   - Drawing lat/lon locations as filled dots using gsn_coordinates
;   - Controlling which contour lines get drawn
;   - Using opacity to emphasize or subdue overlain features
;   - Reversing a color map
;----------------------------------------------------------------------
; This example reads in station data represented by
; 1D arrays, and generates a filled contour plot over a map.
;
; It uses a newer method (NCL V6.4.0 and later) for specifying the
; lat/lon information via special lat1d/lon1d attributes attached
; to the data to be plotted. See station_old_1.ncl if you have an 
; older version of NCL.
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "./function_read_dart_obs.ncl"

begin

  l_create_plot_data = True 
  l_create_plot_data = False
  grid     = "ne11_oQU240"
  grid     = "ne16_oQU240"
  grid     = "ne30pg2_EC30to60E2r2"
  case     = "F20TR_"+grid+"_DART80_compy"
  obsdir   = "/compyfs/zhan391/diag_dart/"+case+"/dart_diagnostics/obs_to_netcdf/"
  outdir   = "/compyfs/zhan391/diag_dart/fig_data"

  obsdate  = "2011-01-18"
  labels   = (/"Jan 18 2011"/)
  TimeUnit = "Hours since 2011-01-01 18:00:00.0"

  fname    = case+".dart.e.eam_obs_seq_final."+obsdate+"*.nc"
  files    = systemfunc ("ls "+obsdir+fname) ; file paths
  print(files)
  nfiles   = dimsizes(files)

  ObsTypeString = "all"; "RADIOSONDE_U_WIND_COMPONENT" ;"all"
  CopyString    = "observation";
  QCString      = "DART quality control"; or "Data QC"
  region        = "0 360 -90 90 0 120000" ; "lon,lat,pressure"
  l_verbose     = False ;True

  string1 = "(a) Conventional observations"
  string2 = "(b) Satellite observations"

;-------------------------------------
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and Y points.
;-------------------------------------
 ;;;conventional: "AIRCRAFT", "ACARS","RADIOSONDE","MARINE"
 ;;;satellite: "GPSRO", "SAT", "AIRS"
 Groups  = (/"AIRCRAFT", "ACARS", "RADIOSONDE","MARINE","GPSRO", "SAT", "AIRS"/)
 Labels  = (/"Aircraft", "ACARS", "Radiosonde", "Other", "GPSRO", "NESDIS", "AIRS"/)
 ngrps   = dimsizes(Groups)
 latplt  = Groups
 lonplt  = Groups
 levplt  = Groups
 timplt  = Groups

 if(l_create_plot_data) then

  do j = 0,ngrps-1

   do i  = 0,nfiles-1,1
    dartobs   = read_dart_obs_netcdf(files(i), ObsTypeString, region, CopyString, QCString, l_verbose)
   ;printVarSummary(dartobs)
    ;;;classify the data based on the type and quality mark;;;;;
    do k = 0,dimsizes(dartobs@oname2)-1
      if(isStrSubset(dartobs@oname2(k),Groups(j)))then 
       indgrp = ind(dartobs@otype.eq.dartobs@oname1(k))
       if(.not.isdefined("lat"))then 
         lat = dartobs@lat(indgrp)
         lon = dartobs@lon(indgrp)
         lev = dartobs@lev(indgrp)
         tim = dartobs@time(indgrp)
       else
         lat1 = lat 
         lon1 = lon
         lev1 = lev
         tim1 = tim
         delete([/lat,lon,lev,tim/])
         lat = array_append_record(lat1,dartobs@lat(indgrp),0)
         lon = array_append_record(lon1,dartobs@lon(indgrp),0)
         lev = array_append_record(lev1,dartobs@lev(indgrp),0)
         tim = array_append_record(tim1,dartobs@time(indgrp),0)
         delete([/lat1,lon1,lev1,tim1/])
       end if 
       delete(indgrp)
      end if 
     end do 
     delete(dartobs)
    end do 

   xx1    = lon+lat
   xq1    = get_unique_values(xx1)
   indx   = get1Dindex(xx1,xq1)
 
   fno  = outdir+"/"+Groups(j)+"_"+grid+".nc"
   system("rm -f " + fno)
   
   print(" ")
   print("output : " + fno)
   print(" ")
   
   lat!0 = "npts"
   lon!0 = "npts"
   lev!0 = "npts"
   tim!0 = "npts"

   flo = addfile(fno,"c")
   flo->time = tim(indx)
   flo->lat  = lat(indx)
   flo->lon  = lon(indx)
   flo->lev  = lev(indx)
   delete([/xx1,xq1,indx/])
   delete([/lat,lon,lev,tim/])

  end do 

 end if 

 ;read data and do plot 
 do j = 0,ngrps-1
  fno  = outdir+"/"+Groups(j)+"_"+grid+".nc"
  flo = addfile(fno,"r")
  latplt@$Groups(j)$ = flo->lat
  lonplt@$Groups(j)$ = flo->lon
  levplt@$Groups(j)$ = flo->lev
  timplt@$Groups(j)$ = flo->time
 end do

;----------------------------------------------------------------------
; Plotting section
;----------------------------------------------------------------------
  plotname = "Fig_dart_obs_"+grid+"_"+obsdate
  wtype = "png"
  wtype@wkWidth = 2500
  wtype@wkHeight = 2500 
  wks = gsn_open_wks(wtype,plotname)       ; send graphics to PNG file

  gsn_define_colormap(wks,"cb_9step")
  markers = (/16,16,16,16,16,16,16/) ;11,8,2/)       ;;line markers
  mkszf   = (/ 1.2, 1.2, 3.0, 2, 0.8, 0.5, 2/)/500.0      ;,0.015,0.015,0.015/)/100   ;;line marker size
  mkthick = (/ 3, 3, 3, 3, 3, 3, 3/)/150.0   ;;line thichness
  lnthick = (/ 3, 3, 3, 3, 3, 3, 3/)/3.0   ;;line thichness
  mkcolor = (/ 2, 6, 28, 37, 6, 28, 2/)
  mkcolor = (/ 2, 6, 28, 37, 45, 3, 77/)
  mkcolor = (/ 2, 6, 28, 37, 22, 3, 51/)
  mkcolor = (/ 55, 3, 22, 37, 2, 28, 6/)
  mkopacy = (/ 2, 1, 2, 2, 1.6, 1.6, 2.0/)/2.0

  ;;;;setup colors etc. 
  lndcolor = "black"

  res                      = True
  res@gsnMaximize          = False   ;True
  res@gsnFrame             = False   ; Want to draw markers later.
  res@gsnDraw              = False

  res@vpXF                 = 0.10
 ;res@vpYF                 = 0.18

  res@tmXTOn               = False
  res@tmYROn               = False
 ;res@tmYLOn               = True
 ;res@tmYLLabelsOn         = True

  FontHeight = 0.018
  res@tiMainFontHeightF    = FontHeight
  res@tmYLLabelFontHeightF = FontHeight
  res@tmXBLabelFontHeightF = FontHeight
  res@gsnStringFontHeightF = FontHeight*1.2

  res@tiMainOffsetYF       = 0.01

  res@tmBorderThicknessF    = 6.0
  res@tmXBMajorThicknessF   = 6.0
  res@tmYLMajorThicknessF   = 6.0
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
  res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0

  mpres                    = res
;---Control appearance of map.
  mpres@gsnMajorLatSpacing     = 30
  mpres@gsnMajorLonSpacing     = 45
  mpres@gsnMinorLatSpacing     = 15
  mpres@gsnMinorLonSpacing     = 15
 ;mpres@mpMinLonF              = 0  ;-180.
 ;mpres@mpMaxLonF              = 360;180.
  mpres@mpCenterLonF           = 180 ;150.0
 ;mpres@mpCenterLatF           = -90.     ; This is necessary to get the correct map
 ;mpres@mpProjection           = "Robinson"       ; choose projection
  mpres@mpLabelsOn             = False
  mpres@mpPerimOn              = False
  mpres@mpFillOn               = False
  mpres@mpGeophysicalLineColor = lndcolor  ; color of cont. outlines
  mpres@mpGeophysicalLineThicknessF = 2.0          ; thickness of outlines
 ;mpres@mpOutlineOn            = True
 ;mpres@mpOutlineDrawOrder     = "PostDraw"
 ;mpres@mpFillDrawOrder        = "Predraw"
 ;mpres@mpGridAndLimbOn        = True
 ;mpres@mpGridMaskMode              = "MaskLand"
 ;mpres@mpGridLatSpacingF      = 90.0
 ;mpres@mpGridLonSpacingF      = 180.0
 ;mpres@mpGridLineThicknessF   = 3.0
 ;mpres@mpGridAndLimbDrawOrder="PreDraw"
 ;mpres@mpGridLineColor        = "Transparent"
 ;mpres@mpGridLineThicknessF   = 1.0
 ;mpres@mpLimbLineColor        = "Black"
 ;mpres@mpLimbLineThicknessF   = 3.0
; mpres@mpInlandWaterFillColor = "LightBlue" ;"Blue"
; mpres@mpOceanFillColor       = "lightskyblue1"
; mpres@mpLandFillColor        = "gray33"

  mpres@tiMainString           = ""
  mpres@gsnLeftString          = string1
  map1                         = gsn_csm_map(wks,mpres)
  mpres@gsnLeftString          = string2
  map2                         = gsn_csm_map(wks,mpres)

  ;---Create the plot; it won't get drawn because gsnFrame was set to False.
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
  gsres@gsMarkerColor = "red" 
  gsres@gsMarkerSizeF = 0.0005        ; Increase marker sizes.

  txres               = True         ; additional info
  txres@txFontHeightF = FontHeight*0.45
  txres@txJust        = "CenterRight"

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Loop through each grouping of markers, and draw them one set at
  ; a time, assigning the proper color and size with gsn_marker.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  dum  = new(ngrps,graphic)
  txt  = new(ngrps,graphic)
  do j = 0,ngrps-1,1
    if (.not.all(ismissing(latplt@$Groups(j)$)))
      gsres@gsMarkerIndex      = markers(j)          ; Use filled dots for markers.
      gsres@gsMarkerColor      = mkcolor(j)
      gsres@gsMarkerSizeF      = mkszf(j) 
      gsres@gsMarkerThicknessF = mkthick(j)
      gsres@gsMarkerOpacityF   = mkopacy(j)
      if(j.le.3)then 
        dum(j) = gsn_add_polymarker(wks,map1,lonplt@$Groups(j)$,latplt@$Groups(j)$,gsres) 
      else
        dum(j) = gsn_add_polymarker(wks,map2,lonplt@$Groups(j)$,latplt@$Groups(j)$,gsres)
      end if 
    end if 
  end do 
  
;  xmkt  = (/0.05,0.165,0.275,0.40,0.55,0.735,0.90/)
;  xtxt  = (/0.14,0.255,0.375,0.475,0.65,0.830,0.975/)
;  ytxt  = (/0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25/)
; 
;  do j = 0,ngrps-1,1
;    gsres@gsMarkerIndex      = markers(j)          ; Use filled dots for markers.
;    gsres@gsMarkerColor      = mkcolor(j)
;    gsres@gsMarkerSizeF      = 0.006
;    gsres@gsMarkerThicknessF = mkthick(j)
;    text  = Labels(j)+" ("+num(.not.ismissing(latplt@$Groups(j)$))+")"
;    print("number of obs: "+ text)
;    txres@txFontColor        = mkcolor(j)
;    if (j.le.3)then
;      gsn_polymarker_ndc(wks,xmkt(j+2),ytxt(j)+0.35,gsres)
;      gsn_text_ndc      (wks,text,xtxt(j+2),ytxt(j)+0.45,txres)
;    else
;      gsn_polymarker_ndc(wks,xmkt(j-1),ytxt(j-1)+ 0.5,gsres)
;      gsn_text_ndc      (wks,text,xtxt(j-1),ytxt(j-1) + 0.5,txres)
;    end if 
;  end do 

  do j = 0,ngrps-1,1
    gsres@gsMarkerIndex      = markers(j)          ; Use filled dots for markers.
    gsres@gsMarkerColor      = mkcolor(j)
    gsres@gsMarkerSizeF      = 0.012
    gsres@gsMarkerThicknessF = mkthick(j)
    text  = Labels(j)+" ("+num(.not.ismissing(latplt@$Groups(j)$))+")"
    print("number of obs: "+ text)
    txres@txFontColor        = mkcolor(j)
    txres@txJust             = "CenterLeft"
    txres@txFontHeightF      = FontHeight*0.70
    if (j.le.3) then 
      gsn_polymarker_ndc(wks,0.80,0.90-0.04*j,gsres)
      gsn_text_ndc      (wks,text,0.82,0.90-0.04*j,txres)
    else
      gsn_polymarker_ndc(wks,0.80,0.60-0.04*j,gsres)
      gsn_text_ndc      (wks,text,0.82,0.60-0.04*j,txres)
    end if 
  end do

 ;;;panel the two figures;;;;;
    pnres                         = True
    pnres@gsnMaximize             = False ;True
    pnres@gsnPanelMainString      = "  "
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.030
;   pnres@gsnPanelBottom          = 0.05
    pnres@gsnPanelXWhiteSpacePercent = 3.0
    pnres@gsnPanelYWhiteSpacePercent = 2.0
    pnres@gsnPanelXF = (/0.10,0.10/)    ; Adjust rightmost plots
    pnres@gsnPanelYF = (/0.95,0.50/)    ; Adjust rightmost plots
    gsn_panel(wks,(/map1,map2/),(/2,1/),pnres)
    system("convert "+plotname+".png" + " -trim " +plotname+".pdf")
    system("rm -f " +plotname+".png")
end 
