;=======================================================================================
; This NCL code calculates radially-averaged tangential and radial wind components
; as well as T anomaly for DCMIP test case #2 (cyclone)
; this code requires the accompanying function set "radialAvg.ncl"
;
; Usage: User should modify "user options" for their particular data set. Currently,
; U, V, T, PS are required as variables.
; If variables are on constant Z surfaces, life is easy.
;
; Grepping for "MODELSPEC" will point to possible areas of the code in need of modification
; for model specific output
;
; Written by Colin Zarzycki (zarzycki@ucar.edu)
; Version 0.1 (6/5/2016) - DCMIP-2016 release
;=======================================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

;=======================================================================================
; User options
;=======================================================================================

 datdir         = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS/sandy/"
 tcnam          = "SANDY"
 timtag         = "2012"
 date           = yyyymmddhh_time(2007, 2017, 6, "integer")
 timplt         = date({2012102118:2012103112})    

;=======================================================================================
; Other settings, required to remain constant for intercomparison
; Generally, don't touch unless experimenting/debugging
;=======================================================================================
 minLat         = 10.               ; max lat to read in (deg)
 maxLat         = 45.              ; min lat to read in (deg)
 minLon         = 260.             ; max lon to read in (deg)
 maxLon         = 310.             ; min lon to read in (deg)

;=======================================================================================
; experiments for plot
;=======================================================================================
 Groups   = (/"ERA5", "E3SMv2_CLIM", "E3SMv2_NDGUV_tau6", "E3SMv2_NDGUVTQ_tau6"/)
 Casdir   = (/"ERA5", "E3SMv2_CLIM", "NDGUV3_tau06",      "NDGUVTQ3_tau06"/)
 Expnam   = (/"ERA5", "E3SM (CLIM)", "E3SM (NDGUV)",      "E3SM (NDGUVTQ)"/)
 ngrps    = dimsizes(Groups)

 do ig = 1,ngrps-1,1
 
  fl0 = datdir+Groups(0)+".2012102118-2012103112.nc"
  f0  = addfile(fl0,"r")

  fl1 = datdir+Groups(ig)+".2012102118-2012103112.nc"
  f1  = addfile(fl1,"r")

 ;=======================================================================================
 ; Get file, coordinate variables, and U, V, PSL
 ;=======================================================================================
  print("Loading data from file...")

  lat   = f0->lat({minLat:maxLat})
  lon   = f0->lon({minLon:maxLon})
  lev   = f0->lev(:)

  nlat  = dimsizes(lat)
  nlon  = dimsizes(lon)
  nlev  = dimsizes(lev)

  lev_p   = 850.0
  U_r     = vinth2p(f0->U_ref(:,:,{minLat:maxLat},{minLon:maxLon}), f0->hyam, f0->hybm, lev_p, \
                    f0->PS_ref(:,{minLat:maxLat},{minLon:maxLon}), 2, f0->P0/100.0, 1, False )
  if(Groups(ig).eq."E3SMv2_CLIM") then 
    U_m   = vinth2p(f1->U(:,:,{minLat:maxLat},{minLon:maxLon}), f1->hyam, f1->hybm, lev_p, \
                    f1->PS(:,{minLat:maxLat},{minLon:maxLon}), 2, f0->P0/100.0, 1, False )
  else 
    U_m   = vinth2p(f1->U_bf_ndg(:,:,{minLat:maxLat},{minLon:maxLon}), f1->hyam, f1->hybm, lev_p, \
                    f1->PS_bf_ndg(:,{minLat:maxLat},{minLon:maxLon}), 2, f0->P0/100.0, 1, False )
  end if 

  U_b     = U_r(:,0,:,:)
  U_b     = U_m(:,0,:,:) - U_r(:,0,:,:)
 
  rad    = 4.*atan(1.)/180.
  clat   = lat            
  clat   = sqrt(cos(rad*clat) )                 ; gw for gaussian grid
  xw     = U_b*conform(U_b, clat, 1)
  copy_VarMeta(U_b, xw)
  
  neof    = 20                                          ; user specified (commonly <=5) 
  eof     = eofunc_n_Wrap(xw, neof, False,0)
  eof_ts  = eofunc_ts_n_Wrap (xw, eof, False,0)
  eof_ts  = dim_standardize_n( eof_ts, 0, 1)      ; normalize
  eof_rg  = eof
  eof_rg  = (/regCoef_n(eof_ts,U_b, 1, 0)/)
  lab_str = "~F33~D~N~~F~U (m s~S~-1~N~)"
  delete([/U_r,U_m,clat,xw,eof,lat,lon,lev/])

 do ie = 0,neof-1,1

  eval    = eof_rg@eval
  pcvar   = cumsum(eof_rg@pcvar,0) 
  uavg    = dim_avg_n_Wrap(U_b,0)
  ueof    = eof_rg(ie,:,:)
 
  figdir         = "./"
  fignam         = tcnam+"_U_eof"+sprinti("%02d",ie)+"_"+Groups(ig)
  wtype          = "png"
  wtype@wkWidth  = 2500
  wtype@wkHeight = 2500
  wks = gsn_open_wks(wtype,figdir+"/"+fignam)
 
  plot     = new(3,"graphic")
  contour0 = new(3,"graphic")
  contour1 = new(3,"graphic")
  contour2 = new(3,"graphic")
    
 ;gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

  res                  = True
  res@gsnDraw          = False
  res@gsnFrame         = False
; res@gsnSpreadColors  = True        ; Span full color map
  res@gsnRightString   = ""
  res@gsnLeftString    = "" 

  res@tmYLLabelsOn           = True
  res@tmYLOn                 = True
  res@tmYROn                 = False
  res@tmXTOn                   = False

  BorderThick = 4.0
  res@tmBorderThicknessF    = BorderThick
  res@tmXBMajorThicknessF   = BorderThick
  res@tmXBMinorThicknessF   = BorderThick*0.5
  res@tmYLMajorThicknessF   = BorderThick
  res@tmYLMinorThicknessF   = BorderThick*0.5
  res@tmXTMajorThicknessF   = BorderThick
  res@tmXTMinorThicknessF   = BorderThick*0.5
  res@tmYRMajorThicknessF   = BorderThick
  res@tmYRMinorThicknessF   = BorderThick*0.5

  FontHeightF = 0.020
  res@tiMainFontHeightF    = FontHeightF
  res@tmYLLabelFontHeightF = FontHeightF
  res@tmXBLabelFontHeightF = FontHeightF
  res@gsnStringFontHeightF = FontHeightF*1.2

  mpres  = res
  res1   = res

;---Control appearance of map.
  mpres@gsnMajorLatSpacing     = 10
  mpres@gsnMajorLonSpacing     = 10
  mpres@gsnMinorLatSpacing     = 5
  mpres@gsnMinorLonSpacing     = 5
 ;mpres@mpProjection           = "Robinson"       ; choose projection
  mpres@mpLabelsOn             = False
  mpres@mpPerimOn              = False
  mpres@mpFillOn               = False
  mpres@mpGeophysicalLineColor = "Grey25"
  mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
 ;mpres@mpOutlineOn            = True
 ;mpres@mpOutlineDrawOrder     = "PostDraw"
 ;mpres@mpFillDrawOrder        = "Predraw"
  mpres@mpOutlineOn            = True
  mpres@mpOutlineBoundarySets  = "National"
  mpres@mpDataBaseVersion      = "MediumRes"   ; necessary for mpDataSetName to be effective
  mpres@mpLandFillColor        = "bisque2"
 ;mpres@mpGridAndLimbOn        = True
 ;mpres@mpGridLatSpacingF      = 90.0
 ;mpres@mpGridLonSpacingF      = 180.0
 ;mpres@mpGridLineThicknessF   = 3.0
 ;mpres@mpGridAndLimbDrawOrder="PreDraw"
 ;mpres@mpMinLonF              = 0  ;-180.
 ;mpres@mpMaxLonF              = 360;180.
  mpres@mpCenterLonF           = 180.0
; mpres@mpCenterLatF           = -90.     ; This is necessary to get the correct map
; mpres@mpOceanFillColor       = "lightskyblue1"
; mpres@mpLandFillColor        = "gray33"

  mpres@mpMaxLatF        = maxLat 
  mpres@mpMinLatF        = minLat 
  mpres@mpMaxLonF        = maxLon 
  mpres@mpMinLonF        = minLon 
  mpres@mpCenterLonF     = 0.5*(minLon+maxLon) ;psminlon
  plot(0)                = gsn_csm_map(wks,mpres)
  plot(1)                = gsn_csm_map(wks,mpres)

  res1@cnFillOn                 = True        ; Turn on contour fill
  res1@cnFillPalette            = "NCV_jaisnd" ;"WhBlGrYeRe" ;"BlAqGrYeOrReVi200"
 ;res1@cnLevelSelectionMode     = "ManualLevels"
 ;res1@cnLevelSpacingF          =  2.0
 ;res1@cnMinLevelValF           =  960 ;989.0
 ;res1@cnMaxLevelValF           =  1040.0
  res1@tiMainString             = "" 
  res1@cnLevelSelectionMode     = "ExplicitLevels"
  res1@cnLevels                 =  fspan(-3,3,11)
  res1@gsnLeftStringOrthogonalPosF = 0.03
  res1@gsnRightString           = "" ;timstr
  res1@cnLinesOn                = True   ; turn off contour lines
  res1@cnLineLabelsOn           = False   ; turn off contour labels
  res1@cnInfoLabelOn            = False
  res1@cnLineLabelDensityF      = 2.0
  res1@cnLineLabelFontHeightF   = FontHeightF
  res1@cnLineLabelBackgroundColor  = -1
  res1@gsnContourLineThicknessesScale = 2.0
  res1@cnLineThicknessF               = 0.8
  res1@cnLineColor                    = "Black"
  res1@gsnContourZeroLineThicknessF   = 2.
  res1@gsnContourNegLineDashPattern   = 1
  res1@lbLabelBarOn             = True   ; Will turn on in panel later
  res1@lbLabelFontHeightF       = FontHeightF*0.6
  res1@lbTitleFontHeightF       = FontHeightF*0.8
 ;res1@pmLabelBarParallelPosF   = 0.50 ;61
 ;res1@pmLabelBarOrthogonalPosF = 0.20
  res1@lbOrientation            = "vertical"  ; vertical label bar
  res1@pmLabelBarHeightF        = 0.60
  res1@pmLabelBarWidthF         = 0.06
  res1@lbTitleString            = lab_str 
  res1@lbTitlePosition          = "Right"         ; title location
  res1@lbTitleDirection         = "Across"        ; letter angle
  res1@lbTitleAngleF            = 90.          ; title angle

  res1@gsnRightString           = Expnam(ig)
  res1@gsnLeftString            = "Mean difference (E3SM - ERA5)"
  contour1(0) = gsn_csm_contour(wks,uavg,res1)

  res1@gsnRightString           = Expnam(ig)
  res1@gsnLeftString            = "EOF analysis (mode"+sprinti("%02d",(ie+1))+")"
  contour1(1) = gsn_csm_contour(wks,ueof,res1)

  overlay(plot(0),contour1(0))
  overlay(plot(1),contour1(1))

; resources for "left" variable
  resL                        = res 
  resL@vpHeightF              = 0.6        ; Changes the aspect ratio
  resL@vpWidthF               = 0.7
 ;resL@vpXF                   = 0.10        ; change start locations
 ;resL@vpYF                   = 0.75        ; the plot
  resL@tmYLOn                 = True
  resL@tmYROn                 = False
  resL@tmXTOn                 = False
  resL@xyLineThicknesses      = 8.                  ; thicker line
  resL@xyDashPatterns         = 0                   ; dashed line for 2nd
  resL@xyLineColors           = "red"                    ; line color
  resL@tiYAxisString          = "Eigenvalues" 
  resL@trYMaxF                = ceil(max(eval))                    ; axis max
  resL@trYMinF                = floor(min(eval))   
  resL@tiMainString           = "" 
  resL@tiXAxisString          = "EOF mode"

; resources for "right" variable
  resR                        = res
  resR@tmYLOn                 = False
  resR@tmYROn                 = True
  resR@tmXTOn                 = False
  resR@xyLineThicknesses      = 8                   ; thicker line
  resR@xyDashPatterns         = 1                   ; dashed line for 2nd
  resR@xyLineColors           = "Blue"                    ; line color
  resR@tiYAxisString          = "Cumulated percent variance" 
  resR@trYMaxF                = 100
  resR@trYMinF                = 0 
  resR@tiMainString           = "" 
 
  xx      = ispan(1,neof,1)
  plot(2) = gsn_csm_xy2(wks,xx,eval,pcvar,resL,resR)
  
  delete([/res,res1,resR,resL/])

  resP                          = True                ; modify the panel plot
  resP@gsnMaximize              = True
  resP@gsnPaperOrientation      = "portrait"
  resP@gsnPanelLabelBar         = False
  resP@gsnPanelMainString       = "" ;Expnam(ig)
  resP@gsnPanelMainFontColor    = "Black"
  resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6 
 ;resP@gsnPanelXWhiteSpacePercent = 4	; set a bit of extra white space between panels in the x and y directions
 ;resP@gsnPanelYWhiteSpacePercent = 4
  resP@gsnPanelYF                 = (/0.95,0.60,0.30/)      ; Adjust bottommost plots.

  gsn_panel(wks,plot(::-1),(/3,1/),resP)             ; now draw as one plot

 end do

end do

end
