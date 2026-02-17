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
 mldir          = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/ml_data/"
 tcnam          = "SANDY"
 timtag         = "2012"
 date           = yyyymmddhh_time(2007, 2017, 6, "integer")
 timplt         = date({2012102118:2012103112})    ; coordinate subscripting

 figdir         = "sub_fig"
 name           = "ndguv_animation_tc"

 wtype          = "png"
 wtype@wkWidth  = 2500
 wtype@wkHeight = 2500
;wks_type                  = "pdf" ; output format, popular options are x11, png, pdf, eps
;wks_type@wkPaperHeightF   = 6.5
;wks_type@wkPaperWidthF    = 11
 wks = gsn_open_wks(wtype,figdir+"/"+name)

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

;=======================================================================================
; Basic information 
;=======================================================================================
 data_on_constant_Z    = False    ; is data already on CONSTANT Z surfaces?
 convert_hybridP_to_Z  = True     ; is data on hybrid pressure levels?
 hasTimeIx             = True     ; does file have time index?
 timeStep              = 0        ; If yes, what time index do you want to plot?
 model                 = "E3SMv2" ; used for mainStr, but also for model specific if statements
 mainStr               = model    ; main string for plot titles
 gravit                = 9.81     ; gravational constant 

 Uname                 = "U"      ; Variable name of zonal wind (typically "U" or "ua")
 Vname                 = "V"      ; Variable name of meridional wind (typically "V" or "ua")
 PSname                = "PS"     ; Variable name of surface pressure (typically "PS","PSL","ps",or "slp")
 ZSname                = "PHIS"   ; Variable name of surface pressure (typically "PHIS")
 Tname                 = "T"      ; Variable name of air temperature (typically "T" or "ta")
 Zname                 = "Z3"     ; Variable name of height (typically "Z" or "za")
 Qname                 = "Q"      ; Variable name of meridional wind (typically "V" or "ua")

;=======================================================================================
; Other settings, required to remain constant for intercomparison
; Generally, don't touch unless experimenting/debugging
;=======================================================================================
 minLat                = 10.               ; max lat to read in (deg)
 maxLat                = 45.              ; min lat to read in (deg)
 minLon                = 260.             ; max lon to read in (deg)
 maxLon                = 310.             ; min lon to read in (deg)
 thetaWindMax          = 60.0             ; max theta wind to plot (m/s)
 radMinMax             = 18.0             ; min/max for radial plot (m/s)
 tAnomMax              = 15.0             ; max for tAnom plots (K)
 offAnomDist           = 5.0              ; offset distance for anom calcs (deg)
 zintmin               = 20.0             ; min height for z interpolation (m)
 zintmax               = 18000.0          ; max height for z interpolation (m)
 nzint                 = 100              ; num points for z interpolation
 Zmidpt                = 2500.0           ; half of Z points BELOW this height, half ABOVE

;=======================================================================================
; experiments for plot
;=======================================================================================
 Groups   = (/"ERA5", "E3SMv2_CLIM", "E3SMv2_NDGUV_tau6", "E3SMv2_NDGUV_tau6"/)
 Casdir   = (/"ERA5", "E3SMv2_CLIM", "NDGUV3_tau06",      "NDGUV3_tau06"/)
 Expnam   = (/"ERA5", "E3SM (CLIM)", "E3SM (NDGUV)",      "E3SM (DeepONet)"/)
 ngrps    = dimsizes(Groups)
 plot     = new(ngrps,"graphic")
 contour0 = new(ngrps,"graphic")
 contour1 = new(ngrps,"graphic")
 contour2 = new(ngrps,"graphic")
 vector0  = new(ngrps,"graphic")
 dum      = new(ngrps,"graphic")

 ;;read land mask;;;;;
 f0       = addfile(datdir+"lmsk.nc","r")
 lmsk     = f0->LANDFRAC(0,{minLat:maxLat},{minLon:maxLon})

 f        = addfile(datdir+"E3SMv2_CLIM.2012102118-2012103112.nc","r")
 time     = f->time
 nt       = dimsizes(time)-1
 PHIS     = f->$ZSname$(0,{minLat:maxLat},{minLon:maxLon})
 ZS       = PHIS
 ZS       = ZS/gravit ; conver to m

 do ti = 0,nt,1

  timeStep = ti  ; If yes, what time index do you want to plot?
  utc_date = cd_calendar(time(timeStep), 0)
  year     = tointeger(utc_date(:,0))    ; Convert to integer for
  month    = tointeger(utc_date(:,1))    ; use sprinti 
  day      = tointeger(utc_date(:,2))
  hour     = tointeger(utc_date(:,3))
  minute   = tointeger(utc_date(:,4))
  second   = utc_date(:,5)
  timstr   = sprinti("%0.2iZ ", hour) + sprinti("%0.2i ", day) + \
             month_abbr(month) + " "  + sprinti("%0.4i", year)
  trkcurt  = cd_calendar(time(timeStep), -3)
  delete([/year,month,day,hour,minute,second/])

  print("timeIndex = "+timeStep + ", date: " +  timstr)
  
  do ig = 0,ngrps-1,1
 
  filename = datdir+Groups(ig)+".2012102118-2012103112.nc"
  f        = addfile(filename,"r")

  if(ig.eq.ngrps-1) then 
    fmlname = mldir + "dn_un_norm.nc"
    f0      = addfile(fmlname,"r")
    f01     = addfile("./time.nc","r")
    mltim   = cd_calendar(f01->time(2312:2559),-3)
    indtml  = ind(trkcurt.eq.mltim)
   ;print(mltim(indtml))
  end if 

;=======================================================================================
; Get file, coordinate variables, and U, V, PSL
;=======================================================================================
  print("Loading data from file...")

  lat   = f->lat({minLat:maxLat})
  lon   = f->lon({minLon:maxLon})
  lev   = f->lev(:)

  nlat  = dimsizes(lat)
  nlon  = dimsizes(lon)
  nlev  = dimsizes(lev)

  if (Groups(ig).eq."ERA5") then 
    Uname1  = Uname+"_ref"
    Vname1  = Vname+"_ref"
    Tname1  = Tname+"_ref"
    Qname1  = Qname+"_ref"
    PSname1 = PSname+"_ref"
  else if (Groups(ig).eq."E3SMv2_CLIM") then
    Uname1  = Uname+""
    Vname1  = Vname+""
    Tname1  = Tname+""
    Qname1  = Qname+""
    PSname1 = PSname+""
  else 
    Uname1  = Uname+"_bf_ndg"
    Vname1  = Vname+"_bf_ndg"
    Tname1  = Tname+"_bf_ndg"
    Qname1  = Qname+"_bf_ndg"
    PSname1 = PSname+"_bf_ndg"
  end if 
  end if 

  U     = f->$Uname1$(timeStep,:,{minLat:maxLat},{minLon:maxLon})
  V     = f->$Vname1$(timeStep,:,{minLat:maxLat},{minLon:maxLon})
  T     = f->$Tname1$(timeStep,:,{minLat:maxLat},{minLon:maxLon})
  Q     = f->$Qname1$(timeStep,:,{minLat:maxLat},{minLon:maxLon})
  PS    = f->$PSname1$(timeStep,{minLat:maxLat},{minLon:maxLon})

  if (ig.eq.0)then 
    U0  = U 
    V0  = V
  end if 

  if (ig.eq.ngrps-1) then
   ;printVarSummary(f0->U)
   ;printVarSummary(f0->V)
    tmpu     = (/f0->U/) 
    tmpu!0   = "time"
    tmpu!1   = "lev"
    tmpu!2   = "lat"
    tmpu!3   = "lon"
    tmpu&lat = f->lat(90:159)
    tmpu&lon = f->lon(250:319)
    tmpv     = tmpu
    tmpv     = (/f0->V/)
    U(10:71,:,:) = (/tmpu(indtml,:,{minLat:maxLat},{minLon:maxLon})/)
    V(10:71,:,:) = (/tmpv(indtml,:,{minLat:maxLat},{minLon:maxLon})/)
    delete([/tmpu,tmpv/])
  end if 

 ;If U and V are not m/s, please convert here
  U@units  = "m/s"
  V@units  = "m/s"
  Q@units  = "kg/kg"
  T@units  = "K"
  T@long_name="Temperature"

;=======================================================================================
; Do Z calculations/interpolations if necessary
;=======================================================================================
  if (convert_hybridP_to_Z) then 
    pm   = T
    pm   = pres_hybrid_ccm (PS,f->P0,f->hyam(::-1),f->hybm(::-1)) 
    dp   = T
    dp   = dpres_hybrid_ccm(PS,f->P0,f->hyai,f->hybi)  ; Pa [kg/(m s2)]   
    Tv   = T
    Tv   = T*(1.+Q*0.61) 
    Z    = cz2ccm(PS,PHIS,Tv,f->P0,f->hyam(::-1),f->hybm(::-1),f->hyai(::-1),f->hybi(::-1))
   ;printMinMax(pm,False)
    Tv@units = "K"
    Z@units  = "gpm"
    pm       = pm / 100.0 ; to hPa
    pm@units = "hPa"
    delete([/pm,dp,Tv/])

  else 
    Z    = f->$Zname$(timeStep,:,{minLat:maxLat},:) 
    Z    = Z / gravit
  end if 
  
  Zlev       = (/100./) 
  Zlev@units = "m"
  Zlev!0     = "lev"
  Zlev&lev   = Zlev

;=======================================================================================
; Interpolate lat/lon variables to constant Z levels
;=======================================================================================
;  print("Interpolating to Z surfaces")
;  U_X     = int2p_n_Wrap(Z,U,Zlev,2,0) 
;  V_X     = int2p_n_Wrap(Z,V,Zlev,2,0) 
;  U_Z     = U_X(0,:,:)
;  V_Z     = V_X(0,:,:)

   lev_p   = 850.0
   U_X     = vinth2p(U, f->hyam, f->hybm, lev_p, PS, 2, f->P0/100.0, 1, False )
   V_X     = vinth2p(V, f->hyam, f->hybm, lev_p, PS, 2, f->P0/100.0, 1, False )
   U_Z     = U_X(0,:,:)
   V_Z     = V_X(0,:,:)
  delete([/U_X,V_X/])

  WIND_Z  = U_Z
  WIND_Z  = sqrt(U_Z^2+V_Z^2)
  WIND_Z@long_name="horizontal wind speed"

  U_Z@_FillValue    = -9999.0 
  V_Z@_FillValue    = -9999.0
  WIND_Z@_FillValue = -9999.0

  print("Plotting Range...")
  printMinMax(U_Z,False)
  printMinMax(V_Z,False)
  printMinMax(WIND_Z,False)

  ;;;mask land 
  ;U_Z     = where(lmsk.eq.1.0, -9999.0, U_Z)
  ;V_Z     = where(lmsk.eq.1.0, -9999.0, V_Z)
  ;WIND_Z  = where(lmsk.eq.1.0, -9999.0, WIND_Z)


;calculate the sea level pressure
  hyam = f->hyam   ; read from a file the mid-layer coef
  hybm = f->hybm   ; read from a file
  p0   = f->P0
  PRES = pres_hybrid_ccm(PS,p0,hyam,hybm)  ; (ntim,klvl,nlat,lmon), [Pa]
  SLP  = pslec(T(nlev-1,:,:), PHIS, PS, PRES(nlev-1,:,:))
  SLP  = SLP/100.0
  SLP@units = "hPa"
  delete([/hyam,hybm,p0,PRES/])
  SLP = where(lmsk.eq.1.0, -9999, SLP)
  SLP@_FillValue = -9999

  PS         = PS /100.0 ; conver to hPa
  PS@units   = "hPa"
  PS         = where(lmsk.eq.1.0, -9999, PS)
  PS@_FillValue = -9999
 ;printMinMax(SLP,False)
 ;printMinMax(PS,False)


;=======================================================================================
; Read storm track
;=======================================================================================
 if (ig.eq.0) then 
  fna    = datdir+tcnam+"_"+Casdir(ig)+"_Track_info_"+ timtag+ ".nc"
  fla    = addfile(fna, "r")
  tmp    = fla->stim
  indxx  = ind(.not.ismissing(tmp))
  if(.not.all(ismissing(indxx))) then
    trktime  = cd_calendar(fla->stim(indxx),-3)
    psminlon = fla->lon(indxx)
    psminlat = fla->lat(indxx)
  end if
  delete([/tmp,indxx/])

 ;print("... PS minimum found at lat: "+psminlat+" lon: "+psminlon)
 ;print("Current track location: " + trkcurt)

end if 


;;;start to plot 
 res                  = True
 res@gsnDraw          = False
 res@gsnFrame         = False
;res@gsnSpreadColors  = True        ; Span full color map
 res@gsnRightString   = ""
 res@gsnLeftString    = "" 

 res@tmYLLabelsOn           = True
 res@tmYLOn                 = True
 res@tmYROn                 = False
 res@tmXTOn                   = False

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

 res0   = res
 res1   = res
 res2   = res
 res_vc = res  
 mpres  = res


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
  plot(ig)               = gsn_csm_map(wks,mpres)


;;set vector;;
 res_vc@vcGlyphStyle              = "LineArrow"
 res_vc@vcLineArrowThicknessF     = 4
 res_vc@vcMinDistanceF            = 0.025
 res_vc@vcRefLengthF              = 0.04
 res_vc@vcRefAnnoOn               = True 
 res_vc@vcRefMagnitudeF           = 20
 res_vc@vcRefAnnoString1          = res_vc@vcRefMagnitudeF+" m s~S~-1~N~"
 res_vc@vcRefAnnoSide             = "Top"
 res_vc@vcRefAnnoString2On        = False
 res_vc@vcRefAnnoPerimOn          = False
 res_vc@vcRefAnnoOrthogonalPosF   = -0.10 ;999
 res_vc@vcRefAnnoParallelPosF     = 0.98
 res_vc@vcRefAnnoFontHeightF      = FontHeightF*0.6
 res_vc@vcRefAnnoBackgroundColor  = "White"
 res_vc@vcVectorDrawOrder         = "PostDraw"
 vector0(ig)  = gsn_csm_vector(wks,U_Z(:,:),V_Z(:,:),res_vc)
 overlay(plot(ig),vector0(ig))

 res1@cnFillOn                 = True        ; Turn on contour fill
 res1@cnFillPalette            = "WhiteBlueGreenYellowRed" ;"precip3_16lev" ;"WhBlGrYeRe" ;"NCV_jaisnd" ;"WhBlGrYeRe" ;"BlAqGrYeOrReVi200"
 res1@cnLevelSelectionMode     = "ManualLevels"
 res1@cnLevelSpacingF          =  0.5
 res1@cnMinLevelValF           =  0 
 res1@cnMaxLevelValF           =  30.0
 res1@tiMainString             = Expnam(ig) 
 res1@gsnLeftString            = "Wind at 850 hPa (m s~S~-1~N~)" 
 res1@gsnLeftStringOrthogonalPosF = 0.03
 res1@gsnRightString           = "" ;timstr
 res1@cnLinesOn                = False   ; turn off contour lines
 res1@cnLineLabelsOn           = False   ; turn off contour labels
 res1@cnInfoLabelOn            = False
 res1@cnLineLabelDensityF      = 2.0
 res1@cnLineLabelFontHeightF   = FontHeightF
 res1@cnLineLabelBackgroundColor     = -1
 res1@gsnContourLineThicknessesScale = 2.0
 res1@cnLineThicknessF               = 0.8
 res1@cnLineColor                    = "Black"
 res1@gsnContourZeroLineThicknessF   = 2.
 res1@gsnContourNegLineDashPattern   = 1

 res1@lbLabelBarOn             = False   ; Will turn on in panel later
 res1@lbLabelFontHeightF       = FontHeightF*0.6
 res1@lbTitleFontHeightF       = FontHeightF*0.8
 res1@pmLabelBarParallelPosF   = 0.50 ;61
 res1@pmLabelBarOrthogonalPosF = 0.20
 res1@pmLabelBarHeightF        = 0.12 ;0.60
 res1@pmLabelBarWidthF         = 0.65 ;0.06
 contour1(ig) = gsn_csm_contour(wks,WIND_Z(:,:),res1)
 overlay(plot(ig),contour1(ig))

 res2@cnFillOn                       = False        ; Turn on contour fill
 res2@cnLinesOn                      = False    ; turn off contour lines
 res2@cnLineLabelsOn                 = False   ; turn off contour labels
 res2@cnInfoLabelOn                  = False
 res2@cnLineLabelDensityF            = 2.0
 res2@cnLevelSpacingF                = 4.0
 res2@cnMinLevelValF                 = 960 
 res2@cnMaxLevelValF                 = 980.0
 res2@cnLineLabelFontHeightF         = FontHeightF
 res2@cnLineLabelBackgroundColor     = -1
 res2@gsnContourLineThicknessesScale = 1.0
 res2@cnLineThicknessF               = 3.0
 res2@cnLineColor                    = "Black"
 contour2(ig) = gsn_csm_contour(wks,WIND_Z(:,:),res2)
 overlay(plot(ig),contour2(ig))

 gres                 = True
 gres@gsMarkerIndex   = 16                  ; filled dot
 gres@gsMarkerColor   = "Magenta"
 gres@gsMarkerSizeF   = 0.007        ; Increase marker sizes.
 indxx =  ind(trktime.le.trkcurt)
 if(.not.all(ismissing(indxx)))then 
   dum(ig) = gsn_add_polymarker(wks,plot(ig),psminlon(indxx),psminlat(indxx),gres)
 end if 
 delete([/indxx/]) ;trktime,psminlon,psminlat/])

 delete([/res,res0,res1,res_vc/])

 end do 

 resP                          = True                ; modify the panel plot
 resP@gsnMaximize              = True
 resP@gsnPaperOrientation      = "portrait"
 resP@gsnPanelLabelBar         = True
 resP@pmLabelBarParallelPosF   =  0.01
 resP@pmLabelBarOrthogonalPosF = -0.02
 resP@lbLabelFontHeightF       = FontHeightF*0.4
 resP@lbTitleFontHeightF       = FontHeightF*0.45
 resP@lbTitleString            = "Wind at 850 hPa (m s~S~-1~N~)" 
 resP@lbOrientation            = "horizontal"  ; vertical label bar
 resP@lbBoxLinesOn             = False
 resP@pmLabelBarHeightF        = 0.04
 resP@pmLabelBarWidthF         = 0.95
 resP@gsnPanelMainString       = "Current time: " + timstr
 resP@gsnPanelMainFontColor    = "Black"
 resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6 
 resP@gsnPanelYWhiteSpacePercent = 0.
 resP@gsnPanelXWhiteSpacePercent = 5.

 gsn_panel(wks,plot,(/1,ngrps/),resP)             ; now draw as one plot

end do

end
