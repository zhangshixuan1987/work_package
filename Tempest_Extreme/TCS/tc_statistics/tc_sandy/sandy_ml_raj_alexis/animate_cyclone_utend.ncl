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
 mldir          = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/ml_data/raj/"
 tcnam          = "SANDY"
 timtag         = "2012"
 date           = yyyymmddhh_time(2007, 2017, 6, "integer")
 timplt         = date({2012102118:2012103112})    ; coordinate subscripting

 figdir         = "sub_fig"
 name           = "utend_animation_tc"

 wtype                     = "png"
 wtype@wkWidth             = 2500
 wtype@wkHeight            = 2500
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

 vname                 = "U"      ; Variable name of zonal wind (typically "U" or "ua")
 PSname                = "PS"     ; Variable name of surface pressure (typically "PS","PSL","ps",or "slp")

;=======================================================================================
; Other settings, required to remain constant for intercomparison
; Generally, don't touch unless experimenting/debugging
;=======================================================================================
 minLat                = 10.               ; max lat to read in (deg)
 maxLat                = 45.              ; min lat to read in (deg)
 minLon                = 260.             ; max lon to read in (deg)
 maxLon                = 310.             ; min lon to read in (deg)

;=======================================================================================
; experiments for plot
;=======================================================================================
 Groups   = (/"E3SMv2_NDGUV_tau6", "E3SMv2_NDGUV_tau6",            "E3SMv2_NDGUV_tau6"/)
 Casdir   = (/"NDGUV3_tau06",      "NDGUV3_tau06",                 "NDGUV3_tau06"/)
 Expnam   = (/"E3SM (NDGUV)",      "E3SM (DeepONet with Nudge)",   "E3SM (DeepONet with ERA5)"/)
 Expnam   = (/"E3SM (NDGUV)",      "DeepONet w/t Nudge",   "DeepONet w/t ERA5"/)
 ngrps    = dimsizes(Groups)
 plot     = new(ngrps,"graphic")
 contour0 = new(ngrps,"graphic")
 contour1 = new(ngrps,"graphic")
 contour2 = new(ngrps,"graphic")
 dum      = new(ngrps,"graphic")

 f        = addfile(datdir+Groups(0)+".ngdtend.2012102118-2012103112.nc","r")
 time     = f->time
 nt       = dimsizes(time)-1

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
 
  filename = datdir+Groups(ig)+".ngdtend.2012102118-2012103112.nc"
  f        = addfile(filename,"r")

  filename = datdir+Groups(ig)+".2012102118-2012103112.nc"
  fx       = addfile(filename,"r")

  if(ig.eq.ngrps-2) then 
    fmlname = mldir + "dn_nudge.nc"
    f0      = addfile(fmlname,"r")
    f01     = addfile("./time.nc","r")
    mltim   = cd_calendar(f01->time(2312:2559),-3)
    indtml  = ind(trkcurt.eq.mltim)
   ;print(mltim(indtml))
  end if 

  if(ig.eq.ngrps-1) then
    fmlname = mldir + "dn_nudge_era.nc"
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

  Uname1  = "Nudge_"+vname
  PSname1 = PSname+"_bf_ndg"

  U     = f->$Uname1$(timeStep,:,{minLat:maxLat},{minLon:maxLon})
  PS    = fx->$PSname1$(timeStep,{minLat:maxLat},{minLon:maxLon})

  if (ig.ge.ngrps-2) then
   ;printVarSummary(f0->U)
   ;printVarSummary(f0->V)
    tmpu     = (/f0->U/) 
    tmpu!0   = "time"
    tmpu!1   = "lev"
    tmpu!2   = "lat"
    tmpu!3   = "lon"
    tmpu&lat = f->lat(90:159)
    tmpu&lon = f->lon(250:319)
    U(:,:,:) = (/tmpu(indtml,:,{minLat:maxLat},{minLon:maxLon})/)
    delete([/tmpu/])
  end if 

 ;If U and V are not m/s, please convert here
  U = U * 1e4 
  U@units  = "10~S~-4~N~ m s~S~-2~N~"

;=======================================================================================
; Interpolate lat/lon variables to constant Z levels
;=======================================================================================
 lev_p   = 500.0
 U_X     = vinth2p(U, f->hyam, f->hybm, lev_p, PS, 2, f->P0/100.0, 1, False )
 U_Z     = U_X(0,:,:)
 delete([/U_X/])

 U_Z@_FillValue    = -9999.0 
;U_Z  = where(lmsk.eq.1.0, -9999.0, U_Z)

 rad  = 4.0*atan(1.0)/180.0
 gw   = cos(U_Z&lat*rad)
 
 if (ig.eq.0) then
   U_0 = U_Z
 end if 

 if(ig.gt.0) then 
   pcor = sprintf("%.2f", pattern_cor(U_Z, U_0,gw, 1))
 else
   pcor = "" 
 end if 
 
;=======================================================================================
; Read storm track
;=======================================================================================
 if (ig.eq.0) then
  fna    = datdir+tcnam+"_ERA5_Track_info_"+ timtag+ ".nc"
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

 BorderThick = 3.0
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

  res1@cnFillOn                      = True        ; Turn on contour fill
  res1@cnFillPalette                  = "MPL_RdBu" 
  res1@tiMainString                   = "" ;Expnam(ig) 
  res1@cnLinesOn                      = False   ; turn off contour lines
  res1@cnLineLabelsOn                 = False   ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@cnLineLabelDensityF            = 2.0
  res1@cnLineLabelFontHeightF         = FontHeightF
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


  res1@gsnRightStringOrthogonalPosF = 0.03
  res1@gsnRightString               = pcor ;timstr
  res1@gsnLeftStringOrthogonalPosF  = 0.03
  res1@gsnLeftString                = Expnam(ig) ;timstr

  res1@cnLevelSelectionMode        = "ManualLevels"
  res1@cnLevelSpacingF             =  0.02
  res1@cnMinLevelValF              =  -2
  res1@cnMaxLevelValF              =   2

  contour1(ig) = gsn_csm_contour(wks,U_Z(:,:),res1)
  overlay(plot(ig),  contour1(ig))

  delete([/res,res0,res1/])

  gres                 = True
  gres@gsMarkerIndex   = 16                  ; filled dot
  gres@gsMarkerColor   = "Black"
  gres@gsMarkerSizeF   = 0.008        ; Increase marker sizes.
  indxx =  ind(trktime.le.trkcurt)
  if(.not.all(ismissing(indxx)))then
    dum(ig) = gsn_add_polymarker(wks,plot(ig),psminlon(indxx),psminlat(indxx),gres)
  end if
  delete([/indxx/]) ;trktime,psminlon,psminlat/])

 end do 

 resP                            = True                ; modify the panel plot
 resP@gsnMaximize                = True
 resP@gsnPaperOrientation        = "portrait"
 resP@gsnPanelLabelBar           = True
 resP@pmLabelBarParallelPosF     =  0.001
 resP@pmLabelBarOrthogonalPosF   =  -0.04
 resP@lbLabelFontHeightF         = FontHeightF*0.6
 resP@lbTitleFontHeightF         = FontHeightF*0.65
 resP@lbTitleString              = "Corrective tendency of zonal wind at "+lev_p+" hPa ("+U@units+")" 
 resP@lbOrientation              = "horizontal"  ; vertical label bar
 resP@lbBoxLinesOn               = False
 resP@pmLabelBarHeightF          = 0.06
 resP@pmLabelBarWidthF           = 0.95
 resP@gsnPanelMainString         = "Current time: " + timstr  ;+ " (correlation: " + pcor+")"
 resP@gsnPanelMainFontColor      = "Black"
 resP@gsnPanelMainFontHeightF    = FontHeightF * 0.8
 resP@gsnPanelYWhiteSpacePercent = 0.
 resP@gsnPanelXWhiteSpacePercent = 5.

 gsn_panel(wks,plot,(/1,3/),resP)             ; now draw as one plot

end do

end
