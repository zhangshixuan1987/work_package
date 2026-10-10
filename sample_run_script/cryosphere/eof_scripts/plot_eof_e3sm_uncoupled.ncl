;==============================================================
; This is a revised scripts based on eof_4_640.ncl available at 
; https://www.ncl.ucar.edu/Applications/Scripts/eof_4_640.ncl
;==============================================================
; Match the AAO pattern at: 
;   http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.loading.shtml
;
; Use method: 
;   http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/aao.shtml
;==============================================================
; Data source was E3SM  geopotential height.
;==============================================================
; This script uses functions eofunc_n_Wrap and 
; eofunc_ts_n_Wrap, which were added in NCL V6.4.0
;==============================================================
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

 model = "e3sm"
 ncase = 1
 cases = new(ncase,string)
 icase = 0
 cases(icase)  = "20200910.alpha5_0.F2010.ne30pg2_r05_oECv3.compy"

 do icase  = 0,ncase-1,1

  case  = cases(icase)

  plot_Label = "eof_"+model+"_"+case

; ==============================================================
; User defined parameters that specify region of globe and
; ==============================================================
  latS   = -90.
  latN   = -20. 
  plev   = 700

  yrStrt = 2
  yrLast = 6
  nyear  = 5

  var    = "Z3"
  title  = str_upper(var)+": "+plev+"hPa "

  neof   = 3         ; First 3 EOF modes only, if you change neof = 5 you will get 5 eof modes. 

  optEOF = True       
  optETS = False

; ==============================================================
; Open the file: Read only the user specified period and level
; ==============================================================
  files  = systemfunc ("ls ./data/"+case+"/FV*.nc") ; file paths
  f      = addfiles (files, "r")
    
  P0mb   = 1000.               ; reference pressure [mb]
  intyp  = 1                   ; 1=linear, 2=log, 3=log-log
  kxtrp  = True                ; True=extrapolate
 ;varflg = 0                   ; all other variable
 ;varflg = 1                   ; temperature is variable
  varflg = -1                  ; geopotential height is variable
  
  nlev  = dimsizes(f[0]->hyam)                ; number of vertical levels
  printVarSummary(f[:]->T(:,nlev-1,:,:))


  tbot  = f[:]->T(:,nlev-1,:,:) 
  hgtvint = vinth2p_ecmwf(f[:]->$var$,f[0]->hyam,f[0]->hybm,\
                          plev,f[:]->PS,intyp,P0mb, \
                          1,kxtrp,varflg,tbot,f[:]->PHIS)

  x     = hgtvint(:,{plev},:,:)

; ==============================================================
; compute climatology and Anomalies
; ==============================================================
  xClm   = clmMonTLL(x)                             ; (12,lat,lon)
  printVarSummary(xClm)

  xAnom  = calcMonAnomTLL ( x, xClm)                ; (time, lat,lon)
  printVarSummary(xAnom)         
  printMinMax(xAnom, True)

; =================================================================
; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
; =================================================================
  rad    = get_d2r("float")
  clat   = xAnom&lat            
  clat   = sqrt( cos(rad*clat) )                 ; gw for gaussian grid
  printVarSummary(clat)

; =================================================================
; weight all observations 
; =================================================================
  xw     = xAnom*conform(xAnom, clat, 1)
  copy_VarMeta(x, xw)
  xw@long_name = "Wgt: "+x@long_name

; =================================================================
; Reorder (lat,lon,time) the *weighted* input data
; Compute EOFs & Standardize time series
; =================================================================
  eof    = eofunc_n_Wrap(xw, neof, optEOF, 0)      
  eof    = -1*eof                                ; *special* match sign of CPC

  eof_ts = eofunc_ts_n_Wrap (xw, eof, optETS, 0)

  printVarSummary( eof )                         ; examine EOF variables
  printVarSummary( eof_ts )

  eof_ts = dim_standardize_n( eof_ts, 0, 1)      ; normalize

; =================================================================
; Regress
; =================================================================

  eof_regres = eof                               ; create an array w meta data
  do ne=0,neof-1
     eof_regres(ne,:,:) = (/ regCoef_n(eof_ts(ne,:), xAnom, 0, 0) /)
  end do

; =================================================================
; Extract the YYYYMM from the time coordinate 
; associated with eof_ts [same as x&time] 
; =================================================================

  yyyymm = cd_calendar(eof_ts&time,-1)  
  yrfrac = yyyymm_to_yyyyfrac(yyyymm, 0.0); not used here

;============================================================
; PLOTS
;============================================================
  wks = gsn_open_wks("pdf",plot_Label)         ; send graphics to PNG file

  plot = new(neof,graphic)                ; create graphic array
                                          ; only needed if paneling
; EOF patterns

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnPolar             = "SH"
 
  res@mpFillOn             = False        ; turn off map fill
  res@mpMaxLatF            = latN
  res@mpCenterLonF         = 180

  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "BlueWhiteOrangeRed"  
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's

                                          ; set symmetric plot min/max
  ;symMinMaxPlt(eof_regres, 16, False, res)       ; contributed.ncl
  res@cnLevelSelectionMode = "ManualLevels"     ; manually set the contour levels with the following 3 resources
  res@cnMinLevelValF  = -40.                    ; set the minimum contour level
  res@cnMaxLevelValF  =  40.                    ; set the maximum contour level
  res@cnLevelSpacingF =   5.                    ; *special* match CPC

; panel plot only resources

  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@gsnPaperOrientation = "portrait"   ; force portrait

  resP@gsnPanelMainString  = title+": "+yrStrt+"-"+yrLast

;*******************************************
; first plot
;*******************************************
  do n=0,neof-1
     res@gsnLeftString  = "EOF "+(n+1)
     res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
     plot(n)=gsn_csm_contour_map_polar(wks,eof_regres(n,:,:),res)
  end do
  gsn_panel(wks,plot,(/neof,1/),resP)     ; now draw as one plot

;*******************************************
; second plot
;*******************************************
; EOF time series  [bar form]

  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling               

; these four rtsources allow the user to stretch the plot size, and
; decide exactly where on the page to draw it.

  rts@vpHeightF = 0.40        ; Changes the aspect ratio
  rts@vpWidthF  = 0.85
  rts@vpXF      = 0.10        ; change start locations
  rts@vpYF      = 0.75        ; the plot


  rts@tiYAxisString = "Standardized"          ; y-axis label      

  rts@gsnYRefLine           = 0.              ; reference line   
  rts@gsnXYBarChart         = True            ; create bar chart 
  rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
  rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

; panel plot only resources
  rtsP                      = True            ; modify the panel plot
  rtsP@gsnMaximize          = True            ; large format
  rtsP@gsnPanelMainString   = title+": "+yrStrt+"-"+yrLast

; create individual plots
  do n=0,neof-1
     rts@gsnLeftString  = "EOF "+(n+1)
     rts@gsnRightString = sprintf("%5.1f", eof_regres@pcvar(n)) +"%"
     plot(n) = gsn_csm_xy (wks,yrfrac,eof_ts(n,:),rts)
  end do
  gsn_panel(wks,plot,(/neof,1/),rtsP)     ; now draw as one plot

end do 
