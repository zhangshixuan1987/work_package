;--------------------------------------------------------------------
; DKRZ NCL example: compute_daily_climatology_and_draw_timeseries.ncl
;
; Description:  compute the daily climatology; raw and smooth
;--------------------------------------------------------------------
begin

 mesh   = (/"ne30pg2", "arneppg2", "tcnapg2"/)
 grid   = (/"NE30pg2","ARNEPpg2 RRM","TCNApg2 RRM"/)
 nmesh  = dimsizes(mesh)
 region = (/"GLB","ARNEP", "TCNA", "ETCNA"/)
 lats   = (/ -90,    10,   10,   10/)
 late   = (/  90,    60,   45,   55/)
 lons   = (/-360,  -155, -100, -120/) + 360.0
 lone   = (/   0,   -95,  -50,  -50/) + 360.0
 nreg   = dimsizes(region)

 do im = 0,nmesh-1,1

  mdir   = "/pscratch/sd/z/zhan391/darpa_scratch/SCREAM-IELM-"+mesh(im)+"_20100101_0120/run/"
  varnam = "TWS";"HCSOI" ;"TWS"

  ymdStrt = 20010201 ;-- start yyyymmdd
  ymdLast = 20200101 ;-- last
  yrStrt  = ymdStrt/10000
  yrLast  = ymdLast/10000
 
  fils   = systemfunc("ls -1 " + mdir + "*.elm.h0.*.nc")
  f      = addfiles(fils,"r")
 
  ymd   = cd_calendar(f[:]->time, -2)   ;-- yyyymmdd
  iStrt = ind(ymd .eq. ymdStrt) ;-- index start
  iLast = ind(ymd .eq. ymdLast) ;-- index last
  delete([/ymd/])
 
  time     = f[:]->time(iStrt:iLast)   ;-- time:units = "hours since"
  ymd      = cd_calendar(time, -2)     ;-- yyyymmdd
  TIME     = cd_calendar(time, 0)      ;-- type float
  year     = toint(TIME(:,0))
  month    = toint(TIME(:,1))
  day      = toint(TIME(:,2))
  date_str = sprinti("%0.4i-",year)+sprinti("%0.2i-",month)+sprinti("%0.2i",day)
  timp     = ispan(1,dimsizes(time),1)

  tws  = f[:]->$varnam$(iStrt:iLast,:)  ;total water storage
  area = f[0]->area
  lmsk = f[0]->landfrac
  lat  = f[0]->lat
  lon  = f[0]->lon

 do ig = 0,nreg-1,1

  indx = ind(lat.ge.lats(ig).and.lat.le.late(ig).and.lon.ge.lons(ig).and.lon.le.lone(ig))
  data = dim_avg_wgt_n(tws(:,indx),area(indx),1,1)

 ;printVarSummary(data)
 ;printMinMax(data,False)
  varmin = min(data) ;9000.0 
  varmax = max(data) ;9200.0 
  varunt = tws@units
  lstr   = tws@long_name + " (" +varnam+")"

 ;-- open workstation
  wks_type = "png"
  wks_type@wkWidth  = 1500
  wks_type@wkHeight = 1500
  fignam   = "fig_monthly_" + mesh(im)+"_"+varnam+"_" + region(ig) +"_tseries"
  wks = gsn_open_wks(wks_type,fignam)  ;-- open a workstation

 ;-- set resources
  res = True
  res@vpHeightF = 0.4    ;-- viewport height
  res@vpWidthF  = 0.8     ;-- viewport width
  res@vpXF      = 0.12        ;-- viewport x-position

  res@tmYMajorGrid                = True          ; implement y grid
  res@tmYMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmYMajorGridLineDashPattern = 2             ; select short dash lines
  res@tmXMajorGrid                = True          ; implement y grid
  res@tmXMajorGridThicknessF      = 1.0           ; 2.0 is default
  res@tmXMajorGridLineDashPattern = 2             ; select short dash lines

  res@tmXTOn  = False       ; have tick marks

  res@gsnStringFontHeightF = 0.018		  ; set the gsnLeft/Center/RightString font height
  res@tmXBLabelFontHeightF = 0.014      ;-- decrease x-axis label font size
  res@tmYLLabelFontHeightF = 0.014      ;-- decrease x-axis label font size
  res@tiXAxisFontHeightF   = 0.016        ;-- x-axis label font size
  res@tiYAxisFontHeightF   = 0.016        ;-- y-axis label font size

  res@gsnLeftString      = grid(im)+" : evolution of " + lstr      ; add the gsn titles
  res@gsnCenterString    = "" 
  res@gsnRightString     = yrStrt+"-"+yrLast

  res@trXMinF = min(timp)   ;-- x-axis min value
  res@trXMaxF = max(timp)   ;-- x-axis max value
  res@trYMinF = varmin    ;-- y-axis min value
  res@trYMaxF = varmax    ;-- y-axis max value

  res@tmXBMode = "Manual"               ;-- don't use default settings for tick marks
  res@tmXBTickSpacingF     = 30             ;-- x major tick mark spacing
  res@tmXBMinorPerMajor    = 2             ;-- set number of minor tick marks
  res@tmXBLabelDeltaF      = 0.15            ;-- move x-axis labels down

  res@tiMainString       = "" 
  res@tiXAxisString      = "Time (months)"          ;-- x-axis string
  res@tiYAxisString      = varnam+" ("+varunt+")"          ;-- x-axis string

  res@xyLineThicknesses = (/2.0, 2.0/)  ;-- make 2nd lines thicker
  res@xyLineColors      = (/"blue","red"/)   ;-- change line color
  res@xyMonoDashPattern = True          ;-- all solid
  res@xyMarkLineModes   = "Lines"         ;-- line mode

 ;-- create the plot
  plot = gsn_csm_xy(wks,timp,data,res)
  delete([/indx,data/])
  end do 
  delete([/fils,f,lat,lon,res,time,ymd,TIME,year,month,day,date_str,timp,tws,area,lmsk/])

 end do 

end
