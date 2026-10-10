  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Purpose: A tempolate to make a eof plot for AAO patterns 
;          Data is vertically interpolated to pressure levels
;          1st panel: EOFs for the first three modes
;          2nd panel: time series plot
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  l_ensemble_mean = False 
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/reg_analysis/"
  obs_name   = (/"ERA5", "ERA5"/)
  exp_name   = (/"ERA5", "v2.LR.historical", "v1.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 Coupled", "v1 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"
  year       = ispan(1979,2014,1)*1.0
  stats_slev = 0.05

  seasons    = (/"ANN", "DJF", "JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  varList    = (/"TS",     "TREFHT"/)
  varUnit    = (/"~S~o~N~C","~S~o~N~C"/)
  nvar       = dimsizes(varList)
  
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       = -50.
  lonW       =   0.
  lonE       = 360.

  ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

   do k = 0,nexps-1

   var0  = varList(0)
   fl01  = systemfunc("ls -1 "+data_dir+"/"+obs_name(0)+"_anom_" + var0+ "_MON_en00_"+time_tag+".nc")
   f01   = addfiles(fl01,"r")
   ListSetType (f01, "join")

   fl11  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + var0+ "_MON_en*_"+time_tag+".nc")
   f11   = addfiles(fl11,"r")
   ListSetType (f11, "join")

   var1  = varList(1)
   fl02  = systemfunc("ls -1 "+data_dir+"/"+obs_name(1)+"_anom_" + var1+ "_MON_en00_"+time_tag+".nc")
   f02   = addfiles(fl02,"r")
   ListSetType (f02, "join")

   fl12  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + var1+ "_MON_en*_"+time_tag+".nc")
   f12   = addfiles(fl12,"r")
   ListSetType (f12, "join")
  
   vnm   = var0+"_o"
   sst0  = f01[:]->$vnm$(0,:,{latS:latN},{lonW:lonE})
   sst1  = f11[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
   sst1!0 = "ens"
   sst1!1 = "time"
   sst1!2 = "lat"
   sst1!3 = "lon"

   vnm   = var1+"_o"
   t20   = f02[:]->$vnm$(0,:,{latS:latN},{lonW:lonE})
   t21   = f12[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
   t21!0 = "ens"
   t21!1 = "time"
   t21!2 = "lat"
   t21!3 = "lon"

   sst0@_FillValue  = -99999
   sst1@_FillValue  = -99999
   t20@_FillValue   = -99999
   t21@_FillValue   = -99999

   sst0 = sst0 - 273.15
   sst1 = sst1 - 273.15
   t20  = t20 - 273.15
   t21  = t21 - 273.15

   print("working on season "+ seasons(i))
   if (seasons(i).eq."ANN") then
     var01 = month_to_annual(sst0(time|:,lat|:,lon|:),1)
     var11 = month_to_annual(sst1(time|:,ens|:,lat|:,lon|:),1)
     var02 = month_to_annual(t20(time|:,lat|:,lon|:),1)
     var12 = month_to_annual(t21(time|:,ens|:,lat|:,lon|:),1)
   else
     var01 = month_to_season(sst0(time|:,lat|:,lon|:),seasons(i))
     var11 = month_to_season(sst1(time|:,ens|:,lat|:,lon|:),seasons(i))
     var02 = month_to_season(t20(time|:,lat|:,lon|:),seasons(i))
     var12 = month_to_season(t21(time|:,ens|:,lat|:,lon|:),seasons(i))
   end if

   vclm0 = dim_avg_n_Wrap(var11,0)
   vclm1 = dim_avg_n_Wrap(var12,0)
   sstd  = var11
   t2d   = var12
   sstd  = sstd - conform(var11,vclm0,(/1,2,3/))
   t2d   = t2d  - conform(var12,vclm1,(/1,2,3/)) 

   printMinMax(vclm0,False)
   printMinMax(vclm1,False)
   printMinMax(sstd,False)
   printMinMax(t2d,False)

   delete([/sst0,sst1,t20,t21,var01,var11,var02,var12,vclm0,vclm1/])
   delete([/fl01,fl11,fl02,fl12/])

   sstd!0 = "year"
   sstd!1 = "ens"
   sstd!2 = "lat"
   sstd!3 = "lon"

   t2d!0 = "year"
   t2d!1 = "ens"
   t2d!2 = "lat"
   t2d!3 = "lon"

   nlat  = dimsizes(t2d(0,0,:,0))
   nlon  = dimsizes(t2d(0,0,0,:))

   if(l_ensemble_mean) then 
    ;calculate ensemble mean then calculate trend;;;  
     tmp1    = dim_avg_n_Wrap(sstd,1)
     tmp2    = dim_avg_n_Wrap(t2d,1)
     rcsst   = regCoef(year,tmp1(lat|:,lon|:,year|:))
     rct2    = regCoef(year,tmp2(lat|:,lon|:,year|:))
     rcsst@_FillValue  = sstd@_FillValue
     rct2@_FillValue   = sstd@_FillValue
     df      = rcsst@nptxy-2   ; degrees of freedom
     tval    = rcsst@tval      ; t-statistic
     b       = tval
     b       = 0.5 
     pb      = betainc(df/(df+tval^2),df/2.0,b)
     prsst   = reshape(pb,(/nlat,nlon/))
     prsst@_FillValue = sstd@_FillValue
     prsst   = where(prsst.le.stats_slev,prsst,prsst@_FillValue)

     df      = rct2@nptxy-2   ; degrees of freedom
     tval    = rct2@tval      ; t-statistic
     b       = tval          
     b       = 0.5        
     pb      = betainc(df/(df+tval^2),df/2.0,b)
     prt2    = reshape(pb,(/nlat,nlon/))
     prt2@_FillValue = sstd@_FillValue
     prt2    = where(prt2.le.stats_slev,prt2,prt2@_FillValue)
     delete([/tmp1,tmp2,pb,df,tval,b/])
   else 
    ;do regression first then caluclate ensemble mean trend 
     nens    = dimsizes(sstd(0,:,0,0))
     tmp1    = regCoef(year,sstd(ens|:,lat|:,lon|:,year|:))
     rcsst   = dim_avg_n_Wrap(tmp1,0)
     df      = dim_avg_n(reshape(rcsst@nptxy,(/nens,nlat,nlon/)),0)-2  
     tval    = dim_avg_n(reshape(rcsst@tval,(/nens,nlat,nlon/)),0) 
     b       = todouble(tval)
     b       = 0.5
     prsst   = betainc(df/(df+tval^2),df/2.0,b)
     prsst@_FillValue  = sstd@_FillValue
     prsst   = where(prsst.le.stats_slev,prsst,prsst@_FillValue)
     delete([/tmp1,df,tval,b/])

     nens    = dimsizes(t2d(0,:,0,0))
     tmp2    = regCoef(year,t2d(ens|:,lat|:,lon|:,year|:))
     rct2    = dim_avg_n_Wrap(tmp2,0)
     df      = dim_avg_n(reshape(rct2@nptxy,(/nens,nlat,nlon/)),0)-2
     tval    = dim_avg_n(reshape(rct2@tval,(/nens,nlat,nlon/)),0) 
     b       = todouble(tval)
     b       = 0.5
     prt2    = betainc(df/(df+tval^2),df/2.0,b)
     prt2@_FillValue = sstd@_FillValue
     prt2    = where(prt2.le.stats_slev,prt2,prt2@_FillValue)
     delete([/tmp2,df,tval,b/])

   end if 

   copy_VarCoords(sstd(0,0,:,:), rcsst)                ; copy lat,lon coords
   copy_VarCoords(t2d(0,0,:,:) , rct2)                ; copy lat,lon coords
   copy_VarCoords(sstd(0,0,:,:), prsst)                ; copy lat,lon coords
   copy_VarCoords(t2d(0,0,:,:) , prt2)                ; copy lat,lon coords

   prsst    = (1.0 - prsst)* 100.0
   prt2     = (1.0 - prt2) * 100.0
   prt2     = where(ismissing(rcsst),prt2,prt2@_FillValue)

   delete([/sstd,t2d/])

   nyrs  = dimsizes(year)
   rcsst = rcsst * nyrs
   rct2  = rct2  * nyrs
   rct2            = where(ismissing(rcsst),rct2,rct2@_FillValue)
   rcsst@long_name = "Trend"
   rcsst@units     = varUnit(0) ;+" year~S~-1~N~"
   rct2@long_name  = "Trend"
   rct2@units      = varUnit(1) ;+" year~S~-1~N~"

  ;printMinMax(rcsst,False)
  ;printMinMax(rct2,False)

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
    wtype           = "png"
    if(wtype.eq."png") then
      lthick          = 2.0
      BorderThick     = 2.0
      wtype@wkWidth   = 1200
      wtype@wkHeight  = 1200
    else
      lthick          = 1.0
      BorderThick     = 1.0
    end if

    fe              = varList(0)+"_"+varList(1)+"_trend_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)
    lndcolor        = "Grey25"
    plot            = new((/nexps/),graphic)
    dum1            = new((/nexps/),graphic)
    dum2            = new((/nexps/),graphic)
    dum3            = new((/nexps/),graphic)
    color1          = read_colormap_file("MPL_RdBu")
    color2          = read_colormap_file("MPL_PuOr")
    color1          = color1(::-1,:)
    color2          = color2(::-1,:)
   end if 

   main_str = "Mean bias"
   exp_str  = exp_strs(k)
   left_str = exp_strs(k)

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@gsnPolar       = "SH"

   res@mpFillOn                    = False        ; turn off map fill
   res@mpMaxLatF                   = latN
  ;res@mpCenterLonF                = 180

   res@tmXTOn         = False
   res@tmYROn         = False

   FontHeight        = 0.022
   res@tiMainFontThicknessF          = 3.0
   res@tiMainFontHeightF             = FontHeight*1.2
   res@tmYLLabelFontHeightF          = FontHeight
   res@tmXBLabelFontHeightF          = FontHeight
   res@tiXAxisFontHeightF            = FontHeight
   res@tiYAxisFontHeightF            = FontHeight
   res@gsnStringFontHeightF          = FontHeight
   res@cnLineLabelFontHeightF        = FontHeight
   res@lbLabelFontHeightF            = FontHeight*0.8
   res@lbTitleFontHeightF            = FontHeight*0.8
   res@gsnPolarLabelFontHeightF      = FontHeight

   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

  ;---Control appearance of map.
  ;res@mpCenterLonF                = 180
  ;res@mpGeophysicalLineColor      = lndcolor  ; color of cont. outlines
  ;res@mpGeophysicalLineThicknessF = lthick          ; thickness of outlines
  ;res@mpOutlineOn                 = True
  ;res@mpOutlineDrawOrder          = "PostDraw"
  ;res@mpFillDrawOrder             = "Predraw"
  ;res@mpGridAndLimbOn             = True
  ;res@mpGridLatSpacingF           = 90.0
  ;res@mpGridLonSpacingF           = 180.0
  ;res@mpGridLineThicknessF        = 3.0
  ;res@mpGridAndLimbDrawOrder      ="PreDraw"
  ;res@mpOceanFillColor            = "lightskyblue1"
  ;res@mpLandFillColor             = "gray33"

  ;*******************************************
  ; first plot
  ;*******************************************
   res0                       = res
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@gsnCenterString       = ""
   res0@tiMainString          = ""
   res0@tiMainOffsetYF        = 0.005

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarWidthF         = 0.62
   res0@pmLabelBarHeightF        = 0.15
   res0@lbLabelFontHeightF       = FontHeight*0.7
   res0@lbTitleFontHeightF       = FontHeight*0.75
   res0@cnFillOn                 = True     ; color plot desired
  ;res0@cnFillMode               = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor    = "Transparent" ;"White"
   res0@cnLevelSelectionMode     = "ExplicitLevels"
   res0@cnLinesOn                = False        ; True is default
   res0@cnInfoLabelOn            = False
   res0@cnLineLabelsOn           = False    ; turn off contour labels
   res0@cnLineColor              = "Black"                   ; color of second contours
   res0@cnLineThicknessF         = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res0@gsnRightStringOrthogonalPosF = 0.03
   res0@gsnLeftStringOrthogonalPosF  = 0.03
   res0@cnFillPalette         = color1
   res0@tiMainString          = exp_str
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@cnLevels             := fspan(-2.0,2.0,21) 
   res0@lbTitleString         = "SAT change over 1979-2014" 
   res0@pmLabelBarOrthogonalPosF = 0.10
   res0@lbBoxLinesOn             = False
   plot(k)  = gsn_csm_contour_map_polar(wks,rct2,res0)

   res1 = True
   res1@lbLabelBarOn             = True
  ;res1@lbOrientation            = "Vertical"
  ;res1@lbTitlePosition          = "Right"     ; title location
  ;res1@lbTitleDirection         = "Across"    ; letter angle
  ;res1@lbTitleAngleF            = 90.         ; title angle
   res1@lbBoxLineColor           = "black"
  ;res1@lbBoxLinesOn             = False
  ;res1@pmLabelBarParallelPosF   = 0.61
   res1@pmLabelBarWidthF         = 0.62
   res1@pmLabelBarHeightF        = 0.15
   res1@lbLabelFontHeightF       = FontHeight*0.7
   res1@lbTitleFontHeightF       = FontHeight*0.75
   res1@cnLevelSelectionMode     = "ExplicitLevels"
   res1@cnLinesOn                = False        ; True is default
   res1@cnInfoLabelOn            = False
   res1@cnLineLabelsOn           = False    ; turn off contour labels
   res1@cnFillOn                 = True     ; color plot desired
   res1@gsnTickMarksOn           = False     ; no tickmarks
   res1@gsnDraw                  = False     ; don't draw
   res1@gsnFrame                 = False     ; don't advance frame
   res1@tiMainString             = ""
   res1@gsnLeftString            = ""
   res1@gsnRightString           = ""
   res1@cnFillPalette            = color2
   res1@cnLevels                := fspan(-1.2,1.2,25) 
   res1@lbTitleString            = "TS change over 1979-2014" 
   res1@lbBoxLinesOn             = False
   res1@pmLabelBarOrthogonalPosF = 0.02
   dum1(k)   = gsn_csm_contour(wks,rcsst,res1)
   overlay(plot(k),dum1(k))

   rescn2 = True
   rescn2@lbLabelBarOn    = False
   rescn2@gsnTickMarksOn  = False     ; no tickmarks
   rescn2@gsnDraw         = False     ; don't draw
   rescn2@gsnFrame        = False     ; don't advance frame
   rescn2@tiMainString    = ""
   rescn2@gsnLeftString   = ""
   rescn2@gsnRightString  = ""
   rescn2@cnMinLevelValF  = (1.0 - stats_slev)*100.0
   rescn2@cnMaxLevelValF  = 100.0
  ;rescn2@cnLevelSpacingF = 1.
   rescn2@cnMonoFillPattern = True
   rescn2@cnMonoFillColor = True
   rescn2@cnFillOn        = True      ; color fill
   rescn2@cnLinesOn       = False
   rescn2@cnLineLabelsOn  = False        ; True is default
   rescn2@cnInfoLabelOn   = False                ; turn off info label
   rescn2@cnFillMode      = "AreaFill" ; raster fill significance patterns
   rescn2@cnFillPattern   = 6
   rescn2@cnFillDotSizeF  = 0.002
   rescn2@cnFillColor     = "black"
   rescn2@cnFillOpacityF  = 0.5
   dum2(k)   = gsn_csm_contour(wks,prsst,rescn2) ; contours are at 95,100
   dum3(k)   = gsn_csm_contour(wks,prt2,rescn2) ; contours are at 95,100
   overlay(plot(k),dum2(k))
   overlay(plot(k),dum3(k))

   delete([/prsst,prt2,rcsst,rct2/])

  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012
   pres1@gsnPanelLabelBar           = False                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = "" ;main_str + " of " + vars + " (" + vunt +")"
   pres1@lbTitleFontHeightF         = FontHeight*0.6
   pres1@lbLabelFontHeightF         = FontHeight*0.55 ;  font height
   pres1@lbBoxLineColor             = "black"
   pres1@lbBoxLinesOn               = True
  ;pres1@pmLabelBarParallelPosF     = 0.61
   pres1@pmLabelBarOrthogonalPosF   = 0.01
   pres1@pmLabelBarWidthF           = 0.95
   pres1@pmLabelBarHeightF          = 0.08
  ;pres1@gsnPanelBottom             = 0.45       ; space for label bar
  ;pres1@gsnPanelTop                = 0.85       ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 2
   pres1@gsnPanelXWhiteSpacePercent = 2
   gsn_panel(wks,plot,(/1,nexps/),pres1)          ; create first panel
  
 end do 

end

