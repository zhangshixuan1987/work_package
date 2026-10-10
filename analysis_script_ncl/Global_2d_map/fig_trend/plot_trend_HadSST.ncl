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

  l_plot_box = False
  flmsk      = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/lnd_msk/landmsk.nc"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim"
  exp_name   = (/"HadSST"/)
  exp_strs   = (/"SST Obs."/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"

  seasons    = (/"ANN", "DJF", "JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11
  stats_slev = 0.05

  varList    = (/"SST"/)
  varUnit    = (/"~S~o~N~C"/)
  nvar       = dimsizes(varList)
  
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  10.
  lonW       =   0.
  lonE       = 360.

 ;sub regions for analysis 
 ;western equatorial Pacific (Niño4 region, 16E-150W, 5S-5N)
 ;western and eastern equatorial Pacific (90W-150W,5-5N) 
 ;Southern Ocean 40-60S, 0-360
  subreg    = (/"Nino34"/) 
  slats     = (/      -5/)   
  slatn     = (/       5/) 
  slonw     = (/     270/)
  slone     = (/     160/)
  regcols   = (/"Black" /)
  reglthk   = (/3       /)
  nsreg     = dimsizes(subreg)

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

   do k = 0,nexps-1

    do j = 0,nvar-1,1

    vars   = varList(j)
    vunt   = varUnit(j)
    print("working on " + vars)

    fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fm0  = addfile(flmsk,"r")
    f1   = addfiles(fl1,"r")
    ListSetType (f1, "join")

    vnm   = vars+"_o"
    vd1   = f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
    lmsk  = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})
    vd1!0 = "ens"
    vd1!1 = "time"
    vd1!2 = "lat"
    vd1!3 = "lon"

    lms1 = conform(vd1,lmsk,(/2,3/))
    vd1  = where(lms1.le.0,vd1,vd1@_FillValue)
    delete([/lmsk,lms1/])
    
    if(vars.eq."PSL") then
     vd1 = vd1/100.0
     vd1@units = "hPa"
     vunt    = "hPa"
    end if 

    print("working on season "+ seasons(i))
    if (seasons(i).eq."ANN") then 
      vdat = vd1(time|:,ens|:,lat|:,lon|:)
      var1 = month_to_annual(vdat,1)
    else 
      vdat = vd1(time|:,ens|:,lat|:,lon|:)
      var1 = month_to_season(vdat,seasons(i))
    end if 

    year = ispan(0,dimsizes(var1(:,0,0,0))-1,1)
    rc   = regCoef_n(year,var1,0,0)
    tval = onedtond(rc@tval, dimsizes(rc))
    df   = onedtond(rc@nptxy, dimsizes(rc)) - 2
    tva1 = dim_avg_n_Wrap(tval,0)
    df1  = dim_avg_n_Wrap(df,0) 
    b    = tva1    ; b must be same size as tval (and df)
    b    = 0.5
    prob = betainc(df1/(df1+tva1^2),df1/2.0,b)       ; prob(nlat,nlon)
    
    vtrend           = var1(0,0,:,:)
    vtprob           = var1(0,0,:,:)
    vtrend@long_name = "Trend"
    vtrend@units     = var1@units+"/year"
    vtprob@long_name = "probability"

    vtrend = (/dim_avg_n(rc,0)*dimsizes(year)/)
    vtprob = where(prob.le.stats_slev,prob,prob@_FillValue)
    delete([/year,rc,tval,df,b,prob,tva1,df1/])

   ;;ensemble mean     
   ;vdm  = dim_avg_n_Wrap(var1,1)
   ;year = ispan(0,dimsizes(vdm(:,0,0))-1,1)
   ;vtrend  = vdm(0,:,:)
   ;vtprob  = vdm(0,:,:)
   ;vtrend@long_name = "Trend"
   ;vtrend@units     = vdm@units+"/year"
   ;vtprob@long_name = "probability"

   ;printVarSummary(vtrend)
   ;printMinMax(vtrend,False)


    rad  = 4.0*atan(1.0)/180.0
    wgty = cos(f1[0]->lat({latS:latN})*rad)
    glm1 = dim_max_n(dim_max_n(vtrend,1),0) 
    glm2 = dim_min_n(dim_min_n(vtrend,1),0) 
    rstr = "Max / Min: "+ sprintf("%.2f",glm1)+" / "+sprintf("%.2f",glm2)
    print(rstr)
    delete([/fl1,f1,glm1,glm2/])
    delete([/vdat,var1,vd1/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (j.eq.0) then 
    ;;;start to process the plot;;;;
    load "./share_fig_info.ncl"
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

    fe              = "Trend_"+exp_name(k)+"_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)
    lndcolor        = "Grey25"
    plot            = new(nvar,graphic)
    dum0            = new(nvar,graphic)
    dum1            = new((/nvar,nsreg/),graphic)
    cmap            = read_colormap_file("BlueDarkRed18")
    color          := array_append_record(cmap(0:7,:),cmap(10:17,:),0)
    delete(cmap)

   end if 

   main_str = "Trend "
   exp_str  = exp_strs(k)
   left_str = exp_strs(k)
   lab_str  =  main_str + " of " + vars + " (" + vunt +" / "+time_tag+")"

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeight = 0.018
   res@gsnStringFontHeightF = FontHeight*2.5
   res@tmXBLabelFontHeightF = FontHeight
   res@tmYLLabelFontHeightF = FontHeight
   res@tiMainFontHeightF    = FontHeight
   res@tiXAxisFontHeightF   = FontHeight
   res@tiYAxisFontHeightF   = FontHeight

   mpres     = res
   resx      = res
   resx1     = res
   
;---Control appearance of map.
   mpres@gsnMajorLatSpacing     = 20
   mpres@gsnMajorLonSpacing     = 45
   mpres@gsnMinorLatSpacing     = 15
   mpres@gsnMinorLonSpacing     = 15
  ;mpres@mpProjection           = "Robinson"       ; choose projection
   mpres@mpLabelsOn             = False
   mpres@mpPerimOn              = False
   mpres@mpFillOn               = False
   mpres@mpGeophysicalLineColor      = lndcolor  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = lthick          ; thickness of outlines
  ;mpres@mpOutlineOn            = True
  ;mpres@mpOutlineDrawOrder     = "PostDraw"
  ;mpres@mpFillDrawOrder        = "Predraw"
  ;mpres@mpGridAndLimbOn        = True
  ;mpres@mpGridLatSpacingF      = 90.0
  ;mpres@mpGridLonSpacingF      = 180.0
  ;mpres@mpGridLineThicknessF   = 3.0
  ;mpres@mpGridAndLimbDrawOrder ="PreDraw"
 ;  mpres@mpLimitMode           = "LatLon"
   mpres@mpCenterLonF           = 180.0         ; default is 0 [GM]
   mpres@mpMinLatF              = latS
   mpres@mpMaxLatF              = latN
   mpres@mpMinLonF              = lonW
   mpres@mpMaxLonF              = lonE
  ;mpres@mpOceanFillColor       = "lightskyblue1"
  ;mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(j) = gsn_csm_map(wks,mpres)

  ;*******************************************
  ; first plot
  ;*******************************************
   res0                       = res
   res0@gsnAddCyclic          = True
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@gsnCenterString       = ""
   res0@tiMainString          = ""
   res0@tiMainOffsetYF        = 0.005

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = lab_str 
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarOrthogonalPosF = 0.25
   res0@pmLabelBarWidthF         = 0.62
   res0@pmLabelBarHeightF        = 0.15
   res0@lbLabelFontHeightF       = FontHeight*0.7
   res0@lbTitleFontHeightF       = FontHeight*0.75

   res0@cnFillOn              = True     ; color plot desired
  ;res0@cnFillMode            = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor = "White"
   res0@cnFillPalette         = color
   res0@cnLevelSelectionMode  = "ExplicitLevels"
   res0@tiMainString          = "" ;exp_str
   res0@gsnLeftString         = left_str
   res0@gsnRightString        = rstr
   res0@cnLevels             := cnlevs1
   res0@cnLinesOn             = False        ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"                   ; color of second contours
   res0@cnLineThicknessF      = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res0@gsnRightStringOrthogonalPosF = 0.03
   res0@gsnLeftStringOrthogonalPosF  = 0.03

   pltrgd  = gsn_csm_contour(wks,vtrend,res0)
   overlay(plot(j) ,pltrgd)
   delete(vtrend)

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
   rescn2@cnFillPattern   = 17 ;6
   rescn2@cnFillDotSizeF  = 0.004
   rescn2@cnFillColor     = "black"
   rescn2@cnFillOpacityF  = 0.6
   dum0(j)   = gsn_csm_contour(wks,vtprob,rescn2) ; contours are at 95,100
   overlay(plot(j),dum0(j))

   if(l_plot_box) then 
    ;---Draw a box using gsn_add_polygon. See how the edges are curved.
    gnres                  = True
    gnres@gsFillColor      = "transparent"
    gnres@gsEdgesOn        = True
    txres                  = True
    txres@txFontHeightF    = 0.015
    txres@txJust           = "CenterLeft"
    do irg = 0,nsreg-1,1
      gnres@gsEdgeColor      = regcols(irg)
      gnres@gsLineColor      = regcols(irg)
      gnres@gsLineThicknessF = reglthk(irg)
      gnres@gsEdgeThicknessF = reglthk(irg)
      latbox = (/slats(irg),slats(irg),slatn(irg),slatn(irg),slats(irg)/)
      lonbox = (/slonw(irg),slone(irg),slone(irg),slonw(irg),slonw(irg)/)
      dum1(j,irg)= gsn_add_polygon(wks, plot(k), lonbox, latbox, gnres)
    end do
   end if 

  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012
   pres1@gsnPanelLabelBar           = False                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = lab_str 
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
   gsn_panel(wks,plot,(/1,nvar/),pres1)          ; create first panel
   delete([/dum1,plot,res,res0,mpres,pres1/]) 
 end do 

end do 

end

