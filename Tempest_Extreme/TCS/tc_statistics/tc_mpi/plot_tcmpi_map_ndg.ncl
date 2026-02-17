;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;This script is used to plot the wrf simulations vs NEXRAD data ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  LOAD FUNCTIONS AND PROCEDURES
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups  = (/"ERA5",     "CLIM",        "NDGUV",         "NDGUVTQ"/)
  LabList = (/"ERA5",     "E3SM (CLIM)", "E3SM (NDG UV)", "E3SM (NDG UV+TQ)"/)
  timtag  = (/"1979-2017","2007-2017",   "2007-2017",     "2007-2017"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  figkey  = "trkdens_10year"

  Figstr  = (/"(a)", "(b)", "(c)", "(d)", "(e)", "(f)","(g)","(h)"/)
  ystr    = 2008
  yend    = 2017
  nyr     = yend - ystr +1

  Varnam  = "trkdens"
  Varstr  = "Track density"
  Varunt  = "storm trainsts per month"
  vbtrk   = new ((/nbsin/),double)
  sigr     = 0.05

 ;********************************
 ; select sub-regions
 ;********************************
  minLat         = -90.             ; max lat to read in (deg)
  maxLat         =  90.             ; min lat to read in (deg)
  minLon         =   0.             ; max lon to read in (deg)
  maxLon         = 360.             ; min lon to read in (deg)
  cenlon         = 200.0 ;
 
  do i = 0,ngrps-1,1

   Varstr = "vmax"
   Varunt = "m/s"
   expnam = Casdir(i)
   if(expnam.eq."ERA5") then 
     fls    = systemfunc("ls " + datdir + expnam +"/"+"*ens00*.nc") 
   else
     fls    = systemfunc("ls " + datdir + expnam +"/"+"*h0*.nc")
   end if 

   f0     = addfile("./lsm.nc","r")
   lsm    = f0->LANDFRAC

   f      = addfiles (fls, "r")
  ;ListSetType (f1, "join")        ; concatenate (default)
  
   if(expnam.eq."ERA5") then 
     time   = f[:]->time
   else 
     time   = f[:]->time
     time   = time - 1
   end if 
   tutc   = cd_calendar(time, 0)
   month  = tointeger(tutc(:,1)) 
   yymm   = cd_calendar(time, -1)
   ym0    = ystr+ "05"
   ym1    = yend+ "11"
  ;wet season (May to November)
   indx   = ind((month.ge.5.and.month.le.11).and.(yymm.ge.ym0.and.yymm.le.ym1))
   lat    = f[0]->lat
   lon    = f[0]->lon
   ifl    = f[:]->ifl(indx,:,:)
   var    = f[:]->vmax(indx,:,:)
   var@units = Varunt
   var = where(ifl.eq.0,-9999,var)
   var@_FillValue = -9999

  if(i.eq.0) then  
    vobs = var 
  end if 

  if(i.eq.0) then 
    vplt  = dim_avg_n_Wrap(vobs,0)
  else if (i.eq.1) then
    rad  = 4.0*atan(1.0)/180.0
    wgty = cos(var&lat*rad)
    xtmp = var(lat|:,lon|:,time|:)   ; reorder but do it only once [temporary]
    ytmp = vobs(lat|:,lon|:,time|:)
    xAve = dim_avg (xtmp)             ; calculate means at each grid point 
    yAve = dim_avg (ytmp)
    xVar = dim_variance (xtmp)       ; calculate variances
    yVar = dim_variance (ytmp)
    xEqv = equiv_sample_size (xtmp, sigr,0)
    yEqv = equiv_sample_size (ytmp, sigr,0)
    xN   = wgt_areaave (xEqv, wgty, 1., 0)    ; wgty could be gaussian weights 
    yN   = wgt_areaave (yEqv, wgty, 1., 0) 
    iflag= False
    prob = var(0,:,:) 
    prob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False) 
    vplt = var(0,:,:) 
    vplt = dim_avg_n(var - vobs,0)
    delete([/rad,wgty,xAve,xVar,xN, yAve,yVar,yN,xtmp,ytmp/])
  else 
    rad  = 4.0*atan(1.0)/180.0
    wgty = cos(var&lat*rad)
    xtmp = var(lat|:,lon|:,time|:)   ; reorder but do it only once [temporary]
    ytmp = vobs(lat|:,lon|:,time|:)
    xAve = dim_avg (xtmp)             ; calculate means at each grid point
    yAve = dim_avg (ytmp)
    xVar = dim_variance (xtmp)       ; calculate variances
    yVar = dim_variance (ytmp)
    xEqv = equiv_sample_size (xtmp, sigr,0)
    yEqv = equiv_sample_size (ytmp, sigr,0)
    xN   = wgt_areaave (xEqv, wgty, 1., 0)    ; wgty could be gaussian weights
    yN   = wgt_areaave (yEqv, wgty, 1., 0)
    iflag= False
    prob = var(0,:,:)
    prob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)
    vplt = var(0,:,:) 
    vplt = dim_avg_n(var- vobs,0)
    delete([/rad,wgty,xAve,xVar,xN, yAve,yVar,yN,xtmp,ytmp/])
  end if 
  end if 

  delete([/tutc,time,yymm,ym0,ym1,indx,var,month,fls,f,ifl/])

  vplt = mask(vplt, (lsm.lt.1), True)
  if(i.gt.0) then
    prob = where(prob.le.sigr,prob,prob@_FillValue)
    prob = 1.0 - prob
    prob = mask(prob, (lsm.lt.1), True)
  end if 

  if(i.eq.0) then 
    ;;;generate the zoom in information;;;;
    fig_type  = "png" ;;figure types: pdf or png
    fig_dir   = "./"
    fignam    = "fig_tcmpi_tropical_bias_nudging"
    fig_type@wkWidth  = 1200
    fig_type@wkHeight = 1200 
    gridthick         = 2.0
    lthick            = 2.0
    lndcolor   = "Black" ;"Grey50"
    wks        = gsn_open_wks(fig_type, fig_dir + fignam)
    plot      = new (ngrps, graphic)
    countourx = new (ngrps, graphic)
    countoury = new (ngrps, graphic)
    cmap0 = read_colormap_file("precip2_17lev") ;"WhViBlGrYeOrRe")
    cmap1 = read_colormap_file("hotcold_18lev") ;BlueDarkRed18")
    cmap2 = read_colormap_file("CBR_coldhot")
    cnlev0 = (/0.1,1,2,5,10,20,30,40,50,60,70,80/)
    cnlev1 = (/-6,-4,-2,-1,-0.1,0.1,1,2,4,6/)
    cnlev2 = (/-6,-4,-2,-1,-0.1,0.1,1,2,4,6/)
  end if 

  dlon = 30
  dlat = 30

  dlat = 45
  dlon = 90

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

   FontHeightF = 0.020
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF

   mpres     = res
   resx      = res
   resx1     = res

   res@gsnStringFontHeightF  = FontHeightF*2.0
   res@lbLabelFontHeightF    = FontHeightF*1.0
   res@lbTitleFontHeightF    = FontHeightF*1.0

;---Control appearance of map.
   mpres@gsnMajorLatSpacing     = 30
   mpres@gsnMajorLonSpacing     = 45
   mpres@gsnMinorLatSpacing     = 15
   mpres@gsnMinorLonSpacing     = 15
  ;mpres@mpProjection           = "Robinson"       ; choose projection
   mpres@mpLabelsOn             = False
   mpres@mpPerimOn              = False
   mpres@mpFillOn               = False
   mpres@mpGeophysicalLineColor = lndcolor  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = 0.1          ; thickness of outlines
  ;mpres@mpOutlineOn            = True
  ;mpres@mpOutlineDrawOrder     = "PostDraw"
  ;mpres@mpFillDrawOrder        = "Predraw"
  ;mpres@mpGridAndLimbOn        = True
  ;mpres@mpGridLatSpacingF      = 90.0
  ;mpres@mpGridLonSpacingF      = 180.0
  ;mpres@mpGridLineThicknessF   = 1.0
  ;mpres@mpGridAndLimbDrawOrder="PreDraw"
  ;mpres@mpMinLonF              = 0  ;-180.
  ;mpres@mpMaxLonF              = 360;180.
   mpres@mpMinLatF              = -90
   mpres@mpMaxLatF              =  90
   mpres@mpCenterLonF           = 180.0
  mpres@mpLandFillColor  = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

  plot(i)                      = gsn_csm_map(wks,mpres)

  res1                                = res 
  res1@cnFillOn                       = True        ; Turn on contour fill
  res1@cnLevelSelectionMode           = "ExplicitLevels"
 ;res1@cnLevelSelectionMode           = "ManualLevels"
 ;res1@cnMinLevelValF                 =  960 ;989.0
 ;res1@cnMaxLevelValF                 =  1040.0
 ;res1@cnLevelSpacingF                =  3.0
  res1@gsnLeftString                  = ""
  res1@gsnRightString                 = ""

  res1@lbLabelBarOn             = True  ; Labelbar already created in 1st plot  
 ;res1@pmLabelBarParallelPosF   =  0.01
  res1@pmLabelBarOrthogonalPosF = 0.10
  res1@lbLabelFontHeightF       = FontHeightF*0.6
  res1@lbTitleFontHeightF       = FontHeightF*0.6
  res1@lbTitleString            = ""
  res1@lbTitlePosition          = "Top";"Bottom"
  res1@lbOrientation            = "horizontal"  ; vertical label bar
  res1@pmLabelBarHeightF        = 0.10
  res1@pmLabelBarWidthF         = 0.55

  res1@cnLinesOn                      = False  ; turn off contour lines
  res1@cnLineLabelsOn                 = False   ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@cnLineLabelDensityF            = 2.0
  res1@cnLineLabelFontHeightF         = FontHeightF*0.8
  res1@cnLineLabelBackgroundColor     = -1
  res1@gsnContourLineThicknessesScale = 2.0
  res1@cnLineThicknessF               = 1.0
  res1@cnLineColor                    = "Black"
  res1@gsnContourZeroLineThicknessF   = 2.
  res1@gsnContourNegLineDashPattern   = 1

  res1@gsnLeftString                  = ""   ;Expnam(i)
  res1@gsnLeftStringOrthogonalPosF    = 0.04
  res1@gsnRightString                 = ""   ;Varstr+" (" +Varunt + ")"
  res1@gsnRightStringOrthogonalPosF   = 0.04
  res1@gsnCenterString                = Expnam(i)
  res1@gsnCenterStringOrthogonalPosF  = 0.04

  if(i.eq.0) then 
    res1@cnFillPalette              := cmap0
    res1@cnLevels                   := cnlev0
    res1@lbTitleString               = Varstr+" (" +Varunt + ")" 
    countourx(i) = gsn_csm_contour(wks,vplt,res1)
  else if (i.eq.1) then 
    res1@cnFillPalette              := cmap1
    res1@cnLevels                   := cnlev1
    res1@lbTitleString               = "~F33~D~N~~F~"+Varstr+" (" +Varunt + ")" 
    countourx(i) = gsn_csm_contour(wks,vplt,res1)
  else
    res1@cnFillPalette              := cmap2
    res1@cnLevels                   := cnlev2
    res1@lbTitleString               = "~F33~D~N~~F~"+Varstr+" (" +Varunt + ")" 
    countourx(i) = gsn_csm_contour(wks,vplt,res1)
  end if 
  end if 

  overlay(plot(i),countourx(i))

  if (i.gt.0) then
    rescn2 = True
    rescn2@lbLabelBarOn    = False
    rescn2@gsnTickMarksOn  = False     ; no tickmarks
    rescn2@gsnDraw         = False     ; don't draw
    rescn2@gsnFrame        = False     ; don't advance frame
    rescn2@tiMainString    = ""
    rescn2@gsnLeftString   = ""
    rescn2@gsnRightString  = ""
   ;rescn2@cnMinLevelValF  = (1.0 - sigr)
   ;rescn2@cnMaxLevelValF  = 1.0
   ;rescn2@cnLevelSpacingF = 1.
    rescn2@cnMonoFillPattern = True
    rescn2@cnMonoFillColor = True
    rescn2@cnFillOn        = True      ; color fill
    rescn2@cnLinesOn       = False
    rescn2@cnLineLabelsOn  = False        ; True is default
    rescn2@cnInfoLabelOn   = False                ; turn off info label
    rescn2@cnFillMode      = "AreaFill" ; raster fill significance patterns
    rescn2@cnFillScaleF    = 0.4
    rescn2@cnFillPattern   = 17
    rescn2@cnFillDotSizeF  = 0.002
    rescn2@cnFillColor     = "black"
    countoury(i) = gsn_csm_contour(wks,prob,rescn2) ; contours are at 95,100
    overlay(plot(i),countoury(i))
  end if 

  delete([/mpres,res,res1/])
 end do 

  resP                     = True
  resP@gsnMaximize              = True
  resP@gsnPaperOrientation      = "portrait"
  resP@gsnPanelLabelBar         = False
  resP@gsnPanelMainString       = "" 
  resP@gsnPanelMainFontColor    = "Black"
  resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6
  resP@gsnPanelMainPosXF        = 0.85
 ;resP@txPosYF                  = 0.7
  resP@gsnPanelBottom           = 0.3
  resP@gsnPanelYWhiteSpacePercent = 0.
  resP@gsnPanelXWhiteSpacePercent = 3.
  gsn_panel(wks,plot,(/1,3/),resP)

end
