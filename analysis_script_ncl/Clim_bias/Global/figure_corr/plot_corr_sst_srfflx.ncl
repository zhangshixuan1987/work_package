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

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim/"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP"   , "v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"

  statistical_siglev = 0.05

  seasons    = (/"MON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions = (/"Global"/)
  nregs   = dimsizes(Regions)
  latS    = -90.
  latN    =  90.
  lonW    =  0.
  lonE    = 360.

  load "./share_var_info.ncl"
  nvar = dimsizes(varList)

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  do j = 0,nvar-1,1

  vars   = varList(j)
  varunt = untList(j)
  print("working on " + vars)

  do k = 0,nexps-1

  if (exp_name(k).eq."ERA5".or.exp_name(k).eq."v2.LR.amip") then 
    fl0 = systemfunc("ls -1 "+data_dir+"/HadSST_anom_SST_" +seasons(i)+"_en*_"+time_tag+".nc")
    f0  = addfiles(fl0,"r")
    ListSetType (f0, "join")
  else
    fl0 = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_SST_" +seasons(i)+"_en*_"+time_tag+".nc")
    f0  = addfiles(fl0,"r")
    ListSetType (f0, "join")
  end if 

  if(exp_name(k).eq."v2.LR.historical") then 
    fl1 = systemfunc("ls -1 "+data_dir+"/v2.LR.amip_anom_" + vars+ "_" +seasons(i)+"_en*_"+time_tag+".nc")
    f1  = addfiles(fl1,"r")
    ListSetType (f1, "join")
  else 
    fl1 = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_" +seasons(i)+"_en*_"+time_tag+".nc")
    f1  = addfiles(fl1,"r")
    ListSetType (f1, "join")
  end if 

  right_str = "Corr: "+vars +" v.s. SST" 

  print("sst file " + fl0)
  print("surface flux file " + fl1)

  vnm0  = "SST_a"
  vd0   = f0[:]->$vnm0$(:,:,{latS:latN},{lonW:lonE})
  vnm1  = vars+"_a"
  vd1   = f1[:]->$vnm1$(:,:,{latS:latN},{lonW:lonE})

  nen0  = dimsizes(vd0(:,0,0,0))
  nen1  = dimsizes(vd1(:,0,0,0))
  nens  = nen0*nen1
  nlat  = dimsizes(vd0(0,0,:,0))
  nlon  = dimsizes(vd0(0,0,0,:))

  corr  = new((/nens,nlat,nlon/),typeof(vd1))
  prob  = new((/nens,nlat,nlon/),typeof(vd1))
  
  do ii = 0,nen0-1,1
   do jj = 0,nen1-1,1
     ll = ii*nen1 + jj
     xx     = vd0(ii,:,:,:)
     yy     = vd1(jj,:,:,:)
     rc1    = escorc_n(xx, yy, 0, 0) 

     ;;;;;;;t-test;;;;;;;;
     df1    = dimsizes(xx(:,0,0))-2.0
     xb1    = rc1*sqrt(df1/(1-rc1^2)) 
     pb1    = student_t(xb1, df1)
     prob1  = where(pb1.le.statistical_siglev,rc1,rc1@_FillValue)
     corr(ll,:,:) = (/rc1/)
     prob(ll,:,:) = (/tofloat(prob1) /)
     delete([/xx,yy,df1,xb1,pb1,prob1,rc1/])
   end do
  end do 

  vard  = vd1(0,0,:,:)
  vpob  = vd1(0,0,:,:)
  vard  = (/dim_avg_n_Wrap(corr,0)/)
  vpob  = (/dim_avg_n_Wrap(prob,0)/)

  delete([/fl0,f0,fl1,f1,vd0,vd1,corr,prob/])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Start to process plot 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (k.eq.0) then 
    ;;;start to process the plot;;;;
    wtype           = "png"
    if(wtype.eq."png") then
      lthick          = 3.0
      BorderThick     = 3.0
      wtype@wkWidth   = 2500
      wtype@wkHeight  = 2500
    else
      lthick          = 1.0
      BorderThick     = 1.0
    end if

    fe              = vars+"_bias_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)

    lndcolor        = "Grey25"
    plot            = new(nexps,graphic)
    dum             = new(nexps,graphic)
    color           = read_colormap_file("BlueDarkRed18")
    ncols           = dimsizes(color(:,0))
    color(ncols/2:ncols/2,:) = 0
  end if 

  if(k.le.1) then 
    cnlevs  = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5/) ;fspan(-1,1,11)
  else 
    cnlevs  = (/-0.1,-0.07,-0.05,-0.03,-0.01,0.01,0.03,0.05,0.07,0.10/) ;fspan(-1,1,11)
  end if 
  
  title_str = "Corr. coef." ;main_str + " of " + vars + " (" + varunt +")"
  exp_str   = exp_strs(k)
  left_str  = ""

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 1.0
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

   mpres     = res
   resx      = res
   resx1     = res

   if(wtype.eq."png") then
    lthick = 3.0
   else
    lthick = 1.0
   end if

;---Control appearance of map.
   mpres@gsnMajorLatSpacing     = 30
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
   mpres@mpCenterLonF           = 180.         ; default is 0 [GM]
   mpres@mpMinLatF              = latS
   mpres@mpMaxLatF              = latN
   mpres@mpMinLonF              = lonW
   mpres@mpMaxLonF              = lonE
  ;mpres@mpOceanFillColor       = "lightskyblue1"
  ;mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(k) = gsn_csm_map(wks,mpres)

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
   res0@gsnStringFontHeightF  = FontHeightF*1.1
   res0@lbLabelFontHeightF    = FontHeightF*0.6
   res0@lbTitleFontHeightF    = FontHeightF*0.8

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = title_str 
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarOrthogonalPosF = 0.2
  ;res0@pmLabelBarWidthF         = 0.95
  ;res0@pmLabelBarHeightF        = 0.08

   res0@cnFillOn              = True     ; color plot desired
  ;res0@cnFillMode            = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor = "White"
   res0@cnFillPalette         = color
   res0@cnLevelSelectionMode  = "ExplicitLevels"
   res0@tiMainString          = exp_str
   res0@gsnLeftString         = left_str
   res0@gsnRightString        = right_str
   res0@cnLevels             := cnlevs
   res0@cnLinesOn             = False        ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"                   ; color of second contours
   res0@cnLineThicknessF      = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1

   pltrgd  = gsn_csm_contour(wks,vard,res0)
   overlay(plot(k) ,pltrgd)

   rescn2 = True
   rescn2@lbLabelBarOn    = False
   rescn2@gsnTickMarksOn  = False     ; no tickmarks
   rescn2@gsnDraw         = False     ; don't draw
   rescn2@gsnFrame        = False     ; don't advance frame
   rescn2@tiMainString    = ""
   rescn2@gsnLeftString   = ""
   rescn2@gsnRightString  = ""
   rescn2@cnMinLevelValF  = (1.0 - statistical_siglev)*100.0
   rescn2@cnMaxLevelValF  = 100.0
  ;rescn2@cnLevelSpacingF = 1.
   rescn2@cnMonoFillPattern = True
   rescn2@cnMonoFillColor = True
   rescn2@cnFillOn        = True      ; color fill
   rescn2@cnLinesOn       = False
   rescn2@cnLineLabelsOn  = False        ; True is default
   rescn2@cnInfoLabelOn   = False                ; turn off info label
   rescn2@cnFillMode      = "AreaFill" ; raster fill significance patterns
   rescn2@cnFillPattern   = 17
   rescn2@cnFillDotSizeF  = 0.0030
   rescn2@cnFillColor     = "black"
   rescn2@cnFillOpacityF  = 0.5
   dum(k) = gsn_csm_contour(wks,vpob,rescn2) ; contours are at 95,100
   overlay(plot(k),dum(k))
   delete([/vpob,vard/])

  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012
   pres1@gsnPanelLabelBar           = False                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = title_str 
   pres1@lbTitleFontHeightF         = FontHeightF*0.6
   pres1@lbLabelFontHeightF         = FontHeightF*0.55 ;  font height
   pres1@lbBoxLineColor             = "black"
   pres1@lbBoxLinesOn               = True
  ;pres1@pmLabelBarParallelPosF     = 0.61
   pres1@pmLabelBarOrthogonalPosF   = 0.01
   pres1@pmLabelBarWidthF           = 0.95
   pres1@pmLabelBarHeightF          = 0.08
  ;pres1@gsnPanelBottom             = 0.45       ; space for label bar
  ;pres1@gsnPanelTop                = 0.85       ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 8
   pres1@gsnPanelXWhiteSpacePercent = 4
   gsn_panel(wks,plot,(/1,nexps/),pres1)          ; create first panel
  
 end do 

end do 

end

