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

  fig_grp    = "wtmbias"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/enso_eof/" + fig_grp
  exp_name   = (/"HadSST","ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"HadSST","ERA5", "v2 AMIP","v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name

  latS       = -20.  ;negative for southern hemisphere
  latN       =  20.  ;negative for southern hemisphere
  lonW       = 120.0
  lonE       = 300.0

 ;define regions for eof analysis
  mode       = ("ENSO_EP") 
  ieof       = 1 
  vars       = (/"SST"/)
  facmod     = (/1.0 /)
  varunt     = (/"~S~o~N~C / std dev"/)

  time_tag   = (/"1979-2014"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  do k = 0,nexps-1

  fl       = data_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars+"_"+seasons(i)+"_"+time_tag+".nc"
  f        = addfile(fl,"r")
  eof      = f->eof_reg(:,ieof-1,:,:)
  eof_pcv  = f->eof_pcv(:,ieof-1)
  yyyy     = f->time 

  eofm     = dim_avg_n_Wrap(eof,0)

  rad      = 4.0*atan(1.0)/180.0
  wgty     = cos(f->lat*rad)

  if (seasons(i).eq."MON") then 
    main_str = "Monthly mean " + mode + " EOF "+ ieof 
  else
    main_str = seasons(i)+" mean " + mode + " EOF "+ ieof
  end if 
  exp_str  = exp_strs(k)
  left_str = sprintf("%.1f", avg(eof_pcv)) + " ~F34~1~N~~F~ " + \
             sprintf("%.2f", stddev(eof_pcv))

  if(exp_name(k).eq."HadSST") then
    eof0      = eof 
    nens0     = dimsizes(eof(:,0,0))
    mstrx1    = "" ;"corr ens " + glCordd2
    mstrx2    = "" ;"diff ens " + glAvedd2
    mstrx3    = "" ;"rmsd ens " + glRmsdd2
    right_str = mstrx1 + "~C~" + mstrx3  ;mstrx2 + "~C~" + mstrx3
    delete([/mstrx1,mstrx2,mstrx3/])
  else
    nens      = dimsizes(eof(:,0,0))
    tmp1      = new(nens*nens0,double)
    tmp2      = new(nens*nens0,double)
    tmp3      = new(nens*nens0,double)
    do ii = 0,nens0-1,1
     do jj = 0,nens-1,1
       tmp1(jj+ii*nens) = pattern_cor(eof(jj,:,:),eof0(ii,:,:),1.0,0)
       tmp2(jj+ii*nens) = wgt_areaave(eof(jj,:,:)-eof0(ii,:,:),wgty,1.0,0)
       tmp3(jj+ii*nens) = wgt_arearmse(eof(jj,:,:),eof0(ii,:,:),wgty,1.0,0)
     end do 
    end do 
    glCordd2  = sprintf("%3.1f", avg(tmp1)) + " ~F34~1~N~~F~ " + \
                sprintf("%3.2f", stddev(tmp1))
    glAvedd2  = sprintf("%3.1f", avg(tmp2)) + " ~F34~1~N~~F~ " + \
                sprintf("%3.2f", stddev(tmp2))
    glRmsdd2  = sprintf("%3.1f", avg(tmp3)) + " ~F34~1~N~~F~ " + \
                sprintf("%3.2f", stddev(tmp3))
    mstrx1    = "corr   w/r obs: " + glCordd2
    mstrx2    = "diff w/r obs: " + glAvedd2
    mstrx3    = "rmsd w/r obs: " + glRmsdd2
    right_str = mstrx1 + "~C~" + mstrx3  ;mstrx2 + "~C~" + mstrx3
    delete([/tmp1,tmp2,tmp3,glCordd2,glAvedd2,glRmsdd2,mstrx1,mstrx2,mstrx3/])
  end if 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Start to process plot 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (k.eq.0) then 
    wtype           = "png"
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    fignam          = mode+"_EOF_"+vars+"_"+seasons(i)+"_"+time_tag +"_"+fig_grp
    wks             = gsn_open_wks(wtype,fignam)
    plot            = new((/nexps/),graphic)
    maxlev   = 11
    cnlevs1  = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5/)
    cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  end if 

  ;============================================================
  ; PLOTS
  ;============================================================
   res                       = True
   res@gsnDraw               = False        ; don't draw yet
   res@gsnFrame              = False        ; don't advance frame yet
   res@gsnAddCyclic          = False        ; data not cyclic

   FontHeightF = 0.024
   res@tiMainFontHeightF     = FontHeightF*1.5
   res@tmYLLabelFontHeightF  = FontHeightF
   res@tmXBLabelFontHeightF  = FontHeightF
   res@lbLabelFontHeightF    = FontHeightF
   res@gsnStringFontHeightF  = FontHeightF*1.2

   res@cnFillOn              = True         ; turn on color fill
   res@cnLinesOn             = False        ; True is default
   res@cnLineLabelsOn        = False        ; True is default

   res@lbLabelBarOn             = False
   res@lbOrientation            = "Vertical"
   res@lbTitleString            = "" 
   res@lbTitlePosition          = "Right"     ; title location
   res@lbTitleDirection         = "Across"    ; letter angle
   res@lbTitleAngleF            = 90.         ; title angle
   res@lbTitleFontHeightF       = FontHeightF ; font height
   res@lbBoxLineColor           = "black"  
  ;res@lbBoxLinesOn             = False
  ;res@pmLabelBarParallelPosF   = 0.61
   res@pmLabelBarOrthogonalPosF = 0.15

   res@tmXTOn            = False
   res@mpCenterLonF      = 180.         ; default is 0 [GM]
   res@mpMinLatF         = latS
   res@mpMaxLatF         = latN
   res@mpMinLonF         = lonW
   res@mpMaxLonF         = lonE

   res@gsnAddCyclic      = False
   res@mpShapeMode       = "FreeAspect"
   res@vpWidthF          = 0.9 ; Resize to a square
   res@vpHeightF         = 0.6
  ;res@vpXF                 = 0.25        ; change start locations
  ;res@vpYF                 = 0.75        ; the plot

   if(wtype.eq."png") then
    lthick = 3.0
   else
    lthick = 1.0
   end if
  
   res@tmBorderThicknessF       = lthick
   res@tmXBMajorThicknessF      = lthick
   res@tmXBMinorThicknessF      = lthick
   res@tmYLMajorThicknessF      = lthick
   res@tmYLMinorThicknessF      = lthick
   res@tmYRMajorThicknessF      = lthick
   res@tmYRMinorThicknessF      = lthick

   res@cnLevelSelectionMode = "ExplicitLevels"

  ;*******************************************
  ; first plot
  ;*******************************************
   res@mpLabelFontHeightF          = FontHeightF*2.0
   res@mpGridLineThicknessF        = lthick* 0.5
   res@mpLimbLineThicknessF        = lthick* 0.5
   res@mpPerimLineThicknessF       = lthick* 0.5
   res@mpGeophysicalLineThicknessF = lthick

   res@gsnLeftString     = ""
   res@gsnRightString    = "" 
   res@gsnCenterString   = "" 
   res@tiMainString      = ""
   res@tiMainOffsetYF    = 0.005

   res0                      = res
   res0@cnFillPalette        = "NCV_blu_red" ; BlueWhiteOrangeRed"
   res0@cnLevelSelectionMode = "ExplicitLevels"
   res0@tiMainString         = exp_str
   res0@gsnLeftString        = left_str
   res0@gsnRightString       = right_str
   res0@cnLevels            := cnlevs1
   res0@cnLinesOn            = True        ; True is default
   res0@cnMinLevelValF       = min(cnlevs1)
   res0@cnMaxLevelValF       = max(cnlevs1)
   res0@cnLevelSpacingF      = 10.0
   res0@cnLineColor          = "Black"                   ; color of second contours
   res0@cnLineThicknessF     = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1

   plot(k) = gsn_csm_contour_map(wks,eofm,res0)

   delete([/eof,eof_pcv,eofm,yyyy/])    
  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012

   pres1@gsnPanelLabelBar           = True                ; add common colorbar
   pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = ""+vars + " (" + varunt +")" 
                                      ;main_str + " of " + vars + " (" + varunt +")"
   pres1@lbTitlePosition            = "Right"     ; title location
   pres1@lbTitleDirection           = "Across"    ; letter angle
   pres1@lbTitleAngleF              = 90.         ; title angle
   pres1@lbTitleFontHeightF         = FontHeightF*0.3
   pres1@lbLabelFontHeightF         = FontHeightF*0.3 ;  font height
   pres1@lbBoxLineColor             = "black"
   pres1@lbBoxLinesOn               = True
   pres1@pmLabelBarParallelPosF     = -0.01 ;0.61
   pres1@pmLabelBarOrthogonalPosF   = 0.01 ;0.0
   pres1@pmLabelBarWidthF           = 0.03 ;0.95
   pres1@pmLabelBarHeightF          = 0.14 ;0.04

  ;pres1@gsnPanelBottom             = 0.45       ; space for label bar
  ;pres1@gsnPanelTop                = 0.85       ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 8
   pres1@gsnPanelXWhiteSpacePercent = 4
   gsn_panel(wks,plot,(/1,nexps/),pres1)          ; create first panel

end do 

end

