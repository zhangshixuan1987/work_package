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

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/Global_2d_map/data"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP","v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_strs

  mode       = ("SAM")
  ieof       = 1

  varList    = (/"SST","TAUX"/)
  vsrList    = (/"SST","TAUU"/)
  facList    = (/  1.0,  1e2/)
  untList    = (/"~S~o~N~C", "10~S~-2~N~ N m~S~-2~N~"/)
  nvars      = dimsizes(varList)

  time_tag   = (/"1979-2014"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11
  statistical_siglev = 0.05

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"Global"/)
  nregs      = dimsizes(Regions)
  MinLat     = -90.
  MaxLat     =  20.
  MinLon     =  0.
  MaxLon     = 360.

 ;declear the season or annual mean data to be used;;;;;
  do i = 0,nseas-1,1 

   do j = 0,nvars-1,1

    vars       = varList(j)
    vstr       = vsrList(j)
    facmod     = facList(j)
    varunt     = untList(j)

    do k = 0,nexps-1,1

      print("working on season "+ seasons(i))
      fl1   = "../data/"+exp_name(k)+"_EOF_PCS_Z500_"+seasons(i)+"_"+time_tag+".nc"
      fl2   = data_dir+"/"+exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_"+time_tag+".nc"
      f1    = addfile(fl1,"r")
      f2    = addfile(fl2,"r")

      rad   = 4.0*atan(1.0)/180.0
      wgty  = cos(f2->lat*rad)

      vts   = f1->eof_ts(:,ieof-1,:)
      nens  = dimsizes(vts(:,0))

      if(exp_name(k).eq."v2.LR.amip") then 
        varstr    = "d"+vars
        var0      = f2->$varstr$
        var0      = var0 * facmod
        var0@units= varunt
        ntime = dimsizes(var0(:,0,0,0))
        nlat  = dimsizes(var0(0,0,:,0))
        nlon  = dimsizes(var0(0,0,0,:))
        var = new((/ntime,nens,nlat,nlon/),typeof(var0))
        do ie = 0,nens-1,1
          var(:,ie,:,:) = (/var0(:,0,:,:)/)
        end do 
        var!0 = "time"
        var!1 = "ens"
        var!2 = "lat"
        var!3 = "lon"
        var&time = var0&time 
        var&lat  = var0&lat 
        var&lon  = var0&lon
        delete([/var0,ntime,nlat,nlon/])
      else
        varstr    = "d"+vars
        var       = f2->$varstr$
        var       = var * facmod
        var@units = varunt
      end if 

      printVarSummary(vts)
      printVarSummary(var)
      xrc1  = var(0,:,:,:)
      xrgpb = var(0,:,:,:)
      do ie = 0,nens-1,1
        xx     = dtrend_msg_n(vts&time,vts(ie,:),True,False,0)
        yy     = dtrend_msg_n(var&time,var(:,ie,:,:),True,False,0)
        rc1    = regCoef_n(xx, yy, 0, 0)
        ;;;;;;;t-test;;;;;;;;
        drc1   = dimsizes(rc1)
        df1    = rc1@nptxy-2   ; degrees of freedom
        tval1  = rc1@tval      ; t-statistic
        xb1    = tval1
        xb1    = 0.50
        pb1    = betainc(df1/(df1+tval1^2),df1/2.0,xb1)
        prob1  = reshape(pb1,(/drc1(0),drc1(1)/))
        prob1  = where(prob1.le.statistical_siglev,rc1,rc1@_FillValue)
        xrc1(ie,:,:)  = (/rc1/)
        xrgpb(ie,:,:) = (/  tofloat(prob1) /)
        delete([/xx,yy,drc1,df1,tval1,xb1,pb1,prob1,rc1/])
      end do
 
      vard  = dim_avg_n_Wrap(xrc1,0)
      vpob  = dim_avg_n_Wrap(xrgpb,0)
      delete([/xrc1,xrgpb,vts,var/]) 
      printVarSummary(vard)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    if(k.eq.0) then 
      wtype           = "png"
      wtype@wkWidth   = 2500
      wtype@wkHeight  = 2500
      fignam          = "SAM_REG_detrend_"+vars+"_"+seasons(i)+"_"+time_tag
      wks             = gsn_open_wks(wtype,fignam)
      plot      = new((/nexps/),graphic)
      dum       = new((/nexps/),graphic)

    end if 

    ;get the information for contour plot
    load "./share_fig_info.ncl"

    ;============================================================
    ; PLOTS
    ;============================================================
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False

    FontHeightF              = 0.028
    res@tiMainFontHeightF    = FontHeightF
    res@tmYLLabelFontHeightF = FontHeightF
    res@tmXBLabelFontHeightF = FontHeightF
    res@lbLabelFontHeightF   = FontHeightF
    res@gsnStringFontHeightF = FontHeightF*1.2

    res@vpWidthF             = 0.9 ; Resize to a square
    res@vpHeightF            = 0.6

    res@mpShapeMode        = "FreeAspect"
    res@gsnAddCyclic       = False

    ;---Control appearance of map.
    res@gsnMajorLatSpacing      = 20
    res@gsnMajorLonSpacing      = 60
    res@gsnMinorLatSpacing      = 10
    res@gsnMinorLonSpacing      = 10

    res@mpMinLonF               = MinLon
    res@mpMaxLonF               = MaxLon
    res@mpMinLatF               = MinLat
    res@mpMaxLatF               = MaxLat
    res@mpCenterLonF            = 180     ; This is necessary to get the correct map

    res@mpFillOn                = False        ; turn off map fill
    res@mpLabelsOn              = False
    res@mpPerimOn               = True
    res@mpFillOn                = False
    res@mpGeophysicalLineColor  = "Grey25"  ; color of cont. outlines
    res@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
    res@mpOutlineOn              = True
    res@mpOutlineDrawOrder       = "PostDraw"
    res@mpFillDrawOrder          = "Predraw"
    res@mpGridAndLimbOn          = False
    res@mpGridLineThicknessF     = 1.0
;   res@mpOceanFillColor         = "lightskyblue1"
;   res@mpLandFillColor          = "gray33"
    res@tmYLLabelsOn             = True
    res@tmYLOn                   = True

    lthick = 1
    res@mpGridLineThicknessF        = lthick
    res@mpLimbLineThicknessF        = lthick
    res@mpPerimLineThicknessF       = lthick
    res@mpGeophysicalLineThicknessF = lthick

    res@gsnLeftString                  = ""
    res@gsnRightString                 = ""
    res@gsnCenterString                = ""
   ;res@gsnLeftStringParallelPosF      = 0.05
   ;res@gsnRightStringParallelPosF     = 0.98
   ;res@gsnLeftStringOrthogonalPosF    = 0.02
   ;res@gsnRightStringOrthogonalPosF   = 0.02
   ;res@gsnCenterStringOrthogonalPosF  = 0.04

    res@cnFillOn                 = True              ; color plot desired
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
    res@cnMissingValFillColor    = "White"
    res@lbLabelBarOn             = True      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF   = 0.61
    res@pmLabelBarOrthogonalPosF = 0.05
    res@lbOrientation            = "vertical"        ; vertical label bar
    res@pmLabelBarHeightF        = 0.65
    res@pmLabelBarWidthF         = 0.12
    res@lbTitleString            = ""
    res@lbTitlePosition          = "Right"                           ; title location
    res@lbTitleDirection         = "Across"                          ; letter angle
    res@lbTitleAngleF            = 90.                               ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2                              ; font height


   ;;;;difference between E3SM and ERA5
    gsn_define_colormap(wks,"temp_diff_18lev")
    res1=res
   ;res1@cnFillPalette                  = "temp_diff_18lev"
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                      := cnlevs1
    res1@cnFillColors                  := (/2,4,6,8,9,10,12,13,14,16,18,20/)
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@cnLinesOn                      = False             ; turn off contour lines
    res1@cnLineLabelsOn                 = False            ; turn off contour labels
    res1@cnInfoLabelOn                  = False
    res1@gsnLeftString                  = left_strs(k) + " (" +seasons(i) +")"
    res1@gsnRightString                 = "~F33~D~N~~F~"+ vars + " (" + varunt+")"
    plot(k)                             = gsn_csm_contour_map(wks,vard,res1)
    delete([/res,res1/])

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
  end do

  ; panel plot only resources
   pres1                            = True                ; mods desired
   pres1@gsnFrame                   = True               ; save panel until both ready
  ;pres1@gsnPanelBottom             = 0.45    ; space for label bar
  ;pres1@gsnPanelTop                = 0.85     ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 8
   gsn_panel(wks,plot,(/1,3/),pres1)          ; create first panel
 
 end do      

end do 

end

