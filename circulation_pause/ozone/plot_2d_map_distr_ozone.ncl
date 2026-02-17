;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  data_dir  = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output/"
  exp_name  = (/"CLIM_30yr","SSP245_30yr","SSP585_30yr"/)
  labList   = (/"AMIP","AMIP+SSP245","AMIP+SSP585"/)           
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014

  time_tag  = (/"1979-2014"/)
  seasons   = (/"SON"/)
  nseas     = dimsizes(seasons)
  maxlev    = 11

  varList   = (/"TCO"/)
  varStrs   = (/"Total column ozone"/)
  varUnit   = (/"DU"/)
  nvars     = dimsizes(varList)

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  latS      = -90. ;negative for southern hemisphere
  latN      = -60. ;negative for southern hemisphere
  lonW      = 120.0
  lonE      = 300.0

 ;;; flags to control the plot style
  lndcolor = "Grey28"

 do j = 0,nvars-1,1

  varnam = varList(j)
  varunt = varUnit(j)

  do i = 0, nseas-1,1

   do k  = 0,nexps-1,1

     fl1 = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/monthly/"+"*.nc")
     f1  = addfiles(fl1,"r")
    ;ListSetType (f1, "join")

     lat  = f1[0]->lat
     lon  = f1[0]->lon
     rad  = 4.0*atan(1.0)/180.0
     clat = lat
     clat = cos(lat*rad)

     tmp0 = f1[:]->SCO
     tmp0 = tmp0+ f1[:]->TCO
     tmp1 = month_to_seasonN(tmp0,(/seasons(i)/))
     vars = tmp1(0,:,:,:)
     vars&time = ispan (ystr,yend,1)
     varm = dim_avg_n_Wrap(vars(:,:,:),0) 
     mvar1 = wgt_areaave_Wrap(vars({1979:1981},{latS:latN},:),clat({latS:latN}),1.,1)
     mstr1 = "Mean: " + sprintf("%5.2f", avg(mvar1))
     delete([/tmp0,tmp1,vars/])
  ;----------------------------------------------------------------------
  ; Plotting section

   if(k.eq.0) then 
     OUTDir  = "./"
     wtype = "png"
     wtype@wkWidth = 2500
     wtype@wkHeight = 2500
    ;wtype@wkPaperSize     = "A4"
    ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
    ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
     fignam = "fig_2d_map_"+exp_name(k)+"_"+varnam+"_" + seasons(i)
     wks1 = gsn_open_wks(wtype,OUTDir+fignam)
     plot = new(nexps,graphic)
     dum  = new(nexps,graphic)
   end if 

  ;---Resources to share between both plots
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

    res2=res

    res@gsnAddCyclic         = True
   ;res@mpShapeMode          = "FreeAspect"
   ;res@vpWidthF             = 0.9 ; Resize to a square
   ;res@vpHeightF            = 0.6

   ;---Control appearance of map.
   ;res@gsnMajorLatSpacing      = 10
   ;res@gsnMajorLonSpacing      = 30
   ;res@gsnMinorLatSpacing      = 5
   ;res@gsnMinorLonSpacing      = 10

    res@mpMinLonF               = lonW
    res@mpMaxLonF               = lonE
    res@mpMinLatF               = latS
    res@mpMaxLatF               = latN
    res@mpCenterLonF            = 180     ; This is necessary to get the correct map

    res@mpFillOn                = False        ; turn off map fill
    res@mpLabelsOn              = False
    res@mpPerimOn               = True
    res@mpFillOn                = False
    res@mpGeophysicalLineColor  = lndcolor  ; color of cont. outlines
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

    res@mpProjection             = "LambertEqualArea"   ;-- change projection
    res@mpDataBaseVersion        = "MediumRes"      ;-- map resolution
    res@mpGridAndLimbOn          =  True            ;-- plot grid lines
    res@mpGridLineColor          = "grey30"         ;-- set grid line color
    res@mpPerimOn                =  False           ;-- don't draw the box around the plot

    res@gsnLeftString                  = ""
   ;res@gsnLeftStringOrthogonalPosF    =  0.02
   ;res@gsnLeftStringParallelPosF      =  0.05
    res@gsnRightString                 = ""
   ;res@gsnRightStringOrthogonalPosF   = 0.02
   ;res@gsnRightStringParallelPosF     = 0.98
    res@gsnCenterString                = "" 
   ;res@gsnCenterStringOrthogonalPosF  = 0.04

    res@cnFillOn                 = True              ; color plot desired
   ;res@cnFillPalette            = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
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

    res1=res
    res1@cnFillPalette                  = "BlAqGrYeOrReVi200" ;"MPL_rainbow" ;"BlAqGrYeOrReVi200"
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                      := ispan(200,360,20) ;(/0,1,2,4,8,10,12,14,16,18,20,24,28/)
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@tiMainString                   = ""
    res1@tiMainOffsetYF                 = 0.01
    res1@cnLinesOn                      = True 
    res1@cnLineLabelsOn                 = False
    res1@cnInfoLabelOn                  = False
    res1@cnLineLabelDensityF            = 2.0
    res1@cnLineLabelFontHeightF         = 0.018
    res1@cnLineLabelBackgroundColor     = -1
    res1@gsnContourLineThicknessesScale = 1.0
    res1@cnLineThicknessF               = 1.0
    res1@cnLineColor                    = "Black"
    res1@gsnContourZeroLineThicknessF   = 2.
    res1@gsnContourNegLineDashPattern   = 1
    res1@tiMainString                   = labList(k) 
    res1@gsnLeftString                  = varnam+" (" + varunt+")"
    res1@gsnRightString                 = mstr1 
    plot(k) = gsn_csm_contour_map(wks1,varm,res1)

    res2@cnLevelSelectionMode           = "ExplicitLevels"
    res2@cnLevels                      := (/150,220,avg(mvar1)*0.7,avg(mvar1)*0.5/) 
    res2@tiMainString                   = ""
    res2@gsnLeftString                  = ""
    res2@gsnRightString                 = ""
    res2@cnFillOn                       = False 
    res2@cnLinesOn                      = True
    res2@cnLineThicknessF               = 5.0 
    res2@cnInfoLabelOn                  = False
    res2@cnLineColor                    = "White"
    dum(k) = gsn_csm_contour(wks1,varm,res2)
    
    overlay(plot(k),dum(k))
    delete([/res,res1,res2/])

   end do 
   
;---Draw both plots in a panel
    pnres                            = True
    pnres@gsnMaximize                = True
    pnres@gsnPanelMainFontHeightF    = 0.016
    pnres@gsnPanelMainString         = ""
    pnres@gsnPanelMainFontColor      = "Black"
    pnres@gsnPanelBottom             = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5
    pnres@gsnPanelBottom             = 0.4        
    gsn_panel(wks1,plot,(/1,3/),pnres)

    delete([/plot,pnres/])

  end do 

 end do 

end

