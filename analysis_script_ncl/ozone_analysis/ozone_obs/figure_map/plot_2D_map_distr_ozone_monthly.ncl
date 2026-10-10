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

  data_dir   = "../data/ozone_vint/"
  exp_name   = (/"ERA5", "v2.LR.historical","v2.LR.amip"/)
  exp_strs   = (/"ERA5", "v2 Coupled",      "v2 AMIP"/)
  nexps      = dimsizes(exp_name)
  LabList    = exp_strs

  time_tag   = (/"1979-2014"/)
  seasons    = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  varnam     = (/"TCO"/)
  facmod     = (/ 1.0 /)
  varunt     = (/"DU"/)
  
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  MinLat    =  -90.  ;negative for southern hemisphere
  MaxLat    =    0.  ;negative for southern hemisphere
  MinLon    = 120.0
  MaxLon    = 300.0

 ;;; flags to control the plot style
  lndcolor = "Grey28"

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  do k = 0,nexps-1

   fl  = systemfunc("ls " + data_dir+exp_name(k)+"_en*_"+varnam+"_MON_"+time_tag+".nc")
   f   = addfiles(fl,"r")
   ListSetType (f, "join")

   lat  = f[0]->lat
   nens = dimsizes(f[:]->$varnam$(:,0,0,0))
   vdat = f[:]->$varnam$

   do ie = 0,nens,1

   if(ie.eq.0) then
    var0 = dim_avg_n_Wrap(vdat,0)
    varm = clmMonTLL(var0)
    enstr = "mean"
   else
    var0 = vdat(ie-1,:,:,:)
    varm = clmMonTLL(var0)
    enstr = "en" + sprinti("%02d",ie)
   end if
   var = varm(i,:,:)

   delete([/varm,var0/])

   ;printVarSummary(var)

   ;;;;calculate the global mean and correlation 
   latrg = lat({MinLat:MaxLat})
   rad   = 4.*atan(1.)/180.    ;
   cost  = cos(latrg*rad)       ; cosine weights
   mvar1 = wgt_areaave_Wrap (var({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
   mstr1 = "Mean: " + sprintf("%5.2f", mvar1)

  ;----------------------------------------------------------------------
  ; Plotting section

   if(ie.eq.0) then 
     OUTDir  = "./"
     wtype = "png"
     wtype@wkWidth = 2500
     wtype@wkHeight = 2500
    ;wtype@wkPaperSize     = "A4"
    ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
    ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
     fignam = "fig_2d_map_"+Regions+"_"+exp_name(k)+"_"+varnam+"_" + seasons(i)
     wks1 = gsn_open_wks(wtype,OUTDir+fignam)
     plot = new(nens+1,graphic)
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

    res@gsnAddCyclic         = True
   ;res@mpShapeMode          = "FreeAspect"
   ;res@vpWidthF             = 0.9 ; Resize to a square
   ;res@vpHeightF            = 0.6

   ;---Control appearance of map.
   ;res@gsnMajorLatSpacing      = 10
   ;res@gsnMajorLonSpacing      = 30
   ;res@gsnMinorLatSpacing      = 5
   ;res@gsnMinorLonSpacing      = 10

;    res@mpMinLonF               = MinLon
;    res@mpMaxLonF               = MaxLon
;    res@mpMinLatF               = MinLat
;    res@mpMaxLatF               = MaxLat
;    res@mpCenterLonF            = 180     ; This is necessary to get the correct map
;
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

;    res@mpProjection             = "LambertEqualArea"   ;-- change projection
;    res@mpDataBaseVersion        = "MediumRes"      ;-- map resolution
;    res@mpGridAndLimbOn          =  True            ;-- plot grid lines
;    res@mpGridLineColor          = "grey30"         ;-- set grid line color
;    res@mpPerimOn                =  False           ;-- don't draw the box around the plot

    res@gsnPolar             = "SH"                  ; specify the hemisphere
    res@mpGridLatSpacingF    = 45.                ; change latitude  line spacing
    res@mpGridLonSpacingF    = 90.                ; change longitude line spacing


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
    res1@cnLevels                      := ispan(220,340,10) ;(/0,1,2,4,8,10,12,14,16,18,20,24,28/)
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
    res1@tiMainString                   = LabList(k) + " ("+enstr+")"
    res1@gsnLeftString                  = varnam+" (" + varunt+")"
    res1@gsnRightString                 = mstr1 

   ;plot(ie) = gsn_csm_contour_map(wks1,var,res1)
    plot(ie) = gsn_csm_contour_map_polar(wks1,var,res1)

    delete([/res,res1/])
   end do 
 
   delete([/fl,f,lat,nens,vdat/])
  
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
    gsn_panel(wks1,plot,(/3,4/),pnres)

    delete([/plot,pnres/])

  end do 

 end do 

end

