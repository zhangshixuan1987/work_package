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

  ModelOutRoot = "../data/"

 ;;;;;;;strings for the data directory;;;;;;;;
  Groups      = (/"ERA5","CLIM","NDGUV","NDGT","NDGQ","NDGUVT","NDGUVTQ"/)
  Expnams     = Groups
  ngrps       = dimsizes(Groups)

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList      = (/"PRECT", "U850","TS"/) ;,"TREFHT"/)
  VstList      = (/"PRECT", "U850","TS"/) ;,"TREFHT"/)
  UnitList     = (/"mm day~S~-1~N~","m s~S~-1","K","K"/)
  PltNums      = (/"a", "b", "c", "c","d","e","f","g","h","i","j","k"/)
  nvars        = dimsizes(VarList)

  do iv = 0,nvars-1,1

    Varname  = VarList(iv)
    Varunt   = UnitList(iv)
    Varstr   = VstList(iv)

    do igp = 0, ngrps-1,1

     fna    = ModelOutRoot+ VarList(iv) + "_" + Groups(igp) + ".nc"
     fla    = addfile(fna, "r")

     vna    = VarList(iv) + "_PDF_" + Groups(igp)
     vpdf   = fla->$vna$

     if(igp.eq.0) then 
       nbin = dimsizes(vpdf@bin_center) 
       xx   = new((/ngrps,nbin/),double)
       yy   = new((/ngrps,nbin/),double)
     end if 

     xx(igp,:) = vpdf@bin_center     
     yy(igp,:) = vpdf

     delete([/vpdf/])

    end do 

    xx@_FillValue = -9999
    yy@_FillValue = -9999
    yy = where(yy.le.1e-6, yy@_FillValue, yy)
    xx = where(yy.le.1e-6, xx@_FillValue, xx)

   ;----------------------------------------------------------------------
   ; Plotting section
    wtype = "pdf"
   ;wtype@wkPaperSize     = "A4"
   ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
 
    if (iv.eq.0) then 
      figdir = "./"
      fignam = "fig_pdf_ndgvar"
      wks    = gsn_open_wks(wtype,figdir+fignam)
      plot   = new((/nvars/),graphic)
    end if 

    gsn_define_colormap(wks,"cb_9step")
    lndcolor = "Grey27"

    LineThick = (/5, 5,  5,  5,  5,  5,  5/) / 2.0
    LineColor = (/2, 3, 14, 28, 37, 45, 58/)
    LinePattn = (/0, 0,  0,  0,  0,  0,  0/)

    print("")
    print("plot figure: "+ fignam)
    print("")


    res                          = True              ; Plot modes desired.
    res@gsnDraw                  = False
    res@gsnFrame                 = False
    res@gsnMaximize              = False              ; Maximize plot

    res@vpWidthF  = 0.9
    res@vpHeightF = 0.4

    FontHeightF = 0.024
    res@tiMainFontThicknessF    = 1.0
    res@tiMainFontHeightF       = FontHeightF*1.0
    res@tmYLLabelFontHeightF    = FontHeightF
    res@tmXBLabelFontHeightF    = FontHeightF
    res@tiXAxisFontHeightF      = FontHeightF*1.0
    res@tiYAxisFontHeightF      = FontHeightF*1.0
    res@gsnStringFontHeightF    = FontHeightF*1.0

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

    res@tmXTOn                = False
    res@tmYROn                = False
    res@tiYAxisString         = "PDF (%)"
    res@tiXAxisString         = Varname + " ("+Varunt+")"

    res@tiMainString          = ""
    res@gsnCenterString       = ""
    res@gsnRightString        = "" 
    res@gsnLeftString         = ""

    if(Varname.eq."PRECT") then 
      res@trYLog                = True
      res@xyComputeYMin         = False
      res@trYMinF               = 1e-6 ;floor(min(yy(0,:)))
      res@trYMaxF               = 1e3  ;ceil(max(yy(0,:)))
      res@trXMinF               = 0.0  ;floor(min(xx(0,:)))
      res@trXMaxF               = ceil(max(xx(0,:))) + 10
    else if (Varname.eq."U850") then
      res@trYLog                = True ;False 
      res@xyComputeYMin         = False
      res@trYMinF               = 1e-6 ;floor(min(yy(0,:)))
      res@trYMaxF               = 1e3  ;ceil(max(yy(0,:)))
      res@trXMinF               = floor(min(xx(0,:)))
      res@trXMaxF               = ceil(max(xx(0,:)))
    else
      res@trYLog                = True ;False
      res@xyComputeYMin         = False
      res@trYMinF               = 1e-6 ;floor(min(yy(0,:)))
      res@trYMaxF               = 1e3  ;ceil(max(yy(0,:)))
      res@trXMinF               = floor(min(xx(0,:))) - 2
      res@trXMaxF               = ceil(max(xx(0,:)))  + 2
    end if 
    end if 

    ;---Bar chart resources
    bres                          = res
    bres@gsnLeftString            = "Probability Distribution of " + Varname 
    bres@gsnRightString           = ""
   ;bres@gsnYRefLine              = 1e-6              ; reference line   
   ;bres@gsnYRefLineColor         = "transparent"
    bres@gsnXYBarChart            = True            ; Turn on bar chart
    bres@xyLineColors             = "transparent" ;"black" ;LineColor
    if(iv.eq.0) then 
      bres@gsnXYBarChartBarWidth  = 500./dimsizes(xx(0,:))
    else if (iv.eq.1) then
      bres@gsnXYBarChartBarWidth  = 30./dimsizes(xx(0,:))
    else 
      bres@gsnXYBarChartBarWidth  = 80./dimsizes(xx(0,:))
    end if 
    end if 
    bres@gsnXYBarChartColors      = "Grey40"        ; Color for bars
    bres@gsnXYBarChartOutlineOnly = False             
    bres@tmYROn                   = False
    plot(iv) = gsn_csm_xy (wks,xx(0,:),yy(0,:),bres)

    xres = res
    xres@xyLineThicknesses        = LineThick
    xres@xyLineColors             = LineColor
    xres@xyMonoDashPattern        = True              ; all solid

    if(iv.eq.0) then 
      xres@pmLegendDisplayMode      = "Always"            ; turn on legend
    else 
      xres@pmLegendDisplayMode      = "Never"            ; turn on legend
    end if 
    xres@pmLegendSide             = "Top"               ; Change location of 

    if(iv.gt.1)then 
      xres@pmLegendParallelPosF   = .15                ; move units right
    else 
      xres@pmLegendParallelPosF   = .85                ; move units right
    end if 
    xres@pmLegendOrthogonalPosF   = -0.75               ; move units down
    xres@pmLegendWidthF           = 0.12                ; Change width and
    xres@pmLegendHeightF          = 0.20                ; height of legend.
    xres@lgPerimOn                = True                ; turn off/on box around
    xres@lgLabelFontHeightF       = FontHeightF*0.6     ; label font height
    xres@xyExplicitLegendLabels   = " "+Expnams(1:)

    tmp_plot = gsn_csm_xy (wks, xx(1:,:), yy(1:,:), xres)

    overlay(plot(iv),tmp_plot)
    delete([/res,xres,bres/])

  end do

  ;---Draw both plots in a panel
  pnres                            = True
  pnres@gsnMaximize                = True
  pnres@gsnPanelMainFontHeightF    = 0.016
  pnres@gsnPanelMainString         = ""
  pnres@gsnPanelMainFontColor      = "Black"
  pnres@gsnPanelBottom             = 0.05
  pnres@gsnPanelYWhiteSpacePercent = 2
  pnres@gsnPanelXWhiteSpacePercent = 2
  pnres@gsnPanelBottom             = 0.4        
  gsn_panel(wks,plot,(/1,nvars/),pnres)
  delete([/pnres/])

end
