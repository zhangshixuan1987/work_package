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

  data_file  = "E3SM_versus_ERA5.nc"
  seasons    = (/"DJF"/)
  nseas      = dimsizes(seasons)

  ;;;read the experiment names and short names;;;
  exp_list     = "experiment_list.txt"
  data         = asciiread(exp_list,-1,"string")
  delim        = ","
  Groups       = str_squeeze(str_get_field(data, 2,delim))
  ngrps        = dimsizes(Groups)
  LabList      = new(ngrps+1,string)
  LabList(0)   = (/"ERA5"/) 
  LabList(1:)  = Groups + " ~F19~:~N~~F~ " + LabList(0)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList      = (/"T","Q","U","RELHUM"/)
  VnmList      = (/"Temperature", "Zonal wind", "Specific humidity", "Relative humidity"/)
  nvars        = dimsizes(VarList)
  color_grp    = 1 ; select a color scale

  ;;;;file names for the plot;;;;
  minz         = 500
  maxz         = 1000
  lats         = -90 
  latn         =  0

  ;;; flags to control the plot style
  plotcntl     = False   ; if True, plot the CNTL origianl field on the leftmost panel
  plot_contour = False  ; if True, plot the contour in the figure  
  plot_leftstr = True   ; if True, the lefadd igure captions as left string 

  plot         = new(nvars*2,graphic)
  dum0         = new(nvars*2,graphic)
  dum1         = new(nvars*2,graphic)
  dum2         = new(nvars*2,graphic)
  dum3         = new(nvars*2,graphic)
  dum4         = new(nvars*2,graphic)

  f = addfile(data_file,"r")

  do itm = 0,nseas-1,1

  season = seasons(itm)

  do igp = 0,ngrps-1,1

  do iv = 0,nvars-1, 1

   Varname = VarList(iv)
   Varstr  = VnmList(iv)

   vobs = "ERA5_"+Varname+"_"+season
   vmod = Groups(igp)+"-ERA5_"+Varname+"_"+season
   var0 = f->$vobs$ 
   var1 = f->$vmod$ 
   var0&lev@units = "hPa"
   var1&lev@units = "hPa"

   printVarSummary(var0)
   printVarSummary(var1)
   Varunt = var0@units
 
   ;;;time and meridionally average;;;
   vav0  = var0 ;dim_avg_n_Wrap(var0,2)
   vadif = var1 ;dim_avg_n_Wrap(var1,2)
   ;;;calculate the global mean;;;
   rad         = 4.0*atan(1.0)/180.0
   clat        = cos(var0&lat*rad)
   vav0@mean   = wgt_areaave(vav0, 1.0, clat, 0)
   vadif@mean  = wgt_areaave(vadif, 1.0, clat, 0)
   vadif@rmse  = wgt_arearmse(vav0,vav0+vadif,1.0, clat, 0)
   vadif@pcor  = pattern_cor(vav0,vav0+vadif,1.0,0)

 ;----------------------------------------------------------------------
 ; Plotting section
  plotFileLabel = "fig_pres_lat_"+Groups(igp)+"_"+Varname+"_"+season
  ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "png"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   
  wks1 = gsn_open_wks(wtype,OUTDir+plotFileLabel)
  load "./var_share_colorbar_preshgt.ncl"


 ;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot

  res@vpWidthF             = 0.8
  res@vpHeightF            = 0.6

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYRMode                = "Automatic"
  res@tiYAxisString           = "Pressure (hPa)"

  res@trYMinF                 = minz    
  res@trYMaxF                 = maxz
  res@trYLog                  = False

  res@tmXBMode                = "Explicit"
  res@tmXBValues              = ispan(toint(lats),toint(latn),15)
  res@tmXBLabels              = (/"90S","75S","60S","45S","30S","15S","0"/)
  res@tmXBMinorValues         = ispan(toint(lats),toint(latn),5)

  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/ 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings

  res@tiMainFontThicknessF    = 3.0
  res@tiMainFontHeightF       = 0.030
  res@tmYLLabelFontHeightF    = 0.028
  res@tmXBLabelFontHeightF    = 0.028
  res@tiXAxisFontHeightF      = 0.030
  res@tiYAxisFontHeightF      = 0.030
  res@gsnStringFontHeightF    = 0.030
  res@cnLineLabelFontHeightF  = 0.028
  res@lbLabelFontHeightF      = 0.028
  res@lbTitleFontHeightF      = 0.028                             ; font height
  
  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
  res@lbLabelStride           = 1                 ; skip every other label
  res@lbLabelFontHeightF      = 0.028
  res@lbPerimOn               = False             ; default has box
  res@lbOrientation           = "vertical"        ; vertical label bar
  res@pmLabelBarHeightF       = 0.6
  res@pmLabelBarWidthF        = 0.12
; res@pmLabelBarParallelPosF  = 0.026
; res@pmLabelBarOrthogonalPosF= 0.15

  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle

  leftstr                     = "" ;"Surface Temperature ("+ts@units+")"
  res@gsnLeftString           = Varstr+" ("+ Varunt +")" ; leftstr ;var_name0    ; long_name is too long!
  res@gsnRightString          = "" ;" ("+Varunt+")";""
 ;res@tiMaingString           = "Annual Mean Response to Present-day climate"

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"

;---Resources for plotting original data
  if(plotcntl.and.igp.eq.1)then
    res1=res
    res1@cnFillPalette                 := ocname
    res1@cnLevelSelectionMode           = "ExplicitLevels"  
    res1@cnLevels                      := orglev
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False  
   ;res1@lbLabelStride                  = 15
    res1@gsnLeftString                  = Varstr+" ("+ Varunt +")" ;Figstr1 + LabList(0) ;Vnmstr + "  ("+Varstr+")"
    res1@gsnLeftStringOrthogonalPosF    = 0.02
   ;res1@gsnLeftStringParallelPosF      = 0.2
    res1@gsnRightString                 = "";Varstr ;Varunt
    res1@gsnRightStringOrthogonalPosF   = 0.02
    ;if(iv.eq.0)then
    ;  res1@tiMainString                 = LabList(0);"CNTL (PD1800)"
    ;else
    ;  res1@tiMainString                 = "";"CNTL (PD1800)"
    ;end if 
    res1@cnLinesOn                      = True ;False
    res1@cnLineLabelsOn                 = False     
    res1@cnInfoLabelOn                  = False
    res1@cnLineLabelDensityF            = 2.0
    res1@cnLineLabelFontHeightF         = 0.018
    res1@cnLineLabelBackgroundColor     = -1
    res1@gsnContourLineThicknessesScale = 1.0
    res1@cnLineThicknessF               = 0.01
    res1@cnLineColor                    = "Black"
    res1@gsnContourZeroLineThicknessF   = 2.
   ;res1@gsnContourNegLineDashPattern   = 1
    res1@lbTitleString                  = "" ;"("+ Varunt+")" 
   ;plot(iv)  = gsn_csm_contour(wks1,vav0,res1)
    plot(iv)  = gsn_csm_pres_hgt(wks1,vav0,res1)
    
  end if 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@cnLevelSelectionMode           = "ExplicitLevels"  
    res2@cnFillColors                  := diffcolors
    res2@cnLevels                      := diffgrp1
    res2@cnLinesOn                      = False 
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@gsnContourZeroLineThicknessF   = 2.
    res2@gsnContourNegLineDashPattern   = 1
    res2@tmYLOn                         = True
   ;res2@tiYAxisString                  = ""
    ;if(iv.eq.0)then
    ;  res2@tiMainString                 = LabList(1);"CNTL (PD1800)"
    ;else
    ;  res2@tiMainString                 = "";"CNTL (PD1800)"
    ;end if
    res2@gsnLeftString                  = "~F33~D~N~~F~"+Varstr+" ("+ Varunt +")"  ;Figstr2 + LabList(1) ;"~F33~D~N~~F~"+ Varstr
    res2@gsnLeftStringOrthogonalPosF    = 0.02
   ;res2@gsnLeftStringParallelPosF      = 0.2 
    res2@gsnRightString                 = "" ;Varstr ;"~F33~D~N~~F~"+ Varstr ;Varunt
    res2@gsnRightStringOrthogonalPosF    = 0.02
    res2@lbTitleString                  = "" ;"("+ Varunt +")"
    ;plot(iv+nvars)                  = gsn_csm_contour(wks1, vadif,res2)
    plot(iv+nvars)                  = gsn_csm_pres_hgt(wks1,vadif,res2)

    polyres                   = True
    polyres@gsMarkerIndex     = 16          ; polymarker style
    polyres@gsMarkerSizeF     = 0.0012          ; polymarker size
    polyres@gsMarkerColor     = "Black" ;"Magenta"         ; polymarker size
   ;dum1(iv)  = gsn_add_polymarker(wks1,plot(iv+nvars),vav0&lat,pre0,polyres)

    if(plot_contour)
     res21 = res2
     res21@gsnLeftString                  = "" 
     res21@gsnRightString                 = ""
     res21@gsnCenterString                = ""
     res21@cnFillOn                       = False
     res21@cnLinesOn                      = True
     res21@cnLineLabelsOn                 = False            ; turn off contour labels
     res21@cnInfoLabelOn                  = False
     res21@cnLineLabelDensityF            = 2.0
     res21@cnLineLabelFontHeightF         = 0.018
     res21@cnLineLabelBackgroundColor     = -1
     res21@gsnContourLineThicknessesScale = 1.0
     res21@cnLineThicknessF               = 1.0
     res21@cnLineColor                    = "Black"
     res21@gsnContourZeroLineThicknessF   = 2.
     res21@gsnContourNegLineDashPattern   = 1
     plot_cont1                  = gsn_csm_contour(wks1,vadif,res21)
    ;plot_cont1                  = gsn_csm_pres_hgt(wks1,vadif,res21)
    ;overlay(plot(iv+nvars),plot_cont1)
     delete([/res21/])
    end if 

   delete([/res/])

   delete([/vav0,var0,var1,vadif/])

;---Draw both plots in a panel
 pnres                         = True
 pnres@gsnFrame                = True
 pnres@gsnMaximize             = False
 pnres@gsnPanelLabelBar        = True
 pnres@gsnPanelMainString      = " "; "Surface Temperature ("+ts@units+")"
 pnres@gsnPanelLabelBar        = False                ; add common colorbar
 pnres@gsnPanelMainFontColor   = "Black"
 pnres@gsnPanelMainFontHeightF = 0.030
;pnres@lbLabelFontHeightF      = 0.018               ; make labels smaller
;pnres@pmLabelBarParallelPosF  = 0.026
;pnres@pmLabelBarOrthogonalPosF= 0.015
;pnres@pmLabelBarHeightF       = 0.1
;pnres@pmLabelBarWidthF        = 0.57
;pnres@gsnPanelBottom          = 0.05
 pnres@gsnPanelYWhiteSpacePercent = 10
 pnres@gsnPanelXWhiteSpacePercent = 5
;pnres@gsnPanelXF = new(nvars*3,float)
 pnres@gsnPanelXF                = (/0.10,0.30,0.50,0.70, \
                                     0.10,0.30,0.50,0.70/)
 gsn_panel(wks1,(/plot(iv+nvars)/),(/1,1/),pnres)

end do

end do ;end of variable loop

end do 
end
