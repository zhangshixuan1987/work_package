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

  l_create_plot_data = False

  ldatalocal = True
  ldatalocal = False

  if (ldatalocal) then
     ModelOutRoot = "./"
  else
     ModelOutRoot = "/pic/projects/E3SM_DA/zhan391/dart_paper_materials/DATA/"
  end if

  ;;;;;;;date strings for the data directory;;;;;;;;
  Casdir     = (/"PDDART_UV_ENS80","PDDART_UVT_ENS80","PDDART_UVTQ_ENS80"/)
  Groups     = (/"EDA-UV", "EDA-UVT","EDA-UVTQ"/)
  ngrps      = dimsizes(Groups)
  LabList    = Groups

  Keys       = (/"Analysis increment", "Ensemble spread"/)
  nkeys      = dimsizes(Keys)

  plot       = new(ngrps*nkeys,graphic)
  PltNums    = (/"(a)", "(b)", "(c)", "(d)", "(e)", "(f)"/)

  FREQ       = "6hr"
  YEAR       = "2011"

;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList    = (/"T"/)
  VstList    = (/"T"/)
  VnmList    = (/"Temperature (T)"/)
  UnitList   = (/"K"/)
  nvars      = dimsizes(VarList)

;;;;;the pressure levels for cross section plot;;;;;;
  P0mb       = 1000.                        ; reference pressure [mb]
  pnew       = asciiread("era5_plev.txt",-1,"float")
  pnew@units = "hPa"
  pnew!0     = "lev"
  pnew&lev   = pnew

 ;;; flags to control the plot style
  plot_contour      = True   ; if True, plot the contour in the figure  

do iv = 0,nvars-1, 1

 Varname = VarList(iv)
 Varunt  = UnitList(iv)
 Varstr  = VstList(iv)
 Vnmstr  = VnmList(iv)

 do igp = 0, ngrps-1,1

  if (l_create_plot_data) then

   ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   input_dir = ModelOutRoot+"/"+Casdir(igp)+"/analysis/"+FREQ+"/"
   File_0    = input_dir + "coord_output_sd.nc"; coordinate information
   File_1    = input_dir + "PS_"+YEAR+"_preassim_mean.nc"
   File_2    = input_dir + "PS_"+YEAR+"_output_mean.nc"
   File_a    = input_dir + Varname+"_"+YEAR+"_preassim_mean.nc" ; first guess ensemble mean 
   File_b    = input_dir + Varname+"_"+YEAR+"_output_mean.nc"   ; analysis ensemble mean
   File_c    = input_dir + Varname+"_"+YEAR+"_preassim_sd.nc"   ; first guess spread  
   File_d    = input_dir + Varname+"_"+YEAR+"_output_sd.nc"     ; analysis spread 

   f0        = addfile(File_0,"r") 
   f1        = addfile(File_1,"r")
   f2        = addfile(File_2,"r")
   fa        = addfile(File_a,"r")
   fb        = addfile(File_b,"r")
   fc        = addfile(File_c,"r")
   fd        = addfile(File_d,"r")

   if(Varname.eq."US")then 
    lat       = f0->slat
   else
    lat       = f0->lat
   end if 

   if(Varname.eq."VS")then                   
    lon       = f0->slon
   else
    lon       = f0->lon
   end if 

   interp     = 2 ; log  
   extrap     = False
   var1       = vinth2p(fa->$Varname$,f0->hyam,f0->hybm,pnew,f1->PS,interp,P0mb,1,extrap)
   var2       = vinth2p(fb->$Varname$,f0->hyam,f0->hybm,pnew,f2->PS,interp,P0mb,1,extrap)
   var3       = vinth2p(fd->$Varname$,f0->hyam,f0->hybm,pnew,f2->PS,interp,P0mb,1,extrap)
   var1!1     = "lev"
   var2!1     = "lev"
   var3!1     = "lev"

   anl_increment  = var1
   anl_increment  = abs(var2 - var1)
   anl_spread     = var3
   delete([/var1,var2,var3/])

   anl_increment&lev@units = "hPa"
   anl_spread&lev@units = "hPa"
 
   ;;calculate the temporal average 
   manl_increment = dim_avg_n_Wrap(anl_increment,0)
   manl_spread    = dim_avg_n_Wrap(anl_spread,0)
   delete([/anl_increment,anl_spread/]) 
 
   if (Varname.eq."Q")then
    manl_increment = manl_increment *1000.0 
    manl_spread    = manl_spread * 1000.0
    manl_increment@units = "g/kg" 
    manl_spread@units = "g/kg"
   end if 
  
   ;;calculate the global mean mentric 
   rad     = 4.*atan(1.)/180.    ;
   cost    = cos(lat*rad)       ; cosine weights
   mean1x  = wgt_areaave(manl_increment,cost,1.0,0)
   mean2x  = wgt_areaave(manl_spread,cost,1.0,0)
   mean1   = sprintf("%5.2f", mean1x)
   mean2   = sprintf("%5.2f", mean2x)
   manl_increment@global_mean = mean1x
   manl_spread@global_mean    = mean2x  
   delete([/mean1x,mean2x,cost,rad/]) 
   
   ;;calculate the zonal mean mentric
   zanl_increment = dim_avg_n_Wrap(manl_increment,2)
   zanl_spread    = dim_avg_n_Wrap(manl_spread,2)
   delete([/manl_increment,manl_spread/])

   printVarSummary(zanl_increment)
   printVarSummary(zanl_spread)
   fno  = Groups(igp)+"_"+Varstr+".nc"
   system("rm -f " + fno)

   print(" ")
   print("output : " + fno)
   print(" ")

   vnm1 = Varstr+"_inc"
   vnm2 = Varstr+"_sd"

   flo = addfile(fno,"c")
   flo->lat    = lat
   flo->lon    = lon
   flo->lev    = pnew
   flo->$vnm1$ = zanl_increment
   flo->$vnm2$ = zanl_spread
   delete([/lat,lon,zanl_increment,zanl_spread,fno/])

   ;;;read the data and do plot;;;
   lat  = flo->lat
   lon  = flo->lon
   lev  = flo->lev
   var1 = flo->$vnm1$
   var2 = flo->$vnm2$

 else
   
   ;;;read the data and do plot;;;
   vnm1 = Varstr+"_inc"
   vnm2 = Varstr+"_sd"
   fno  = Groups(igp)+"_"+Varstr+".nc"
   flo  = addfile(fno,"r")
   lat  = flo->lat
   lon  = flo->lon
   lev  = flo->lev
   var1 = flo->$vnm1$
   var2 = flo->$vnm2$

 end if 

  print(min(var1)+ " " +max(var1))

;----------------------------------------------------------------------
; Plotting section
;----------------------------------------------------------------------

  if(igp.eq.0)then 

   ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
   OUTDir   = "./"
   plotname = OUTDir+"Fig03"

   wtype = "pdf"
  ;wtype@wkWidth = 2500
  ;wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length

   print("")
   print("plot " + plotname)
   print("")

   wks1 = gsn_open_wks(wtype,plotname)

  end if 
  
 ;;;set up the plot information;;;;
  load "./var_share_colorbar_preshgt.ncl"

  Figstr1  = PltNums(igp)+" "+ Keys(0) + " of "+Varstr + " ("+Varunt+")"
  Figstr2  = PltNums(igp+ngrps)+" "+ Keys(1) + " of "+Varstr + " ("+Varunt+")"

;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@gsnAddCyclic         = False
  res@vpHeightF            = 0.45
  res@vpWidthF             = 0.8

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"

  res@tmBorderThicknessF    = 1.0
  res@tmXBMajorThicknessF   = 1.0
  res@tmYLMajorThicknessF   = 1.0
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
  res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0
  res@tiMainFontThicknessF   = 1.0

  FontHeight = 0.024
  res@gsnStringFontHeightF = FontHeight*1.2
  res@tmXBLabelFontHeightF = FontHeight
  res@tmYLLabelFontHeightF = FontHeight
  res@tmYRLabelFontHeightF = FontHeight
  res@tiMainFontHeightF    = FontHeight
  res@tiXAxisFontHeightF   = FontHeight
  res@tiYAxisFontHeightF   = FontHeight
  res@lbLabelFontHeightF   = FontHeight
  res@lbTitleFontHeightF   = FontHeight 
  res@cnLineLabelFontHeightF = FontHeight

  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
  res@lbLabelStride           = 1                 ; skip every other label
  res@lbPerimOn               = False             ; default has box
  res@lbOrientation           = "vertical"        ; vertical label bar
  res@pmLabelBarHeightF       = 0.45
  res@pmLabelBarWidthF        = 0.12
; res@pmLabelBarParallelPosF  = 0.026
; res@pmLabelBarOrthogonalPosF= 0.15
  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle

  res@gsnLeftString          = "" ; leftstr ;var_name0    ; long_name is too long!
  res@gsnRightString         = "" ;" ("+Varunt+")";""
  res@tiMainString           = ""

  res@tmXTOn                 = False
  res@tmYROn                 = False
  res@tmYRMode               = "Automatic"
  res@tmYLOn                 = True  
  res@tmYLLabelsOn           = True

  res@trYMinF                = 10.0              ; less white space
 ;res@trYMaxF                = 1000.0
  res@trYReverse             = True                     ; reverse Y-axis
  res@tiYAxisString          = "Pressure (hPa)"
  res@tmYLMode               = "Explicit"                ; explicit labels
  res@tmYLValues             = (/1000, 850, 700, 500, 300, 200, 100, 10/)
                               ;(/ 950, 700, 500, 300,200,100,50,30,10,5,1/)
  res@tmYLLabels             = ""+res@tmYLValues         ; make strings

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  res2=res
  res2@cnLevelSelectionMode           = "ExplicitLevels"  
  res2@cnFillPalette                  = xocname
 ;res2@cnFillColors                   = colorgroup1
  res2@cnLevels                      := diffgrp1
  res2@cnLinesOn                      = False 
  res2@cnLineLabelsOn                 = False            ; turn off contour labels
  res2@cnInfoLabelOn                  = False
  res2@tiMainString                   = "";"PD300 - CNTL"
  res2@gsnLeftString                  = Figstr1
  res2@gsnLeftStringOrthogonalPosF    = 0.02
 ;res2@gsnLeftStringParallelPosF      = 0.2 
  res2@gsnRightString                 = Groups(igp)
  res2@gsnRightStringOrthogonalPosF   = 0.02
  res2@lbTitleString                  = "" ;"("+ Varunt +")"
  res2@lbLabelStrings                 = sprintf("%5.1f",diffgrp1)
  res2@tmYLLabelsOn          = True
  res2@gsnLeftStringOrthogonalPosF    = 0.02
  res2@gsnRightStringOrthogonalPosF   = 0.04
  plot(igp*nkeys)                     = gsn_csm_contour(wks1,var1,res2)
 ;plot_rgd1                           = gsn_csm_pres_hgt(wks1,var1,res2)

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
    plot_cont1                  = gsn_csm_contour(wks1,var1,res21)
   ;plot_cont1                  = gsn_csm_pres_hgt(wks1,var1,res21)
    overlay(plot(igp*nkeys),plot_cont1)
    delete([/res21/])
  end if 

  res3=res
  res3@cnLevelSelectionMode        = "ExplicitLevels"
  res3@cnFillPalette               = docname
 ;res3@cnFillColors                = colorgroup2
  res3@cnLevels                   := diffgrp2
  res3@tiMainString                = "";"PD300 - CNTL"
  res3@gsnLeftString               = Figstr2
 ;res3@gsnLeftStringOrthogonalPosF = 0.04
 ;res3@gsnLeftStringParallelPosF   = 0.2
  res3@gsnRightString              = Groups(igp)
  res3@lbTitleString               = ""
  res3@lbLabelStrings              = sprintf("%5.1f",diffgrp2)

  plot(igp*nkeys+1)                = gsn_csm_contour(wks1,var2,res3)
 ;plot_rgd2                        = gsn_csm_pres_hgt(wks1,var2,res3)

  if(plot_contour)
   res31 = res3
   res31@gsnLeftString                  = "" 
   res31@gsnRightString                 = ""
   res31@gsnCenterString                = ""
   res31@cnFillOn                       = False
   res31@cnLinesOn                      = True
   res31@cnLineLabelsOn                 = False            ; turn off contour labels
   res31@cnInfoLabelOn                  = False
   res31@cnLineLabelDensityF            = 2.0
   res31@cnLineLabelFontHeightF         = 0.018
   res31@cnLineLabelBackgroundColor     = -1
   res31@gsnContourLineThicknessesScale = 1.0
   res31@cnLineThicknessF               = 1.0
   res31@cnLineColor                    = "Black"
   res31@gsnContourZeroLineThicknessF   = 2.
   res31@gsnContourNegLineDashPattern   = 1
   plot_cont1                  = gsn_csm_contour(wks1,var2,res31)
  ;plot_cont1                  = gsn_csm_pres_hgt(wks1,var2,res31)
   overlay(plot(igp*nkeys+1),plot_cont1)
   delete([/res31/])
  end if
  
  delete([/lat,lon,lev,var1,var2/])
 end do

;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = False
    pnres@gsnPanelLabelBar        = False
    pnres@gsnPanelMainString      = "" 
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.024
    pnres@lbLabelFontHeightF      = .01            ; default is HUGE
    pnres@lbPerimOn               = False             ; default has box
    pnres@lbOrientation           = "vertical"        ; vertical label bar
    pnres@lbTitleString           = "("+Varunt+")" ; "~F33~D~N~~F~"+ "RH"+" ("+Varunt+")"
    pnres@lbTitlePosition         = "Right"
    pnres@lbTitleDirection        = "Across"
    pnres@lbTitleAngleF           = 90
    pnres@lbTitleFontHeightF      = 0.01
    pnres@pmLabelBarHeightF       = 0.22
    pnres@pmLabelBarWidthF        = 0.06
    pnres@pmLabelBarParallelPosF  = -0.005
   ;pnres@pmLabelBarOrthogonalPosF= -0.15
    pnres@gsnPanelYWhiteSpacePercent = 2
    pnres@gsnPanelXWhiteSpacePercent = 5
   ;pnres@gsnPanelXF = (/0.15,0.39,0.63/)    ; Adjust rightmost plots
    gsn_panel(wks1,(/plot/),(/3,2/),pnres)
    delete([pnres])

end do ; end of var loop

end
