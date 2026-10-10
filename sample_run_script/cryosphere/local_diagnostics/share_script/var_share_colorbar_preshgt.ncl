;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp1),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;setup for fill pattern
  FillColor           = "Grey28";"Grey";"Black"
  FillPattern         = 0 ;17 ;6 ;17;11
  FillScaleF          = 0.1
  FillDotSizeF        = 0.0002    ; increase dot size of contours
  FillOpacityF        = 0.2

 ; FillColor           = "White" ;"Lightgrey";"Grey28";"Grey";"Black"
 ; FillPattern         = 0   ;11 ;17 ;6 ;17;11
 ; FillScaleF          = 0.5
 ; FillDotSizeF        = 0.002    ; increase dot size of contours
 ; FillOpacityF        = 0.8

  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot

  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap

  colorbar_without_white := True
  if (plot_ttest_hatrv .or. plot_ttest_hatch) then
     colorbar_without_white := False
  end if

 if (Varname .eq. "T") then
  mskthres = 0.001  ;threshold used to mask out the regions with T < 0.001

  ;ocname := "amwg"  ; color bar name for cntrol model plot
  ;orglev := ispan(180,300,10)
  ocname    := "prcp_3"  ; color bar name for cntrol model plot
  ocname    := "BlAqGrYeOrReVi200"
  orgcolors := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  orglev    := (/200,210,220,230,240,250,260,270,280,300/)

  dcname   := "BlueDarkRed18"
  diffgrp1 := (/-5.0,-3.0, -2.0,-1.0,-0.5, 0.5, 1.0, 2.0, 3.0, 5.0/)
 end if 

 if (Varname .eq. "Q") then

  mskthres = 0.002  ;threshold used to mask out the regions with Q < 0.01

  ;ocname := "amwg"  ; color bar name for cntrol model plot
  ;orglev := (/0.01,0.05,0.1,0.2,0.5,1,2,5,8,10,12,15/)
  ocname    := "prcp_3"  ; color bar name for cntrol model plot
  orgcolors := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  ocname    := "BlAqGrYeOrReVi200"
  orglev    := (/0.01,0.05,0.1,0.2,0.5,1,2,5,8,10/)

  dcname   := "BlueDarkRed18"
  diffgrp1 := (/-2.0,-1.0, -0.5,-0.2,-0.1, 0.1, 0.2, 0.5, 1.0, 2.0/)

 end if 


 if (Varname .eq. "U") then
  mskthres = 0.001  ;threshold used to mask out the regions with T < 0.001

  ;ocname := "amwg"  ; color bar name for cntrol model plot
  ;orglev := ispan(180,300,10)
  ocname    := "cmocean_curl" ;"GreenMagenta16"  ; color bar name for cntrol model plot
  orgcolors := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  orglev    := (/-40,-30,-25.0,-20.0,-15,-10,-5,-2,2,5,10,15,20,25,30,40/)

  dcname   := "BlueDarkRed18"
  diffgrp1 := (/-5.0,-3.0, -2.0,-1.0,-0.5, 0.5, 1.0, 2.0, 3.0, 5.0/)

 end if


