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

  colorbar_without_white := False

 if (Varname .eq. "U") then
  mskthres     = 0.001  ;threshold used to mask out the regions with T < 0.001
  ocname      := "prcp_3"  ; color bar name for cntrol model plot
  orgcolors   := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  orglev      := (/-30,-20,-10,-5,0,5,10,20,30/)
  colorgroup1 := (/15,26,30,40,45,50,70,80,90,102,2/)
  diffgrp1    := (/-4.0,-2.0,-1.0,-0.5,-0.2,0.2,0.5,1.0,2.0,3.0/)
  ndiff       := dimsizes(diffgrp1)
 end if


 if (Varname .eq. "T") then
  mskthres     = 0.001  ;threshold used to mask out the regions with T < 0.001
  ocname      := "prcp_3"  ; color bar name for cntrol model plot
  orgcolors   := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  orglev      := (/210,220,230,240,250,260,270,280,290/)
  colorgroup1 := (/15,26,30,40,45,50,70,80,90,102,2/)
  diffgrp1    := (/-3.0,-2.0,-1.0,-0.5,-0.2,0.2,0.5,1.0,2.0,3.0/)
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "CLOUD") then
  mskthres     = 0.02 ;threshold used to mask out the regions with CLOUD < 0.01
  ocname      := "prcp_3"  ; color bar name for cntrol model plot
  orgcolors   := (/3,4,5,7,10,12,13,14,17,18/) ;fspan(10,90,81)
  orglev      := (/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4/)
  colorgroup1 := (/15,26,30,40,45,50,70,80,90,102,2/)
  diffgrp1    := (/-0.2,-0.15,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.15,0.2/)
  ndiff       := dimsizes(diffgrp1)
 end if


 if (Varname .eq. "RELHUM") then
  mskthres     = 0.05   ;threshold used to mask out the regions with RH < 0.05
  ocname      := "prcp_3"  ; color bar name for cntrol model plot
  orgcolors   := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  orglev      := (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  colorgroup1 := (/15,26,30,40,45,50,70,80,90,102,2/)
  diffgrp1    := (/-0.3,-0.2,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.2,0.3/)*100.0
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "Q") then
  mskthres     = 0.05   ;threshold used to mask out the regions with RH < 0.05
  ocname      := "prcp_3"  ; color bar name for cntrol model plot
  orgcolors   := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
  orglev      := (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  colorgroup1 := (/15,26,30,40,45,50,70,80,90,102,2/)
  diffgrp1    := (/-0.3,-0.2,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.2,0.3/)*10.0
  ndiff       := dimsizes(diffgrp1)
 end if 
