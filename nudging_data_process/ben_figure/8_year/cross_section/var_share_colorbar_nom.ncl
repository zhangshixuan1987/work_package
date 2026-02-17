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

  colorbar_without_white := False

;  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
;  ocname      := "prcp_3"  ; color bar name for cntrol model plot
;  orgcolors   := (/3,4,5,7,10,11,12,13,14,17,18/) ;fspan(10,90,81)
; ;olorgroup2  := (/15,26,30,40,45,50,70,80,90,102,2/)
;  colorgroup1 := (/15,26,40,45,50,70,80,102,2/)

  gsn_define_colormap(wks1,"CBR_coldhot")
  colorgroup1 := (/3,4,5,6,0,8,9,10,11/)

 if (Varname .eq. "U") then
  mskthres     = 0.001  ;threshold used to mask out the regions with T < 0.001
  orglev      := (/-30,-20,-10,-5,0,5,10,20,30/)
 ;diffgrp1    := (/-3.0,-2.0,-1.0,-0.5,0.5,1.0,2.0,3.0/)
  diffgrp1    := (/-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6/)
  ndiff       := dimsizes(diffgrp1)
 end if


 if (Varname .eq. "T") then
  mskthres     = 0.001  ;threshold used to mask out the regions with T < 0.001
  orglev      := (/210,220,230,240,250,260,270,280,290/)
 ;diffgrp1    := (/-3.0,-2.0,-1.0,-0.5,0.5,1.0,2.0,3.0/)
  diffgrp1    := (/-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6/)*2.0
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "RELHUM") then
  mskthres     = 0.05   ;threshold used to mask out the regions with RH < 0.05
  orglev      := (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  diffgrp1    := (/-20,-10,-5.0,-1,1,5,10,20/)
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "Q") then
  mskthres     = 0.05   ;threshold used to mask out the regions with RH < 0.05
  orglev      := (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
  diffgrp1    := (/-0.4,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.4/)
  ndiff       := dimsizes(diffgrp1)
 end if 
