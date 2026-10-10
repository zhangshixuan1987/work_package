;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp1),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;;setup for fill pattern
  ;FillColor           = "Grey75" ;"Grey65";"Grey28";"Grey";"Black"
  ;FillPattern         = 0 ;17 ;6 ;17;11
  ;FillScaleF          = 0.1
  ;FillDotSizeF        = 0.0002    ; increase dot size of contours
  ;FillOpacityF        = 1.0 ;0.15

  FillColor           = "Black";"Green3" ;"Magenta"; "Grey28";"Green3";"Magenta";"Grey28";"Grey";"Black"
  FillPattern         = 17 ;6 ;17;11
  FillScaleF          = 0.3 ;0.5
  FillDotSizeF        = 0.003    ; increase dot size of contours
  FillOpacityF        = 0.8

  FillColor1          = "Black";"Green3";"Magenta";"Grey28";"Grey";"Black"
  FillPattern1        = 17 ;6 ;17;11
  FillScaleF1         = 0.5
  FillDotSizeF1       = 0.002    ; increase dot size of contours
  FillOpacityF1       = 0.8

  if (plot_ttest_mskout ) then
     colorbar_without_white := True
  else if (plot_ttest_hatrv .or. plot_ttest_hatch) then
     colorbar_without_white := True
  else
     colorbar_without_white := False
  end if
  end if

 if (Varname .eq. "PSL") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-120.0,-60.0,-40.0,-20.0,-10.0,10.0,20.0,40.0,60.0,120.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "PSL".or.Varname .eq. "PS") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-120.0,-60.0,-40.0,-20.0,-10.0,10.0,20.0,40.0,60.0,120.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "SHFLX") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/10.0
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "LHFLX") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/10.0
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "TREFHT".or.Varname.eq."TS") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-2.0,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,2.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "U10") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-1.0,-0.8,-0.5,-0.2,-0.1,0.1,0.2,0.5,0.8,1.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "PRECT") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-0.5,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.5/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "CLDTOT") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-1.0,-0.8,-0.5,-0.2,-0.1,0.1,0.2,0.5,0.8,1.0/)/10.0
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "TAUX".or. Varname .eq. "TAUY".or.Varname .eq. "TAUXY") then
  dcname   := read_colormap_file("CBR_coldhot")
  ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-2.0,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,2.0/)/10.0
  ndiff    := dimsizes(diffgrp1)
 end if


