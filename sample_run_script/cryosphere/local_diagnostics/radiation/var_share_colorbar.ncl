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

  FillColor           = "Green3" ;"Magenta"; "Grey28";"Green3";"Magenta";"Grey28";"Grey";"Black"
  FillPattern         = 17 ;6 ;17;11
  FillScaleF          = 0.3 ;0.5
  FillDotSizeF        = 0.003    ; increase dot size of contours
  FillOpacityF        = 0.8

  FillColor1          = "Green3";"Magenta";"Grey28";"Grey";"Black"
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

 if (Varname .eq. "TH7001000") then
  minvar := 8
  maxvar := 24
  maxLev := 17
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := fspan(2,22,11)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/8,58,88,118, 138,168,198,248/)
     diffgrp1    := (/-0.8,-0.4,-0.2,0,0.2,0.4,0.8/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-0.8,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.8/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "CLDLOW") then
  minvar := 0
  maxvar := 0.8
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/) ;fspan(10,90,81)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     ;colorgroup1 := (/8,58,88,118,138,168,198,248/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-0.16,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.16/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-0.12,-0.08,-0.03,-0.01,0.01,0.03,0.08,0.12/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.01,11)

  dcname   := "CBR_coldhot"
  diffgrp1 := (/-0.5,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.5/)
  difflabs := sprintf("%3.2f",diffgrp1)
  ndiff       := dimsizes(diffgrp1)

 end if


 if (Varname .eq. "CLDHGH") then
  minvar := 0
  maxvar := 0.8
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.01,0.05,0.08,0.12,0.16,0.20,0.24,0.30,0.4,0.5/) ;fspan(10,90,81)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     ;colorgroup1 := (/8,58,88,118,138,168,198,248/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-0.16,-0.08,-0.04,-0.02,0,0.02,0.04,0.08,0.16/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-0.12,-0.08,-0.03,-0.01,0.01,0.03,0.08,0.12/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.008,11)

  dcname   := "CBR_coldhot"
  diffgrp1 := (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  difflabs := sprintf("%3.2f",diffgrp1)
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "CLDTOT") then
  ocname   := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev   := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ;fspan(10,90,81)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.01,11)

  dcname   := "CBR_coldhot"
  diffgrp1 := (/-0.20,-0.16,-0.12,-0.08,-0.04,0.04,0.08,0.12,0.16,0.20/) 
  ndiff1   := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "TGCLDLWP") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 140/);

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     ;colorgroup1 := (/8,58,88,118,138,168,198,248/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-20,-10,-5,-1,0,1,5,10,20/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-40,-20,-10, -5,5,10,20,40/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,3,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "TGCLDIWP") then
  minvar := 0
  maxvar := 40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0, 2, 5, 8, 10, 15, 20, 30, 40, 50,60/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-8,-4,-2,-1,0,1,2,4,8/)
  else
     colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,0,11,13,15,17,19/)
     diffgrp1    := (/-8,-4,-2,-1,1,2,4,8/)
  end if

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-40.0,-30.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,30.0,40.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "NUMLIQ") then
  minvar := 0
  maxvar := 40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0, 5, 10, 20, 30, 40, 50, 60, 80, 100, 120/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-80,-40,-10,-5,0,5,10,40,80/)
  else
     colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,0,11,13,15,17,19/)
     diffgrp1    := (/-80,-40,-20,-10,10,20,40,80/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "NUMICE") then
  minvar := 0
  maxvar := 40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/10, 20, 50, 80, 100, 120, 150, 200, 250, 300, 350/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-80,-40,-20,-10,0,10,20,40,80/)
  else
     colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,0,11,13,15,17,19/)
     diffgrp1    := (/-80,-40,-20,-10,10,20,40,80/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "CCN3") then
  minvar := 0
  maxvar := 40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-50,-30,-10,-5,0,5,10,30,50/)
  else
     colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/2,4,6,8,10,0,11,13,15,17,19/)
     diffgrp1    := (/-30,-20,-10,-5,-1,1,5,10,20,30/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "LWCF") then 
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-8, -4, 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 50/)

  xcname   :="precip3_16lev";"WhiteBlueGreenYellowRed" ;"MPL_YlGn"; "precip_11lev"
  diffgrp0 := fspan(0,0.4,11)
 
  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-8.0,-4.0,-2.0,-1.0,-0.5,0.5,1.0,2.0,4.0,8.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "SWCF") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;cmap = read_colormap_file("WhBlGrYeRe")
  ;ocname := cmap(::-1,:)
  orglev := (/-80,-75,-60,-50,-40,-30,-20,-10,-5,0,5/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)
 
  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-15.0,-10.0,-5.0,-1.0,1.0,5.0,10.0,15.0,20.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "NETCF") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;cmap = read_colormap_file("WhBlGrYeRe")
  ;ocname := cmap(::-1,:)
  orglev := (/-60,-50,-40,-30,-20,-10,-5,-1,0,1,5,10,20,30/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     ;colorgroup1 := (/8,58,88,118,138,168,198,248/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-8,-4,-2, -1,1,2,4,8/)
     diffgrp1    := (/-12,-8,-4,-2,0,2,4,8,12/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-8,-4,-2, -1,1,2,4,8/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-16.0,-12.0,-8.0,-5.0,-1.0,1.0,5.0,8.0,12.0,16.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "PRECT") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.3,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-1.0,-0.8,-0.5,-0.2,-0.1,0.1,0.2,0.5,0.8,1.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "PRECC") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.0,0.1,0.2,0.5,0.8,1.0,2.0,3.0,4.0,6.0,8.0/)
  orglev := (/0.01,0.05,0.08,0.1,0.2,0.5,0.8,1.0,1.5,2.0/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-1.0,-0.8,-0.5,-0.2,-0.1,0.1,0.2,0.5,0.8,1.0/)
  ndiff    := dimsizes(diffgrp1)

 end if 

 if (Varname .eq. "PRECL") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.1,0.2,0.3,0.4,0.5,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0/)
  orglev := (/0.01,0.1,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,8.0/)
  orglev := (/0.1,0.2,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.2,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-1.0,-0.8,-0.5,-0.2,-0.1,0.1,0.2,0.5,0.8,1.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "OMEGA500") then
  minvar := -40
  maxvar :=  40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-80,-60,-40,-20,-10,0,10,20,40,60,80/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_PRGn") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
    ;colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/12,28,44,56,73,87,103,119/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-12,-8,-4,-1,0,1,4,8,12/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
    ;colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,-2,2,5,10,20/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "U10") then
  minvar :=  0
  maxvar :=  30
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0,1,2,3,4,5,6,7,8,10,12/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_PRGn") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
    ;colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/12,28,44,56,73,87,103,119/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-2.0,-1.0,-0.5,-0.1,0,0.1,0.5,1.0,2.0/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
    ;colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-10,-5,-2,-1,1,2,5,10/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "SHFLX") then
  minvar :=  0
  maxvar :=  30
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0,5,10,20,30,40,50,60,70,80,100/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_PRGn") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
    ;colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/12,28,44,56,73,87,103,119/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-8,-5,-3,-1,0,1,3,5,8/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
    ;colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-10,-5,-2,-1,1,2,5,10/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "LHFLX") then
  minvar :=  0
  maxvar := 120
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/10,30,50,70,90,110,120,130,140,150,180/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_PRGn") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
    ;colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/12,28,44,56,73,87,103,119/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/-12,-8,-4,-1,0,1,4,8,12/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
    ;colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,-2,2,5,10,20/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "FLUT" .or. Varname .eq. "FLNTOA") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/120,130,140,150,160,170,180,190,200,210,220/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.4,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-12.0,-8.0,-4.0,-2.0,-1.0,1.0,2.0,4.0,8.0,12.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "FLNS") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/30,32,34,36,38,40,42,44,46,48,50/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-16.0,-12.0,-8.0,-4.0,4.0,8.0,12.0,16.0,20.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname.eq."FSNT" .or. Varname .eq. "FSNTOA") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/30,40,50,60,70,80,90,100,110,120,150/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,2,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-16.0,-12.0,-8.0,-4.0,4.0,8.0,12.0,16.0,20.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "FSUTOA") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/60,70,80,90,100,110,120,130,140,150,180/)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-15.0,-10.0,-5.0,-1.0,1.0,5.0,10.0,15.0,20.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "UVS") then

  if(refpres.eq.1000.)then
   refvec1 = 10. 
   refvec2 = 1.
   difflevs := (/-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0/)
  else if (refpres.eq.850.)
   refvec1 = 10.
   refvec2 = 1.
   difflevs := (/-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0/)
  else if (refpres.eq.700.)
   refvec1 = 20.
   refvec2 = 1.
   difflevs := (/-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0/)
  else if (refpres.eq.500.)
   refvec1 = 25.
   refvec2 = 1.5 
   difflevs := (/-2.0,-1.0,-0.5,-0.25,0,0.25,0.5,1.0,2.0/)
  else if (refpres.eq.200.)
   refvec1 = 30.
   refvec2 = 2.
   difflevs := (/-2.0,-1.0,-0.5,-0.25,0,0.25,0.5,1.0,2.0/)
  else
   refvec1 = 10.
   refvec2 = 1.
  end if 
  end if 
  end if 
  end if 
  end if 

  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := ispan(190,240,5) ;(/-10,0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90/)
  ocname  := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot

  dcname   := "BlueDarkRed18"
  difflevs := (/-2.0,-1.0,-0.5,-0.2,0,0.2,0.5,1.0,2.0/)

 end if

 if (Varname .eq. "AODVIS" .or. Varname .eq. "AODALL")then

  minvar := 0
  maxvar := 1
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname  := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := fspan(0,0.2,11) ;(/0.1,0.2,0.5,0.8,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0,15.0/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_BrBG") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/10,25,40,55, 77,92,107,122/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/ -1.0, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0/)/10.0
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     colorgroup1 := (/10,25,40,55,0,77,92,107,122/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/ -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0/)/10.0
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  
  ndiff       := dimsizes(diffgrp1)
  
 end if

 if (Varname .eq. "BURDEN2") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname  := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := fspan(0,4,11) ;(/0.1,0.2,0.5,0.8,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0,15.0/)
  orglev := (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,2.0,4.0/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_BrBG") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/10,25,40,55, 77,92,107,122/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/ -1.0, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     colorgroup1 := (/10,25,40,55,0,77,92,107,122/)
     diffgrp1    := (/ -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if


 if (Varname .eq. "BURDEN1")then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname  := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,2.0,4.0/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_BrBG") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/10,25,40,55, 77,92,107,122/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/ -1.0, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/ -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "BURDEN3") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname  := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.01,0.1,0.2,0.4,0.6,0.8,1.0,2.0,4.0,8.0,12.0/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  gsn_define_colormap(wks1, "BlueDarkRed18") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     ;colorgroup1 := (/8,58,88,118,138,168,198,248/)
     colorgroup1 := (/2,4,6,8,10,11,13,15,17,19/)
     diffgrp1    := (/ -1.0, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/ -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if
