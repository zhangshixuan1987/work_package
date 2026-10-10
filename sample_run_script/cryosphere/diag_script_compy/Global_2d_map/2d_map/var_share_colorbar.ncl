;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp1),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;regions for the regional map plot (pacific ocean region)
  mpMinLatF    =-90  ; 20  ;-180.
  mpMaxLatF    = 90  ; 60  ;180.
  mpMinLonF    =  0  ;100  ;-180.
  mpMaxLonF    = 360 ;260  ;180.
  mpCenterLonF = 240 ;180  ;150.0

 ;;;subregions used by the analyses;;;;;;;;;;;;;;;;
 ;Low stratocumulus clouds
 ;Northeast Pacific (NEP) 15–35N, 120–140W
 ;Southeast Pacific (SEP) 10–30S, 75–95W
 ;Southeast Atlantic (SEA) 10–30S, 10W–10E
 ;Southeast Indian Ocean (SEI) 20–40S, 90–110E
 ;Stormtrack
 ;North Pacific (NP) 40–60N, 165–185E
 ;Northeast Atlantic (NEA) 10–30N, 25–45W
 ;North Atlantic (NA) 45–65N, 30–50W
 ;High Clouds
 ;Southern Ocean (SO) 45–65S, 130-170E; 180W–180E
  Subreg       = (/ "NEP",  "SEP",  "SEA",  "SEI",  "NP", "NEA",  "NA",   "SO"/)
  minlat       = (/  15.0,  -10.0,  -10.0,  -20.0,  40.0,  10.0,  45.0,   -45.0/)
  maxlat       = (/  35.0,  -30.0,  -30.0,  -40.0,  60.0,  30.0,  65.0,   -65.0/)
  minlon       = (/-140.0,  -95.0,  -10.0,   90.0, 165.0, -45.0, -50.0,   130.0/)
  maxlon       = (/-120.0,  -75.0,   10.0,  110.0, 185.0, -25.0, -30.0,   170.0/)
  nsubreg      = dimsizes(Subreg)

 ;;setup for fill pattern
  ;FillColor           = "Grey75" ;"Grey65";"Grey28";"Grey";"Black"
  ;FillPattern         = 0 ;17 ;6 ;17;11
  ;FillScaleF          = 0.1
  ;FillDotSizeF        = 0.0002    ; increase dot size of contours
  ;FillOpacityF        = 1.0 ;0.15

  FillColor           = "Grey28";"Grey";"Black"
  FillPattern         = 17 ;6 ;17;11
  FillScaleF          = 0.5
  FillDotSizeF        = 0.002    ; increase dot size of contours
  FillOpacityF        = 0.8

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
  orglev := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ;fspan(10,90,81)

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
  orglev := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ;fspan(10,90,81)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/8,58,88,118,138,168,198,248/)
     diffgrp1    := (/-0.12,-0.08,-0.03,0,0.03,0.08,0.12/) 
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-0.12,-0.08,-0.03,-0.01,0.01,0.03,0.08,0.12/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "CLDHGH") then
  minvar := 0
  maxvar := 0.8
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ;fspan(10,90,81)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/8,58,88,118, 138,168,198,248/)
     diffgrp1    := (/-0.07,-0.05,-0.03,0.,0.03,0.05,0.07/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-0.07,-0.05,-0.03,-0.01,0.01,0.03,0.05,0.07/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "CLDTOT") then
  minvar := 0
  maxvar := 0.8
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ;fspan(10,90,81)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "NCV_blu_red") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/8,58,88,118, 138,168,198,248/)
     diffgrp1    := (/-0.5,-0.3,-0.1,0,0.1,0.3,0.5/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/8,58,88,118,0,138,168,198,248/)
     diffgrp1    := (/-0.4,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.4/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "TGCLDLWP") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 100, 120, 140/);

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,0,5,10,20/)
  else
     colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,-2,2,5,10,20/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "TGCLDIWP") then
  minvar := 0
  maxvar := 40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/1, 2, 5, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-8,-4,-1,0,1,4,8/)
  else
     colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     diffgrp1    := (/-8,-4,-2,-1,1,2,4,8/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "LWCF") then 
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-10,0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,0,5,10,20/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/15,26,33,42,0,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,-2,2,5,10,20/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "SWCF") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ;ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  cmap = read_colormap_file("WhBlGrYeRe")
  ocname := cmap(::-1,:)

  orglev := (/-120,-90,-80,-75,-60,-50,-40,-30,-20,-10,-5,0,10,30/)

  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  dcname := "BlueDarkRed18"
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-40,-20,-10,0,10,20,40/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-40,-20,-10,-5,5,10,20,40/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  colorgroup1 =colorgroup1(::-1)
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "FSUTOA") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;cmap = read_colormap_file("WhBlGrYeRe")
  ;ocname := cmap(::-1,:)

  orglev := (/0,20,40,60,80,100,120,140,160,180,200/)

  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  dcname := "BlueDarkRed18"
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-40,-20,-10,0,10,20,40/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-40,-20,-10,-5,5,10,20,40/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  colorgroup1 =colorgroup1(::-1)
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "FLUT") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;cmap = read_colormap_file("WhBlGrYeRe")
  ;ocname := cmap(::-1,:)

  orglev := (/100,120,140,160,180,200,220,240,260,280,300/)

  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  dcname := "BlueDarkRed18"
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,0,5,10,20/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,-2,2,5,10,20/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  colorgroup1 =colorgroup1(::-1)
  ndiff       := dimsizes(diffgrp1)
 end if


 if (Varname .eq. "NETCF") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-120,-90,-80,-75,-60,-50,-40,-30,-20,-10,-5,0,10,30/)

  gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
  dcname := "BlueDarkRed18"
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/15,26,33,42,75,85,102,2/)
     diffgrp1    := (/-40,-20,-10,0,10,20,40/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-40,-20,-10,-5,5,10,20,40/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "PRECT") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.1,0.2,0.5,0.8,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0,15.0/)

  gsn_define_colormap(wks1, "MPL_BrBG") ; choose colormap
  dcname := "BlueDarkRed18"
  if (colorbar_without_white) then ; do not include white in colorbar
     colorgroup1 := (/10,25,40,55, 77,92,107,122/)
     diffgrp1    := (/ -2.0, -1.0, -0.5, 0, 0.5, 1.0, 2.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     colorgroup1 := (/10,25,40,55,0,77,92,107,122/)
     diffgrp1    := (/ -2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "PRECC") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.1,0.2,0.5,0.8,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0,15.0/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_BrBG") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     diffgrp1    := (/ -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6/)
     colorgroup1 := (/10,25,40,55, 77,92,107,122/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     diffgrp1    := (/ -0.6, -0.4, -0.2, -0.1, 0.1,0.2, 0.4, 0.6/)
     colorgroup1 := (/10,25,40,55,0,77,92,107,122/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if 

 if (Varname .eq. "PRECL") then
  minvar := 0
  maxvar := 4
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.1,0.2,0.3,0.4,0.5,0.8,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_BrBG") ; choose colormap
  if (colorbar_without_white) then ; do not include white in colorbar
     diffgrp1    := (/ -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6/)
     colorgroup1 := (/10,25,40,55, 77,92,107,122/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  else
     diffgrp1    := (/ -0.6, -0.4, -0.2, -0.1, 0.1,0.2, 0.4, 0.6/)
     colorgroup1 := (/10,25,40,55,0,77,92,107,122/)
     diffgrp2    := (/-40,-20,-10,-5,5,10,20,40/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "OMEGA500") then
  minvar := -40
  maxvar :=  40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-50,-40,-30,-10,0,10,20,30,40,50/)

  dcname := "BlueDarkRed18"
  gsn_define_colormap(wks1, "MPL_PRGn") ; choose colormap
  colorbar_without_white := True
  if (colorbar_without_white) then ; do not include white in colorbar
    ;colorgroup1 := (/15,26,33,42,75,85,102,2/)
     colorgroup1 := (/12,28,44,56,73,87,103,119/)
     diffgrp1    := (/-4,-2,-1,0,1,2,4/)
     diffgrp2    := (/-50,-20,-10,0,10,20,50/)
  else
    ;colorgroup1 := (/15,26,33,42,50,75,85,102,2/)
     colorgroup1 := (/15,30,40,45,50,75,85,102,2/)
     diffgrp1    := (/-20,-10,-5,-2,2,5,10,20/)
     diffgrp2    := (/-50,-20,-10,-5,5,10,20,50/)
  end if
  ndiff       := dimsizes(diffgrp1)
 end if
