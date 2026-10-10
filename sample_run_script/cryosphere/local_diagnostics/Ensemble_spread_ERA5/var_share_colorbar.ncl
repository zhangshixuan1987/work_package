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

 if (Varname .eq. "CLDLOW") then
  minvar := 0
  maxvar := 0.8
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/) ;fspan(10,90,81)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.2,11)

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
  orglev := (/0.01,0.05,0.08,0.12,0.16,0.24,0.32,0.36,0.40,0.44/) ;fspan(10,90,81)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := "CBR_coldhot"
  diffgrp1 := (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  difflabs := sprintf("%3.2f",diffgrp1)
  ndiff       := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "CLDTOT") then
  ocname   := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev   := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/) ;fspan(10,90,81)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.2,11)

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
  orglev := (/5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80/);

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,3,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "PSL") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/965,970,975,985,990,995,1000,1005,1010,1015,1020/);

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/5
  ndiff    := dimsizes(diffgrp1)
 end if


 if (Varname .eq. "PS") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/945,955,965,975,985,995,1005,1015,1025/);

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/5
  ndiff    := dimsizes(diffgrp1)
 end if



 if (Varname .eq. "TGCLDIWP") then
  minvar := 0
  maxvar := 40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0, 2, 5, 8, 10, 15, 20, 30, 40, 50,60/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-40.0,-30.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,30.0,40.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "LWCF") then 
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0, 4, 8, 12, 16, 20, 24, 26, 28, 30, 32/)

  xcname   :="precip3_16lev";"WhiteBlueGreenYellowRed" ;"MPL_YlGn"; "precip_11lev"
  diffgrp0 := fspan(0,0.4,11)
 
  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-9.0,-7.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,7.0,9.0/)
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
  orglev := (/-80,-70,-65,-60,-55,-50,-40,-30,-20,-10,-5,0/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1.2,11)
 
  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
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
  orglev := (/-60,-50,-40,-30,-25,-20,-15,-10,-5,0,5/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-16.0,-12.0,-8.0,-5.0,-1.0,1.0,5.0,8.0,12.0,16.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "RESTOM") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;cmap = read_colormap_file("WhBlGrYeRe")
  ;ocname := cmap(::-1,:)
  ;orglev := (/-60,-50,-40,-30,-25,-20,-15,-10,-5,0,5/)
  orglev := (/-120,-80,-40,-10,0,10,20,40,80,120/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-50.0,-20.0,-10.0,-5.0,-2.0,2.0,5.0,10.0,20.0,50.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "QNET") then
  minvar := -100
  maxvar := 20
  maxLev := 13
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;cmap = read_colormap_file("WhBlGrYeRe")
  ;ocname := cmap(::-1,:)
  orglev := (/-60,-40,-20,-10,-5, 0, 5, 10, 20, 40,60/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-50.0,-20.0,-10.0,-5.0,-2.0,2.0,5.0,10.0,20.0,50.0/)
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
  diffgrp0 := fspan(0,0.1,11)

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
  orglev := (/0.01,0.1,0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.4/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

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

 if (Varname .eq. "FLUT" .or. Varname .eq. "FLNTOA") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/140,160,180,200,210,220,230,240,250,260/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.4,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-10.0,-8.0,-6.0,-4.0,-2.0,2.0,4.0,6.0,8.0,10.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "FLNS") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/36,38,40,42,44,46,48,52,56,60,64,68/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.4,11)

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
  orglev := (/40,60,80,100,120,140,180,220,260,300,340/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-16.0,-12.0,-8.0,-4.0,4.0,8.0,12.0,16.0,20.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "FSNTOAC") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/40,60,80,100,140,180,220,260,280,300,320/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-16.0,-12.0,-8.0,-4.0,4.0,8.0,12.0,16.0,20.0/)/2
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname.eq."FSDS") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/40,60,80,100,120,140,160,180,200,220,240/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-16.0,-12.0,-8.0,-4.0,4.0,8.0,12.0,16.0,20.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname.eq."FLDS") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/40,80,120,160,200,240,280,320,360,400/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

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
  orglev := (/40,50,60,70,80,85,90,95,100,110,120/)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-15.0,-10.0,-5.0,-1.0,1.0,5.0,10.0,15.0,20.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "FSUTOAC") then
  minvar := 0
  maxvar := 40
  maxLev := 11
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/20,30,40,50,60,70,80,90,100,110,120/)

  dcname   := read_colormap_file("CBR_coldhot")
 ;dcname   := dcname(::-1,:)
  diffgrp1 := (/-16.0,-12.0,-8.0,-4.0,-2.0,2.0,4.0,8.0,12.0,16.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "TREFHT") then
  minvar := -40
  maxvar :=  40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/228,236,244,252,260,268,272,276,280,290/)

  xcname   := "MPL_YlGn"; "precip_11lev"
  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-2.0,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,2.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname.eq. "TS") then
  minvar := -40
  maxvar :=  40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/220,228,236,242,250,258, 266, 272, 280,288/)

  xcname   := "MPL_YlGn"; "precip_11lev"
  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-2.0,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,2.0/)
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "OMEGA500") then
  minvar := -40
  maxvar :=  40
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-80,-60,-40,-20,-10,0,10,20,40,60,80/)

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,3,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/5.0
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "U10") then
  minvar :=  0
  maxvar :=  30
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0,1,2,3,4,5,6,7,8,10,12/)

  xcname   := "MPL_YlGn"; "precip_11lev"
  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "SHFLX") then
  minvar :=  0
  maxvar :=  30
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0,5,10,20,30,40,50,60,70,80,100/)

  xcname   := "precip3_16lev"; "precip_11lev"
  ;xcname   := "MPL_YlGn"; "precip_11lev"
  diffgrp0 := fspan(0,3,11)


  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/5.0
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "LHFLX") then
  minvar :=  0
  maxvar := 120
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-15.0,0.0,15.0,30.0,45.0,60.0,75.0,90.0,105.0,120.0/)

  xcname   := "MPL_YlGn"; "precip_11lev"
  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,3,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-60.0,-40.0,-20.0,-10.0,-5.0,5.0,10.0,20.0,40.0,60.0/)/5.0
  ndiff    := dimsizes(diffgrp1)

 end if

 if (Varname .eq. "TMQ") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/1, 2, 4, 6, 8 ,10 , 12, 16, 20, 30, 40/);

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,0.1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-5.0,-4.0,-3.0,-2.0,-1.0,1.0,2.0,3.0,4.0,5.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "SHFLX") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/0, 1, 2, 5, 8, 10, 20, 30, 40, 50, 60/);

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,1,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-12.0,-8.0,-4.0,-2.0,-1.0,1.0,2.0,4.0,8.0,12.0/)
  ndiff    := dimsizes(diffgrp1)
 end if

 if (Varname .eq. "LHFLX") then
  minvar := 0
  maxvar := 80
  maxLev := 9
  mnmxint:= nice_mnmxintvl(minvar,maxvar,maxLev, False)
  ocname := "WhBlGrYeRe" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  ;ocname := "NCV_banded" ;"amwg"  ; color bar name for cntrol model plot
  orglev := (/-5, 0, 5, 10, 20, 40, 60, 80, 100, 120/);

  xcname   := "precip3_16lev"; "precip_11lev"
  diffgrp0 := fspan(0,2,11)

  dcname   := read_colormap_file("CBR_coldhot")
  dcname   := dcname(::-1,:)
  diffgrp1 := (/-20.0,-10.0,-5.0,-1.0,1.0,5.0,10.0,20.0/)
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
