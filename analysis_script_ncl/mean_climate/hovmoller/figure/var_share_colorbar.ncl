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

 FillColor           = "Grey28";"Grey";"Black"
 FillPattern         = 17 ;6 ;17;11
 FillScaleF          = 0.5
 FillDotSizeF        = 0.002    ; increase dot size of contours
 FillOpacityF        = 0.8


 ;gsn_define_colormap(wks1,"StepSeq25")
 ;gsn_merge_colormaps(wks1,"StepSeq25","cb_9step")
 ;colors  = (/55,17,2,96/)
  gsn_define_colormap(wks1,"cb_9step")
  colors  = (/1,2,14,28,3/)

  markers = (/16,4,7,8,2,11/)       ;;line markers
  mkszf   = (/0.8,0.8,0.8,0.8,0.8/)/2.0   ;;line marker size
  lnthick = (/3,3,3,3,3/)*1.5   ;;line thichness
  dashind = (/0,0,0,0,0/)       ;;line dash pattern

  ;;set up for fig size
  vpHeightF  = 0.8
  vpWidthF   = 0.4
  vpXF       = 0.25
  FontHeight = 0.024

 if (Varname .eq. "PRECT") then

  gsn_define_colormap (wks1,"percent_11lev")
  colorgroup1 := (/2,3,4,5,6,7,8,9,10,11,12,13/)

  xmin1    := 0.
  xmax1    := 29.

  if (irg.eq.0)then
   ymin1   :=  0
   ymax1   := 16.
   yint1   :=  2
   ylevs   := (/0,0.5,1,2,3,4,5,6,8,10,12/)
  end if

  if (irg.eq.1)then
   ymin1   := 0
   ymax1   := 16.
   yint1   := 2
   ylevs   := (/0,0.5,1,2,3,4,5,6,8,12,16/)
  end if

  if (irg.eq.2)then
   ymin1   := 0
   ymax1   := 18.
   yint1   := 2
   ylevs   := (/0,1,2,3,4,5,6,8,10,12,14/)
  end if

  amJust               = "TopRight"    ; Use bottom right corner of box
  amParallelPosF       =   0.48             ; Move legend to right
  amOrthogonalPosF     =  -0.48             ; Move legend down.
 end if

 if (Varname .eq. "FLUT") then

  gsn_define_colormap(wks1, "WhBlGrYeRe")
  colorgroup1 := (/5,10,20,30,40,50,60,70,80,90,100/)

  xmin1   := 0.
  xmax1   := 29.

  if (irg.eq.0)then
   ymin1   := 180
   ymax1   := 260.
   yint1   := 8
   yref    := (/130,210/)
   ylevs   := (/200,208,216,224,228,232,236,240,244,248/)
  end if 

  if (irg.eq.1)then
   ymin1   := 200
   ymax1   := 254.
   yint1   := 6
   yref    := (/180,240/)
   ylevs   := (/200,208,212,216,220,224,228,232,236,242/)
  end if

  if (irg.eq.2)then
   ymin1   := 192
   ymax1   := 292.
   yint1   := 10
   yref    := (/80,180/)
   ylevs   := (/210,216,222,228,232,248,260,272,284,296/)
  end if
 
  amJust               = "TopRight"    ; Use bottom right corner of box
  amParallelPosF       =   0.48             ; Move legend to right
  amOrthogonalPosF     =  -0.48             ; Move legend down.
 end if

