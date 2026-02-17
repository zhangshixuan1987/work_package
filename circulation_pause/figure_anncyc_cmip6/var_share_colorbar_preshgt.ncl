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

  if (varnam.eq."T") then 
     cnlev1  := fspan(200,320,13)
     cnlev2  := fspan(-2.0,2.0,11)
  else if (varnam.eq."U") then 
     cnlev1  := fspan(-30.0,30.0,13)
     cnlev2  := fspan(-2.4,2.4,13)
  else
    print("please specify the countour levels for variable " + vars)
  end if
  end if

