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

 cmap       = read_colormap_file("BlAqGrYeOrReVi200")
 ocname    := cmap
 dcname     = "ViBlGrWhYeOrRe"
 gsn_define_colormap(wks1,dcname)
 diffcolors = (/15,26,33,42,44, 50, 65, 75,85,102,2/)
 delete(cmap)

 if (Varname .eq. "T") then
     orglev     := (/200,210,220,230,240,250,260,270,280,285,290,295/)
     diffgrp1   := (/-6.0,-4.0,-2.0,-1.0,-0.5,0.5,1.0,2.0,4.0,6.0/)
     diffgrp2   := (/-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04/)
 end if 

 if (Varname .eq. "U") then
     orglev    := (/-20,-15,-10,-5,0,5,10,15,20,25,30,40/)
     diffgrp1   := (/-6.0,-4.0,-2.0,-1.0,-0.5,0.5,1.0,2.0,4.0,6.0/)
 end if

 if (Varname .eq. "CLOUD".or.Varname .eq. "CLOUDFRAC_CLUBB") then
     orglev    := (/0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8/)
     diffgrp1   := (/-6.0,-4.0,-2.0,-1.0,-0.5,0.5,1.0,2.0,4.0,6.0/)/10.0
     diffgrp2    := (/-20,-10,-5,-1,1,5,10,20/)
 end if 

 if (Varname .eq. "RELHUM") then
     orglev    := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)
     diffgrp1  := (/-30.,-20.,-10.,-5.,-1.,-0.1,0.1,1.,5.,10.,20.,30./)
     diffgrp2  := (/-2,-1,-0.5,-0.1,0.1,0.5,1,2/)
 end if

 if (Varname .eq. "Q".or.Varname.eq."rvm") then
     orglev    := (/0.01,0.05,0.1,0.2,0.5,1,2,5,8,10/)
     diffgrp1  := (/-1.,-0.5,-0.2,-0.1,-0.05,-0.02,0.02,0.05,0.1,0.2,0.5,1./)
     diffgrp2  := (/-2,-1,-0.5,-0.1,0.1,0.5,1.0,2.0/)
 end if

