;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;sub region for sythentic analysis
  lats  = -25.         ;21. -2  ;-30
  late  =  -5.         ;36. -2  ;-10
  lons  = -95 + 360. ;-137. -8.0 + 360
  lone  = -75 + 360. ;-122. -8.0 + 360
  ;;draw box for the region
  yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
  xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line

 ;;setup for fill pattern
  FillColor           = "Grey75" ;"Grey65";"Grey28";"Grey";"Black"
  FillPattern         = 0 ;17 ;6 ;17;11
  FillScaleF          = 0.1
  FillDotSizeF        = 0.0002    ; increase dot size of contours
  FillOpacityF        = 1.0 ;0.15

  colorbar_without_white := False

  ;;;;color groups (see fifth column in var_info.txt);;;;

  if(color_grp(iv).eq.1)then 

    cmap       = read_colormap_file("nice_gfdl")
    ocname    := cmap(60:180,:)
    dcname     = "ViBlGrWhYeOrRe"
    orgcolors := (/3,4,5,7,10,11,12,13,14,17/)
    orgcolors  = orgcolors(::-1)
 
    dcname     = "BlueDarkRed18" ;"NCV_blu_red"
    gsn_define_colormap(wks1, dcname)
    diffcolors = (/2,3,5,6,0,8,9,10,12/)  
    diffcolors = (/2,4,6,8,0,12,14,16,18/)
    delete(cmap)

  else if(color_grp(iv).eq.2)then

    cmap       = read_colormap_file("BlAqGrYeOrReVi200")
    ocname    := cmap
    ocname    := ocname(::-1,:)
    dcname     = "ViBlGrWhYeOrRe"
    gsn_define_colormap(wks1,dcname)
    diffcolors = (/15,26,33,42,50,75,85,102,2/)
    delete(cmap)

  else if(color_grp(iv).eq.3)then

    cmap       = read_colormap_file("BlAqGrYeOrReVi200")
    ocname    := cmap
    dcname     = "ViBlGrWhYeOrRe"
    gsn_define_colormap(wks1,dcname)
    diffcolors = (/15,26,33,42,50,75,85,102,2/)
    delete(cmap)

  else if(color_grp(iv).eq.4)then

    ocname    := "precip3_16lev" ;"percent_11lev" ;"prcp_1" ;"NCV_banded" ;"amwg"  ; color bar name for cntrol model plota
    orgcolors := (/3,4,5,7,10,11,12,13,14,17/)
    orgcolors  = orgcolors(::-1)

    dcname     = "precip4_diff_19lev" ;"MPL_BrBG"
    gsn_define_colormap(wks1,dcname) ; choose colormap
    diffcolors = (/2,7,9,10,0,13,14,16,21/)   
  else

    ocname    := "t2m_29lev"
    orgcolors := (/3,4,5,7,10,11,12,13,14,17/)
    orgcolors  = orgcolors(::-1)

    gsn_define_colormap(wks1, "ViBlGrWhYeOrRe") ; choose colormap
    diffcolors = (/15,26,33,42,50,75,85,102,2/)

  end if 
  end if 
  end if 
  end if 


 if (Varname .eq. "CLDTOT") then
   orglev    := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)
   diffgrp   := (/-0.20,-0.10,-0.05,-0.02,0.02,0.05,0.10,0.20/)
   plot_contour1 = False
 end if

 if (Varname .eq. "CLDLOW") then
   orglev    := (/0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95/)
   diffgrp   := (/-0.20,-0.10,-0.05,-0.02,0.02,0.05,0.10,0.20/)
   plot_contour1 = False
 end if

 if (Varname .eq. "CLDHGH") then
   orglev    := (/0.01,0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.7/)
   diffgrp   := (/-0.20,-0.10,-0.05,-0.02,0.02,0.05,0.10,0.20/)
   plot_contour1 = False
 end if

 if (Varname .eq. "SWCF") then
   orglev    := (/-100,-90,-80,-60,-50,-40,-30,-20,-10,-5,-1/)
   diffgrp   := (/-20,-10,-5,-2,2,5,10,20/)*1.0
   plot_contour1 = False
 end if

 if (Varname .eq. "LWCF") then
   orglev    := (/-5, 0, 5, 10, 15, 20, 25, 30, 40, 50, 60/)
   diffgrp   := (/-15,-10,-5,-2,2,5,10,15/)*1.0
   plot_contour1 = False
 end if

 if (Varname .eq. "NETCF") then
   orglev    := (/-60,-50,-40,-30,-20,-25,-10,-5, 0, 5, 10/)
   diffgrp   := (/-15,-10,-5,-2,2,5,10,15/)*1.0
   plot_contour1 = False
 end if 

 if (Varname .eq. "FSNTOA") then
   orglev    := (/10,30,60,90,120,150,180,210,240,270,300,330,360,390,420/)
   diffgrp   := (/-20,-10,-5,-2,2,5,10,20/)*1.0
   plot_contour1 = False
 end if

 if (Varname .eq. "FLUT") then
   orglev    := (/100,115,130,145,160,175,190,205,220,235,250,265,280,295,310/)
   diffgrp   := (/-15,-10,-5,-2,2,5,10,15/)*1.0
   plot_contour1 = False
 end if

 if (Varname .eq. "TGCLDLWP") then
   orglev    := (/0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120/);
   diffgrp   := (/-20,-15,-10,-5,5,10,15,20/)*1.0
   plot_contour1 = False
 end if 

 if (Varname .eq. "TGCLDIWP") then
   orglev    := (/0,0.1,0.5,1, 2, 5, 8, 10, 15, 20, 30, 40/)
   diffgrp   := (/-10,-5,-2,-1.0,1.0,2,5,10/)*1.0
   plot_contour1 = False
 end if

 if (Varname .eq. "PRECC") then
   orglev    := (/0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0/)
   diffgrp   := (/ -1.0, -0.5, -0.2, -0.1, 0.1,0.2, 0.5, 1.0/)
   plot_contour1 = False
 end if

 if (Varname .eq. "PRECL") then
   orglev    := (/0.01,0.05,0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0/)
   diffgrp   := (/ -1.0, -0.5, -0.2, -0.1, 0.1,0.2, 0.5, 1.0/)
   plot_contour1 = False
 end if

 if (Varname .eq. "PRECT") then
   orglev    := (/0.01,0.05,0.1,0.5,1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0/)
   diffgrp   := (/ -1.0, -0.5, -0.2, -0.1, 0.1,0.2, 0.5, 1.0/)
   plot_contour1 = False
 end if

 if (Varname .eq. "TS") then
   orglev    := (/-30,-20,-10,-5,0,5,10,15,20,25,30/)*1.0
   diffgrp   := (/ -0.6, -0.4, -0.2, -0.1, 0.1,0.2, 0.4, 0.6/)
   plot_contour1 = False
 end if

