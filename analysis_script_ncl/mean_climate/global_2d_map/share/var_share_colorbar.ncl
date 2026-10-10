;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp1),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;;colors for the original plot
  omap0      = read_colormap_file("prcp_3")
  orgcolors := (/omap0(3,:),omap0(4,:),omap0(5,:),omap0(6,:), omap0(7,:), \
                 omap0(10,:),omap0(11,:),omap0(12,:),omap0(13,:), \
                 omap0(14,:),omap0(15,:), omap0(16,:), omap0(17,:)/)

  ;colors for the difference plot 
  gsn_merge_colormaps(wks1,"temp_diff_18lev","WhViBlGrYeOrReWh")

  ;colorbar setup for difference plots
 ;cncolors   = (/72,62,55,30,5,7,9,10, 0, 120,119,14,16, 104,91,85,77 /)
 ;cncolors   = (/72,62,55,5,7,10,120,0,0,14,16,91,85,77/)
  cncolors   = (/72,55,5,7,10,0,120,14,16,85,77/)

  if (Varname .eq. "FLUT") then
    orglev    := (/180,190,200,210,220,230,240,250,260,270,280/)
    cnlevels  := (/-30,-25,-20,-15,-10,-5,5,10,15,20,25,30/)
  else if (Varname .eq. "CLDTOT") then
    orglev    := (/0.2,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95/) 
    cnlevels  := (/-0.40,-0.30,-0.20,-0.15,-0.10,-0.05,0.05,0.10,0.15,0.20,0.30,0.40/)
  else if (Varname .eq. "NETCF") then
    orgcolors  = orgcolors(::-1,:) 
    orglev    := (/-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10/)
    cnlevels  := (/-50,-40,-30,-20,-10,-5,5,10,20,30,40,50/)
  else if (Varname .eq. "SWCF") then
    orgcolors  = orgcolors(::-1,:)
    orglev    := (/-140,-120,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10/)
    cnlevels  := (/-50,-40,-30,-20,-10,-5,5,10,20,30,40,50/)
  else if (Varname .eq. "LWCF") then
    orglev    := (/-2,0,2,5,10,15,20,25,30,40,60,80/)
    cnlevels  := (/-25,-20,-15,-10,-5,-2,2,5,10,15,20,25/)
  else if (Varname .eq. "PRECT") then
    orglev    := (/0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0/)
    cnlevels   = (/-2.0,-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5,2.0/)
  else if (Varname .eq. "PRECST") then
    orglev    := (/0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0/)
    cnlevels   = (/-1.5,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,1.5/)
  else if (Varname .eq. "PRECL") then
    orglev    := (/0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,8.0,10.0/)
    cnlevels   = (/-4.0,-3.0,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,3.0,4.0/)
  else if (Varname .eq. "SICONC") then
    orglev    := (/1,2,5,10,20,30,40,50,60,70,80,90/)*1.0
    cnlevels   = (/-1.5,-1.0,-0.5,-0.2,-0.1,0.1,0.2,0.5,1.0,1.5/)*10.0
  else
    print("please specify color bar ...")
    exit
  end if 
  end if 
  end if
  end if
  end if
  end if 
  end if 
  end if 
  end if 
