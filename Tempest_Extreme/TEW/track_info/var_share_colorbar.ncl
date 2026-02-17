;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp1),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;sub region for sythentic analysis
  lats  =  0.         ;21. -2  ;-30
  late  = 50.         ;36. -2  ;-10
  lons  = -100 + 360. ;-137. -8.0 + 360
  lone  = -40 + 360.  ;-122. -8.0 + 360

  ;;draw box for the region
  yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
  xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line

  ;;color for the map boundary line
   lndcolor = "Black" ;"Grey50"

  ;colorbar setup for difference plots
  cncolors  := (/72,62,55,5,7,10,0,120,14,16,91,85,77/)
  cnlevels  := (/0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25/) ;1.5,2.0/)
