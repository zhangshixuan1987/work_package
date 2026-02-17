;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;this script is used to difine the absolute difference levels (diffgrp1),;;; 
;;;;;;;;the relative difference levels (diffgrp2), then choose color bar for plot;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 ;sub region for sythentic analysis
  lats  = -25.       ;21. -2  ;-30
  late  =  -5.       ;36. -2  ;-10
  lons  = -95 + 360. ;-137. -8.0 + 360
  lone  = -75 + 360. ;-122. -8.0 + 360

  ;;draw box for the region
  yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
  xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line

  ;;color for the map boundary line
   lndcolor = "Black" ;"Grey50"

  ;colorbar setup for difference plots
  cncolors  := (/72,62,55,5,7,10,0,120,14,16,91,85,77/)
  cnlevels  := ispan(10,360,35) ;(/0,1,2,4,6,8,10,12,14,16,18,20,30,40/)
