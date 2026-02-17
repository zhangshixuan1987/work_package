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

  ;gsn_define_colormap(wks1,"StepSeq25")
  ;gsn_merge_colormaps(wks1,"StepSeq25","cb_9step")
  ;colors  = (/55,17,2,96,45,59/)
  gsn_define_colormap(wks1,"cb_9step")
  colors  = (/1,2,6,28/)

  marker  = (/4,16,11,9/)       ;;line markers
  mksizf  = (/1.5,1.0,1.0,1.2/)/80.0 
  mkthik  = (/3,3,3,3/)
  mkopcy  = (/0.8,1.0,1.0,1.0/)
  lnpatt  = (/0,0,0,0/)
  lnthik  = (/6,4,4,4/)/2   ;;line thichness
  lnopcy  = (/0.8,1.0,1.0,1.0/)

  ;;color for the map boundary line
  lndcolor = "Black" ;"Grey50"

  cmap = read_colormap_file("percent_11lev")
  cmap(0,:) = 0 
 ;colorbar setup for difference plots
  cnlevels  := ispan(100,1000,100) ;(/0,1,2,4,6,8,10,12,14,16,18,20,30,40/)
