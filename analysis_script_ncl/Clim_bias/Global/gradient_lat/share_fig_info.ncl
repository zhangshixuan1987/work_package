 ;sub region for sythentic analysis
  lats  = -40. ;21. -2  ;-30
  late  = -25. ;36. -2  ;-10
  lons  = 20   ;-137. -8.0 + 360
  lone  = 75   ;-122. -8.0 + 360
 ;draw box for the region
  yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
  xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line

 
 maxlev   = 11
  if (vars.eq."TREFHT" .or. vars.eq."TS".or. vars.eq."SST") then
     cnlevs2  = fspan(-3,3,21)
  else if (vars.eq."PSL".or.vars.eq."PS") then 
     cnlevs1  = (/-200.,-150.,-100.,-50.,-10.,0.,10.,50.,100.,150.,200./)
     cnlevs2  = (/-12.,-8.,-4.,-2.,-0.5,0.5,2.,4.,8.,12./)
  else 
    print("please specify the countour levels for variable " + vars) 
  end if 
  end if 
