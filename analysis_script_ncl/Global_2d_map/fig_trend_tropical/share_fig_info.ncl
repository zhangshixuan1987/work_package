  maxlev   = 11
  if (vars.eq."TREFHT" .or. vars.eq."TS".or. vars.eq."SST") then
     cnlevs1  = fspan(-0.2,0.2,9) 
     cnlevs2  = fspan(-5,5,11)
  else if (vars.eq."PSL") then 
     cnlevs1  = fspan(-0.2,0.2,9) 
     cnlevs2  = fspan(-4,4,11)
  else 
    print("please specify the countour levels for variable " + vars) 
  end if 
  end if 
