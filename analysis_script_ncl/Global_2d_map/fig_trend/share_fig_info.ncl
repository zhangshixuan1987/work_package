  maxlev   = 11
  if (vars.eq."TREFHT" .or. vars.eq."TS".or. vars.eq."SST") then
     cnlevs1  = fspan(-1.5,1.5,11) 
     cnlevs2  = fspan(-5,5,11)
  else if (vars.eq."PSL") then 
     cnlevs1  = fspan(-2,2,11) 
     cnlevs2  = fspan(-5,5,11)
  else 
    print("please specify the countour levels for variable " + vars) 
  end if 
  end if 
