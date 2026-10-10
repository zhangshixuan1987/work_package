  maxlev   = 11
  if (vars.eq."SST".or.vars.eq."TREFHT" .or. vars.eq."TS") then
     cnlevs1  = (/-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5/)
  else if (vars.eq."TAUX".or.vars.eq."TAUY") then
     cnlevs1  = (/-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05/)
     cnlevs1  = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.4,0.8,1.2/)
  else 
    print("please specify the countour levels for variable " + vars) 
  end if 
  end if 
