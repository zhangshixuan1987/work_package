  maxlev   = 11
  if (vars.eq."TREFHT" .or. vars.eq."TS".or. vars.eq."SST") then
     cnlevs1  = fspan(0.0,4.0,11) 
     cnlevs2  = fspan(0.0,4.0,11) 
  else if (vars.eq."PSL") then 
     cnlevs1  = fspan(0.0,20.0,11)  
     cnlevs2  = fspan(0.0,20.0,11)
  else if (vars.eq."TMQ") then
     cnlevs1  = fspan(0.0,20.0,11)
     cnlevs2  = fspan(0.0,20.0,11)
  else if (vars.eq."TGCLDLWP".or.vars.eq."TGCLDIWP") then
     cnlevs1  = fspan(0.0,2.0,11)
     cnlevs2  = fspan(0.0,2.0,11)
  else if (vars.eq."PRECT") then
     cnlevs1  = fspan(0.0,1.0,11)
     cnlevs2  = fspan(0.0,1.0,11)
  else if (vars.eq."PRECST") then
     cnlevs1  = fspan(0.0,1.0,11)
     cnlevs2  = fspan(0.0,1.0,11)
  else if (vars.eq."TAUX".or.vars.eq."TAUY") then
     cnlevs1  = fspan(0.0,1.0,11)
     cnlevs2  = fspan(0.0,1.0,11)
  else if (vars.eq."ALBEDO" .or. vars.eq."ALBEDO_SRF" ) then
     cnlevs1  = fspan(0.0,1.0,11)
     cnlevs2  = fspan(0.0,1.0,11)
  else if (vars.eq."U10") then
     cnlevs1  = fspan(0.0,1.0,11)
     cnlevs2  = fspan(0.0,1.0,11)
  else if (vars.eq."CLDTOT".or.vars.eq."CLDLOW") then
     cnlevs1  = fspan(0.0,1.0,11)
     cnlevs2  = fspan(0.0,1.0,11)
  else if (vars.eq."LWCF" .or. vars.eq."SWCF" .or. vars.eq."NETCF") then
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else if (vars.eq."LWCF_SRF" .or. vars.eq."SWCF_SRF" .or. vars.eq."NETCF_SRF") then
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else if (vars.eq."SHFLX") then
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else if (vars.eq."NET_FLUX_SRF") then 
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else if (vars.eq."LHFLX" .or. vars.eq."RESSURF") then
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else if (vars.eq."FLNS" .or. vars.eq."FSNS" .or. vars.eq."NetRad_SRF") then
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else if (vars.eq."FLNSC" .or. vars.eq."FSNSC" .or. vars.eq."NetRad_SRF_clr") then
     cnlevs1  = fspan(0.0,10.0,11)
     cnlevs2  = fspan(0.0,10.0,11)
  else 
    print("please specify the countour levels for variable " + vars) 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if
  end if 
  end if 
  end if 
  end if
  end if
  end if
  end if
  end if
  end if 
  end if 
  end if 

