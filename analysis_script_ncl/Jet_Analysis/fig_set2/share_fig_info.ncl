  maxlev   = 11
  if (vars.eq."TREFHT" .or. vars.eq."TS") then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)
     cnlevs2  = (/-0.2,-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1,0.2/)
  else if (vars.eq."PSL") then 
     cnlevs1  = (/-120.,-80.,-40.,-20.,-10.,0.,10.,20.,40.,80.,120./)
     cnlevs2  = (/-60.,-40.,-20.,-10.,-1.,1.,10.,20.,40.,60./)
  else if (vars.eq."TMQ") then
     cnlevs1  = (/-0.4,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.4/)
     cnlevs2  = (/-0.2,-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1,0.2/)
  else if (vars.eq."TGCLDLWP".or.vars.eq."TGCLDIWP") then
     cnlevs1  = (/-2.0,-1.0,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1.0,2.0/)
     cnlevs2  = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  else if (vars.eq."PRECT".or.vars.eq."PRECST") then
     cnlevs1  = (/-0.12,-0.08,-0.04,-0.02,-0.01,0.0,0.01,0.02,0.04,0.08,0.12/)
     cnlevs2  = (/-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.04,0.06,0.08/)
  else if (vars.eq."TAUX".or.vars.eq."TAUY") then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)*10.0
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)*10.0
  else if (vars.eq."ALBEDO" .or. vars.eq."ALBEDO_SRF" ) then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."U10") then
     cnlevs1  = (/-0.4,-0.3,-0.2,-0.1,-0.01,0.0,0.01,0.1,0.2,0.3,0.4/)
     cnlevs2  = (/-0.2,-0.1,-0.05,-0.02,-0.01,0.01,0.02,0.05,0.1,0.2/)
  else if (vars.eq."CLDTOT") then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."LWCF" .or. vars.eq."SWCF" .or. vars.eq."NETCF") then
     cnlevs1  = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.4,0.8,1.2/)
     cnlevs2  = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  else if (vars.eq."LWCF_SRF" .or. vars.eq."SWCF_SRF" .or. vars.eq."NETCF_SRF") then
     cnlevs1  = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.4,0.8,1.2/)
     cnlevs2  = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  else if (vars.eq."SHFLX") then
     cnlevs1  = (/-2.4,-1.2,-0.6,-0.3,-0.1,0.0,0.1,0.3,0.6,1.2,2.4/)
     cnlevs2  = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  else if (vars.eq."LHFLX" .or. vars.eq."RESSURF") then
     cnlevs1  = (/-2.4,-1.2,-0.6,-0.3,-0.1,0.0,0.1,0.3,0.6,1.2,2.4/)
     cnlevs2  = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  else if (vars.eq."FLNS" .or. vars.eq."FSNS" .or. vars.eq."NetRad_SRF") then
     cnlevs1  = (/-2.4,-1.2,-0.6,-0.3,-0.1,0.0,0.1,0.3,0.6,1.2,2.4/)
     cnlevs2  = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  else if (vars.eq."FLNSC" .or. vars.eq."FSNSC" .or. vars.eq."NetRad_SRF_clr") then
     cnlevs1  = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.4,0.8,1.2/)
     cnlevs2  = (/-0.8,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2,0.4,0.8/)
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

