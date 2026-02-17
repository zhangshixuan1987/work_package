  maxlev   = 11
  if (vars.eq."SST".or.vars.eq."TREFHT" .or. vars.eq."TS") then
     cnlevs1  = (/-0.4,-0.3,-0.2,-0.1,-0.01,0.0,0.01,0.1,0.2,0.3,0.4/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."PSL") then 
     cnlevs1  = (/-200.,-150.,-100.,-50.,-10.,0.,10.,50.,100.,150.,200./)
     cnlevs2  = (/-80.,-40.,-20.,-10.,-5.,5.,10.,20.,40.,80./)
  else if (vars.eq."TMQ") then
     cnlevs1  = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.4,0.8,1.2/)
     cnlevs2  = (/-0.6,-0.4,-0.2,-0.1,-0.01,0.01,0.1,0.2,0.4,0.6/)
  else if (vars.eq."TGCLDLWP".or.vars.eq."TGCLDIWP") then
     cnlevs1  = (/-4.0,-3.0,-2.0,-1.0,-0.1,0.0,0.1,1.0,2.0,3.0,4.0/)
     cnlevs2  = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/)
  else if (vars.eq."PRECT".or.vars.eq."PRECST") then
     cnlevs1  = (/-0.4,-0.3,-0.2,-0.1,-0.01,0.0,0.01,0.1,0.2,0.3,0.4/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."TAUX".or.vars.eq."TAUY") then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)*10.0
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)*10.0
  else if (vars.eq."ALBEDO" .or. vars.eq."ALBEDO_SRF" ) then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."U10") then
     cnlevs1  = (/-0.8,-0.5,-0.3,-0.1,-0.01,0.0,0.01,0.1,0.3,0.5,0.8/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."CLDTOT") then
     cnlevs1  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.0,0.01,0.05,0.1,0.2,0.3/)
     cnlevs2  = (/-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3/)
  else if (vars.eq."LWCF" .or. vars.eq."SWCF" .or. vars.eq."NETCF") then
     cnlevs1  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.0,0.1,0.5,1.0,2.0,4.0/)
     cnlevs2  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,4.0/)
  else if (vars.eq."LWCF_SRF" .or. vars.eq."SWCF_SRF" .or. vars.eq."NETCF_SRF") then
     cnlevs1  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.0,0.1,0.5,1.0,2.0,4.0/)
     cnlevs2  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,4.0/)
  else if (vars.eq."SHFLX") then
     cnlevs1  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.0,0.1,0.5,1.0,2.0,4.0/)
     cnlevs2  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,4.0/)
  else if (vars.eq."LHFLX" .or. vars.eq."RESSURF") then
     cnlevs1  = (/-8.0,-4.0,-2.0,-1.0,-0.2,0.0,0.2,1.0,2.0,4.0,8.0/)
     cnlevs2  = (/-8.0,-4.0,-2.0,-1.0,-0.2,0.2,1.0,2.0,4.0,8.0/)
  else if (vars.eq."FLNS" .or. vars.eq."FSNS" .or. vars.eq."NetRad_SRF") then
     cnlevs1  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.0,0.1,0.5,1.0,2.0,4.0/)
     cnlevs2  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,4.0/)
  else if (vars.eq."FLNSC" .or. vars.eq."FSNSC" .or. vars.eq."NetRad_SRF_clr") then
     cnlevs1  = (/-2.0,-1.0,-0.5,-0.1,-0.01,0.0,0.01,0.1,0.5,1.0,2.0/)
     cnlevs2  = (/-2.0,-1.0,-0.5,-0.1,-0.01,0.01,0.1,0.5,1.0,2.0/)
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

