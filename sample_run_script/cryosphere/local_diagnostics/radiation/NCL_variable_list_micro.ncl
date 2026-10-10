
 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvar3     = 8
  VarList3  = new(nvar3, "string")
  VstList3  = new(nvar3, "string")
  VnmList3  = new(nvar3, "string")
  UnitList3 = new(nvar3, "string")
  
  ivar3 = 0 
  VarList3(ivar3) = "NUMLIQ"
  VstList3(ivar3) = "NUMLIQ"
  VnmList3(ivar3) = "Integrated liq num"
  UnitList3(ivar3)= "x10~S~6~N~ kg~S~-1~N~"

  ivar3 = ivar3+1
  VarList3(ivar3) = "NUMICE"
  VstList3(ivar3) = "NUMICE"
  VnmList3(ivar3) = "Integrated ice num"
  UnitList3(ivar3)= "x10~S~3~N~ kg~S~-1~N~"

  ivar3 = ivar3+1
  VarList3(ivar3) = "CCN3"
  VstList3(ivar3) = "CCN3"
  VnmList3(ivar3) = "CCN3"
  UnitList3(ivar3)= "# cm~S~-3~N~"

  ivar3 = ivar3+1
  VarList3(ivar3) = "OMEGA500"
  VstList3(ivar3) = "OMEGA500"
  VnmList3(ivar3) = "OMEGA500"
  UnitList3(ivar3)= "hPa day~S~-1~N~"

  ivar3 = ivar3+1
  VarList3(ivar3) = "TH7001000"
  VstList3(ivar3) = "TH7001000"
  VnmList3(ivar3) = "Lower tropospheric stability"
  UnitList3(ivar3)= "K"

  ivar3 = ivar3+1
  VarList3(ivar3) = "U10"
  VstList3(ivar3) = "U10"
  VnmList3(ivar3) = "10 m wind speed"
  UnitList3(ivar3)= "m s~S~-1~N~"

  ivar3 = ivar3+1
  VarList3(ivar3) = "LHFLX"
  VstList3(ivar3) = "LHFLX"
  VnmList3(ivar3) = "Latent heat flux"
  UnitList3(ivar3)= "W m~S~-2~N~"

  ivar3 = ivar3+1
  VarList3(ivar3) = "SHFLX"
  VstList3(ivar3) = "SHFLX"
  VnmList3(ivar3) = "Sensible heat flux"
  UnitList3(ivar3)= "W m~S~-2~N~"

