
 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvar3     = 8
  VarList4  = new(nvar3, "string")
  VstList4  = new(nvar3, "string")
  VnmList4  = new(nvar3, "string")
  UnitList4 = new(nvar3, "string")

  ivar4 = 0 
  VarList4(ivar4) = "LWCF"
  VstList4(ivar4) = "LWCF"
  VnmList4(ivar4) = "Longwave cloud radiative effect"
  UnitList4(ivar4)= "W m~S~-2~N~"
 
  ivar4 = ivar4+1
  VarList4(ivar4) = "SWCF"
  VstList4(ivar4) = "SWCF"
  VnmList4(ivar4) = "Shortwave cloud radiative effect"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "NETCF"
  VstList4(ivar4) = "NETCF"
  VnmList4(ivar4) = "Cloud radiative effect"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FLNS"
  VstList4(ivar4) = "FLNS"
  VnmList4(ivar4) = "Net longwave flux at surface"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "SOLIN"
  VstList4(ivar4) = "SOLIN"
  VnmList4(ivar4) = "Solar insolation"
  UnitList4(ivar4)= "W m~S~-2~N~"


  ivar4 = ivar4+1
  VarList4(ivar4) = "FLUT"
  VstList4(ivar4) = "FLUT"
  VnmList4(ivar4) = "Upwelling longwave flux at top of model"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSNTOA"
  VstList4(ivar4) = "FSNTOA"
  VnmList4(ivar4) = "Downwelling solar flux at top of atmosphere"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSUTOA"
  VstList4(ivar4) = "FSUTOA"
  VnmList4(ivar4) = "Upwelling solar flux at top of atmosphere"
  UnitList4(ivar4)= "W m~S~-2~N~"
