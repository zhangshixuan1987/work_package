
 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvar4     = 14
  VarList4  = new(nvar4, "string")
  VstList4  = new(nvar4, "string")
  VnmList4  = new(nvar4, "string")
  UnitList4 = new(nvar4, "string")

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
  VarList4(ivar4) = "FLUT"
  VstList4(ivar4) = "FLUT"
  VnmList4(ivar4) = "TOA outgoing longwave flux"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSNTOA"
  VstList4(ivar4) = "FSNTOA"
  VnmList4(ivar4) = "Net solar flux at top of atmosphere"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSNTOAC"
  VstList4(ivar4) = "FSNTOAC"
  VnmList4(ivar4) = "Clear-sky net solar flux at top of atmosphere"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSUTOA"
  VstList4(ivar4) = "FSUTOA"
  VnmList4(ivar4) = "Upwelling solar flux at top of atmosphere"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSUTOAC"
  VstList4(ivar4) = "FSUTOAC"
  VnmList4(ivar4) = "Clear-sky upwelling solar flux at top of atmosphere"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSDS"
  VstList4(ivar4) = "FSDS"
  VnmList4(ivar4) = "Downwelling solar flux at surface"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FSNS"
  VstList4(ivar4) = "FSNS"
  VnmList4(ivar4) = "Net solar flux at surface"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FLNS"
  VstList4(ivar4) = "FLNS"
  VnmList4(ivar4) = "Net longwave flux at surface"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "FLDS"
  VstList4(ivar4) = "FLDS"
  VnmList4(ivar4) = "Downwelling longwave flux at surface"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "RESTOM"
  VstList4(ivar4) = "RESTOM"
  VnmList4(ivar4) = "Top of model residual energy balance"
  UnitList4(ivar4)= "W m~S~-2~N~"

  ivar4 = ivar4+1
  VarList4(ivar4) = "QNET"
  VstList4(ivar4) = "QNET"
  VnmList4(ivar4) = "Net heat fluxes at surface"
  UnitList4(ivar4)= "W m~S~-2~N~"
