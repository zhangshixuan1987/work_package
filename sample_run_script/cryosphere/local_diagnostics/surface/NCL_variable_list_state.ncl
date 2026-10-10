
 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvar3     = 11
  VarList3  = new(nvar3, "string")
  VstList3  = new(nvar3, "string")
  VnmList3  = new(nvar3, "string")
  UnitList3 = new(nvar3, "string")
  
  ivar3 = 0 
  VarList3(ivar3) = "TREFHT"
  VstList3(ivar3) = "TREFHT"
  VnmList3(ivar3) = "2m air temperature"
  UnitList3(ivar3)= "K"

  ivar3 = ivar3+1
  VarList3(ivar3) = "TS"
  VstList3(ivar3) = "TS"
  VnmList3(ivar3) = "Surface temperature"
  UnitList3(ivar3)= "K"

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
  VarList3(ivar3) = "QREFHT"
  VstList3(ivar3) = "QREFHT"
  VnmList3(ivar3) = "2 m humidity"
  UnitList3(ivar3)= "g kg~S~-1~N~"

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

  ivar3 = ivar3+1
  VarList3(ivar3) = "TMQ"
  VstList3(ivar3) = "TMQ"
  VnmList3(ivar3) = "Total precipitable water"
  UnitList3(ivar3)= "kg m~S~2~N~"


  ivar3 = ivar3+1
  VarList3(ivar3) = "PS"
  VstList3(ivar3) = "PS"
  VnmList3(ivar3) = "Surface pressure"
  UnitList3(ivar3)= "hPa"

  ivar3 = ivar3+1
  VarList3(ivar3) = "PSL"
  VstList3(ivar3) = "PSL"
  VnmList3(ivar3) = "Sea level pressure"
  UnitList3(ivar3)= "hPa"

