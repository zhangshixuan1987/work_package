
 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvar3     = 4
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
  VarList3(ivar3) = "PS"
  VstList3(ivar3) = "PS"
  VnmList3(ivar3) = "Surface pressure"
  UnitList3(ivar3)= "hPa"

  ivar3 = ivar3+1
  VarList3(ivar3) = "PSL"
  VstList3(ivar3) = "PSL"
  VnmList3(ivar3) = "Sea level pressure"
  UnitList3(ivar3)= "hPa"

