
 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  nvar1     = 8
  VarList1  = new(nvar1, "string")
  VstList1  = new(nvar1, "string")
  VnmList1  = new(nvar1, "string")
  UnitList1 = new(nvar1, "string")
  ERA5Var1= new(nvar1, "string")
  
  ivar1 = 0 
  VarList1(ivar1) = "CLDTOT"
  VstList1(ivar1) = "CLDTOT"
  VnmList1(ivar1) = "Total cloud cover"
  UnitList1(ivar1)= "fraction"
  ERA5Var1(ivar1) = "tcc"

  ivar1 = ivar1+1
  VarList1(ivar1) = "CLDLOW"
  VstList1(ivar1) = "CLDLOW"
  VnmList1(ivar1) = "Low cloud cover"
  UnitList1(ivar1)= "fraction"
  ERA5Var1(ivar1) = "lcc"

  ivar1 = ivar1+1
  VarList1(ivar1) = "CLDHGH"
  VstList1(ivar1) = "CLDHGH"
  VnmList1(ivar1) = "High cloud cover"
  UnitList1(ivar1)= "fraction"
  ERA5Var1(ivar1) = "hcc"

  ivar1 = ivar1+1
  VarList1(ivar1) = "TGCLDLWP"
  VstList1(ivar1) = "TGCLDLWP"
  VnmList1(ivar1) = "Liquid water path"
  UnitList1(ivar1)= "g m~S~-2~N~"
  ERA5Var1(ivar1) = "tclw"

  ivar1 = ivar1+1
  VarList1(ivar1) = "TGCLDIWP"
  VstList1(ivar1) = "TGCLDIWP"
  VnmList1(ivar1) = "Ice water path"
  UnitList1(ivar1)= "g m~S~-2~N~"
  ERA5Var1(ivar1) = "tciw"

  ivar1 = ivar1+1
  VarList1(ivar1) = "PRECC"
  VstList1(ivar1) = "PRECC"
  VnmList1(ivar1) = "Convective precipitation"
  UnitList1(ivar1)= "mm day~S~-1~N~"
  ERA5Var1(ivar1) = "mcpr"

  ivar1 = ivar1+1
  VarList1(ivar1) = "PRECL"
  VstList1(ivar1) = "PRECL"
  VnmList1(ivar1) = "Stratiform precipitation"
  UnitList1(ivar1)= "mm day~S~-1~N~"
  ERA5Var1(ivar1) = "mlspr"

  ivar1 = ivar1+1
  VarList1(ivar1) = "PRECT"
  VstList1(ivar1) = "PRECT"
  VnmList1(ivar1) = "Total precipitation"
  UnitList1(ivar1)= "mm day~S~-1~N~"
  ERA5Var1(ivar1) = "mtpr"

