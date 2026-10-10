;;;;;;;;;;the experiment names;;;;;;;;;;;;;
  nexps          = 7
  Exp_long_names  = new(nexps,string)
  Exp_short_names = new(nexps,string)

  iexp            = 0
  Exp_long_names(iexp)  = "v1_CTRL"
  Exp_short_names(iexp) = "Exp A"

  iexp            = iexp+1
  Exp_long_names(iexp)  = "v1_All_Shorter"
  Exp_short_names(iexp) = "Exp B"

  iexp            = iexp+1
  Exp_long_names(iexp)  = "v1_MacMic_Shorter"
  Exp_short_names(iexp) = "Exp C"

  iexp            = iexp+1
  Exp_long_names(iexp)  = "v1_All_Except_MacMic_Shorter"
  Exp_short_names(iexp) = "Exp D"

  iexp            = iexp+1
  Exp_long_names(iexp)  = "v1_CPL+DeepCu_Shorter"
  Exp_short_names(iexp) = "Exp F"

  iexp            = iexp+1
  Exp_long_names(iexp)  = "v1_Dribble"
  Exp_short_names(iexp) = "Exp G"

  iexp            = iexp+1
  Exp_long_names(iexp)  = "v1_CPL+DeepCu+Tau_Shorter"
  Exp_short_names(iexp) = "Exp F_tau"

