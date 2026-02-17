
;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "pdf"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+plotFileLabel)

;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "pdf"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+plotFileLabel)
  gsn_define_colormap(wks,"amwg")
;;;;;;;;;;;;;;;;;;;;;;;;
  gsn_define_colormap(wks,"seaice_2")
  gsn_define_colormap(wks,"StepSeq25")
  colors  = (/1,17,3,7/); (/"palegreen4","orange4","navy"/)
  colors  = (/1,17,3,12/); (/"palegreen4","orange4","navy"/)
  barcol  = (/14,25/)
  markers = (/16,12,11,8/)
  mkszf   = (/0.8,0.8,0.8,0.8,0.8/)
  dashind = (/0,0,0,0/)
  lnthick = (/10.0,10.0,10.0,1.0,1.0,1.0/) ;*2.0

