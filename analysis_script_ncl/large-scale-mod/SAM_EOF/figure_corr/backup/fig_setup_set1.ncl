
;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "pdf"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+plotFileLabel)

;;;;;;;;;;;;;;;;;;;;;;;;
 ;gsn_define_colormap(wks,"t2m_29lev")
 ;colors   = (/26,2,6,31/)
  gsn_define_colormap(wks,"srip_reanalysis")
  colors   = (/1,1,14,16,2/)
  BarChartFillScale = 0.3
  markers  = (/16,16,4,9,6/)
  mkszf    = (/0.6,0.6,0.8,0.8,0.8/)*0.03
  dashind  = (/0,2,0,0,0/)
  lnthick  = (/1.5,1.5,3.0,3.0,3.0/)*1.2
  llth     = (/3.0,3.0,3.0,3.0,3.0/)

