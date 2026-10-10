
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
  colors   = (/16,2,12,14,7,10/)
  BarChartFillScale = 0.3
  markers  = (/9,6,16,16,16,16/)
  mkszf    = (/0.8,0.8,0.6,0.6,0.6,0.6/)*0.03
  dashind  = (/0,0,1,1,1,1/)
  lnthick  = (/3.0,3.0,3.0,3.0,3.0,3.0/)*0.6
  llth     = (/3.0,3.0,3.0,3.0,3.0/)

