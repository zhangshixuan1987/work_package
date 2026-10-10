
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
  colors   = (/20,14,18,12/)
 ;gsn_define_colormap(wks,"cb_9step")
 ;colors   = (/3,10,43,51/)

  BarChartFillScale = 0.3

  barWidth = (/0.8,0.8,0.8,0.8/)/1.5
  barclr   = (/4,-1,-1,-1/)
  markers  = (/11,8,7,6/)
  mkszf    = (/0.8,0.8,0.8,0.8/)
  dashind  = (/0,0,0,0/)
  lnthick  = (/5.0,5.0,5.0,5.0/)
  llth     = (/3.0,3.0,3.0,3.0/)
  bopc     = (/0.5,0.8,0.8,0.8/)
  pattns   = (/0,7,0,0,0,0/)

