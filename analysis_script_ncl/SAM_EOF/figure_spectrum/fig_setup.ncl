
;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "pdf"
  bthick = 1.0
  if(wtype.eq."png") then 
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    bthick = 3.0
  end if
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+plotFileLabel)

;;;;;;;;;;;;;;;;;;;;;;;;
  gsn_define_colormap(wks,"seaice_2")
  gsn_define_colormap(wks,"StepSeq25")
  gsn_define_colormap(wks,"t2m_29lev")
 ;colors  = (/1,3,12,17/); (/"palegreen4","orange4","navy"/)
  colors  = (/1,27,10,25,18/);
  barcol  = (/14,25/)
  markers = (/16,12,11,8/)
  mkszf   = (/0.8,0.8,0.8,0.8/)
  dashind = (/0,0,0,0/)
  lnthick = (/1.0,2.0,1.5,1.5/)

