;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ldatalocal = True
  ldatalocal = False

  ObsRoot = "/pic/projects/uq_climate/zhan391/climate_analysis/obs_data_20140804/"

  if (ldatalocal) then
     ModelOutRoot = "./E3SMV1_SIMULATION/"
  else
     ModelOutRoot = "/pic/projects/uq_climate/zhan391/e3sm_cam4/run_amwg_10yr/climo/"
  end if
  
;;;;;;;;;weight files for ESMF regrid function;;;;;
  wgt_dir  = "./"
  wgt_file = wgt_dir+"SE_NE30_to_FV_192x288.bilinear.nc"


;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Global"/)
  nregs        = dimsizes(Regions)

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir       = (/"ql0_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_F_intel_constance_DT1800_fullmodel_rad", "ql1_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_F_intel_constance_DT1800_fullmodel_rad", "ql1_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_B3_F_intel_constance_DT1800_fullmodel_rad"/)
  Groups       = (/"ql0_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_F_intel_constance_DT1800_fullmodel_rad", "ql1_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_F_intel_constance_DT1800_fullmodel_rad", "ql1_fmin0_1.9x2.5_gx1v6_E3SM_MASTER_B3_F_intel_constance_DT1800_fullmodel_rad"/)
  ngrps        = dimsizes(Groups)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  labels       = (/"Baseline model","Revised splitting","Revised splitting & closure"/)
  nlabls       = dimsizes(labels)

;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Varname      = (/"CLDTOT"/)
  nvars        = dimsizes(Varname)  
  Varunt       = ("Total Cloud (%)")
  nvars        = dimsizes(Varname)
  Plttim       = (/"ANN"/)
  nptim        = dimsizes(Plttim)
  mainstr      = "Total Cloud Cover"

  plotFileLabel ="PD_zonal_mean_total_cloud"

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Filnob1 = systemfunc("cd "+ModelOutRoot+Casdir(0)+"/"+"; ls *"+Plttim(0)+"_mean*")
  Filob1  = addfile(ModelOutRoot+Casdir(0)+"/"+ Filnob1,"r")
  Filnam0 = systemfunc("cd "+ModelOutRoot+Casdir(1)+"/"+"; ls *"+Plttim(0)+"_mean*")
  File0   = addfile(ModelOutRoot+Casdir(1)+"/"+Filnam0,"r")
  Filnam1 = systemfunc("cd "+ModelOutRoot+Casdir(2)+"/"+"; ls *"+Plttim(0)+"_mean*")
  File1   = addfile(ModelOutRoot+Casdir(2)+"/"+Filnam1,"r")
 ;print(Filnam0)
 ;print(Filnam1)
 ;print(Filob)

  ;;read the observations;;;
  latob1     = Filob1->lat
  lonob1     = Filob1->lon
  nlatob1    = dimsizes(latob1)
  nlonob1    = dimsizes(lonob1)
  pcpob1     = Filob1->CLDTOT
  pcpob1     = pcpob1*100.0

  lat       = File0->lat
  lon       = File0->lon
  nlat      = dimsizes(lat)
  nlon      = dimsizes(lon)
  pcp1800   = File0->CLDTOT
  pcp0300   = File1->CLDTOT
  pcp1800   = pcp1800*100.0
  pcp0300   = pcp0300*100.0

;;;processe the zonal mean data;;;;
  plotobs1  = new((/1,nlatob1/),typeof(pcpob1)) 
  plot1800  = new((/1,nlat/),typeof(pcp1800)) 
  plot0300  = new((/1,nlat/),typeof(pcp0300)) 
 
  plotobs1(0,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(pcpob1,0),1)
  plot1800(0,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(pcp1800,0),1)
  plot0300(0,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(pcp0300,0),1)

  printVarSummary(plotobs1)
  printVarSummary(plot1800)
  printVarSummary(plot0300)

;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "png"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+"fig_1d_pcp_"+plotFileLabel)
  gsn_define_colormap(wks,"amwg")
;;;;;;;;;;;;;;;;;;;;;;;;
  gsn_define_colormap(wks,"seaice_2")
  gsn_define_colormap(wks,"StepSeq25")
  colors  =(/18,4,9/); (/"palegreen4","orange4","navy"/) 
  markers = (/12,11,8/)
  mkszf   = (/0.8,0.8,0.8,0.8/)
  dashind = (/0,0,2/) 
  lnthick = (/5.0,5.0,5.0,5.0,5.0/)
  plot    = new(3,graphic)                      

  res               = True           ; plot mods desired
  res@vpXF          = 0.17
  res@vpYF          = 0.8
  res@vpWidthF      = 0.8                      ; set width of plot
  res@vpHeightF     = 0.6                      ; set height of plot

  res@gsnDraw  = False
  res@gsnFrame = False

  res@tmYROn = False            ; Turn off right tickmarks.
  res@tmXTOn = False            ; Turn off top tickmarks
  FontHeight = 0.028
  res@gsnStringFontHeightF  = FontHeight*0.8
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight
  res@tiYAxisFontHeightF    = FontHeight
  res@tiMainFontThicknessF  = 8.0
  res@tiYAxisFontThicknessF = 5.0
  res@tiXAxisFontThicknessF = 5.0

  res@trXMinF               = -90 ;min(latrgd)
  res@trXMaxF               =  90 ;max(latrgd)
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = ispan(-90,90,30)
  res@tmXBLabels            = "" + res@tmXBValues
  res@tmXBMinorValues       = ispan(-90,90,10)
  
  res@gsnYRefLine           = -10.0             ; create a reference line   
  res@gsnYRefLineDashPattern= 2             
  res@gsnYRefLineColor      = "Black"
  res@gsnYRefLineThicknessF = 3.0
  
  res@gsnLeftString         ="" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor= "Black"
  res@gsnRightString         = "";"(PD)"
  res@gsnRightStringFontColor= "Black"
;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res0                      = res
  res0@tiMainString         = mainstr
  res0@tiXAxisString        = "Latitude (~S~o~N~N)"
  res0@tiYAxisString        = Varunt
  res0@trYMinF              = 0.0   ;floor(min(pltobs)) ;min(ynflx)-0.1 ;You could also just use
  res0@trYMaxF              = 100.0   ;ceil(max(pltobs)) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
  res0@xyMarkLineModes      = "Lines"
  res0@xyLineThicknessF     = 8.0            ; make thicker
  res0@xyDashPatterns       = dashind(0)
  res0@xyMarkerSizes        = 0.012
  res0@xyLineColors         = colors(0)
  res0@xyMarkers            = markers(0)
  res0@xyMarkerColors       = colors(0)
  plot(0)                   = gsn_csm_xy(wks,latob1,plotobs1,res0)

;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res1                      = res
  res1@tiMainString         = mainstr
  res1@tiXAxisString        = "Latitude (~S~o~N~N)"
  res1@tiYAxisString        = Varunt
  res1@trYMinF              = 0.0 ;floor(min(pltdat(:,:,0:nlat-1))) ;min(ynflx)-0.1 ;You could also just use
  res1@trYMaxF              = 100.0  ;ceil(max(pltdat(:,:,0:nlat-1))) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
  res1@xyMarkLineModes      = "Lines"
  res1@xyLineThicknessF     = 8.0            ; make thicker
  res1@xyDashPatterns       = dashind(1)
  res1@xyMarkerSizes        = 0.012
  res1@xyLineColors         = colors(1)
  res1@xyMarkers            = markers(1)
  res1@xyMarkerColors       = colors(1)
  plot(1)                   = gsn_csm_xy(wks,lat,plot1800,res1)

;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res2                      = res
  res2@tiMainString         = mainstr
  res2@tiXAxisString        = "Latitude (~S~o~N~N)"
  res2@tiYAxisString        = "W m~S~-2~N~";Varname(iv)+ " ("+var@units+")"
  res2@trYMinF              = 0.0 ;floor(min(pltdat(:,:,0:nlat-1))) ;min(ynflx)-0.1 ;You could also just use
  res2@trYMaxF              = 100.0  ;ceil(max(pltdat(:,:,0:nlat-1))) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
  res2@xyMarkLineModes      = "Lines"
  res2@xyLineThicknessF     = 8.0            ; make thicker
  res2@xyDashPatterns       = dashind(2)
  res2@xyMarkerSizes        = 0.012
  res2@xyLineColors         = colors(2)
  res2@xyMarkers            = markers(2)
  res2@xyMarkerColors       = colors(2)
  plot(2)                   = gsn_csm_xy(wks,lat,plot0300,res2)
 
  overlay(plot(0),plot(1))
  ;overlay(plot(0),plot(2))

;Attach a legend
;;manual legend

;  label1=(/"  SWCF", "  LWCF", "  NETCF"/)
;  dum1  = new(dimsizes(label1),graphic)
;
;  txx = new( dimsizes(label1),"float" )
;  txx = (/ -70, -70, -35/)
;  txy = (/ -70,  30, -20/)
;  txres = True
;  txres@txFontHeightF = FontHeight*0.7
;  txres@txJust = "CenterRight"
;  do iy= 0,dimsizes(label1)-1
;     txres@txFontColor = colors(iy)
;     txt = "~F22~"+label1(iy)
;;    txt = txt +" ("+sprintf("%3.1f",slope(iy))+")~F~"
;     dum1(iy) = gsn_add_text(wks,plot(0), txt, txx(iy),txy(iy),txres)
;  end do


   lgres                    = True
   lgres@lgItemType         = "Lines"
   lgres@lgLineLabelFontThicknessF= 3    
  ;lgres@lgLineDashSegLenF  = 0.4
   lgres@vpWidthF           = 0.30           ; width of legend (NDC)
   lgres@vpHeightF          = 0.18           ; height of legend (NDC)
   lgres@lgLabelFontHeightF = .25            ; legend label font thickness
   lgres@lgMonoDashIndex    = False
   lgres@lgMonoLineThickness= False
   lgres@lgPerimOn          = False
   ;lgres@lgRightMarginF     = 0.5
   ;lgres@lgLeftMarginF      = 0.5
   ;lgres@lgTopMarginF       = 0.5
   ;lgres@lgBottomtMarginF   = 0.5

   lgres@lgPerimThicknessF  = 2.0 

   lgres1 = lgres
   lgres1@lgLineColors       = colors ;(/"Black","Black","Black"/)
   lgres1@lgDashIndexes      = (/0,0,1,1/)
   lgres1@lgLineThicknesses  = (/8.0,8.0,8.0,8.0/)
   label1  = (/"  "+labels(0),"  "+labels(1)/)
   legend1 = gsn_create_legend (wks,dimsizes(label1), label1,lgres1)

;;;;;;;;;;;;
  amres                  = True
  amres@amJust           = "BottomLeft"    ; Use bottom right corner of box
  amres@amParallelPosF   =  -0.45             ; Move legend to right
  amres@amOrthogonalPosF =   0.45              ; Move legend down.
  annoid0 = gsn_add_annotation(plot(0),legend1,amres)  ; add legend to the top figure

  draw(plot(0))
  frame(wks)
;
;  panres = True                                         ; panel resource list
;  panres@gsnMaximize                = True
;  panres@gsnPaperMargin             = 0.1
;  panres@gsnPanelBottom             = 0.05     ; move bottom up so size is 0.3
;  panres@gsnPanelYWhiteSpacePercent = 0.8       ;set spacing vertically between 2 panels
;  gsn_panel(wks,plot,(/1,1/),panres)
;  frame(wks)
end

