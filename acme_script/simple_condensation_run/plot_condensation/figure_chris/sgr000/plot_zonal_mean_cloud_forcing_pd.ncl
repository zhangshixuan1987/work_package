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
  labels       = (/"Original","B1+Revised splitting","B3+Revised splitting"/)
  nlabls       = dimsizes(labels)

;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Varname      = (/"LWCF","SWCF"/)
  nvars        = dimsizes(Varname)  
  Varunt       = ("Cloud forcing (W m~S~-2~N~)")
  nvars        = dimsizes(Varname)
  Plttim       = (/"ANN"/)
  nptim        = dimsizes(Plttim)

  plotFileLabel ="PD_zonal_mean_cloud_forcing"

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Filnob  = systemfunc("cd "+ModelOutRoot+Casdir(0)+"/"+"; ls *"+Plttim(0)+"_mean*")
  Filob   = addfile(ModelOutRoot+Casdir(0)+"/"+ Filnob,"r")
  Filnam0 = systemfunc("cd "+ModelOutRoot+Casdir(1)+"/"+"; ls *"+Plttim(0)+"_mean*")
  File0   = addfile(ModelOutRoot+Casdir(1)+"/"+Filnam0,"r")
  Filnam1 = systemfunc("cd "+ModelOutRoot+Casdir(2)+"/"+"; ls *"+Plttim(0)+"_mean*")
  File1   = addfile(ModelOutRoot+Casdir(2)+"/"+Filnam1,"r")
 ;print(Filnam0)
 ;print(Filnam1)
 ;print(Filob)

  ;;read the observations;;;
  latob     = Filob->lat
  lonob     = Filob->lon
  nlatob    = dimsizes(latob)
  nlonob    = dimsizes(lonob)
  swcfob    = Filob->SWCF
  lwcfob    = Filob->LWCF
  necfob    = swcfob
  necfob    = swcfob+lwcfob

  lat       = File0->lat
  lon       = File0->lon
  nlat      = dimsizes(lat)
  nlon      = dimsizes(lon)
  swcf1800  = File0->SWCF
  lwcf1800  = File0->LWCF
  necf1800  = swcf1800+lwcf1800
  swcf0300  = File1->SWCF
  lwcf0300  = File1->LWCF
  necf0300  = swcf0300
  necf0300  = swcf0300+lwcf0300

;;;processe the zonal mean data;;;;
  plotobs   = new((/3,nlatob/),typeof(swcfob)) 
  plot1800  = new((/3,nlat/),typeof(swcf1800)) 
  plot0300  = new((/3,nlat/),typeof(swcf0300)) 
 
  plotobs(0,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(swcfob,0),1)
  plotobs(1,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(lwcfob,0),1)
  plotobs(2,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(necfob,0),1)

  plot1800(0,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(swcf1800,0),1)
  plot1800(1,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(lwcf1800,0),1)
  plot1800(2,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(necf1800,0),1)

  plot0300(0,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(swcf0300,0),1)
  plot0300(1,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(lwcf0300,0),1)
  plot0300(2,:) = dim_avg_n_Wrap(dim_avg_n_Wrap(necf0300,0),1)

  printVarSummary(plotobs)
  printVarSummary(plot1800)
  printVarSummary(plot0300)

;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "png"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+"fig_1d_cre_"+plotFileLabel)
;;;;;;;;;;;;;;;;;;;;;;;;
  colors  = (/"palegreen4","orange4","navy","palegreen4","orange4","navy"/)
  markers = (/16,11,8,16,11,8/)
  mkszf   = (/0.8,0.8,0.8,0.8/)
  dashind = (/0,0,0,1,1,1/) 
  lnthick = (/5.0,5.0,5.0,5.0,5.0,5.0/)
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
  FontHeight = 0.030
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
  
  res@gsnYRefLine           = 0.0             ; create a reference line   
  res@gsnYRefLineDashPattern= 2             
  res@gsnYRefLineColor      = "Black"
  res@gsnYRefLineThicknessF = 3.0
  
  res@gsnLeftString         ="" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor= "Black"
  res@gsnRightString         = "";"(PD)"
  res@gsnRightStringFontColor= "Black"
;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res0                      = res
  res0@tiMainString         = "Cloud Redative Forcing" ;"Maximum Surface Wind"
  res0@tiXAxisString        = "Latitude [deg N]"
  res0@tiYAxisString        = "W m~S~-2~N~";Varname(iv)+ " ("+var@units+")"
  res0@trYMinF              = -110.0 ;floor(min(pltobs)) ;min(ynflx)-0.1 ;You could also just use
  res0@trYMaxF              = 60.0   ;ceil(max(pltobs)) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
  res0@xyMarkLineModes      = "Lines"
  res0@xyLineThicknessF     = 8.0            ; make thicker
  res0@xyDashPatterns       = dashind(0:nvars)
  res0@xyMarkerSizes        = 0.012
  res0@xyLineColors         = colors(0:nvars)
  res0@xyMarkers            = markers(0:nvars)
  res0@xyMarkerColors       = colors(0:nvars)
  plot(0)                   = gsn_csm_xy(wks,latob,plotobs,res0)


;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res1                      = res
  res1@tiMainString         = "Cloud Redative Effect (CRE)" ;"Maximum Surface Wind"
  res1@tiXAxisString        = "Latitude [deg N]"
  res1@tiYAxisString        = "W m~S~-2~N~";Varname(iv)+ " ("+var@units+")"
  res1@trYMinF              = -110.0 ;floor(min(pltdat(:,:,0:nlat-1))) ;min(ynflx)-0.1 ;You could also just use
  res1@trYMaxF              =  80.0  ;ceil(max(pltdat(:,:,0:nlat-1))) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
  res1@xyMarkLineModes      = "Lines"
  res1@xyLineThicknessF     = 3.0            ; make thicker
  res1@xyDashPatterns       = dashind(0:nvars)
  res1@xyMarkerSizes        = 0.012
  res1@xyLineColors         = colors(0:nvars)
  res1@xyMarkers            = markers(0:nvars)
  res1@xyMarkerColors       = colors(0:nvars)
  plot(1)                   = gsn_csm_xy(wks,lat,plot1800,res1)
;;;;;;;first figure for PD origional values;;;;;;;;;;;;;;;;;;;;
  res2                      = res
  res2@tiMainString         = "Cloud Redative Effect (CRE)" ;"Maximum Surface Wind"
  res2@tiXAxisString        = "Latitude [deg N]"
  res2@tiYAxisString        = "W m~S~-2~N~";Varname(iv)+ " ("+var@units+")"
  res2@trYMinF              = -110.0 ;floor(min(pltdat(:,:,0:nlat-1))) ;min(ynflx)-0.1 ;You could also just use
  res2@trYMaxF              = 60.0  ;ceil(max(pltdat(:,:,0:nlat-1))) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
  res2@xyMarkLineModes      = "Lines"
  res2@xyLineThicknessF     = 3.0            ; make thicker
  res2@xyDashPatterns       = dashind(nvars+1:)
  res2@xyMarkerSizes        = 0.012
  res2@xyLineColors         = colors(nvars+1:)
  res2@xyMarkers            = markers(nvars+1:)
  res2@xyMarkerColors       = colors(nvars+1:)
  plot(2)                   = gsn_csm_xy(wks,lat,plot0300,res2)
 
  overlay(plot(0),plot(1))
  overlay(plot(0),plot(2))

;Attach a legend
;;manual legend

  label1=(/"  SWCF", "  LWCF", "  NETCF"/)
  dum1  = new(dimsizes(label1),graphic)

  txx = new( dimsizes(label1),"float" )
  txx = (/ -60, -70, -35/)
  txy = (/ -80,  30, -20/)
  txres = True
  txres@txFontHeightF = FontHeight*0.7
  txres@txJust = "CenterRight"
  do iy= 0,dimsizes(label1)-1
     txres@txFontColor = res1@xyLineColors(iy)
     txt = "~F22~"+label1(iy)
;    txt = txt +" ("+sprintf("%3.1f",slope(iy))+")~F~"
     dum1(iy) = gsn_add_text(wks,plot(0), txt, txx(iy),txy(iy),txres)
  end do


   lgres                    = True
   lgres@lgItemType         = "Lines"
   lgres@lgLineLabelFontThicknessF= 3    
   lgres@vpWidthF           = 0.32           ; width of legend (NDC)
   lgres@vpHeightF          = 0.16           ; height of legend (NDC)
   lgres@lgLabelFontHeightF = .18            ; legend label font thickness
   lgres@lgMonoDashIndex    = False
   lgres@lgMonoLineThickness= False
   lgres@lgPerimOn          = False
   lgres@lgPerimThicknessF  = 2.0 

   lgres1 = lgres
   lgres1@lgLineColors       = (/"Black","Black","Black"/)
   lgres1@lgDashIndexes      = (/0,0,1/)
   lgres1@lgLineThicknesses  = (/8.0,3.0,3.0/)
   label1  = (/"  "+labels(0),"  "+labels(1),"  "+labels(2)/)
   legend1 = gsn_create_legend (wks,dimsizes(label1), label1,lgres1)

;;;;;;;;;;;;
  amres                  = True
  amres@amJust           = "BottomRight"    ; Use bottom right corner of box
  amres@amParallelPosF   =   0.48             ; Move legend to right
  amres@amOrthogonalPosF =   0.5              ; Move legend down.
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

