;*************************************************
; histo_4.ncl
;
; Concepts illustrated:
;   - Drawing a histogram using exact min/max of data
;   - Generating dummy data using "rand"
;   - Drawing X axis labels at a 45 degree angle
;
;************************************************
;
; This file is loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************

load "./function_read_dart_diag.ncl"

begin

  l_create_plot_data = True

  datadir = "/compyfs/zhan391/v3_dart_cda_scratch/diag_dart/"
  outdir  = "/compyfs/zhan391/v3_dart_cda_scratch/diag_dart/fig_data/"

  period     = "2011011500-2011012012" 
  Casename   = (/"ne11_oQU240","ne16_oQU240","ne30pg2_EC30to60E2r2"/)
  Expname    = (/"ne11(2.8deg)","ne16(1.9deg)","ne30(1deg)"/)
  ncases     = dimsizes(Casename)

  regions = (/"SH","NH","NAmer","Tropics"/)
  regstrs = (/"Southern Hemisphere","Northern Hemisphere","North America","Tropics"/)
  nregs      = dimsizes(regions)

  Metrics    = (/"bias", "RMSE"/)
  nmecs      = dimsizes(Metrics)
  
  ;;;read the variable information;;
  fvar         = "./variable.info" 
  data         = asciiread(fvar,-1,"string")
  ObsType      = str_squeeze(str_get_field(data, 1,","))
  VarName      = str_squeeze(str_get_field(data, 3,","))
  VarUnit      = str_squeeze(str_get_field(data, 4,","))
  xmin1        = stringtofloat(str_get_field(data, 5,","))
  xmax1        = stringtofloat(str_get_field(data, 6,","))
  xmin2        = stringtofloat(str_get_field(data, 7,","))
  xmax2        = stringtofloat(str_get_field(data, 8,","))
  fignum1      = str_squeeze(str_get_field(data, 9,","))
  fignum2      = str_squeeze(str_get_field(data, 10,","))
  nvars        = dimsizes(VarName)

  do irg = 0,nregs-1,1

  do iv = 0,nvars-1,1

   ObsTypeString = ObsType(iv)
   
   do icas = 0,ncases-1,1

    casestr      = "F20TR_"+Casename(icas)+"_DART80_compy"
    ModelOutRoot = datadir+casestr+"/dart_diagnostics/obs_diag/"
    fname        = ModelOutRoot+casestr+".dart.e.eam_obs_diag_output."+period+".nc"
    region       = regstrs(irg)
    type         = "guess"
  
    CopyString1 = "bias" 
    bias = read_dart_diag_netcdf(fname, CopyString1, ObsTypeString, region, type)
    
    CopyString2 = "rmse"
    rmse  = read_dart_diag_netcdf(fname, CopyString2, ObsTypeString, region, type)

    bias@_FillValue = -999.0
    rmse@_FillValue = -999.0

    mbias = dim_avg_n_Wrap(bias,0)
    mrmse = dim_avg_n_Wrap(rmse,0)
   
    mbias@_FillValue = bias@_FillValue
    mrmse@_FillValue = rmse@_FillValue

    if (str_squeeze(ObsTypeString) .eq. "RADIOSONDE_SPECIFIC_HUMIDITY") then
      print(mbias)
      print(mrmse)   
    end if 

    pltlev       = rmse@lev
    pltlev@units = rmse@lev_unit
    pltlev@range = rmse@levp
    pltlev!0     = "lev"

    if (icas.eq.0) then
      pltdat   = new((/2,ncases,dimsizes(pltlev)/),typeof(mrmse))
      pltdat@_FillValue = -999.0
    end if
   
    pltdat(0,icas,:)   = (/mbias(:)/)
    pltdat(1,icas,:)   = (/mrmse(:)/)
   end do

   pltdat!0  = "ncases"
   pltdat!1  = "nvar"
   pltdat!2  = "lev"
   pltdat&lev = pltlev

   ;;save the plot data;;;;;
   ff  = outdir+ObsTypeString+"_profile.nc"
   system("rm -f " + ff)
 
   print(" ")
   print("output : " + ff)
   print(" ")
 
   if (l_create_plot_data) then 
     flo = addfile(ff,"c")
     flo->plev            = pltlev
     flo->plev@units      = pltlev@units 
     flo->$ObsTypeString$ = pltdat
   end if 

   fno    = addfile(ff,"r")
   pltdat = fno->$ObsTypeString$
   pltlev = fno->plev
   mdat   = dim_avg_n(pltdat,2)

  ;----------------------------------------------------------------------
  ; Plotting section ;;;;

  if (iv.eq.0) then 
    plot = new((/2,nvars/),graphic)
    plot_name = "Fig_rmse_bias_prof_" + regions(irg) 
    print(" ")
    print("plot figure : " + plot_name+".pdf")
    print(" ")
    wtype = "pdf"
    wks = gsn_open_wks(wtype,plot_name)             ; send graphics to PNG file
    gsn_define_colormap(wks,"cb_9step")
    colors  = (/14,28,2/)
    markers = (/16, 4, 3, 6, 16/)
    mkszf   = (/0.7, 0.9, 0.8, 0.8/)*2.0
    dashind = (/0,0,0,0,0/)       ;;line dash pattern
    lnthick = (/3,3,3,3,3/)/1.5   ;;line thichness
    mkthk   = mkszf
    mkthk   = 1.0
  end if 
  print(" ")
  print("plot variable: " + ObsTypeString)
  print(" ")

;---Resources to share between both plots
  res          = True                   ; plot mods desired
  res@vpHeightF = 0.8
  res@vpWidthF =  0.5
  res@gsnDraw  = False
  res@gsnFrame = False

  FontHeight = 0.035
  res@gsnStringFontHeightF = FontHeight*1.1
  res@tmXBLabelFontHeightF = FontHeight
  res@tmYLLabelFontHeightF = FontHeight
  res@tmYRLabelFontHeightF = FontHeight
  res@tiMainFontHeightF    = FontHeight*1.1
  res@tiXAxisFontHeightF   = FontHeight
  res@tiYAxisFontHeightF   = FontHeight

  res@tmBorderThicknessF    = 1.5
  res@tmXBMajorThicknessF   = 1.5
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0

  res@trYReverse        = True                       ; reverse Y-axis
  res@tmYLMode          = "Explicit"                ; explicit labels
  res@tmYLValues        = (/ 1000,850,700,600,500,400,300,200,100,1/)
  res@tmYLLabels        = ""+res@tmYLValues         ; make strings
  res@tmXTOn            = False
  res@trYMaxF           = 1000.0
  res@trYMinF           = 100 
  res@tmYROn            = False
  
  res@gsnXRefLine             = (/0./)      ; four X reference lines
  res@gsnXRefLineThicknesses  = 1.0
  res@gsnXRefLineDashPatterns = 2
  res@gsnXRefLineColors       = "Black"

  res@xyMarkLineMode         = "MarkLines"
  res@xyLineColors           = colors
  res@xyMarkers              = markers
  res@xyMarkerColors         = colors
  res@xyMarkerSizes          = 0.015*mkszf
  res@xyMarkerThicknesses    = lnthick
  res@xyDashPatterns         = dashind
  res@xyLineThicknesses      = lnthick

  res@tiMainString           = ""
  res@tiMainOffsetYF         = 0.02
  res@tiMainFontThicknessF   = 3.0
  res@gsnRightString         = "" 
  res@gsnLeftString          = "" 
  res@tiYAxisString          = ""
  res@tiXAxisString          = ""

  res1 = res
  if(iv.eq.0) then
    res1@tmYLLabelsOn        = True      
    res1@tiYAxisString       = "Pressure (hPa)"
    res1@gsnLeftStringOrthogonalPosF = 0.0
  else
    res1@tmYLLabelsOn        = False
    res1@tiYAxisString       = ""
    res1@gsnLeftStringOrthogonalPosF = 0.02
  end if

  res1@tiXAxisString         = "("+VarUnit(iv)+")"
  res1@gsnLeftString         = fignum1(iv)+ " "+Metrics(0)+" in " + VarName(iv)
  res1@trXMaxF              := xmax1(iv)
  res1@trXMinF              := xmin1(iv)
;  res1@tmXBMode              = "Explicit"                ; explicit labels
;  res1@tmXBValues            = fspan(res1@trXMinF, res1@trXMaxF,7)
;  res1@tmXBMinorValues       = fspan(res1@trXMinF, res1@trXMaxF,13)
;  res1@tmXBLabels            = sprintf("%.1f",res1@tmXBValues)
  plot(0,iv)  = gsn_csm_xy (wks,pltdat(0,:,:),pltlev,res1)   ; create plot
 
  res2 = res
  if(iv.eq.0) then 
    res2@tmYLLabelsOn        = True 
    res2@tiYAxisString       = "Pressure (hPa)"
    res2@gsnLeftStringOrthogonalPosF = 0.0
  else
    res2@tmYLLabelsOn        = False
    res2@tiYAxisString       = ""
    res2@gsnLeftStringOrthogonalPosF = 0.02
  end if 
  res2@tiXAxisString         = "("+VarUnit(iv)+")"
  res2@gsnLeftString         = fignum2(iv)+ " "+Metrics(1)+" in " + VarName(iv)
  res2@trXMaxF              := xmax2(iv)
  res2@trXMinF              := xmin2(iv)
;  res2@tmXBMode              = "Explicit"                ; explicit labels
;  res2@tmXBValues            = fspan(res2@trXMinF, res2@trXMaxF,7)
;  res2@tmXBMinorValues       = fspan(res2@trXMinF, res2@trXMaxF,13)
;  res2@tmXBLabels            = sprintf("%.1f",res2@tmXBValues)
  plot(1,iv)  = gsn_csm_xy (wks,pltdat(1,:,:),pltlev,res2)   ; create plot

;Attach a legend
  lgres                    = True
  lgres@lgItemType         = "MarkLines"
  lgres@lgLineLabelFontThicknessF= 1.0
 ;lgres@lgLineDashSegLenF  = 0.01
  lgres@lgBoxMinorExtentF  = 0.25
  lgres@lgBoxMajorExtentF  = 0.1
  lgres@lgMonoDashIndex    = False
  lgres@lgMonoLineThickness= False
  lgres@lgMonoMarkerThickness = False
  lgres@lgPerimOn          = False ;True
 ;lgres@lgRightMarginF     = 0.5
 ;lgres@lgLeftMarginF      = 0.5
 ;lgres@lgTopMarginF       = 0.5
 ;lgres@lgBottomtMarginF   = 0.5
  lgres@lgPerimThicknessF  = 1.0

  lgres1 = lgres
  lgres1@vpWidthF           = 0.28           ; width of legend (NDC)
  lgres1@vpHeightF          = 0.28           ; height of legend (NDC)
  lgres1@lgLabelFontHeightF = FontHeight * 10.0
  lgres1@lgLineColors       = colors
  lgres1@lgDashIndexes      = dashind
  lgres1@lgLineThicknesses  = lnthick
  lgres1@lgMarkerColors     = colors
  lgres1@lgMarkerIndexes    = markers
  lgres1@lgMarkerThicknesses= mkthk

  ;Labels for the legend
  labls = Expname
  do ii= 0,ncases-1
   labls(ii) = "  " + Expname(ii) ;+ " ("+sprintf("%5.2f",mdat(0,ii))+")"
  end do
  legend1 = gsn_create_legend (wks,dimsizes(labls), labls,lgres1)

  labls = Expname
  do ii= 0,ncases-1
   labls(ii) = "  " + Expname(ii) ;+ " ("+sprintf("%5.2f",mdat(1,ii))+")"
  end do
  legend2 = gsn_create_legend (wks,dimsizes(labls), labls,lgres1)

  amres1                  = True
  amres1@amJust           = "TopRight"
  amres1@amParallelPosF   =   0.40
  amres1@amOrthogonalPosF =  -0.45
  if(iv.eq.1)then 
   annoid1 = gsn_add_annotation(plot(0,iv),legend1,amres1)  ; add legend to the top figure
   annoid2 = gsn_add_annotation(plot(1,iv),legend2,amres1)  ; add legend to the top figure
  end if 

 end do

;---Draw both plots in a panel
  pnres                         = True
  pnres@gsnMaximize             = False
  pnres@gsnPanelMainString      = ""
  pnres@gsnPanelMainFontColor   = "Black"
  pnres@gsnPanelMainFontHeightF = 0.020
  pnres@gsnPanelLabelBar        = False
  pnres@gsnPanelBottom          = 0.05
  pnres@gsnPanelYWhiteSpacePercent = 2
  pnres@gsnPanelXF = (/0.15,0.35,0.55,0.75/)  
  gsn_panel(wks,(/ndtooned(plot)/),(/2,4/),pnres)
  delete([/pnres/])

end do 

end
