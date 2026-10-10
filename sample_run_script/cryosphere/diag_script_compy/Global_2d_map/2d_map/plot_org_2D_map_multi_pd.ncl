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

  if (ldatalocal) then
     ObsDatRoot   = "./ceres_ebaf/"
     ModelOutRoot = "./"
  else
     AnlDatRoot   = "/compyfs/zhan391/e3sm_dart/acme_init/Observations/"
                    ;/compyfs/zhan391/e3sm_dart/acme_init/era5_reanalysis/"
     ObsDatRoot   = "/compyfs/zhan391/e3sm_dart/acme_init/Observations/"
     ModelOutRoot = "/compyfs/zhan391/e3sm_dart/post/"
  end if

;;;;;;;;;weight files for ESMF regrid function;;;;;
  wgt_dir  = "./"
  wgt_file = wgt_dir+"SE_NE30_to_FV_192x288.bilinear.nc"

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"SEP"/)
  nregs        = dimsizes(Regions)
;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir       = (/"CERES_EBAF","ERA5","PDNDG_ENS1_F20TRC5-AMIP_f19_d025","PDNDGUVT_ENS1_F20TRC5-AMIP_f19_d025","PDDAUV_ENS80_F20TRC5-AMIP_f19_d025_201101","PDDA_ENS80_F20TRC5-AMIP_f19_d025_201101","PDDAUVTQ_ENS80_F20TRC5-AMIP_f19_d025_201101"/)
  Groups       = (/"CERES-EBAF","ERA5","NDG_ERA5UV","NDG_ERA5UVT","DARTEDA-UV","DARTEDA-UVT","DARTEDA-UVTQ"/)
  ngrps        = dimsizes(Groups)

  Plttim       = (/"CERES_EBAF-TOA_Ed4.1_f1.9x2.5_201101.nc","era5_an_sig_levl_f1.9x2.5_201101.nc","*mean.cam.h0.2011-01.nc*","*mean.cam.h0.2011-01.nc*","*mean.cam.h0.2011-01.nc*","*mean.cam.h0.2011-01.nc*","*mean.cam.h0.2011-01.nc*"/)
  nptim        = dimsizes(Plttim)

  LabList      = Groups 

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList      = (/"CLDTOT","SWCF","LWCF","NETCF","FSUTOA","FLUT"/)
  VstList      = (/"CLDTOT","SW CRE", "LW CRE", "CRE", "FSUTOA", "FLUT"/)
  UnitList     = (/"fraction","W m~S~-2~N~","W m~S~-2~N~","W m~S~-2~N~","W m~S~-2~N~","W m~S~-2~N~"/)
  PltNums      = (/"", "", "", "", "", "", "", "", "", "", ""/)
  nvars        = dimsizes(VarList)

 ;;; flags to control the plot style
  plotcntl          = True   ; if True, plot the CNTL origianl field on the leftmost panel

  plotrdif          = False  ; True  ; if True, plot the relative difference with respect to CNTL

  plot_ttest_mskout = False  ; if True, mask the non-significant region as missing values

  plot_ttest_hatch  = False  ; if True, plot the significant region as hatch region

  plot_ttest_hatrv  = False  ; if True, plot the non-significant region as hatch region

  plot_contour      = False  ; True   ; if True, plot the contour in the figure

  plot_leftstr      = True   ; if True, the lefadd igure captions as left string

  plot_lndmsk       = False ;
 ;;;these flags control the zonal mean plot
  plot_zonal_mean   = False  ; whether or not plot zonal mean
  add_fill_line     = False  ; Turn off (on if True) fill color below the zonal mean line
  add_std_shading   = False ; If true, the grey shading is added for uncertainty of zonal mean
  add_std_line      = False  ; If true, the dash red lines are added for upper and bottom bounds
                            ; of uncertainty for zonal mean

 ;;color for the map boundary line
  lndcolor = "Black" ;"Grey50"

 ;;flags to perform significant test;;;;;;;;
  l_studentttest    = False ; True if perform student-t test
  l_pairedttest     = False  ; True if perform parired-t test
  l_kstest          = False  ; True if perform K-S test, see https://www.ncl.ucar.edu/Document/Functions/Built-in/kolsm2_n.shtml

 plot         = new(ngrps,graphic)

 do iv = 0,nvars-1, 1

  Varname = VarList(iv)
  Varunt  = UnitList(iv)
  Varstr  = VstList(iv)
  Figstr  = PltNums(iv)

  plotFileLabel = Varname+"_PD"

 ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "png"
  wtype@wkWidth = 2500
  wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_org_"+plotFileLabel)

 do igp = 1, ngrps-1,1

  refcase = Casdir (0)
  tescase = Casdir (igp)
  
;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Filnam0 = systemfunc("cd "+ObsDatRoot+refcase+"/REGRID/"+"; ls -1 *"+Plttim(0))
  File0   = addfile(ObsDatRoot+refcase+"/REGRID/"+Filnam0,"r")
  if(igp.eq.1)then
    Filnam1 = systemfunc("cd "+AnlDatRoot+tescase+"/monthly_mean/"+"; ls -1 *"+Plttim(igp))
    File1   = addfile(AnlDatRoot+tescase+"/monthly_mean/"+Filnam1,"r")  
  else
    Filnam1 = systemfunc("cd "+ModelOutRoot+tescase+"/"+"; ls -1 *"+Plttim(igp))
    File1   = addfile(ModelOutRoot+tescase+"/"+Filnam1,"r")
  end if 
  ;ListSetType (File0, "join")
  ;ListSetType (File1, "join")

  lat       = File0->lat
  lon       = File0->lon

  if(Varname.eq."CLDTOT")then
   obsval = File0->cldarea_total_daynight_mon
   obsval = obsval/100.0 ; convert to fraction
  else if(Varname.eq."SWCF")
   obsval = -File0->toa_sw_all_mon
   obsval = -File0->toa_sw_all_mon + File0->toa_sw_clr_c_mon
  else if(Varname.eq."LWCF")
   obsval = - File0->toa_lw_all_mon
   obsval = - File0->toa_lw_all_mon + File0->toa_lw_clr_c_mon
  else if(Varname.eq."NETCF")
   obsval = -File0->toa_sw_all_mon
   obsval = -File0->toa_sw_all_mon + File0->toa_sw_clr_c_mon - File0->toa_lw_all_mon + File0->toa_lw_clr_c_mon
  else if(Varname.eq."FSUTOA")
   obsval = File0->toa_sw_all_mon 
  else if(Varname.eq."FLUT")
   obsval = File0->toa_lw_all_mon 
  end if 
  end if
  end if
  end if
  end if
  end if 

  if(igp.eq.1)then

  if(Varname.eq."CLDTOT")then
   anlval = File1->tcc
   anlval = anlval*1.0 ;/100.0 ; convert to fraction
  else if(Varname.eq."SWCF")
   anlval = File1->tsr;-File1->mtnswrfcs
   anlval = File1->tsr - File1->tsrc ;-File1->mtnswrfcs + File1->mtnswrf
   anlval = anlval/86400.0
  else if(Varname.eq."LWCF")
   anlval = File1->ttr ;-File1->mtnlwrfcs
   anlval = File1->ttr - File1->ttrc ;-File1->mtnlwrfcs + File1->mtnlwrf
   anlval = anlval/86400.0
  else if(Varname.eq."NETCF")
   anlval = File1->ttr ;-File1->mtnswrfcs
   anlval = File1->ttr - File1->ttrc  + File1->tsr - File1->tsrc;-File1->mtnswrfcs + File1->mtnswrf - File1->mtnlwrfcs + File1->mtnlwrf
  anlval = anlval/86400.0
  else if(Varname.eq."FSUTOA")
   anlval = File1->tisr-File1->tsr 
   anlval = anlval/86400.0
  else if(Varname.eq."FLUT")
   anlval = -File1->ttr 
   anlval = anlval/86400.0
  end if
  end if
  end if
  end if
  end if
  end if
  
  end if 

  mobsval = dim_avg_n_Wrap(obsval,0) 
  if (igp.eq.1)then 
   opxcre0300 = obsval
   owxcre0300 = anlval 
  else
   opxcre0300 = obsval
   if(Varname.eq."PRECT") then
    owxcre0300 = File1->PRECC
    owxcre0300 = File1->PRECC + File1->PRECL  
   else if(Varname.eq."NETCF") then
    owxcre0300 = File1->SWCF
    owxcre0300 = File1->LWCF + File1->SWCF
   else
    owxcre0300 = File1->$Varname$
   end if
   end if 
  end if 
  copy_VarCoords(opxcre0300,owxcre0300)

  printVarSummary(opxcre0300)
  printVarSummary(owxcre0300)

  if (Varname .eq. "TS") then
   opxcre0300 = opxcre0300 - 273.15
   owxcre0300 = owxcre0300 - 273.15
  end if

  if (Varname .eq. "OMEGA500") then
   opxcre0300 = opxcre0300*24.0*3600.0/100.0
   owxcre0300 = owxcre0300*24.0*3600.0/100.0
  end if

  if (Varname .eq. "TGCLDLWP" .or. Varname .eq. "TGCLDIWP") then
   opxcre0300 = opxcre0300*1000.0  ;; g/m2
   owxcre0300 = owxcre0300*1000.0
  end if

  if (Varname .eq. "PRECT" .or. Varname .eq. "PRECC" .or. Varname .eq. "PRECL") then
   opxcre0300 = opxcre0300*1000.0 *24.0 *3600.0
   owxcre0300 = owxcre0300*1000.0 *24.0 *3600.0 ;; mm/day
  end if

  printVarSummary(opxcre0300)

  wxcre0300  = owxcre0300 ; dim_avg_n_Wrap(owxcre0300,1)
  pxcre0300  = opxcre0300 ; wxcre0300
  pxcre0300  = opxcre0300 ; conform(wxcre0300,mobsval,(/1,2/)) ;dim_avg_n_Wrap(opxcre0300,1)

  printVarSummary(pxcre0300)
  print(min(pxcre0300)+ " "+max(pxcre0300))
  print(min(wxcre0300)+ " "+max(wxcre0300))

  pncre0300  = dim_avg_n_Wrap(pxcre0300,0)
  wncre0300  = dim_avg_n_Wrap(wxcre0300,0)
  printVarSummary(opxcre0300)
  delete([/opxcre0300,owxcre0300/])
  ;;;mask out the data over land
  ;pncre0300 = where(lndmsk.gt.0.0,32766,pncre0300)
  ;wncre0300 = where(lndmsk.gt.0.0,32766,wncre0300)
  ;pncre0300@_FillValue = 32766
  ;wncre0300@_FillValue = 32766
  ;print(num(ismissing(pncre0300)))
 
 ;;calculate the differences
  wncredf1  = wncre0300
  wncredf1  = wncre0300 - pncre0300

  if (plotrdif)then
   mskxx     = wncre0300
   mskxx     = where(wncre0300.ne.0.0,wncre0300,32766)
   mskxx@_FillValue = 32766
   wncredf1  = wncre0300
   wncredf1  = (wncre0300 - pncre0300)*100/mskxx
   wncredf1  = where(ismissing(wncredf1),0.0,wncredf1)
  end if 

  ;;;;find out the regions for the average;;;;;
  rad     = 4.*atan(1.)/180.    ;
  cost    = lat
  cost    = cos(lat*rad)       ; cosine weights
 ;mean1x  = wgt_areaave_Wrap( pncre0300({lat1:lat2},{lon1:lon2}), cost({lat1:lat2}), 1.0, 0 )
 ;mean2x  = wgt_areaave_Wrap( wncredf1({lat1:lat2},{lon1:lon2}), cost({lat1:lat2}), 1.0, 0 )
  meanob  = wgt_areaave_Wrap( mobsval, cost, 1.0, 0 )
  mean1x  = wgt_areaave_Wrap( wncre0300, cost, 1.0, 0 )
  mean2x  = wgt_areaave_Wrap( wncredf1,  cost, 1.0, 0 )
  mean1   = sprintf("%5.2f", mean1x)
  mean2   = sprintf("%5.2f", mean2x)
  rmse2x  = wgt_arearmse(wncre0300,pncre0300, cost, 1.0, 0 )
  rmse2   = sprintf("%5.2f", rmse2x)

  sigr = 0.05                        ; critical sig lvl for r

  if (l_studentttest) then
    rad  = 4.*atan(1.)/180.    ;
    wgty = cos(lat*rad)
    xtmp = pxcre0300(lat|:,lon|:,ncl_join|:)
    ytmp = wxcre0300(lat|:,lon|:,ncl_join|:)
    xAve = dim_avg (xtmp)              ; calculate means at each grid point
    yAve = dim_avg (ytmp)
    xVar = dim_variance (xtmp)         ; calculate variances
    yVar = dim_variance (ytmp)
    xEqv    = equiv_sample_size (xtmp, sigr,0)
    yEqv    = equiv_sample_size (ytmp, sigr,0)
    xN      = dim_avg_wgt_n(xEqv, wgty, 1, 0)    ; wgty could be gaussian weights
    yN      = dim_avg_wgt_n(yEqv, wgty, 1, 0)
    iflag   = False                        ; population variance similar
    prob    = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)
    delete([/xtmp,ytmp,xAve,yAve,xVar,yVar,xEqv,yEqv,xN,yN/])
  end if

  if (l_pairedttest) then
    xtmp     = pxcre0300(lat|:,lon|:,ncl_join|:)
    ytmp     = wxcre0300(lat|:,lon|:,ncl_join|:)
    nexpmembers = dimsizes(xtmp(0,0,:))
    vdiff    = ytmp - xtmp
    print(min(vdiff)+ " " +max(vdiff))
    diff_avg = dim_avg(vdiff)
    diff_std = dim_stddev(vdiff)
    diff_std@_FillValue = -999.
    print(min(diff_std)+ " " +max(diff_std))
    print(min(abs(diff_avg))+ " " +max(abs(diff_avg)))
    allzeros = diff_std.eq.0. ;(diff_avg.le.1e-3) .and. (diff_std.eq.0.)
    diff_std_tmp = where( .not.allzeros, diff_std, diff_std@_FillValue )
    printVarSummary(diff_std_tmp)

    tval   = where( .not.allzeros, diff_avg*sqrt(1.*nexpmembers)/diff_std_tmp, 0.)
    df     = nexpmembers - 1
    paramA = conform(tval, df/2.0, 0) ; see example 2 at https://www.ncl.ucar.edu/Document/Functions/Built-in/betainc.shtml
    paramB = conform(tval, 0.5, 0)     ; see example 2 at https://www.ncl.ucar.edu/Document/Functions/Built-in/betainc.shtml
    ;printVarSummary(tval)
    ;print(min(paramA) + " "+max(paramA))
    ;print(min(paramB) + " "+max(paramB))
    prob = where( .not.allzeros, betainc( df/(df+tval^2), paramA, paramB), 2. )
   ;prob = prob*0.5 ; one-sided
    delete([/xtmp,ytmp,nexpmembers,vdiff,diff_avg,diff_std,allzeros,diff_std_tmp,tval,df,paramA,paramB/])
  end if

  if (l_kstest) then
    xtmp     = pxcre0300 ;(lat|:,lon|:,ncl_join|:)
    ytmp     = wxcre0300 ;(lat|:,lon|:,ncl_join|:)
    prob     = kolsm2_n(pxcre0300, wxcre0300,0)
    delete([/xtmp,ytmp/])
  end if

  if(.not.l_pairedttest) then
    prob = (/wncredf1/)
    prob = (/sigr/)
  end if 

  v1mnts1 = wncredf1 ;smth9_Wrap(wncredf2, 0.50,  0.25, True)
  v1mnts1 = where(prob.le.sigr,wncredf1, 32766)
  if(plot_ttest_hatrv) then
    v1mnts1 = wncredf1
    v1mnts1 = where(prob.gt.sigr,wncredf1, 32766)
  end if
  v1mnts1@_FillValue = 32766
  delete([/prob/])

;;;processe the zonal mean data;;;;
 if (plot_zonal_mean) then
  nlat = dimsizes(lat)
  plotdiff      = new((/3,nlat/),typeof(pxcre0300))
  vardiff       = wxcre0300 - pxcre0300
  zmean         = dim_avg_n_Wrap(vardiff,2)
  zrmsd         = sqrt(dim_variance_n_Wrap(zmean,0))
  plotdiff(0,:) = dim_avg_n_Wrap(zmean,0)
  plotdiff(1,:) = plotdiff(0,:) - zrmsd
  plotdiff(2,:) = plotdiff(0,:) + zrmsd
 end if

 if(plot_lndmsk)then
  wncredf1  = where(lndmsk.eq.1,32766,wncredf1)
  pncre0300 = where(lndmsk.eq.1,32766,pncre0300)
  wncre0300 = where(lndmsk.eq.1,32766,wncre0300)
  wncredf1@_FillValue = 32766
  pncre0300@_FillValue = 32766
  wncre0300@_FillValue = 32766
 end if 


  load "./var_share_colorbar.ncl"

;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False
    res@tiMainFontHeightF    = 0.018
    res@tmYLLabelFontHeightF = 0.022
    res@tmXBLabelFontHeightF = 0.022

    mpres                    = res

    res@lbLabelFontHeightF   = 0.018
    res@gsnStringFontHeightF = 0.028
    res@cnFillOn             = True              ; color plot desired
   ;res@cnFillPalette        = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn            = False             ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour labels
   ;res@cnFillMode           = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor = "White"

    res@lbLabelBarOn            = True      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF  = 0.55
    res@pmLabelBarOrthogonalPosF= 0.15
    res@lbOrientation           = "Horizontal";"vertical"        ; vertical label bar
    res@pmLabelBarHeightF       = 0.2 ;0.65
    res@pmLabelBarWidthF        = 0.6 ;0.1
    ;res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
    ;res@lbTitlePosition         = "Right"                           ; title location
    ;res@lbTitleDirection        = "Across"                          ; letter angle
    ;res@lbTitleAngleF           = 90.                               ; title angle
    ;res@lbTitleFontHeightF      = 0.022                              ; font height


    leftstr                  = "" ;"Surface Temperature ("+ts@units+")"
    res@gsnLeftString        = "" ; leftstr ;var_name0    ; long_name is too long!
    res@gsnRightString       = "" ;" ("+Varunt+")";""
 ;  res@tiMaingString        = "Annual Mean Response to Present-day climate"

;---Control appearance of map.
    mpres@gsnMajorLatSpacing     = 30 ;20
    mpres@gsnMajorLonSpacing     = 45 ;30
    mpres@gsnMinorLatSpacing     = 15 ;10
    mpres@gsnMinorLonSpacing     = 15 ;15
   ;mpres@mpProjection           = "Robinson"       ; choose projection
    mpres@mpLabelsOn             = False
    mpres@mpPerimOn              = False
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor = lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 1.0          ; thickness of outlines
   ;mpres@mpOutlineOn            = True
   ;mpres@mpOutlineDrawOrder     = "PostDraw"
   ;mpres@mpFillDrawOrder        = "Predraw"
   ;mpres@mpGridAndLimbOn        = True
   ;mpres@mpGridLatSpacingF      = 90.0
   ;mpres@mpGridLonSpacingF      = 180.0
   ;mpres@mpGridLineThicknessF   = 3.0
   ;mpres@mpGridAndLimbDrawOrder="PreDraw"

    ;mpres@mpMinLatF              = mpMinLatF
    ;mpres@mpMaxLatF              = mpMaxLatF
    ;mpres@mpMinLonF              = mpMinLonF
    ;mpres@mpMaxLonF              = mpMaxLonF
    ;mpres@mpCenterLonF           = mpCenterLonF
    mpres@mpCenterLonF           = 240 ;150.0
;   mpres@mpCenterLatF           = -90.     ; This is necessary to get the correct map
;   mpres@mpOceanFillColor       = "lightskyblue1"
;   mpres@mpLandFillColor        = "gray33"

;---Resources for plotting original data
    if (igp.eq.1)then
     if(plotcntl)then
      mpres@tmYLLabelsOn           = True
      mpres@tmYLOn                 = True
     else
      mpres@tmYLLabelsOn           = False
      mpres@tmYLOn                 = False
     end if 
    else
     mpres@tmYLLabelsOn           = False
     mpres@tmYLOn                 = False
    end if 
    plot(igp) = gsn_csm_map(wks1,mpres)

   if (plotcntl.and.igp.eq.1)then
    mpres@tmYLLabelsOn           = True
    mpres@tmYLOn                 = True
    plot(0) = gsn_csm_map(wks1,mpres)
    res1=res
    res1@gsnAddCyclic                    = True
    res1@cnFillPalette                   = ocname
    ;res1@cnLevelSelectionMode            = "ManualLevels"
    ;res1@cnMinLevelValF                  = mnmxint(0)
    ;res1@cnMaxLevelValF                  = mnmxint(1)
    ;res1@cnLevelSpacingF                 = mnmxint(2)
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                       = orglev
   ;res1@cnFillColors                   = orgcolors
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False  
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = Groups (0)
    res1@gsnLeftString                  = Varstr + "  ("+ Varunt + ")"
    res1@gsnLeftStringOrthogonalPosF    = 0.02
   ;res1@gsnLeftStringParallelPosF      = 0.05
    res1@gsnRightString                 = "Mean: " + mean1
    res1@gsnRightStringOrthogonalPosF   = 0.02
   ;res1@gsnRightStringParallelPosF     = 0.9
    res1@cnLinesOn                      = plot_contour        
    res1@cnLineLabelsOn                 = False     
    res1@cnInfoLabelOn                  = False
    res1@cnLineLabelDensityF            = 2.0
    res1@cnLineLabelFontHeightF         = 0.018
    res1@cnLineLabelBackgroundColor     = -1
    res1@gsnContourLineThicknessesScale = 1.0
    res1@cnLineThicknessF               = 1.0
    res1@cnLineColor                    = "Black"
    res1@gsnContourZeroLineThicknessF   = 2.
    res1@gsnContourNegLineDashPattern   = 1
    plot_org1  = gsn_csm_contour(wks1,wncre0300,res1)
   
    resb                   = True ; polyline mods desired
    resb@gsLineColor       = "black" ; color of lines
    resb@gsLineThicknessF  = 5.0 ; thickness of lines
    txres                  =  True
    txres@txFontHeightF    = 0.014
    txres@txFontColor      = resb@gsLineColor
    txres@txFontThicknessF = 10 ;resb@gsLineThicknessF
    txres@txJust           = "CenterCenter"
    ;do irg = 0,nsubreg-1
    ; if (irg.eq.0)then
    ;  reg1=new((/4,nsubreg/),graphic)
    ;  reg2=new((/nsubreg/),graphic)
    ; end if 
    ; yptsb = (/maxlat(irg), minlat(irg), minlat(irg), maxlat(irg), maxlat(irg)/) ;define lat points for each line
    ; xptsb = (/minlon(irg), minlon(irg), maxlon(irg), maxlon(irg), minlon(irg)/) ; define lon points for each line
    ; do i = 0,3
    ;  resb@gsLineColor = "Black" ; color of lines
    ;  reg1(i,irg)=gsn_add_polyline(wks1,plot_org1,xptsb(i:i+1),yptsb(i:i+1),resb) 
    ; end do
    ; reg2(irg) = gsn_add_text(wks1,plot_org1,Subreg(irg),0.5*(maxlon(irg)+minlon(irg)),0.5*(maxlat(irg)+minlat(irg)),txres)
    ;end do 
    overlay(plot(0),plot_org1)
    delete([/res1/])
    end if 
 
    res2=res
    res2@gsnAddCyclic                   = True
    ;res2@cnLevelSelectionMode          = "ManualLevels"
    ;res2@cnMinLevelValF                = mnmxint(0)
    ;res2@cnMaxLevelValF                = mnmxint(1)
    ;res2@cnLevelSpacingF               = mnmxint(2)
     res2@cnLevelSelectionMode           := "ExplicitLevels"
     res2@cnLevels                       := diffgrp1
     res2@cnFillPalette                  := dcname
     res2@lbBoxSeparatorLinesOn           = True
   ;res2@cnFillColors                   = orgcolors
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@gsnAddCyclic                   = True
    res2@tiMainString                   = Groups(igp)
    res2@gsnLeftString                  = Varstr + "  ("+ Varunt + ")"
    res2@gsnLeftStringOrthogonalPosF    = 0.02
   ;res2@gsnLeftStringParallelPosF      = 0.05
    if (plotrdif)then
      res2@gsnLeftString                = Varstr + "  (%)"
      res2@gsnRightString               = "Global mean: " + mean2+"; "+ "RMSE: "+ rmse2
    else
      res2@gsnLeftString                = Varstr + "  ("+ Varunt + ")"
      res2@gsnRightString               = "Global mean: " + mean2+"; "+ "RMSE: "+ rmse2
    end if
    res2@gsnRightStringOrthogonalPosF   = 0.02
   ;res2@gsnRightStringParallelPosF     = 0.9
    res2@cnLinesOn                      = plot_contour
    res2@cnLineLabelsOn                 = False
    res2@cnInfoLabelOn                  = False
    res2@cnLineLabelDensityF            = 2.0
    res2@cnLineLabelFontHeightF         = 0.018
    res2@cnLineLabelBackgroundColor     = -1
    res2@gsnContourLineThicknessesScale = 1.0
    res2@cnLineThicknessF               = 1.0
    res2@cnLineColor                    = "Black"
    res2@gsnContourZeroLineThicknessF   = 2.
    res2@gsnContourNegLineDashPattern   = 1
    plot_rgd1  = gsn_csm_contour(wks1,wncredf1,res2)

    resc                   = True ; polyline mods desired
    resc@gsLineColor       = "black" ; color of lines
    resc@gsLineThicknessF  = 5.0 ; thickness of lines
    txres1                  =  True
    txres1@txFontHeightF    = 0.014
    txres1@txFontColor      = resb@gsLineColor
    txres1@txFontThicknessF = 10 ;resb@gsLineThicknessF
    txres1@txJust           = "CenterCenter"
    ;do irg = 0,nsubreg-1
    ; if (irg.eq.0)then
    ;  reg3=new((/4,nsubreg/),graphic)
    ;  reg4=new((/nsubreg/),graphic)
    ; end if
    ; yptsc = (/maxlat(irg), minlat(irg), minlat(irg), maxlat(irg), maxlat(irg)/) ;define lat points for each line
    ; xptsc = (/minlon(irg), minlon(irg), maxlon(irg), maxlon(irg), minlon(irg)/) ; define lon points for each line
    ; do i = 0,3
    ;  resc@gsLineColor = "Black" ; color of lines
    ;  reg3(i,irg)=gsn_add_polyline(wks1,plot_rgd1,xptsc(i:i+1),yptsc(i:i+1),resc)
    ; end do
    ; reg4(irg) = gsn_add_text(wks1,plot_rgd1,Subreg(irg),0.5*(maxlon(irg)+minlon(irg)),0.5*(maxlat(irg)+minlat(irg)),txres)
    ;end do
    overlay(plot(igp),plot_rgd1)

    res3=res
    res3@cnFillOn              = True     ; color plot desired
    res3@lbLabelBarOn          = False    ; turn off contour labels
    res3@cnLinesOn             = False    ; turn off contour lines
    res3@cnLineLabelsOn        = False    ; turn off contour labels
    res3@cnMissingValFillColor = -1
    res3@cnMonoFillColor       = True     ; default color is fground [black]
    res3@cnMonoFillPattern     = True     ; want multiple patterns
    res3@cnMonoFillScale       = True     ; want patterns w/ diff densities
    res3@cnFillColor           = FillColor
    res3@cnFillPattern         = FillPattern
    res3@cnFillScaleF          = FillScaleF
    res3@cnFillDotSizeF        = FillDotSizeF
    res3@cnFillOpacityF        = FillOpacityF
    res3@tiMainString          = "" ;label
    res3@cnInfoLabelOn         = False
    str=unique_string("tmp1")
    plot@str   = gsn_csm_contour(wks1,v1mnts1,res3)
    ;overlay(plot(igp),plot@str)
    delete([/res3/])

    delete([/res,resc/])
    delete([/mpres/])

   end do ; end of experiment group loop


;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = True
    pnres@gsnPanelLabelBar        = False
    pnres@gsnPanelMainString      = "";Varstr; "Surface Temperature ("+ts@units+")"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@lbOrientation           = "vertical"        ; vertical label bar
    pnres@pmLabelBarHeightF       = 0.14
    pnres@pmLabelBarWidthF        = 0.03
    pnres@gsnPanelMainFontHeightF = 0.030
    pnres@lbLabelFontHeightF      = 0.008               ; make labels smaller
    pnres@pmLabelBarParallelPosF  = -0.005
    pnres@pmLabelBarOrthogonalPosF= 0.005
;   pnres@pmLabelBarHeightF       = 0.1
;   pnres@pmLabelBarWidthF        = 0.57
;   pnres@lbTitleString           = "(" + Varunt + ")";Varstr+ " ("+Varunt+")";
;   pnres@lbTitlePosition         = "Right"                           ; title location
;   pnres@lbTitleDirection        = "Across"                          ; letter angle
;   pnres@lbTitleAngleF           = 90.                               ; title angle
;   pnres@lbTitleFontHeightF      = 0.012                              ; font height
;   pnres@gsnPanelBottom          = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 2
    gsn_panel(wks1,(/plot(:)/),(/1,ngrps/),pnres)
    delete([/pnres/])

 end do ; end of var loop

end

