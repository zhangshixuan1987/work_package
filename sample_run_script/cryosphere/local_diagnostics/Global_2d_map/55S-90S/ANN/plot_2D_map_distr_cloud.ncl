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
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Polar_SH"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    = -55.
  MinLon    =   0.
  MaxLon    = 360.

  ;;;the time for the plot;;;
  Plot_time = (/"ANN"/)
  nplttim   = dimsizes(Plot_time)
  
  ;;;flags to determine which ERA5 datasets will be compared;;;; 
  l_era5_single_level   = True
  l_era5_pressure_level = False
  l_era5_surface_level  = False
  l_plot_era5_spread    = True ;; determine if we plot era5 ensemble spread
  l_plot_reference      = True ;; determine if we plot a reference climo difference between observations

  load "./NCL_experiment_information.ncl"

 ;;;;;;;date strings for the data directory;;;;;;;;
  Expnames  = Exp_short_names
  Groups    = Exp_long_names
  ngrps     = dimsizes(Groups)

  LabList     = Exp_short_names
  LabList(0)  = Exp_short_names(0)
  LabList(1:) = Exp_short_names(1:) + " ~F19~:~N~~F~ " + Exp_short_names(0)
  plot        = new(ngrps+2,graphic)

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  load "./NCL_variable_list_cloud.ncl"
  VarList      = VarList1
  VstList      = VstList1
  VnmList      = VnmList1
  UnitList     = UnitList1
  nvars        = dimsizes(VarList)

 ;;; flags to control the plot style
  plotcntl          = True   ; if True, plot the CNTL origianl field on the leftmost panel

  plot_ttest_mskout = False  ; if True, mask the non-significant region as missing values

  plot_ttest_hatch  = False ; if True, plot the significant region as hatch region

  plot_ttest_hatrv  = False  ; if True, plot the non-significant region as hatch region

  plot_contour      = False   ; if True, plot the contour in the figure

  plot_leftstr      = True   ; if True, the lefadd igure captions as left string

  plot_msktrop      = False   ; if True, mask out the tropics region [-20,20]

 ;;color for the map boundary line
  lndcolor = "Grey28"

do ison = 0,nplttim-1

 SEASON_STRING = Plot_time(ison)


do iv = 0,nvars-1, 1

 Varname = VarList(iv)
 Varunt  = UnitList(iv)
 Varstr  = VstList(iv)

 plotFileLabel = Varname+""

 ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
 OUTDir  = "./"
 wtype = "png"
 wtype@wkWidth = 2500
 wtype@wkHeight = 2500 
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks1 = gsn_open_wks(wtype,OUTDir+"fig_"+SEASON_STRING+"_2d_map_"+plotFileLabel)

 ;;;set up the plot information;;;;
 load "./var_share_colorbar.ncl"

 do igp = 1, ngrps-1
 
  FERA5_HRES    = addfile(ERA5_HRES,"r")
  l_plot_era5_spread = ERA5_ENS.ne.""

  if(l_plot_era5_spread)then 
    FERA5_ENS = addfile(ERA5_ENS,"r")
  end if 

  print(Exp_data_files(igp))
  FE3SM      = addfile(Exp_data_files(igp),"r")
  
  if (l_plot_reference) then
   tmpfil               = systemfunc ("ls -1 "+ERA5_CLIM+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FERA5_CLIM           = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+CERES_EBAF_TOA+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FCERES_EBAF_SRF_CLIM = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+CERES_EBAF_SFC+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FCERES_EBAF_TOA_CLIM = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+GPCP_v2p3+"/GPCP*"+SEASON_STRING+"*climo.nc") ; file paths
   FGPCP_CLIM           = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+ISCCPCOSP+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FISCCPCOSP           = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+OAFlux+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FOAFlux              = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+HadISST+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FHadISST             = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+SSMI+"/SSMI_"+SEASON_STRING+"*climo.nc") ; file paths
   FSSMI                = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+SSMI+"/SSMI_SEAICE_"+SEASON_STRING+"*climo.nc") ; file paths
   FSSMIICE             = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+CRU+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FCRU                 = addfile(tmpfil,"r")
   tmpfil               = systemfunc ("ls -1 "+MERRA2+"/*"+SEASON_STRING+"*climo.nc") ; file paths
   FMERRA2              = addfile(tmpfil,"r")
  end if 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 time   = FE3SM->time
 lat    = FE3SM->lat
 lon    = FE3SM->lon
 
 Missing_Value = FERA5_HRES->mlspr@_FillValue

;;;;construct  the array to save the data;;;;;;;;;
 ntime  = dimsizes(time) ;; 3 years = 36 months
 nlat   = dimsizes(lat)  ;; latitude point
 nlon   = dimsizes(lon)  ;; longitutde point

 Var1   = new((/ntime,nlat,nlon/),double,FERA5_HRES->ttr@_FillValue)
 Var1!0 = "time"
 Var1!1 = "lat"
 Var1!2 = "lon"
 Var1&lat = lat
 Var1&lon = lon

 Var2     = Var1
 Vsd1     = Var1

 Vrf1   = new((/1,nlat,nlon/),double,FERA5_HRES->ttr@_FillValue)
 Vrf1!0 = "time"
 Vrf1!1 = "lat"
 Vrf1!2 = "lon"
 Vrf1&lat = lat
 Vrf1&lon = lon
 Vrf2   = Vrf1

 ;;directly read monthly mean data ;;;;
 if (Varname.eq."CLDTOT")then
  Var1(:,:,:)    = (/FERA5_HRES->tcc/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->tcc/)
  end if 
  if(l_plot_reference)then
   Vrf1(:,:,:)   = (/linint2 (FERA5_CLIM->lon,FERA5_CLIM->lat,FERA5_CLIM->clt, True,lon,lat, 0)/)
   Vrf2(:,:,:)   = (/linint2 (FISCCPCOSP->lon,FISCCPCOSP->lat,FISCCPCOSP->CLDTOT_ISCCPCOSP,True,lon,lat,0)/)
   Vrf1          = Vrf1/100.0
   Vrf2          = Vrf2/100.0
   RefLabels     = "ERA5 - ISCCPCOSP (CLIM)"
  end if
 end if 

 if (Varname.eq."CLDLOW")then
  Var1(:,:,:)    = (/FERA5_HRES->lcc/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->lcc/)
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   =  Missing_Value
   Vrf2(:,:,:)   =  Missing_Value
   Vrf1          = Vrf1/100.0
   RefLabels     = "ERA5 - ISCCPCOSP (CLIM)"
  end if
 end if

 if (Varname.eq."CLDHGH")then
  Var1(:,:,:)    = (/FERA5_HRES->hcc/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->hcc/)
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   = Missing_Value
   Vrf2(:,:,:)   = Missing_Value
   RefLabels     = "ERA5 - ISCCPCOSP (CLIM)"
  end if
 end if 

 if (Varname.eq."PRECL")then
  Var1(:,:,:)    = (/FERA5_HRES->mlspr/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->mlspr/)
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   = Missing_Value
   Vrf2(:,:,:)   = Missing_Value
   RefLabels     = ""
  end if
 end if

 if (Varname.eq."PRECC")then
  Var1(:,:,:)    = (/FERA5_HRES->mcpr/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->mcpr/)
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   = Missing_Value
   Vrf2(:,:,:)   = Missing_Value
   RefLabels     = ""
  end if
 end if

 if (Varname.eq."PRECT")then
  Var1(:,:,:)    = (/FERA5_HRES->mtpr/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->mtpr/)
  end if
  if(l_plot_reference)then
  ;Vrf1(:,:,:)   = (/area_conserve_remap (FERA5_CLIM->lon,FERA5_CLIM->lat,FERA5_CLIM->pr, lon,lat, False)/)
  ;Vrf2(:,:,:)   = (/area_conserve_remap (FGPCP_CLIM->lon,FGPCP_CLIM->lat,FGPCP_CLIM->PRECT,lon,lat,False)/)
   Vrf1(:,:,:)   = (/linint2 (FERA5_CLIM->lon,FERA5_CLIM->lat,FERA5_CLIM->pr, True,lon,lat, 0)/)
   Vrf2(:,:,:)   = (/linint2 (FGPCP_CLIM->lon,FGPCP_CLIM->lat,FGPCP_CLIM->PRECT,True,lon,lat,0)/)
   Vrf1          = Vrf1* 24.0 * 3600.0 ;convert to mm/day
   RefLabels     = "ERA5 - GPCP_V2.3 (Clim)"
  end if
 end if

 if(Varname.eq."TGCLDLWP")then
  Var1(:,:,:)    = (/FERA5_HRES->tclw/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->tclw/)
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   = Missing_Value
   Vrf2(:,:,:)   = Missing_Value
   RefLabels     = ""
  end if
 end if 

 if(Varname.eq."TGCLDIWP")then
  Var1(:,:,:)    = (/FERA5_HRES->tciw/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = (/FERA5_ENS->tciw/)
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   = Missing_Value
   Vrf2(:,:,:)   = Missing_Value
   RefLabels     = ""
  end if
 end if

 if(Varname.eq."TGCLDCWP")then
  Var1(:,:,:)    = (/FERA5_HRES->tciw+FERA5_HRES->tclw/)
  if (l_plot_era5_spread)then
   Vsd1(:,:,:)   = Missing_Value
  end if
  if(l_plot_reference)then
   Vrf1(:,:,:)   = Missing_Value ;(/linint2 (FERA5_CLIM->lon,FERA5_CLIM->lat,FERA5_CLIM->pr, True,lon,lat, 0)/)
   Vrf2(:,:,:)   = Missing_Value ;(/linint2 (FGPCP_CLIM->lon,FGPCP_CLIM->lat,FGPCP_CLIM->PRECT,True,lon,lat,0)/)
   RefLabels     = ""
  end if
 end if

 ;;;;read E3SM data;;;;;;;;;
  if (Varname.eq."PRECT")then
   Var2   = (/FE3SM->PRECC  + FE3SM->PRECL/)
  else
   Var2   = (/FE3SM->$Varname$/)
  end if 
  

  Var1@_FillValue = Missing_Value
  Var2@_FillValue = Missing_Value
  Vsd1@_FillValue = Missing_Value
  Vrf1@_FillValue = Missing_Value
  Vrf2@_FillValue = Missing_Value

 if(Varname.eq."PRECT".or.Varname.eq."PRECL".or.Varname.eq."PRECC")then
   Var1   = Var1*24*3600.0
   Vsd1   = Vsd1*24*3600.0
   Var2   = Var2*24*3600.0*1000.0
 end if 

 if(Varname.eq."TGCLDLWP".or.Varname.eq."TGCLDIWP".or.Varname.eq."TGCLDCWP")then
   Var1   = Var1*1000.0
   Vsd1   = Vsd1*1000.0 
   Var2   = Var2*1000.0
 end if

 ;;;;;; calculate the  3-year mean climo;;;;
 if (SEASON_STRING .eq. "ANN")then 
  mpivar  = dim_avg_n_Wrap(Var1,0) 
  mpdvar  = dim_avg_n_Wrap(Var2,0)
  spivar  = dim_avg_n_Wrap(Vsd1,0)
 else
  xx1 = month_to_seasonN(Var1,SEASON_STRING)
  xx2 = month_to_seasonN(Var2,SEASON_STRING)
  xx3 = month_to_seasonN(Vsd1,SEASON_STRING)
  mpivar  = dim_avg_n_Wrap(xx1(0,:,:,:),0)
  mpdvar  = dim_avg_n_Wrap(xx2(0,:,:,:),0)
  spivar  = dim_avg_n_Wrap(xx3(0,:,:,:),0) 
  delete([/xx1,xx2,xx3/])  
 end if 

 rpivar  = dim_avg_n_Wrap(Vrf1,0)
 rpdvar  = dim_avg_n_Wrap(Vrf2,0)

 mpdmpi  = mpivar
 mpdmpi  = mpdvar - mpivar
 rpdmpi  = rpivar
 rpdmpi  = rpdvar - rpivar

 printVarSummary(mpdmpi)
 printVarSummary(rpdmpi)

 Figstr  = "" ;PltNums(iv*(ngrps-1)+ igp-1)

 if (plot_ttest_mskout.or.plot_ttest_hatch.or.plot_ttest_hatrv)then
   ;;;do a statistical test on the difference;;;;;;
   ;sigr  = 0.05                        ; critical sig lvl for r
   ;rad   = 4.*atan(1.)/180.    ;
   ;wgty  = cos(lat*rad)
   ;xtmp  = pivar(lat|:,lon|:,nens|:)
   ;ytmp  = pdvar(lat|:,lon|:,nens|:)
   ;xAve  = dim_avg (xtmp)              ; calculate means at each grid point
   ;yAve  = dim_avg (ytmp)
   ;xVar  = dim_variance (xtmp)         ; calculate variances
   ; yVar  = dim_variance (ytmp)
   ;xEqv  = equiv_sample_size (xtmp, sigr,0)
   ;yEqv  = equiv_sample_size (ytmp, sigr,0)
   ;iflag = False                        ; population variance similar
   ;prob  = ttest(xAve,xVar,xEqv, yAve,yVar,yEqv, iflag, False)
   ;delete([/xtmp,ytmp,xAve,yAve,xVar,yVar,xEqv,yEqv/])
   prob   = mpdmpi
   prob   = where(mpdmpi.ge.spdmpi,1.0,0.0)

   ;;;mask out non-significant point with FillValue;;;;;
   v1mnts1 = mpdmpi

   v1mnts1 = where(prob.le.sigr,mpdmpi, mpdmpi@_FillValue)

   if(plot_ttest_hatrv) then
     v1mnts1 = mpdmpi
     v1mnts1 = where(prob.gt.sigr,mpdmpi, mpdmpi@_FillValue)
   end if
   v1mnts1@_FillValue = mpdmpi@_FillValue

   delete([/prob/])

 end if 

  ;;;;calculate the mean and spread for ERA5;;;;;;;
  latrg     = lat({MinLat:MaxLat})
  rad       = 4.*atan(1.)/180.    ;
  cost      = cos(latrg*rad)       ; cosine weights
  mpimean   = wgt_areaave_Wrap (mpivar({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  spimean   = wgt_areaave_Wrap (spivar({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mpistr    = "Mean: " + sprintf("%5.2f", mpimean)
  spistr    = "Mean: " + sprintf("%5.2f", spimean)
  
  ;;;;;;;;;calculate mean bias, RMSE and pattern correlation;;;;;;;
  pdmpim    = wgt_areaave_Wrap(mpdmpi({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  pdmpirmse = wgt_arearmse(mpivar({MinLat:MaxLat},{MinLon:MaxLon}),mpdvar({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  pdmpicorr = pattern_cor(mpivar({MinLat:MaxLat},{MinLon:MaxLon}),mpdvar({MinLat:MaxLat},{MinLon:MaxLon}),cost,0)
  bias      = "Bias: "+sprintf("%5.2f", dim_avg(pdmpim))    ;+ " ~F34~1~N~~F~ " + sprintf("%5.1e", pdmpis)   
  rmse      = "RMSE: "+sprintf("%5.2f", dim_avg(pdmpirmse)) ;+ " ~F34~1~N~~F~ " + sprintf("%5.1e", dim_stddev(pdmpirmse)) 
  pcrr      = "Corr: "+sprintf("%5.2f", dim_avg(pdmpicorr)) ;+ " ~F34~1~N~~F~ " + sprintf("%5.1e", dim_stddev(pdmpicorr)) 

  rpdmpim    = wgt_areaave_Wrap(rpdmpi({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  rpdmpirmse = wgt_arearmse(rpivar({MinLat:MaxLat},{MinLon:MaxLon}),rpdvar({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  rpdmpicorr = pattern_cor(rpivar({MinLat:MaxLat},{MinLon:MaxLon}),rpdvar({MinLat:MaxLat},{MinLon:MaxLon}),cost,0)
  rbias      = "Bias: "+sprintf("%5.2f", dim_avg(rpdmpim))    
  rrmse      = "RMSE: "+sprintf("%5.2f", dim_avg(rpdmpirmse)) 
  rpcrr      = "Corr: "+sprintf("%5.2f", dim_avg(rpdmpicorr)) 

;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False

    FontHeightF = 0.022
    res@tiMainFontHeightF        = FontHeightF*1.2
    res@tmYLLabelFontHeightF     = FontHeightF
    res@tmXBLabelFontHeightF     = FontHeightF

    mpres                    = res

    res@lbLabelFontHeightF   = FontHeightF/1.2
    res@gsnStringFontHeightF = FontHeightF*1.2

    res@cnFillOn                 = True              ; color plot desired
   ;res@cnFillPalette            = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor    = "White"
    res@lbLabelBarOn             = True      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF   = 0.61 
    res@pmLabelBarOrthogonalPosF = 0.15
    res@lbOrientation            = "vertical"        ; vertical label bar
    res@pmLabelBarHeightF        = 0.55
    res@pmLabelBarWidthF         = 0.12
    res@lbTitleString            = Varstr+ " ("+Varunt+")";
    res@lbTitlePosition          = "Right"                           ; title location
    res@lbTitleDirection         = "Across"                          ; letter angle
    res@lbTitleAngleF            = 90.                               ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2                              ; font height


    leftstr                  = "" ;"Surface Temperature ("+ts@units+")"
    res@gsnLeftString        = "" ; leftstr ;var_name0    ; long_name is too long!
    res@gsnRightString       = "" ;" ("+Varunt+")";""
 ;  res@tiMaingString        = "Annual Mean Response to Present-day climate"

    ;---Control appearance of map.
   mpres@gsnPolarLabelFontHeightF = FontHeightF
   mpres@gsnPolar               = "SH"
   mpres@mpFillOn               = False        ; turn off map fill
   mpres@mpMaxLatF              = MaxLat
   mpres@mpCenterLonF            = 0
   mpres@mpLabelsOn              = False
   mpres@mpPerimOn               = True
   mpres@mpFillOn                = False
   mpres@mpGeophysicalLineColor  = lndcolor  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
   mpres@mpOutlineOn              = True
   mpres@mpOutlineDrawOrder       = "PostDraw"
   mpres@mpFillDrawOrder          = "Predraw"
   mpres@mpGridAndLimbOn          = True
   mpres@mpGridLineThicknessF     = 3.0
;  mpres@mpOceanFillColor         = "lightskyblue1"
;  mpres@mpLandFillColor          = "gray33"

    ;if (igp.eq.0) then
    ; mpres@tmYLLabelsOn           = True
    ; mpres@tmYLOn                 = True
    ; plot1(igp) = gsn_csm_map(wks1,mpres)
    ;else if (igp.lt.ngrps-1) then
    ; mpres@tmYLLabelsOn           = False
    ; mpres@tmYLOn                 = False
    ; plot1(igp) = gsn_csm_map(wks1,mpres)
    ;else
    ; mpres@tmYLLabelsOn           = False
    ; mpres@tmYLOn                 = False
    ; mpres@tmYROn                 = True
     plot(igp+2) = gsn_csm_map(wks1,mpres)
    ;end if
    ;end if

    mpres@tmYLLabelsOn           = True
    mpres@tmYLOn                 = True
    plot(igp+2) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
  if(plotcntl)then
    plot(0)                   = gsn_csm_map(wks1,mpres)
    plot(1)                   = gsn_csm_map(wks1,mpres)
    plot(2)                   = gsn_csm_map(wks1,mpres)
    res1=res
    res1@gsnAddCyclic         = True
    res1@cnFillPalette                   = ocname
   ;res1@cnLevelSelectionMode            = "ManualLevels"
   ;res1@cnMinLevelValF                  = mnmxint(0)
   ;res1@cnMaxLevelValF                  = mnmxint(1)
   ;res1@cnLevelSpacingF                 = mnmxint(2)
    res1@cnLevelSelectionMode           := "ExplicitLevels"
    res1@cnLevels                       := orglev
    res1@lbLabelStrings                 := sprintf("%5.2f",orglev)
   ;res1@cnFillColors                   := orgcolors
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = LabList(0) + " (reanalysis)" ;Figstr + LabList0(igp)
    res1@tiMainOffsetYF                 = 0.01
    res1@gsnLeftString                  = ""
    res1@gsnLeftStringOrthogonalPosF    = 0.02
   ;res1@gsnLeftStringParallelPosF      = 0.05
    res1@gsnRightString                 = ""; mpistr ;
    res1@gsnRightStringOrthogonalPosF   = 0.02
    res1@gsnRightStringParallelPosF     = 0.98
    res1@gsnCenterString                = "" ; mpistr
    res1@gsnCenterStringOrthogonalPosF  = 0.04
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
    plot_org1  = gsn_csm_contour(wks1,mpivar,res1)
    overlay(plot(0),plot_org1)
    delete([/res1/])
  
    res1x=res
    res1x@gsnAddCyclic                   = True
    res1x@cnFillPalette                  = xcname
    res1x@cnLevelSelectionMode           = "ExplicitLevels"
    res1x@cnLevels                      := diffgrp0
    res1x@lbBoxLinesOn                   = True
   ;res1x@lbBoxSeparatorLinesOn          = False
   ;res1x@lbLabelStride                  = 15
    res1x@cnLinesOn                      = False             ; turn off contour lines
    res1x@cnLineLabelsOn                 = False            ; turn off contour labels
    res1x@cnInfoLabelOn                  = False
    res1x@tiMainString                   = LabList(0)+" (ensemble spread)" ;Figstr + LabList1(igp) ;"PD300 - CNTL"
    res1x@tiMainOffsetYF                 = 0.01
    res1x@gsnLeftString                  = "" ;Varstr +"" + "  ("+ Varunt + ")" ;Groups(igp) + " ~F19~:~N~~F~ " + Groups(0)
    res1x@gsnLeftStringOrthogonalPosF     = 0.04
   ;res1x@gsnLeftStringParallelPosF      = 0.05
    res1x@gsnRightString                  = "" ;spistr
    res1x@gsnRightStringOrthogonalPosF    = 0.04
   ;res1x@gsnRightStringParallelPosF     = 0.9
    res1x@gsnCenterString                 = "" ;"Global mean: "+ mdfstr
    res1x@gsnCenterStringOrthogonalPosF     = 0.04
   ;res1x@lbLabelStrings                 := sprintf("%5.2f",diffgrp0)
    plot_org2  = gsn_csm_contour(wks1,spivar,res1x)
    overlay(plot(1),plot_org2)
    delete([/res1x/])

    tres                      = True                ; text mods desired
    tres@txFontHeightF        = FontHeightF               ; make smaller
    txidx1 = gsn_create_text(wks1, mpistr, tres)
    txidx2 = gsn_create_text(wks1, spistr, tres)
    amres = True
    amres@amParallelPosF   =  0.5    ; This is the right edge of the plot.
    amres@amOrthogonalPosF =  0.5    ; This is the bottom edge of the plot.
    annoidx1 = gsn_add_annotation(plot(0), txidx1, amres)
    annoidx2 = gsn_add_annotation(plot(1), txidx2, amres)

    if(l_plot_reference)
     res2x=res
     res2x@gsnAddCyclic                   = True
     res2x@cnFillPalette                  = dcname
     res2x@cnLevelSelectionMode           = "ExplicitLevels"
     res2x@cnLevels                       := diffgrp1
     res2x@lbLabelStrings                 := sprintf("%5.2f",diffgrp1)
   ; res2x@cnFillColors                   := colorgroup1
     res2x@lbBoxLinesOn                   = True
    ;res2x@lbBoxSeparatorLinesOn          = False
    ;res2x@lbLabelStride                  = 15
     res2x@cnLinesOn                      = False             ; turn off contour lines
     res2x@cnLineLabelsOn                 = False            ; turn off contour labels
     res2x@cnInfoLabelOn                  = False
     res2x@tiMainString                   = RefLabels 
     res2x@tiMainOffsetYF                 = 0.01
     res2x@gsnLeftString                  = "" ;Varstr +"" + "  ("+ Varunt + ")" ;Groups(igp) + " ~F19~:~N~~F~ " + Groups(0)
     res2x@gsnLeftStringOrthogonalPosF     = 0.04
    ;res2x@gsnLeftStringParallelPosF       = 0.05
     res2x@gsnRightString                  = "" ;spistr
     res2x@gsnRightStringOrthogonalPosF    = 0.04
    ;res2x@gsnRightStringParallelPosF      = 0.9
     res2x@gsnCenterString                 = "" ;"Global mean: "+ mdfstr
     res2x@gsnCenterStringOrthogonalPosF   = 0.04
    ;res2x@lbLabelStrings                 := sprintf("%5.2f",diffgrp0)
     plot_org3  = gsn_csm_contour(wks1,rpdmpi,res2x)
     overlay(plot(2),plot_org3)
     delete([/res2x/])
     
     figstring = ""+rrmse + "~C~"+ "    "+rbias+ "~C~" + "    "+rpcrr
     tres                      = True                ; text mods desired
     tres@txFontHeightF        = FontHeightF         ; make smaller
     txidx3 = gsn_create_text(wks1, figstring, tres)
     amres = True
     amres@amParallelPosF   =  0.5    ; This is the right edge of the plot.
     amres@amOrthogonalPosF =  0.6    ; This is the bottom edge of the plot.
     annoidx3 = gsn_add_annotation(plot(2), txidx3, amres)
     delete([/tres,amres/])
    end if 

   end if 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@gsnAddCyclic                   = True
    res2@cnFillPalette                  = dcname
   ;res2@cnLevelSelectionMode          = "ManualLevels"
   ;res2@cnMinLevelValF                = mnmxint(0)
   ;res2@cnMaxLevelValF                = mnmxint(1)
   ;res2@cnLevelSpacingF               = mnmxint(2)
    res2@cnLevelSelectionMode           = "ExplicitLevels"
    res2@cnLevels                       := diffgrp1
    res2@lbLabelStrings                 := sprintf("%5.2f",diffgrp1)
  ; res2@cnFillColors                   := colorgroup1
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@tiMainString                   = LabList(igp) ;Figstr + LabList1(igp) ;"PD300 - CNTL"
    res2@tiMainOffsetYF                 = 0.01
    res2@gsnLeftString                  = " " ;Varstr +"" + "  ("+ Varunt + ")" ;Groups(igp) + " ~F19~:~N~~F~ " + Groups(0)
    res2@gsnLeftStringOrthogonalPosF     = 0.04
    ;res2@gsnLeftStringParallelPosF      = 0.05 
    res2@gsnRightString                  = " " ;bias+ "~C~" + rmse + "~C~"+pcrr
    res2@gsnRightStringOrthogonalPosF    = 0.04
    res2@gsnRightStringParallelPosF     = 0.95
    res2@gsnCenterString                 = "" ;"Global mean: "+ mdfstr 
    res2@gsnCenterStringOrthogonalPosF     = 0.04
    if(plot_ttest_mskout) then
     plot_rgd1                  = gsn_csm_contour(wks1,v1mnts1,res2)
    else
     plot_rgd1                  = gsn_csm_contour(wks1,mpdmpi,res2)  
    end if 

    if(plot_contour)
     res21 = res2
     res21@gsnLeftString                  = ""
     res21@gsnRightString                 = ""
     res21@gsnCenterString                = ""
     res21@cnFillOn                       = False
     res21@cnLinesOn                      = True
     res21@cnLineLabelsOn                 = False            ; turn off contour labels
     res21@cnInfoLabelOn                  = False
     res21@cnLineLabelDensityF            = 2.0
     res21@cnLineLabelFontHeightF         = 0.018
     res21@cnLineLabelBackgroundColor     = -1
     res21@gsnContourLineThicknessesScale = 1.0
     res21@cnLineThicknessF               = 1.0
     res21@cnLineColor                    = "Black"
     res21@gsnContourZeroLineThicknessF   = 2.
     res21@gsnContourNegLineDashPattern   = 1
     plot_cont1                  = gsn_csm_contour(wks1,mpdmpi,res21)
     overlay(plot_rgd1,plot_cont1)
     delete([/res21/])
    end if
    overlay(plot(igp+2),plot_rgd1)
    delete([/res2/])

    ;---Add test shading to plot
    if(plot_ttest_hatch.or.plot_ttest_hatrv)then
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
    ;wrf_smooth_2d(v1mnts1,3)
     plot@str                = gsn_csm_contour(wks1,v1mnts1,res3)
     overlay(plot(igp+2),plot@str)
     delete([/res3/])
    end if


    figstring = ""+rmse + "~C~"+ "    "+bias+ "~C~" + "    "+pcrr
    tres                      = True                ; text mods desired
    tres@txFontHeightF        = FontHeightF         ; make smaller
    txid1 = gsn_create_text(wks1, figstring, tres)

    amres = True
    amres@amParallelPosF   =  0.5    ; This is the right edge of the plot.
    amres@amOrthogonalPosF =  0.6    ; This is the bottom edge of the plot.
    annoid1 = gsn_add_annotation(plot(igp+2), txid1, amres)
    delete([/res/])
 end do ;end of experiment group loop 

;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = True
    pnres@gsnPanelLabelBar        = False
    pnres@gsnPanelMainFontHeightF = 0.016
    pnres@gsnPanelMainString      = VnmList(iv)+ " ("+SEASON_STRING+")"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@lbOrientation           = "vertical"        ; vertical label bar
    pnres@pmLabelBarHeightF       = 0.12
    pnres@pmLabelBarWidthF        = 0.03
    pnres@lbLabelFontHeightF      = 0.006               ; make labels smaller
    pnres@pmLabelBarParallelPosF  = -0.005
    pnres@pmLabelBarOrthogonalPosF= 0.005
    pnres@lbTitleString           = Varstr +"" + "  ("+ Varunt + ")" 
    pnres@lbTitlePosition         = "Right"                           ; title location
    pnres@lbTitleDirection        = "Across"                          ; letter angle
    pnres@lbTitleAngleF           = 90.                               ; title angle
    pnres@lbTitleFontHeightF      = 0.006                              ; font height
    pnres@gsnPanelBottom          = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5

    pnres@gsnPanelBottom      = 0.4        ; move bottom up from 0.0 to 0.4
    plot1 = new(6,graphic)
    if (l_plot_reference) then
      plot1      = plot
    else
      plot1(0:1) = plot(0:1)
      plot1(3:)  = plot(2:)
    end if 
    gsn_panel(wks1,plot1,(/2,3/),pnres)

   delete([/pnres/])

end do ; end of var loop

end do ; end od season loop

end

