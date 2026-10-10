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
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.

 ;;;;;;;date strings for the data directory;;;;;;;;
  Groups      = (/"ERA5","EAMv2"/)
  ngrps       = dimsizes(Groups)

  npanel      = 6 
  LabList     = new(npanel,string)
  LabList(0)  = Groups(0) + " reanalysis"
  LabList(1)  = "E3SM - ERA5"
  LabList(2)  = "Estimated with fixed PS  (ERA5)"
  LabList(3)  = "Estimated with fixed PS  (E3SM - ERA5)"
  LabList(4)  = "Estimated with fixed T2m (ERA5)"
  LabList(5)  = "Estimated with fixed T2m (E3SM - ERA5)"
  Rstring     = "PSL (hPa)"

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DATA_DIR     = "/lcrc/group/acme/ac.szhang/acme_scratch/data/polar_diag/"
  ERA5_File    = DATA_DIR + "ERA5_analysis_360x180_monthly_1981-2010.nc"
  E3SM_File    = DATA_DIR + "20210727.v2.LR.F20TR_eam.h0_monthly_1981-2010.nc"
  intyp        = 2


  FERA5        = addfile(ERA5_File,"r")
  pres         = FERA5->lev
  zer5         = FERA5->phis
  ter5         = FERA5->t(:,:,:,:)
  qer5         = FERA5->q(:,:,:,:)
  psr5         = FERA5->sp(:,:,:)
  slp5         = FERA5->msl(:,:,:)
  tsr5         = FERA5->t2m(:,:,:)
 
  FE3SM        = addfile(E3SM_File,"r")
  lat          = FE3SM->lat
  psrm         = FE3SM->PS(:,:,:)
  slpm         = FE3SM->PSL(:,:,:)
  tsmd         = FE3SM->TREFHT(:,:,:) 
  zmod         = vinth2p(FE3SM->Z3(:,:,:,:),FE3SM->hyam,FE3SM->hybm,pres,psrm,intyp,FE3SM->P0/100.0,1,False)
  qmod         = vinth2p(FE3SM->Q(:,:,:,:),FE3SM->hyam,FE3SM->hybm,pres,psrm,intyp,FE3SM->P0/100.0,1,False)
  tmod         = vinth2p(FE3SM->T(:,:,:,:),FE3SM->hyam,FE3SM->hybm,pres,psrm,intyp,FE3SM->P0/100.0,1,False)
 
  ;;;;;;;;;derive the sea level pressure;;; 
  dimt  = dimsizes(tmod)
  ntim  = dimt(0)
  klvl  = dimt(1)
  nlat  = dimt(2)
  mlon  = dimt(3)
  ptmp  = psrm 
  ptmp  = tofloat(pres(klvl-1))
  
  ;;;;set 1 ;;;;;
  ter51 = tsr5
  ter51 = tsr51*(1.+qer5*0.61)
  tmod1 = tmod
  tmod1 = tmod1*(1.+qmod*0.61)
  tmp1  = tsr51(time|:,lat|:,lon|:,lev|:)
  tmp2  = tmod1(time|:,lat|:,lon|:,lev|:)
  tmp3  = zer5
  tmp3  = tmp3/9.80655
  zhr51 = hydro (conform(tmp1,pres,3), tmp1, tmp3)
  zhmd1 = hydro (conform(tmp2,pres,3), tmp2, tmp3)
  delete([/ter51,tmod1,tmp1,tmp2,tmp3/])

  ;;;;set 2 ;;;;;
  ter51 = ter5
  ter51 = ter51*(1.+qer5*0.61)
  tmod1 = tmod
  tmod1 = tmod1*(1.+qer5*0.61)
  tmp1  = tsr51(time|:,lat|:,lon|:,lev|:)
  tmp2  = tmod1(time|:,lat|:,lon|:,lev|:)
  tmp3  = zer5
  tmp3  = tmp3/9.80655
  zhr52  = hydro (conform(tmp1,pres,3), tmp1, tmp3)
  zhmd2  = hydro (conform(tmp2,pres,3), tmp2, tmp3)
  delete([/ter51,tmod1,tmp1,tmp2,tmp3/])



  var1   = dim_avg_n_Wrap(zer5,0)
  var2   = dim_avg_n_Wrap(zmod,0)
  var3   = dim_avg_n_Wrap(zhr51,0)
  var4   = dim_avg_n_Wrap(zhmd1,0)
  var5   = dim_avg_n_Wrap(zhr52,0)
  var6   = dim_avg_n_Wrap(zhmd2,0)

  var2   = var2 - var1
  var4   = var4 - var3
  var6   = var6 - var5

  var1   = var1/100.0
  var2   = var2/100.0
  var3   = var3/100.0
  var4   = var4/100.0
  var5   = var5/100.0
  var6   = var6/100.0


  ;;;;calculate the global mean and correlation
  latrg   = lat({MinLat:MaxLat})
  rad     = 4.*atan(1.)/180.    ;
  cost    = cos(latrg*rad)       ; cosine weights
  mvar1   = wgt_areaave_Wrap (var1({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mvar2   = wgt_areaave_Wrap (var2({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mvar3   = wgt_areaave_Wrap (var3({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mvar4   = wgt_areaave_Wrap (var4({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)

  mstr1   = "Mean: " + sprintf("%5.2f", mvar1)
  mstr2   = "Mean: " + sprintf("%5.2f", mvar2)
  mstr3   = "Mean: " + sprintf("%5.2f", mvar3)
  mstr4   = "Mean: " + sprintf("%5.2f", mvar4)

  ;----------------------------------------------------------------------
  ; Plotting section

   OUTDir  = "./"
   wtype = "png"
   wtype@wkWidth = 2500
   wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_map_e3sm_vs_era5_PSL")
   plot        = new(npanel,graphic)


  ;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False

    FontHeightF              = 0.024
    res@tiMainFontHeightF    = FontHeightF*1.2
    res@tmYLLabelFontHeightF = FontHeightF
    res@tmXBLabelFontHeightF = FontHeightF
    mpres                    = res

    res@lbLabelFontHeightF   = FontHeightF*0.6
    res@gsnStringFontHeightF = FontHeightF*1.5

    res@cnFillOn                 = True              ; color plot desired
   ;res@cnFillPalette            = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor    = "White"
    res@lbLabelBarOn             = True      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF   = 0.61
    res@pmLabelBarOrthogonalPosF = 0.05
    res@lbOrientation            = "vertical"        ; vertical label bar
    res@pmLabelBarHeightF        = 0.65
    res@pmLabelBarWidthF         = 0.12
    res@lbTitleString            = ""
    res@lbTitlePosition          = "Right"                           ; title location
    res@lbTitleDirection         = "Across"                          ; letter angle
    res@lbTitleAngleF            = 90.                               ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2                              ; font height

    ;---Control appearance of map.
    mpres@gsnMajorLatSpacing      = 30
    mpres@gsnMajorLonSpacing      = 45
    mpres@gsnMinorLatSpacing      = 15
    mpres@gsnMinorLonSpacing      = 15
    mpres@mpFillOn                = False        ; turn off map fill
    mpres@mpCenterLonF            = 0
    mpres@mpLabelsOn              = False
    mpres@mpPerimOn               = True
    mpres@mpFillOn                = False
    mpres@mpGeophysicalLineColor  = "Black";lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
    mpres@mpOutlineOn              = True
    mpres@mpOutlineDrawOrder       = "PostDraw"
    mpres@mpFillDrawOrder          = "Predraw"
    mpres@mpGridAndLimbOn          = False
    mpres@mpGridLineThicknessF     = 1.0
;   mpres@mpOceanFillColor         = "lightskyblue1"
;   mpres@mpLandFillColor          = "gray33"
    mpres@tmYLLabelsOn             = True
    mpres@tmYLOn                   = True
    do igp = 0,npanel-1
      plot(igp) = gsn_csm_map(wks1,mpres)
    end do

    res@gsnLeftString                  = ""
    res@gsnLeftStringOrthogonalPosF    = 0.02
   ;res@gsnLeftStringParallelPosF      = 0.05
    res@gsnRightString                 = ""
    res@gsnRightStringOrthogonalPosF   = 0.02
    res@gsnRightStringParallelPosF     = 0.98
    res@gsnCenterString                = ""
    res@gsnCenterStringOrthogonalPosF  = 0.04

;---Resources for plotting original data
    res1=res
    res1@gsnAddCyclic                  = True
    res1@cnFillPalette                 = "WhBlGrYeRe"
    res1@cnLevelSelectionMode          = "ExplicitLevels"
    res1@cnLevels                     := (/965,970,975,985,990,995,1000,1005,1010,1015,1020/);
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = ""
    res1@tiMainOffsetYF                 = 0.01
    res1@gsnLeftString                  = LabList(0)
    res1@gsnRightString                 = Rstring
    res1@cnLinesOn                      = True
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

    plot_org1  = gsn_csm_contour(wks1,var1,res1)
    res1@gsnLeftString                  = LabList(2)
    overlay(plot(0),plot_org1)
    plot_org2  = gsn_csm_contour(wks1,var3,res1)
    overlay(plot(2),plot_org2)
    plot_org3  = gsn_csm_contour(wks1,var5,res1)
    overlay(plot(4),plot_org3)

    delete([/res1/])

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                   = True
    res2@cnFillPalette                  = "BlueDarkRed18" ;"NCV_blu_red"
    res2@cnLevelSelectionMode           = "ExplicitLevels"
    res2@cnLevels                      := (/-20,-10,-5,-2,-1,1,2,5,10,20/)
  ; res2@cnFillColors                  := colorgroup1
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@gsnLeftString                  = LabList(1)
    res2@gsnRightString                 = "~F33~D~N~~F~"+ Rstring
    res2@gsnCenterString                = "" ;"Global mean: "+ mdfstr
    res2@gsnCenterStringOrthogonalPosF  = 0.04
    plot_rgd1                           = gsn_csm_contour(wks1,var2,res2)
    overlay(plot(1),plot_rgd1)

    res2@gsnLeftString                  = LabList(3)
    plot_rgd2                           = gsn_csm_contour(wks1,var4,res2)
    overlay(plot(3),plot_rgd2)

    res2@gsnLeftString                  = LabList(5)
    plot_rgd3                           = gsn_csm_contour(wks1,var6,res2)
    overlay(plot(5),plot_rgd3)

;---Draw both plots in a panel
    pnres                            = True
    pnres@gsnMaximize                = True
    pnres@gsnPanelMainFontHeightF    = 0.016
    pnres@gsnPanelMainString         = ""
    pnres@gsnPanelMainFontColor      = "Black"
    pnres@gsnPanelBottom             = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5
    pnres@gsnPanelBottom             = 0.4
    gsn_panel(wks1,plot,(/3,2/),pnres)

   delete([/pnres/])


end 
