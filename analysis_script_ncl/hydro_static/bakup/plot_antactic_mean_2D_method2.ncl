;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
external MAKEIC "./MAKEIC.so"
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
  grav      = 9.80665

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
  DATA_DIR    = "/global/cfs/cdirs/e3sm/zhan391/data/polar_diag/"
  ERA5_File   = DATA_DIR + "ERA5_analysis_360x180_monthly_1981-2010.nc"
  E3SM_File   = DATA_DIR + "20210727.v2.LR.F20TR_eam.h0_monthly_1981-2010.nc"
  intyp       = 2

  f1          = addfile(ERA5_File,"r")
  pfix        = f1->lev
  zsr5        = f1->phis(:,:,:)
  ter5        = f1->t(:,:,:,:)
  psr5        = f1->sp(:,:,:)
  zer5        = f1->z(:,:,:,:)
  zer5        = zer5/grav 

  f2          = addfile(E3SM_File,"r")
  lat         = f2->lat
  zsmd        = f2->PHIS(:,:,:)
  psrm        = f2->PS(:,:,:)
  tmod        = f2->T(:,:,:,:) 
  zmod        = f2->Z3(:,:,:,:)
  pm          = pres_hybrid_ccm(psrm,f2->P0,f2->hyam,f2->hybm)
  pi          = pres_hybrid_ccm(psrm,f2->P0,f2->hyai,f2->hybi)

  dimx = dimsizes(tmod)
  ntime = dimx(0)
  nlev  = dimx(1)
  nlat  = dimx(2)
  nlon  = dimx(3)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;interpolate to ERA5 pressure levels
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  tbot   = tmod(nlev-1,:,:)
  intyp  = 2                             ; 1=linear, 2=log, 3=log-log
  kxtrp  = True                          ; True=extrapolate
  varflg = 1                             ; temperature is variable
  tpm    = vinth2p_ecmwf(tmtmp,f2->hyam,f2->hybm,pfix,psam,intyp,f2->P0/100.0, \
               1,kxtrp,varflg,tbot,zsmd)
  varflg = -1
  zpm    = vinth2p_ecmwf(zmtmp,f2->hyam,f2->hybm,pfix,psam,intyp,f2->P0/100.0, \
               1,kxtrp,varflg,tbot,zsmd)
  delete([/tmtmp,pmtmp,pitmp,zsold,zsnew,psold,psnew/])

  printVarSummary(tpm)
  printVarSummary(zpm)
  ptop   = min(pfix)
  ztop   = 0
  dpmod  = dpres_plevel(pfix, psam, ptop, 0)  ; dp(30)
  dper5  = dpres_plevel(pfix, psr5, ptop, 0)  ; dp(30)
  dzmod  = dz_height(zpm,  zsr5(0,:,:), ztop, 0)
  dzer5  = dz_height(zer5, zsr5(0,:,:), ztop, 0)

  ;;open netcdf file to save the data for plots;;;;;
  setfileoption("nc", "Format",  "NetCDF4")
  out_file_nam1 = "Hydrostatic_analysis_out_ERA5.nc"
  out_file_nam2 = "Hydrostatic_analysis_out_ERA5.nc"

  system("rm " + out_file_name)
  fout = addfile(out_file_name,"c")
  fout->ter5  = ter5
  fout->zer5  = zer5
  fout->dzer5 = zsr5
  fout->dper5 = zsr5
  fout->plev  = pfix
  fout->zser5 = zsr5
  fout->pser5 = psr5

  fout->zeam  = zer5
  fout->team  = ter5
  fout->zseam = (/zsmd/)
  fout->pseam = psr5
  fout->dzeam = zsr5
  fout->dpeam = zsr5

  dp = P/R * /T * g * dz 

  do i  = 0,0 ; ntime -1,1

    zmtmp = todouble(zmod(i,:,:,:))
    tmtmp = todouble(tmod(i,:,:,:))
    pmtmp = todouble(pm(i,:,:,:))
    pitmp = todouble(pi(i,:,:,:))
    zsold = todouble(zsmd(i,:,:))
    zsnew = todouble(zsr5(i,:,:))
    psold = todouble(psrm(i,:,:))
    psnew = todouble(psrm(i,:,:))

    ;----------
    ; Adjust Ps
    ;----------
    print ("     Adjusting PS")
    MAKEIC::psadj(nlev ,nlev+1, nlat, nlon, tmtmp, \
                  pmtmp , pitmp , zsold , zsnew ,psold , \
                  psnew)
    psam  = psrm(i,:,:)
    psam  = tofloat(psnew)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;interpolate to ERA5 pressure levels 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    tbot   = tmtmp(nlev-1,:,:)
    intyp  = 2                             ; 1=linear, 2=log, 3=log-log
    kxtrp  = True                          ; True=extrapolate
    varflg = 1                             ; temperature is variable
    tpm    = vinth2p_ecmwf(tmtmp,f2->hyam,f2->hybm,pfix,psam,intyp,f2->P0/100.0, \
                 1,kxtrp,varflg,tbot,zsnew)
    varflg = -1   
    zpm    = vinth2p_ecmwf(zmtmp,f2->hyam,f2->hybm,pfix,psam,intyp,f2->P0/100.0, \
                 1,kxtrp,varflg,tbot,zsnew)
    delete([/tmtmp,pmtmp,pitmp,zsold,zsnew,psold,psnew/])

    printVarSummary(tpm)
    printVarSummary(zpm)
    ptop   = min(pfix)
    ztop   = 0
    dpmod  = dpres_plevel(pfix, psam, ptop, 0)  ; dp(30)
    dper5  = dpres_plevel(pfix, psr5, ptop, 0)  ; dp(30)
    dzmod  = dz_height(zpm,  zsr5(0,:,:), ztop, 0)  
    dzer5  = dz_height(zer5, zsr5(0,:,:), ztop, 0)  

    print(min(dpmod)+ " " +max(dpmod))
    print(min(dper5)+ " " +max(dper5))
    print(min(dzmod)+ " " +max(dzmod))
    print(min(dzer5)+ " " +max(dzer5))
    fout->dzer5(i,:,:) = (/dzer5/)
    fout->dzeam(i,:,:) = (/dzmod/)
    fout->dper5(i,:,:) = (/dper5/)
    fout->dpeam(i,:,:) = (/dpmod/)
    fout->zeam(i,:,:)  = (/zpm/)
    fout->team(i,:,:)  = (/tpm/)
    fout->pseam(i,:,:) = (/psam/)

  end do 

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
