  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Purpose: A tempolate to make a eof plot for AAO patterns 
;          Data is vertically interpolated to pressure levels
;          1st panel: EOFs for the first three modes
;          2nd panel: time series plot
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  flmsk      = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/lnd_msk/landmsk.nc"
  data_dir   = "./"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"E3SM AMIP", "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  time_tag   = "2001-2014"

  seasons    = (/"ANN"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  varList    = (/"dTs","dFQn","dSAF","dCRE","dnSAF","dFLW","dFLWC","dLH","dSH"/)
  vstList    = (/"~F33~D~N~~F~TS","~F33~D~N~~F~Qn","~F33~D~N~~F~SAF(SW)","~F33~D~N~~F~CRE", \
                 "~F33~D~N~~F~no_SAF(SW)","~F33~D~N~~F~FLW(down)","~F33~D~N~~F~FLWC(down,clear-sky)", \
                 "-~F33~D~N~~F~(LH)","-~F33~D~N~~F~(SH)"/)
  varUnit    = (/"~S~o~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C",\
                 "W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C"/)
  nvar       = dimsizes(varList)

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  20.
  lonW       =   0.
  lonE       = 360.

  solat1     = -40.0
  solat2     = -60.0

  ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

   do j = 0,nvar-1,1

    vars   = varList(j)
    vunt   = varUnit(j)
    print("working on " + vars)


    fm0  = addfile(flmsk,"r")
    lmsk = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})

    if(j.eq.0) then 
      left_str = exp_strs(1) + " " + vstList(j)
      fl0  = data_dir+exp_name(1)+"_energy_MON_"+time_tag+".nc"
      f0   = addfile(fl0,"r")
    else
      left_str = exp_strs(0) + " " + vstList(j)
      fl0  = data_dir+exp_name(0)+"_energy_MON_"+time_tag+".nc"
      f0   = addfile(fl0,"r")
    end if 

    if(vars.eq."dSUM") then
      vdat = f0->dSAF1(:,:,{latS:latN},{lonW:lonE})
      vdat = vdat + f0->dSAF2(:,:,{latS:latN},{lonW:lonE})
      do i = 3,nvar-2,1
       vnm  = varList(i)
       vdat = vdat + f0->$vnm$(:,:,{latS:latN},{lonW:lonE})
      end do
    else if(vars.eq."dSAF") then
      vdat = f0->dSAF1(:,:,{latS:latN},{lonW:lonE})
      vdat = vdat + f0->dSAF2(:,:,{latS:latN},{lonW:lonE})
    else
      vnm  = varList(j)
      vdat = f0->$vnm$(:,:,{latS:latN},{lonW:lonE})
    end if
    end if

    vdat!0 = "ens"
    vdat!1 = "time"
    vdat!2 = "lat"
    vdat!3 = "lon"  

   print("working on season "+ seasons(i))
   if (seasons(i).eq."ANN") then 
     tmp = vdat(time|:,ens|:,lat|:,lon|:)
     va1 = month_to_annual(tmp,1)
   else 
     tmp = vdat(time|:,ens|:,lat|:,lon|:)
     va1 = month_to_season(tmp,seasons(i))
   end if 
   delete([/tmp,vdat/])
 
   vde  = dim_avg_n_Wrap(va1,1) 
   vdm  = dim_avg_n_Wrap(vde,0)
   vsd  = dim_stddev_n_Wrap(vde,0)
   
   ;;calculate the mean over southern ocean region 
   
   rad  = 4.0*atan(1.0)/180.0
   wgty = cos(f0->lat({solat1:solat2})*rad)
   tmp1 = wgt_areaave(vdm({solat1:solat2},:),wgty,1.0,0)
   tmp2 = wgt_areaave(vsd({solat1:solat2},:),wgty,1.0,0)
   glm1 = sprintf("%3.1f", tmp1) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", tmp2)
   right_str = "40-60S:  "+ glm1
  
   delete([/fl0,f0/])
   delete([/tmp1,tmp2,glm1/])
   delete([/va1,vde/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (j.eq.0) then 
    wtype           = "png"
    if(wtype.eq."png") then
      lthick          = 2.0
      BorderThick     = 2.0
      wtype@wkWidth   = 1200
      wtype@wkHeight  = 1200
    else
      lthick          = 1.0
      BorderThick     = 1.0
    end if

    fe              = "Fig_bias_std_map_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)

    lndcolor        = "Grey25"
    plot            = new(nvar,graphic)
    color2          = read_colormap_file("cmocean_balance")
    ncols           = dimsizes(color2(:,0))
    color2(ncols/2-5:ncols/2+5,:) = 0

   end if 

  ;;;start to process the plot;;;;
  load "./share_fig_info.ncl"

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeight = 0.018
   res@gsnStringFontHeightF = FontHeight*2.5
   res@tmXBLabelFontHeightF = FontHeight
   res@tmYLLabelFontHeightF = FontHeight
   res@tiMainFontHeightF    = FontHeight
   res@tiXAxisFontHeightF   = FontHeight
   res@tiYAxisFontHeightF   = FontHeight

   mpres     = res
   resx      = res
   resx1     = res
   
;---Control appearance of map.
   mpres@gsnMajorLatSpacing     = 20
   mpres@gsnMajorLonSpacing     = 45
   mpres@gsnMinorLatSpacing     = 15
   mpres@gsnMinorLonSpacing     = 15
  ;mpres@mpProjection           = "Robinson"       ; choose projection
   mpres@mpLabelsOn             = False
   mpres@mpPerimOn              = False
   mpres@mpFillOn               = False
   mpres@mpGeophysicalLineColor      = lndcolor  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = lthick          ; thickness of outlines
  ;mpres@mpOutlineOn            = True
  ;mpres@mpOutlineDrawOrder     = "PostDraw"
  ;mpres@mpFillDrawOrder        = "Predraw"
  ;mpres@mpGridAndLimbOn        = True
  ;mpres@mpGridLatSpacingF      = 90.0
  ;mpres@mpGridLonSpacingF      = 180.0
  ;mpres@mpGridLineThicknessF   = 3.0
  ;mpres@mpGridAndLimbDrawOrder ="PreDraw"
 ;  mpres@mpLimitMode           = "LatLon"
   mpres@mpCenterLonF           = 180.0         ; default is 0 [GM]
   mpres@mpMinLatF              = latS
   mpres@mpMaxLatF              = latN
   mpres@mpMinLonF              = lonW
   mpres@mpMaxLonF              = lonE
  ;mpres@mpOceanFillColor       = "lightskyblue1"
  ;mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(j) = gsn_csm_map(wks,mpres)

  ;*******************************************
  ; first plot
  ;*******************************************
   res0                       = res
   res0@gsnAddCyclic          = True
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@gsnCenterString       = ""
   res0@tiMainString          = ""
   res0@tiMainOffsetYF        = 0.005

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = "" ;main_str + " of " + vars + " (" + vunt +")"
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarOrthogonalPosF = 0.25
   res0@pmLabelBarWidthF         = 0.62
   res0@pmLabelBarHeightF        = 0.15
   res0@lbLabelFontHeightF       = FontHeight*0.7
   res0@lbTitleFontHeightF       = FontHeight*0.75

   res0@cnFillOn              = True     ; color plot desired
  ;res0@cnFillMode            = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor = "White"
   res0@cnFillPalette         = color2
   res0@cnLevelSelectionMode  = "ExplicitLevels"
   res0@tiMainString          = "" ;exp_str
   res0@gsnLeftString         = left_str
   res0@gsnRightString        = right_str
   res0@cnLevels             := cnlevs2
   res0@cnLinesOn             = False        ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"                   ; color of second contours
   res0@cnLineThicknessF      = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res0@gsnRightStringOrthogonalPosF = 0.03
   res0@gsnLeftStringOrthogonalPosF  = 0.03

   pltrgd  = gsn_csm_contour(wks,vsd,res0)
   overlay(plot(j) ,pltrgd)
   delete(vsd)

  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012
   pres1@gsnPanelLabelBar           = False                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = "" ;main_str + " of " + vars + " (" + vunt +")"
   pres1@lbTitleFontHeightF         = FontHeight*0.6
   pres1@lbLabelFontHeightF         = FontHeight*0.55 ;  font height
   pres1@lbBoxLineColor             = "black"
   pres1@lbBoxLinesOn               = True
  ;pres1@pmLabelBarParallelPosF     = 0.61
   pres1@pmLabelBarOrthogonalPosF   = 0.01
   pres1@pmLabelBarWidthF           = 0.95
   pres1@pmLabelBarHeightF          = 0.08
  ;pres1@gsnPanelBottom             = 0.45       ; space for label bar
  ;pres1@gsnPanelTop                = 0.85       ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 2
   pres1@gsnPanelXWhiteSpacePercent = 2
   gsn_panel(wks,plot,(/nvar,1/),pres1)          ; create first panel
  
 end do 

end

