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
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim_rad"
  obs_name   = "ERA5"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"E3SM AMIP", "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "2001-2014"

  ylabels    = (/"Jan","Feb","Mar","Apr","May","Jun",\
                 "Jul","Aug","Sep","Oct","Nov","Dec"/)
  ystring    = "Zonally Averaged Annual Cycle"
  maxlev     = 11

  varList    = (/"TS", "TREFHT", "PSL"/)
  varUnit    = (/"~S~o~N~C", "~S~o~N~C", "hPa"/)
  nvar       = dimsizes(varList)

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  10.
  lonW       =   0.
  lonE       = 360.

 ;sub regions for analysis 
 ;western equatorial Pacific (Niño4 region, 16E-150W, 5S-5N)
 ;western and eastern equatorial Pacific (90W-150W,5-5N) 
 ;Southern Ocean 40-60S, 0-360
  subreg    = (/"Nino34"/) 
  slats     = (/      -5/)   
  slatn     = (/       5/) 
  slonw     = (/     270/)
  slone     = (/     160/)
  regcols   = (/"Black" /)
  reglthk   = (/3       /)
  nsreg     = dimsizes(subreg)

  do j = 0,nvar-1,1

  vars   = varList(j)
  vunt   = varUnit(j)
  print("working on " + vars)

   do k = 0,nexps-1

    fl0  = systemfunc("ls -1 "+data_dir+"/"+obs_name+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fm0  = addfile(flmsk,"r")
    f0   = addfiles(fl0,"r")    
    f1   = addfiles(fl1,"r")   
    ListSetType (f0, "join")    
    ListSetType (f1, "join")

    vnm  = vars+"_o"
    vd0  = f0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
    vd1  = f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
    lmsk = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})

    printMinMax(vd0,False)
    printMinMax(vd1,False)

    dim0 = dimsizes(vd0)
    dim1 = dimsizes(vd1)
     
    nens = dim0(0)*dim1(0) 
    ntim = dim1(1)
    nlat = dim1(2)
    nlon = dim1(3)

    vdat = new((/nens,ntim,nlat,nlon/),typeof(vd0))
    vdat!0 = "ens"
    vdat!1 = "time"
    vdat!2 = "lat"
    vdat!3 = "lon"
    vdat&time = vd1&time 
    vdat&lat  = vd1&lat
    vdat&lon  = vd1&lon
    do ii = 0,dim0(0)-1,1
      do jj = 0,dim1(0)-1,1
        vdat(ii*dim1(0)+jj,:,:,:) = vd1(jj,:,:,:) - vd0(ii,:,:,:)
      end do 
    end do 
    vdat@_FillValue = -9999.0

    lms0 = conform(vdat,lmsk,(/2,3/))
    vdat = where(lms0.le.0,vdat,vdat@_FillValue)
    delete([/dim0,dim1,nens,ntim,nlat,nlon,lmsk,lms0/])

    if(vars.eq."PSL") then
     vdat = vdat/100.0
     vdat@units = "hPa"
     vunt    = "hPa"
    end if 

    if(vars.eq."TMQ") then
     vdat = vdat/1000.0
     vdat@units = "kg m~S~-2~N~"
     vunt    = "kg m~S~-2~N~"
    end if

    if(vars.eq."TAUX".or.vars.eq."TAUY") then
     vdat = vdat*100.0
     vdat@units = "10~S~-2~N~ N m~S~-2~N~"
     vunt    = "10~S~-2~N~ N m~S~-2~N~"
    end if

    vmn = dim_avg_n_Wrap(vdat,0) 
    vsd = dim_stddev_n_Wrap(vdat,0)
    vmn1 = clmMonTLL(vmn)
    vsd1 = clmMonTLL(vsd)
    vdm  = dim_avg_n_Wrap(vmn1,2)
    vds  = dim_avg_n_Wrap(vsd1,2)

   ;rad  = 4.0*atan(1.0)/180.0
   ;wgty = cos(f1[0]->lat({latS:latN})*rad)
   ;glm1 = dim_max_n(dim_max_n(vdm,1),0) 
   ;glm2 = dim_min_n(dim_min_n(vdm,1),0) 
   ;glm3 = wgt_areaave(vdm,1.0,wgty,0)
   ;rstr = "Max / Min / Mean:  "+ sprintf("%.2f",glm1)+" / "+sprintf("%.2f",glm2)+" / "+sprintf("%.2f",glm3)
   ;print(rstr)
   rstr = "ens. avg and stddev"
   ;delete([/glm1,glm2,glm3/])
    delete([/fl0,f0,fl1,f1/])
    delete([/vdat,vmn,vsd,vmn1,vsd1,vd0,vd1/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
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

    fe              = vars+"_sub_bias_vs_"+obs_name+"_anacyc"
    wks             = gsn_open_wks(wtype,fe)
    lndcolor        = "Grey25"
    plot            = new(nexps,graphic)
    dum             = new(nexps,graphic)
    color           = read_colormap_file("BlueDarkRed18") ;BlAqGrYeOrReVi200") ;BlGrYeOrReVi200")

   end if 

   main_str = "mean bias"
   exp_str  = exp_strs(k)
   left_str = exp_strs(k)
   lab_str  =  main_str + " of " + vars + " (" + vunt +")"

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
   res@vpWidthF       = 0.6
   res@vpHeightF      = 0.6

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
   res@gsnStringFontHeightF  = FontHeight*1.1
   res@tmXBLabelFontHeightF  = FontHeight
   res@tmYLLabelFontHeightF  = FontHeight
   res@tiMainFontHeightF     = FontHeight
   res@tiXAxisFontHeightF    = FontHeight*1.1
   res@tiYAxisFontHeightF    = FontHeight*1.1
   res@lbLabelFontHeightF    = FontHeight*0.95
   res@lbTitleFontHeightF    = FontHeight*1.05
   
   res@tiXAxisString    = "Latitude"
   res@tiXAxisOffsetYF  = 0.11

   res@tmXBMode         = "Explicit"
   res@tmXBValues       = ispan(-90,0,15)
   res@tmXBMinorValues  = ispan(-90,0,5)
   res@tmXBLabels       = (/"90S","75S","60S","45S","30S","15S","0"/)
   res@tiYAxisString    = "Month"
   res@tmYLMode         = "Explicit"      
   res@tmYLValues       = ispan(0,11,1)
   res@tmYLLabels       = ylabels
   res@tmYLLabelJust    = "CenterLeft"
  ;*******************************************
  ; first plot
  ;*******************************************
   res0                       = res
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@gsnCenterString       = ""
   res0@tiMainString          = ""
   res0@tiMainOffsetYF        = 0.005

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = lab_str 
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarOrthogonalPosF = 0.10
   res0@pmLabelBarWidthF         = 0.62
   res0@pmLabelBarHeightF        = 0.08
   res0@lbTitlePosition          = "Bottom"

   res0@cnFillOn              = True     ; color plot desired
  ;res0@cnFillMode            = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor = "White"
   res0@cnFillPalette         = color
   res0@cnLevelSelectionMode  = "ExplicitLevels"
   res0@tiMainString          = "" ;exp_str
   res0@gsnLeftString         = left_str
   res0@gsnRightString        = rstr
   res0@cnLevels             := cnlevs2
   res0@cnLinesOn             = False        ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"               ; color of second contours
   res0@cnLineThicknessF      = lthick           ; line thickness
   res0@gsnContourZeroLineThicknessF = lthick*2.0  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res0@gsnRightStringOrthogonalPosF = 0.03
   res0@gsnLeftStringOrthogonalPosF  = 0.03

   plot(k)  = gsn_csm_contour(wks,vdm,res0)

   res1 = res0
   res1@cnFillOn                     = False     ; color plot desired
   res1@tiMainString                 = "" ;exp_str
   res1@gsnLeftString                = ""
   res1@gsnRightString               = ""
   res1@cnLevelSelectionMode         = "ManualLevels"	
   res1@cnMinLevelValF               = minlev
   res1@cnMaxLevelValF               = maxlev
   res1@cnLevelSpacingF              = intlev
  ;res1@cnLevels                    := cnlevs1
   res1@cnLinesOn                    = True 
   res1@cnInfoLabelOn                = False
   res1@cnLineLabelsOn               = True    ; turn off contour labels
   res1@cnLineLabelFontHeightF       = FontHeight*1.1
   res1@cnConstFLabelPerimOn         = False 
  ;res1@cnLineLabelBackgroundColor   = "Transparent"
   res1@cnLineLabelFormat            = "0@;*.2f"                      ; Only one value after decimal point
   res1@cnLineLabelDensityF          = 2.0
   res1@cnLineColor                  = "Black"                   ; color of second contours
   res1@cnLineThicknessF             = lthick*1.5                ; line thickness
   res1@gsnContourZeroLineThicknessF = lthick*2.0  
   res1@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res1@gsnRightStringOrthogonalPosF = 0.03
   res1@gsnLeftStringOrthogonalPosF  = 0.03

   dum(k)  = gsn_csm_contour(wks,vds,res1)
   overlay(plot(k),dum(k))

   delete([/vdm,vds/])

  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012
   pres1@gsnPanelLabelBar           = False                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = lab_str 
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
   gsn_panel(wks,plot,(/1,nexps/),pres1)          ; create first panel
   delete([/dum,plot,res,res0,pres1/]) 

end do 

end

