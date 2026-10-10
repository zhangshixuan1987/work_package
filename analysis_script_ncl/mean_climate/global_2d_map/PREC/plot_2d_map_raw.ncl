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

  l_create_plot_data = False
  l_create_plot_data = True
  plot_contour       = False   ; if True, plot the contour in the figure
  plot_leftstr       = True    ; if True, the lefadd igure captions as left string

  data_dir = "/global/cfs/cdirs/e3sm/zhan391/data/CVDP_RGD/"
  latS     = -45
  latN     = -90
  lonW     = 0
  lonE     = 360

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Obsnam   = (/"GPCP", "ERA5"/) 
  Casdir   = (/"Observation","CMIP6_MME","v2.LR.historical"/)
  Groups   = (/"OBS","CMIP6-MMM","E3SMv2"/)
  ngrps    = dimsizes(Groups)
  LabList  = Groups
  plot     = new (ngrps, graphic)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"PRECT", "PRECST"/)
  VstList  = (/"PRECT", "PRECST"/)
  UnitList = (/"mm day~S~-1~N~", "mm day~S~-1~N~"/)
  nvars    = dimsizes(VarList)
  Figstr   = (/"a", "b", "c", "d", "e", "f", "g", "h"/)
  Figstr   = "" 

  do iv = 0,nvars-1, 1

   Varname = VarList(iv)
   Varunt  = UnitList(iv)
   Varstr  = VstList(iv)

   if(Obsnam(iv).eq."CERES-OAFlux") then
     yst  = 200101
     yed  = 201412
   else
     yst  = 197901
     yed  = 201412
   end if

   do igp = 0, ngrps-1,1

   if (Groups(igp).eq."OBS") then
     Input = data_dir + "/" + Obsnam(iv) + "/"
   else if (Casdir(igp).eq."CMIP6_MME") then
     Input = data_dir + Casdir(igp)+"/historical/*/"
   else
     Input = data_dir + Casdir(igp)+"/"
   end if
   end if 

   if (Groups(igp).eq."OBS") then
     if(Varname.eq."FLNTC") then
       flm   = systemfunc("ls "+Input+"*en00*.FLNTOAC.*.nc")
     else
       flm   = systemfunc("ls "+Input+"*en00*."+Varname+".*.nc")
     end if
   else if (Groups(igp).eq."CMIP6_MME") then
     flm       = systemfunc("ls "+Input+"*en00*."+Varname+".*.nc")
   else
     flm       = systemfunc("ls "+Input+"*en*."+Varname+".*.nc")
   end if
   end if 

   fm        = addfiles(flm, "r")
   ListSetType (fm, "join")

   timm      = cd_calendar(fm[0]->time,1)
   indm      = ind(timm(:).ge.yst.and.timm(:).le.yed)
   vmdat     = fm[:]->$Varname$(:,indm,{latS:latN},{lonW:lonE})
   
   if(Varname.eq."TGCLDCWP".or.Varname.eq."TGCLDLWP".or.Varname.eq."TGCLDIWP".or.Varname.eq."TMQ") then
     vmdat   = where(vmdat.lt.0.0,0.0,vmdat)
   end if

   vplt  = dim_avg_n_Wrap(dim_avg_n_Wrap(vmdat,1),0)

   delete([/flm,fm,timm,indm,vmdat/])

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(igp.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    plotname = "fig_2d_map_ANN_raw_"+Varname
    wks1 = gsn_open_wks(wtype,plotname)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")
   end if 
   load "../share/var_share_colorbar.ncl"

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 0.5
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeightF = 0.024
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF

   mpres     = res
   resx      = res
   resx1     = res

   res@gsnStringFontHeightF  = FontHeightF*1.2
   res@lbLabelFontHeightF    = FontHeightF*0.6
   res@lbTitleFontHeightF    = FontHeightF*0.8

   res@cnFillOn              = True     ; color plot desired
  ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = False   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.50 ;61
   res@pmLabelBarOrthogonalPosF = 0.20
  ;res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.12 ;0.60
   res@pmLabelBarWidthF         = 0.65 ;0.06
  ;res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
  ;res@lbTitlePosition          = "Right"         ; title location
  ;res@lbTitleDirection         = "Across"        ; letter angle
  ;res@lbTitleAngleF            = 90.          ; title angle

   res@gsnCenterString    = ""
   res@gsnLeftString      = ""
   res@gsnRightString     = "" 
   res@tiMainString       = "" 

;---Control appearance of map.
   mpres@mpMaxLatF                   = latS
   mpres@mpCenterLonF                = 0.0
   mpres@gsnPolar                    = "SH"  
   mpres@mpGridLatSpacingF           = 15.  
   mpres@mpGridLonSpacingF           = 45.  
   mpres@gsnPolarLabelFontHeightF    = FontHeightF 
   mpres@gsnPolarLabelDistance       = 1.08           ; default is 1.04
   mpres@mpLabelsOn                  = False
   mpres@mpPerimOn                   = False
   mpres@mpFillOn                    = False
   mpres@mpOutlineOn                 = True 
   mpres@mpGeophysicalLineColor      = "gray42"
   mpres@mpGeophysicalLineThicknessF = 0.01
   mpres@mpGridLineThicknessF        = 0.01
   mpres@mpOceanFillColor            = "lightskyblue1"
   mpres@mpLandFillColor             = "gray33"
   mpres@mpGridAndLimbOn             = True
   mpres@mpLimbLineThicknessF        = 0.01
   mpres@mpGridLineDashPattern       = 2 
   mpres@mpGridLineThicknessF        = 0.01
   mpres@mpLimbLineThicknessF        = 0.01
   mpres@mpLimbLineDashPattern       = 2

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False
   plot(igp)              = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnAddCyclic                   = True
   res1@cnLevelSelectionMode           = "ExplicitLevels"  
   res1@cnLevels                       = orglev
   res1@cnFillColors                   = orgcolors
   res1@lbBoxLinesOn                   = True
  ;res1@lbBoxSeparatorLinesOn          = False  
  ;res1@lbLabelStride                  = 15
   res1@gsnLeftString                  = "" ;Figstr(igp) + "  " +LabList(igp) 
   res1@gsnLeftStringOrthogonalPosF    = 0.08
  ;res1@gsnLeftStringParallelPosF      = 0.2
   res1@gsnRightString                 = "" ;Varstr+ " ("+Varunt+")";
   res1@gsnRightStringOrthogonalPosF   = 0.055
   res1@gsnCenterString                = Figstr(igp) + "  " +LabList(igp)
   res1@gsnCenterStringOrthogonalPosF  = 0.1
   res1@cnLinesOn                      = plot_contour  
   res1@cnLineLabelsOn                 = False  
   res1@cnInfoLabelOn                  = False
   res1@cnLineLabelDensityF            = 2.0
   res1@cnLineLabelFontHeightF         = 0.018
   res1@cnLineLabelBackgroundColor     = -1
   res1@gsnContourLineThicknessesScale = 1.0
   res1@cnLineThicknessF               = 0.001
   res1@cnLineColor                    = "Black"
   res1@gsnContourZeroLineThicknessF   = 0.2
   res1@gsnContourNegLineDashPattern   = 0
   plot_org1  = gsn_csm_contour(wks1,vplt,res1)
   overlay(plot(igp) ,plot_org1)

   delete([/res,res1,mpres/])

  end do ; end of experiment group loop

;---Draw both plots in a panel
;***************************************
; panel first set of plots
;***************************************
  pres1                          = True
  pres1@gsnFrame                 = True      ; don't advance frame yet
  pres1@gsnMaximize              = False
  pres1@gsnPanelLabelBar         = True 
  pres1@pmLabelBarParallelPosF   = -0.010
  pres1@pmLabelBarOrthogonalPosF = -0.03
  pres1@lbLabelFontHeightF       = FontHeightF * 0.6
  pres1@lbOrientation            = "Horizontal"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.07
  pres1@pmLabelBarWidthF         = 0.75
  pres1@lbTitleFontHeightF       = FontHeightF * 0.6
  pres1@lbTitleString            = Varstr+ " ("+Varunt+")";
  pres1@lbTitlePosition          = "Top"         ; title location
 ;pres1@lbTitleDirection         = "Across"        ; letter angle
 ;pres1@lbTitleAngleF            = 90.          ; title angle
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = FontHeightF*1.1

  gsn_panel(wks1,(/plot(0:)/),(/1,ngrps/),pres1)

 end do

end
