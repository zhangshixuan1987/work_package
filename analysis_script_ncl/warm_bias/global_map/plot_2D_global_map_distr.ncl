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
  ldatalocal = True

  l_create_plot_data = True

  if (ldatalocal) then
    ModelOutRoot = "../../data/1-year/"
  else
    ModelOutRoot = "/pic/projects/uq_climate/zhan391/Nudging_analysis/F20TRC5-CMIP6_NUDG_ISSUE/amwg_clim/se_climo/"
  end if

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs  = dimsizes(Regions)

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups    = (/"v1_CTRL", "v1_All_Shorter", "v1_CTRL+DRIB+RadCPL","v1_All_Shorter+DRIB+RadCPL"/)
  ngrps     = dimsizes(Groups)
  Plttim    = (/"360x180_ANN"/)
  nptim     = dimsizes(Plttim)
  plot      = new (6, graphic)

  LabList      = new(ngrps+1,string)
  LabList(0)   = Groups(0)
  LabList(1)   = Groups(1) + " ~F19~:~N~~F~ " + Groups(0) 
  LabList(2)   = Groups(3) + " ~F19~:~N~~F~ " + Groups(2)
  LabList(3)   = Groups(2) + " ~F19~:~N~~F~ " + Groups(0)
  LabList(4)   = Groups(1) + " ~F19~:~N~~F~ " + Groups(3)

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fname     = "./var_info.txt"
  indat     = asciiread(fname,-1,"string")
  VarList   = str_squeeze(str_get_field(indat, 1,","))
  VstList   = str_squeeze(str_get_field(indat, 2,","))
  UnitList  = str_squeeze(str_get_field(indat, 3,","))
  Factors   = stringtofloat(str_get_field(indat, 4,","))
  print(Factors)
  color_grp = stringtointeger(str_get_field(indat, 5,","))
  nvars     = dimsizes(VarList)
   
 ;;; flags to control the plot style
  plotcntl          = True   ; if True, plot the CNTL origianl field on the leftmost panel
  plot_contour      = False  ; if True, plot the contour in the figure

 ;;color for the map boundary line 
 lndcolor = "Black" ;"Grey50"

 do iv = 0,nvars-1, 1

  Varname = VarList(iv)
  Varunt  = UnitList(iv)
  Varstr  = VstList(iv)
  Varfac  = Factors(iv)

 do igp = 0, ngrps-1,1

  if (igp.eq.0)then 
    refcase = Groups (0)
    tescase = Groups (1)
  else if(igp.eq.1) then
    refcase = Groups (2)
    tescase = Groups (3)
  else if(igp.eq.2) then 
    refcase = Groups (0)
    tescase = Groups (2)
  else
    refcase = Groups (3)
    tescase = Groups (1)
  end if 
  end if
  end if 

  ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if(l_create_plot_data) then

     fn0 = ModelOutRoot+refcase+"_"+Plttim(0)+".nc"
     fn1 = ModelOutRoot+tescase+"_"+Plttim(0)+".nc"
     File0  = addfile(fn0,"r")
     File1  = addfile(fn1,"r")

     print(" ")
     print("input : " + fn0)
     print("input : " + fn1)
     print(" ")

     lat    = File0->lat
     lon    = File0->lon

     if(Varname.eq."PRECT") then
        pxcre0300 = File0->PRECC(:,:,:)
        wxcre0300 = File1->PRECC(:,:,:)
        pxcre0300 = File0->PRECC(:,:,:) + File0->PRECL(:,:,:)
        wxcre0300 = File1->PRECC(:,:,:) + File1->PRECL(:,:,:)  
     else if(Varname.eq."NETCF") then
        pxcre0300 = File0->LWCF(:,:,:)
        wxcre0300 = File1->SWCF(:,:,:)
        pxcre0300 = File0->LWCF(:,:,:) + File0->SWCF(:,:,:)
        wxcre0300 = File1->LWCF(:,:,:) + File1->SWCF(:,:,:)
     else
        pxcre0300 = File0->$Varname$(:,:,:)
        wxcre0300 = File1->$Varname$(:,:,:)
     end if
     end if

     if(Varname.eq."TS") then
       pxcre0300 = pxcre0300+Varfac
       wxcre0300 = wxcre0300+Varfac
     else
       pxcre0300 = pxcre0300*Varfac
       wxcre0300 = wxcre0300*Varfac
     end if 
  
     pncre0300 = dim_avg_n_Wrap(pxcre0300,0)
     wncre0300 = dim_avg_n_Wrap(wxcre0300,0)
  
     wncredf1  = wncre0300
     wncredf1  = wncre0300 - pncre0300
     wncredf2  = wncre0300

     pncre0300@long_name = Varstr
     pncre0300@units     = Varunt
     pncre0300@info      = "Reference"
     pncre0300@exp_name  = LabList(igp) 

     wncredf1@long_name  = Varstr
     wncredf1@units      = Varunt
     wncredf1@info      = "Absolute difference"
     wncredf1@exp_name   = LabList(igp)

     wncredf2@long_name  = Varstr
     wncredf2@units      = Varunt
     wncredf2@info      = "Relative difference"
     wncredf2@exp_name   = LabList(igp)

     mskxx  = pncre0300
     mskxx  = where(pncre0300.ne.0.0,pncre0300,32766)
     mskxx@_FillValue = 32766
     wncredf2  = (wncre0300 - mskxx)*100/mskxx
     wncredf2  = where(ismissing(wncredf2),0.0,wncredf2)
     if (Varname .eq. "NETCF") then
       wncredf2 = where(abs(wncredf1).le.2.0,0.0,wncredf2)
     end if  

     ;;;calculate the global mean;;;
     rad      = 4.0*atan(1.0)/180.0
     clat     = cos(lat*rad)
     ogav     = wgt_areaave(pncre0300, clat, 1.0, 0)
     dgav     = wgt_areaave(wncredf1, clat, 1.0, 0)
     rmse     = wgt_arearmse(pncre0300,wncre0300,clat,1.0,0)
     corr     = pattern_cor(pncre0300,wncre0300,tofloat(clat),0)

     pncre0300@mean = wgt_areaave(pncre0300, clat, 1.0, 0)
     wncredf1@mean  = wgt_areaave(wncredf1, clat, 1.0, 0)
     wncredf1@rmse  = wgt_arearmse(pncre0300,wncre0300,clat,1.0,0)
     wncredf1@pcor  = pattern_cor(pncre0300,wncre0300,tofloat(clat),0)
     wncredf2@mean  = wgt_areaave(wncredf2, clat, 1.0, 0)
   
     Rstr1   = "Global  mean:   " + sprintf("%5.2f",pncre0300@mean)
     Rstr2   = "Global  mean:   " + sprintf("%5.2f",wncredf1@mean)
     Rstr2   = Rstr2 + "~C~" + "Global  RMSE: "+ sprintf("%5.2f",wncredf1@rmse)
     Rstr2   = Rstr2 + "~C~" + "Pattern CORR: "+ sprintf("%5.2f",wncredf1@pcor)

     fno = VarList(iv) + ".nc"

     if(igp.eq.0) then
      system("rm -f " + fno)
      flo = addfile(fno,"c")
     else
      flo = addfile(fno,"w")
     end if 

     print(" ")
     print("output : " + fno)
     print(" ")

     vna = VarList(iv) + "_refo_gp" + sprinti("%02d",igp) 
     vnb = VarList(iv) + "_adif_gp" + sprinti("%02d",igp)
     vnc = VarList(iv) + "_rdif_gp" + sprinti("%02d",igp)

     flo->$vna$=pncre0300
     flo->$vnb$=wncredf1
     flo->$vnc$=wncredf2

     delete(flo)
   
  else ;; not l_use_plot_data

     fna = VarList(iv) + ".nc"
     fla = addfile(fna, "r")

     vna = VarList(iv) + "_refo_gp" + sprinti("%02d",igp)
     vnb = VarList(iv) + "_adif_gp" + sprinti("%02d",igp)
     vnc = VarList(iv) + "_rdif_gp" + sprinti("%02d",igp)

     pncre0300=fla->$vna$
     wncredf1=fla->$vnb$
     wncredf2=fla->$vnc$

     Rstr1   = "Global  mean:   " + sprintf("%5.2f",pncre0300@mean)
     Rstr2   = "Global  mean:   " + sprintf("%5.2f",wncredf1@mean)
     Rstr2   = Rstr2 + "~C~" + "Global  RMSE: "+ sprintf("%5.2f",wncredf1@rmse)
     Rstr2   = Rstr2 + "~C~" + "Pattern CORR: "+ sprintf("%5.2f",wncredf1@pcor)

  end if

;----------------------------------------------------------------------
; Plotting section

  ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir   = "./"
  plotname = OUTDir+"fig_2d_map_"+Varname

  system("mkdir -p " + OUTDir)

  wtype = "pdf"
 ;wtype@wkPaperSize  = "A4"
 ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length

  if(igp.eq.0)then 
    wks1 = gsn_open_wks(wtype,plotname)
  end if 

  print("plotname : " + plotname + ".pdf" )

  ;;;set up the plot information;;;;
  load "./var_share_colorbar.ncl"

  ;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot

   res@tmXTOn         = False
   res@tmYROn         = False

   FontHeightF                       = 0.018
   res@tiMainFontHeightF             = FontHeightF*0.95
   res@tmYLLabelFontHeightF          = FontHeightF
   res@tmXBLabelFontHeightF          = FontHeightF
   res@tiXAxisFontHeightF            = FontHeightF
   res@tiYAxisFontHeightF            = FontHeightF
   res@gsnStringFontHeightF          = FontHeightF*1.6

   res@tiMainFontThicknessF  = 1.0
   res@tmBorderThicknessF    = 1.2
   res@tmXBMajorThicknessF   = 1.2
   res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
   res@tmYLMajorThicknessF   = 1.2
   res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0

   mpres     = res

   res@cnFillOn              = True     ; color plot desired
  ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = True   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.5 ;61
   res@pmLabelBarOrthogonalPosF = 0.05
   res@lbOrientation            = "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.63
   res@pmLabelBarWidthF         = 0.1
   res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
   res@lbTitlePosition          = "Right"         ; title location
   res@lbTitleDirection         = "Across"        ; letter angle
   res@lbTitleAngleF            = 90.          ; title angle
   res@lbLabelFontHeightF       = FontHeightF*0.8
   res@lbTitleFontHeightF       = FontHeightF*0.8


   res@gsnLeftString      = "" 
   res@gsnRightString     = "" 
   res@tiMainString       = ""

;---Control appearance of map.
    mpres@gsnMajorLatSpacing     = 30
    mpres@gsnMajorLonSpacing     = 45
    mpres@gsnMinorLatSpacing     = 15
    mpres@gsnMinorLonSpacing     = 15
    mpres@mpLabelsOn             = False
    mpres@mpPerimOn              = False
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor      = lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 1.0    ; thickness of outlines
    mpres@mpMinLonF           = 0  ;-180.
    mpres@mpMaxLonF           = 360;180.
    mpres@mpCenterLonF        = 180 ;150.0
;   mpres@mpCenterLatF        = -90.  ; This is necessary to get the correct map
;   mpres@mpOceanFillColor    = "lightskyblue1"
    mpres@mpLandFillColor     = "gray33"

    mpres@tmYLLabelsOn    = True
    mpres@tmYLOn          = True

    plot(igp+2) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
  if(plotcntl.and.igp.eq.0)then
    plot(0) = gsn_csm_map(wks1,mpres)
    res1=res
    res1@cnFillPalette                 := ocname 
    res1@cnLevelSelectionMode           = "ExplicitLevels"  
    res1@cnLevels                      := orglev
    res1@gsnCenterString                = ""
    res1@gsnCenterStringOrthogonalPosF  = 0.04
   ;res1@gsnCenterStringParallelPosF    = 0.2
    res1@gsnLeftString                  = Varstr+ " ("+Varunt+")" 
    res1@gsnLeftStringOrthogonalPosF    = 0.04
   ;res1@gsnLeftStringParallelPosF      = 0.2
    res1@gsnRightString                 = Rstr1
    res1@gsnRightStringOrthogonalPosF   = 0.04
    res1@tiMainString                   = LabList(0);
    res1@lbBoxLinesOn                   = True
    if(max(abs(orglev)).gt.1.0)then          
      res1@lbLabelStrings              := orglev ;sprintf("%5.1f",orglev)
    else
      res1@lbLabelStrings              := sprintf("%5.2f",orglev)
    end if 
    res1@cnLinesOn                      = plot_contour1 ;True 
    res1@cnLineLabelsOn                 = False  
    res1@cnInfoLabelOn                  = False
    res1@cnLineLabelDensityF            = 2.0
    res1@cnLineLabelFontHeightF         = 0.018
    res1@cnLineLabelBackgroundColor     = -1
    res1@cnSmoothingOn                  = True
    res1@gsnContourLineThicknessesScale = 1.0
    res1@cnLineThicknessF               = 0.001
    res1@cnLineColor                    = "Black"
    res1@gsnContourZeroLineThicknessF   = 0.2
    res1@gsnContourNegLineDashPattern   = 0
    plot_org1  = gsn_csm_contour(wks1,pncre0300,res1)
    overlay(plot(0),plot_org1)
  end if 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@cnLevelSelectionMode          = "ExplicitLevels"  
   ;res2@cnFillPalette                 = dcname 
    res2@cnFillColors                 := diffcolors
    res2@cnLevels                     := diffgrp
    res2@cnLinesOn                     = False    ; turn off contour lines
    res2@cnLineLabelsOn                = False   ; turn off contour labels
    res2@cnInfoLabelOn                 = False
    res2@tiMainString                  = LabList(igp+1)
    if(max(abs(diffgrp)).gt.1.0)then
      res2@lbLabelStrings              = sprintf("%5.1f",diffgrp)
    else
      res2@lbLabelStrings              = sprintf("%5.2f",diffgrp)
    end if
    res2@gsnCenterString               = ""
    res2@gsnCenterStringOrthogonalPosF = 0.04
   ;res2@gsnCenterStringParallelPosF   = 0.2
    res2@gsnLeftString                 = "~F33~D~N~~F~"+ Varstr + "  ("+Varunt+")"
    res2@gsnLeftStringOrthogonalPosF   = 0.04
   ;res2@gsnLeftStringParallelPosF     = 0.2 
    res2@gsnRightString                = Rstr2 ;"Global mean: "+ mean2 ;"unit: "+Varunt
    res2@gsnRightStringOrthogonalPosF  = 0.04
    plot_rgd1      = gsn_csm_contour(wks1,wncredf1,res2)  
    overlay(plot(igp+2),plot_rgd1)
  
   if(plot_contour)
     res21 =    res2
     res21@gsnLeftString            = ""
     res21@gsnRightString           = ""
     res21@gsnCenterString          = ""
     res21@cnFillOn                 = False
     res21@cnLinesOn                = True
     res21@cnLineLabelsOn           = False   ; turn off contour labels
     res21@cnInfoLabelOn            = False
     res21@cnLineLabelDensityF      = 2.0
     res21@cnLineLabelFontHeightF         = FontHeightF
     res21@cnLineLabelBackgroundColor     = -1
     res21@gsnContourLineThicknessesScale = 1.0
     res21@cnLineThicknessF         = 1.0
     res21@cnLineColor              = "Black"
     res21@gsnContourZeroLineThicknessF   = 2.
     res21@gsnContourNegLineDashPattern   = 1
     plot_cont1      = gsn_csm_contour(wks1,wncredf1,res21)
     overlay(plot(igp+2),plot_cont1)
     delete([/res21/])
   end if

   delete([/res,res2,mpres/])

  end do 

;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = False ;True
    pnres@gsnPanelMainString      = "  ";Varstr; "Surface Temperature ("+ts@units+")"
    pnres@gsnPanelLabelBar        = False                ; add common colorbar
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.030
    pnres@gsnPanelXWhiteSpacePercent = 3.0
    pnres@gsnPanelYWhiteSpacePercent = 10.0
    if(plotcntl)then
      gsn_panel(wks1,(/plot(:)/),(/3,2/),pnres)
    else
      gsn_panel(wks1,(/plot(1:)/),(/1,ngrps-1/),pnres)
    end if

 end do ; end of iv loop

end
