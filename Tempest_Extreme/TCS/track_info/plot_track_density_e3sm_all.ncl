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

 ;;; flags to control the plot style
  l_plot_latlon_box = False ; if True, add lat-lon box for SEP region
  datdir    = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TEW/track_info"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList   = (/"ERA5", "E3SM CLIM","E3SM SSP245", "E3SM SSP585"/)
  Groups    = (/"ERA5", "CLIM","SSP245", "SSP585"/) 
  sset      = "set1"
  ngrps     = dimsizes(Groups)
  plot      = new (ngrps, graphic)
  ystr      = 1979
  yend      = 2014
  timtag    = ystr+"-"+yend 
  nyr       = yend - ystr +1

  do i = 0,ngrps-1,1

   fna  = datdir+"/"+Groups(i)+"_"+sset+"_TEW_statistics.nc"
   fla  = addfile(fna, "r")
   lat  = fla->lat
   lon  = fla->lon
   time = fla->time 

   Varnam = "trkdens"
   Varstr = "Track density"
   Varunt = "storm trainsts per month"
   figkey = Varnam 

   ;for NH calculate May-November; 
   ;for SH calculate November-April.
   tim1 = new(nyr*7,integer)
   tim2 = new(nyr*6,integer)
   do ii = 0,nyr-1 
     kk1 = 0 
     kk2 = 0
     do jj = 1, 12
       if(jj.ge.5.and.jj.le.11)then 
         tim1(ii*7+kk1) = toint(sprinti("%04d",ystr + ii) + sprinti("%02d",jj))
         kk1 = kk1 + 1
       end if 
       if(.not.(jj.ge.5.and.jj.lt.11))then
         tim2(ii*6+kk2) = toint(sprinti("%04d",ystr + ii) + sprinti("%02d",jj))
         kk2 = kk2 + 1
       end if
     end do
   end do 

   indx1 = get1Dindex(time,tim1) 
   indx2 = get1Dindex(time,tim2)
   var0 = dim_avg_n_Wrap(fla->$Varnam$(indx1,:,:),0)
   var1 = dim_avg_n_Wrap(fla->$Varnam$(indx2,:,:),0)
   var  = var0
   var({5:90},:)  = var0({5:90},:)
   var({-90:5},:) = var1({-90:5},:) 
   delete([/indx1,indx2,var0,var1,tim1,tim2/])
   delete([/fna,fla,lat,lon,time/])

   print(" ")
   print("........................... ")
   print("igp : " + i) 
   print("........................... ")

   printMinMax(var,False)

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(i.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "track_" + figkey 
    wks1 = gsn_open_wks(wtype,plotname)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")
   end if 

   ;;;set up the plot information;;;;
   Figstr = "" 
   lndcolor = "Black" ;"Grey50"
  ;colorbar setup for difference plots
   cncolors  := (/72,62,55,5,7,10,0,120,14,16,91,85,77/)
  ;cnlevels  := (/0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,8.0,10.0,15.0/)
   cnlevels  := (/1.0,1.5,2.0,2.5,3.0,4.0,6.0,8.0,10.0,15.0,20,30/)

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

   FontHeightF = 0.018
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF

   mpres     = res
   resx      = res
   resx1     = res

   res@gsnStringFontHeightF  = FontHeightF*3.0
   res@lbLabelFontHeightF    = FontHeightF*0.9
   res@lbTitleFontHeightF    = FontHeightF*0.8

   res@cnFillOn              = True     ; color plot desired
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = True   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.50 ;61
   res@pmLabelBarOrthogonalPosF = 0.30
  ;res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.16 ;0.60
   res@pmLabelBarWidthF         = 0.60 ;0.06
   res@lbTitleString            = Varstr+ " ("+Varunt+")";
   res@lbTitlePosition          = "Top" ;"Right"         ; title location
  ;res@lbTitleDirection         = "Across"        ; letter angle
  ;res@lbTitleAngleF            = 90.          ; title angle

   res@gsnCenterString    = ""
   res@gsnLeftString      = ""
   res@gsnRightString     = "" 
   res@tiMainString       = "" 

;---Control appearance of map.
    mpres@gsnMajorLatSpacing     = 30
    mpres@gsnMajorLonSpacing     = 45
    mpres@gsnMinorLatSpacing     = 15
    mpres@gsnMinorLonSpacing     = 15
   ;mpres@mpProjection           = "Robinson"       ; choose projection
    mpres@mpLabelsOn             = False
    mpres@mpPerimOn              = False
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor = lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 0.1          ; thickness of outlines
   ;mpres@mpOutlineOn            = True
   ;mpres@mpOutlineDrawOrder     = "PostDraw"
   ;mpres@mpFillDrawOrder        = "Predraw"
   ;mpres@mpGridAndLimbOn        = True
   ;mpres@mpGridLatSpacingF      = 90.0
   ;mpres@mpGridLonSpacingF      = 180.0
   ;mpres@mpGridLineThicknessF   = 1.0
   ;mpres@mpGridAndLimbDrawOrder="PreDraw"
   ;mpres@mpMinLonF              = 0  ;-180.
   ;mpres@mpMaxLonF              = 360;180.
    mpres@mpMinLatF              = -60  
    mpres@mpMaxLatF              =  60 
    mpres@mpCenterLonF           = 180.0
;   mpres@mpCenterLatF           = -90.     ; This is necessary to get the correct map
;   mpres@mpOceanFillColor       = "lightskyblue1"
;   mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(i) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnAddCyclic                = True
   res1@cnLevelSelectionMode        = "ExplicitLevels"  
  ;res1@cnFillPalette               = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
   res1@cnFillPalette               = "wh-bl-gr-ye-re";"precip2_17lev"
   res1@cnLevels                   := cnlevels
  ;res1@cnFillColors               := cncolors
   res1@lbBoxLinesOn                = True
  ;res1@lbBoxSeparatorLinesOn       = False  
  ;res1@lbLabelStride               = 15
   res1@gsnLeftString               = Figstr + "  " +LabList(i)
   res1@gsnLeftStringOrthogonalPosF = 0.04
  ;res1@gsnLeftStringParallelPosF   = 0.2
   res1@gsnRightString              = "TC track density" 
   res1@gsnRightStringOrthogonalPosF= 0.055
   res1@cnLinesOn                   = False ;plot_contour  
   res1@cnLineLabelsOn              = False  
   res1@cnInfoLabelOn               = False
   res1@cnLineLabelDensityF         = 2.0
   res1@cnLineLabelFontHeightF      = 0.018
   res1@cnLineLabelBackgroundColor  = -1
   res1@gsnContourLineThicknessesScale = 1.0
   res1@cnLineThicknessF            = 0.001
   res1@cnLineColor                 = "Black"
   res1@gsnContourZeroLineThicknessF= 0.2
   res1@gsnContourNegLineDashPattern= 0
  ;res1@sfXArray     = lon
  ;res1@sfYArray     = lat
   plot_org1  = gsn_csm_contour(wks1,var,res1)
   overlay(plot(i) ,plot_org1)

   ;;;;;add a box for selected region;;;;
   
   if(l_plot_latlon_box) then 
      regx = new((/4/),graphic)
      resb                  = True ; polyline mods desired
      resb@gsLineColor      = "Magenta" ; color of lines
      resb@gsLineThicknessF = 2.5 ; thickness of lines
      do i = 0,3
       regx(i)=gsn_add_polyline(wks1,plot_rgd1,xptsa(i:i+1),yptsa(i:i+1),resb)
      end do
   end if

 delete([/res,res1,mpres/])

end do ; end of experiment group loop

;---Draw both plots in a panel
;***************************************
; panel first set of plots
;***************************************
  pres1                          = True
  pres1@gsnFrame                 = True      ; don't advance frame yet
  pres1@gsnMaximize              = False
  pres1@gsnPanelLabelBar         = False
  pres1@pmLabelBarParallelPosF   = -0.010
  pres1@pmLabelBarOrthogonalPosF = 0.0
  pres1@lbLabelFontHeightF       = FontHeightF * 0.45
  pres1@lbOrientation            = "vertical"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.21
  pres1@pmLabelBarWidthF         = 0.05
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = 0.028
  gsn_panel(wks1,(/plot/),(/2,2/),pres1)

end
