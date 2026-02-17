;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TC Basin definitions.(Knutson et. al.,2020,BAMS)
; Basin   = (/"GLB","NA","NEP","NWP","NI","SI","SWP","SA"/)
; Blats   = (/  -90,   0,    0,    0,   0, -90,  -90,  -90,  -90/)
; Blatn   = (/   90,  90,   90,   90,  90,   0,    0,    0,    0/)
; Blonw   = (/    0, 265,  180,  100,  30,  20,  135,  295,    0/)
; Blone   = (/  360, 360,  265,  180, 100, 135,  295,  360,   20/)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  datdir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS/track_info"
  basin  = "NA" 
  latr   = (/0,60/)
  lonr   = (/265,360/)
  month  = (/5, 6, 7, 8, 9, 10, 11/)
  season = "MJJASON" ;NH, calculate May-November;   
  ystr   = 2008
  yend   = 2017
  sset   = "set1"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList= (/"IBTrACS","ERA5", "E3SM (CLIM)","E3SM (NDG-UV)", "E3SM(NDG-UVTQ)", "E3SM(NDG-UVTQ)"/)
  Groups = (/"IBTrACS","ERA5", "CLIM",       "NDGUV",         "NDGUVT",         "NDGUVTQ"/) 
  ngrps  = dimsizes(Groups)

  plot      = new (ngrps, graphic)
  nyr       = yend - ystr +1

  do i = 0,ngrps-1,1

   if(i.le.1) then 
     timtag = "1979-2017"
   else 
     timtag = "2007-2017"
   end if 

   fna  = datdir+"/"+Groups(i)+"_"+sset+"_Track_statistics_"+timtag+".nc"
   fla  = addfile(fna, "r")
   lat  = fla->lat
   lon  = fla->lon
   time = fla->time 

   Varnam = "trkdens"
   Varstr = "Track density"
   Varunt = "storm trainsts per month"

   tim1 = new((/nyr*dimsizes(month)/),integer)
   kk1  = 0 
   do yy = ystr, yend 
    do mm = 0,dimsizes(month)-1,1
       tim1(kk1) = toint(sprinti("%04d",yy) + sprinti("%02d",month(mm)))
       kk1 = kk1 + 1
     end do
   end do 
   indx1 = get1Dindex(time,tim1) 
   var   = dim_avg_n_Wrap(fla->$Varnam$(indx1,:,:),0)
   delete([/indx1,tim1/])
   delete([/fna,fla,lat,lon,time/])

   printMinMax(var,False)

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(i.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "track_" + Varnam + "_" + basin
    wks1 = gsn_open_wks(wtype,plotname)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")
   end if 

   ;;;set up the plot information;;;;
   Figstr = "" 
   cncolors  := (/72,62,55,5,7,10,0,120,14,16,91,85,77/)
   cnlevels  := (/0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6/) ;1.5,2.0/)

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

   res@gsnStringFontHeightF  = FontHeightF*1.5
   res@lbLabelFontHeightF    = FontHeightF*0.9
   res@lbTitleFontHeightF    = FontHeightF*0.8

   res@cnFillOn              = True     ; color plot desired
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = False   ; Will turn on in panel later
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
    mpres@mpLimitMode            = "LatLon"
    mpres@mpMinLonF              = lonr(0)
    mpres@mpMaxLonF              = lonr(1)
    mpres@mpMinLatF              = latr(0)
    mpres@mpMaxLatF              = latr(1)
    mpres@mpCenterLonF           = 0.5*(lonr(0)+lonr(1))
    mpres@mpCenterLatF           = 0.5*(latr(0)+latr(1)) ; This is necessary to get the correct map

    mpres@gsnMajorLatSpacing     = 10
    mpres@gsnMajorLonSpacing     = 10
    mpres@gsnMinorLatSpacing     = 5
    mpres@gsnMinorLonSpacing     = 5

    mpres@mpDataBaseVersion      = "MediumRes"         ; better resolution
    mpres@mpOutlineBoundarySets  = "AllBoundaries"     ; more outlines
    mpres@mpDataSetName          = "Earth..4"     
    mpres@mpOutlineBoundarySets  = "GeophysicalAndUSStates"; turn on states

    mpres@mpPerimOn              = True
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor = "Black"  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 1.0          ; thickness of outlines
    mpres@mpOutlineOn            = True
;   mpres@mpOceanFillColor       = "lightskyblue1"
;   mpres@mpLandFillColor        = "gray33"

    mpres@tmYLLabelsOn     = True
    mpres@tmXBLabelsOn     = True
    mpres@tmYLOn           = True
    mpres@tmYROn           = False
    mpres@tmXBOn           = True

   plot(i) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnLeftString               = Figstr + "  " +LabList(i)
   res1@gsnLeftStringOrthogonalPosF = 0.02
  ;res1@gsnLeftStringParallelPosF   = 0.2
   res1@gsnRightString              = "Mean: " +sprintf("%.2f",avg(var))
   res1@gsnRightStringOrthogonalPosF= 0.02

   res1@cnLevelSelectionMode        = "ExplicitLevels"  
  ;res1@cnFillPalette               = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
   res1@cnFillPalette               = "wh-bl-gr-ye-re";"precip2_17lev"
   res1@cnLevels                   := cnlevels
  ;res1@cnFillColors               := cncolors
   res1@lbBoxLinesOn                = True
  ;res1@lbBoxSeparatorLinesOn       = False  
  ;res1@lbLabelStride               = 15
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
   plot_org1  = gsn_csm_contour(wks1,var,res1)
   overlay(plot(i) ,plot_org1)

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
  pres1@pmLabelBarParallelPosF   = 0.0
  pres1@pmLabelBarOrthogonalPosF = -0.02
  pres1@lbLabelFontHeightF       = FontHeightF * 0.6 
  pres1@lbOrientation            = "Horizontal" ;"vertical"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.05
  pres1@pmLabelBarWidthF         = 0.75
  pres1@lbTitleString            = Varstr+ " ("+Varunt+")";
  pres1@lbTitlePosition          = "Top" ;"Right"         ; title location
  pres1@lbTitleFontHeightF       = FontHeightF * 0.6
 ;pres1@lbTitleDirection         = "Across"        ; letter angle
 ;pres1@lbTitleAngleF            = 90.          ; title angle
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = 0.028
  gsn_panel(wks1,(/plot/),(/3,2/),pres1)

end
