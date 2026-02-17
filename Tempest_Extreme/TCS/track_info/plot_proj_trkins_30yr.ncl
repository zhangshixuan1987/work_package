;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This script used to plot e3sm diagnostics for TCs;;;;;;
;TC Basin definitions.(Knutson et. al.,2020,BAMS)
;North Atlantic: 0-90N, 265-360E
;Northeast Pacific: 0-90N, 180-265E
;Northwest Pacific: 0-90N, 100-180E
;North Indian: 0-90N, 30-100E
;South Indian: 90S-0, 20-135E
;Southwest Pacific: 90S-0, 135-295E
;South Atlantic(SA1+SA2): 90S-0, South America (295-360E) to Africa (0-20E)
;The North Atlantic-northeast Pacific boundary is on a
;diagonal tracing a path through Mexico and Central America.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ;TC Basin definitions.(Knutson et. al.,2020,BAMS)
  l_plot_latlon_box = False ; if True, add lat-lon box for region
  Basin = (/"NA","NEP","NWP","NI","SI","SWP","SA1","SA2"/)
  Blats = (/   0,    0,    0,   0, -90,  -90,  -90,   90/)
  Blatn = (/  90,   90,   90,  90,   0,    0,    0,    0/)
  Blonw = (/ 265,  180,  100,  30,  20,  135,  295,    0/)
  Blone = (/ 360,  265,  180, 100, 135,  295,  360,   20/)
  nbsin = dimsizes(Basin)

  datdir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS/track_info"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList   = (/"Present-day (PD)", "SSP245 - PD", "SSP585 - PD"/)
  Groups    = (/"CLIM",             "SSP245",      "SSP585"/)
  timtag    = (/"1979-2014",        "1979-2014",   "1979-2014"/)
  ngrps     = dimsizes(Groups)
  sset      = "set1"
  figkey    = "trkins_30year_proj"

  Figstr    = (/"(a)", "(b)", "(c)"/)
  ystr      = 1981
  yend      = 2010
  nyr       = yend - ystr +1

  do i = 0,ngrps-1,1

   fna  = datdir+"/"+Groups(i)+"_"+sset+"_Track_statistics_"+timtag(i)+".nc"
   fla  = addfile(fna, "r")
   lat  = fla->lat
   lon  = fla->lon
   time = fla->time 

   Varnam = "tcint"
   Varstr = "Track intensity"
   Varunt = "m/s"

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
   var0 = dim_max_n_Wrap(fla->$Varnam$(indx1,:,:),0)
   var1 = dim_max_n_Wrap(fla->$Varnam$(indx2,:,:),0)
   var  = var0
   var({5:90},:)  = var0({5:90},:)
   var({-90:5},:) = var1({-90:5},:) 
   delete([/indx1,indx2,var0,var1,tim1,tim2/])
   delete([/fna,fla,lat,lon,time/])

   print(" ")
   print("........................... ")
   print("igp : " + i) 
   print("........................... ")

   if(i.eq.0) then
     vref = var 
   end if 

   if(i.gt.0) then
     var = var - vref
   end if

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
    wks1     = gsn_open_wks(wtype,plotname)
    plot     = new (ngrps, graphic)
    regx     = new ((/ngrps,nbsin,4/), graphic) 

    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")

   end if 

   ;;;set up the plot information;;;;
   lndcolor   = "Black" ;"Grey50"
   cnlevels  := (/0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25/) * 30.0 
   cnlevel2  := (/-0.4,-0.3,-0.2,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.2,0.3,0.4/) * 30.0

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

   res@gsnStringFontHeightF  = FontHeightF*2.0
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
   res@lbTitleString            = ""
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
   mpres@mpMinLatF              = -90  
   mpres@mpMaxLatF              =  90 
   mpres@mpCenterLonF           = 180.0

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(i) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnAddCyclic                = True
   res1@cnLevelSelectionMode        = "ExplicitLevels"  
  ;res1@cnFillPalette               = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
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

   if(i.le.0) then 
     res1@gsnLeftString               = Figstr(i) + "  " +LabList(i)
     res1@gsnLeftStringOrthogonalPosF = 0.04
     res1@gsnRightString              = Varstr 
     res1@gsnRightStringOrthogonalPosF= 0.04
     res1@cnFillPalette               = "wh-bl-gr-ye-re";"precip2_17lev"
     res1@cnLevels                   := cnlevels
     res1@lbTitleString               = "("+Varunt+")";
     plot_org1  = gsn_csm_contour(wks1,var,res1)
     overlay(plot(i) ,plot_org1)
   else
     res1@gsnLeftString               = Figstr(i) + "  " +LabList(i)
     res1@gsnLeftStringOrthogonalPosF = 0.04
     res1@gsnRightString              = Varstr+" response"
     res1@gsnRightStringOrthogonalPosF= 0.04
     cmap = read_colormap_file("BlueDarkRed18")
     cmap(8:9,:) = 0
     res1@cnFillPalette              := cmap 
     res1@cnLevels                   := cnlevel2
     res1@lbTitleString               = "("+Varunt+")";
     plot_org1  = gsn_csm_contour(wks1,var,res1)
     overlay(plot(i) ,plot_org1)
   end if 

   ;;;;;add a box for selected region;;;;
   if(l_plot_latlon_box) then 
    do ii = 0,nbsin-1,1
      lats  = Blats(ii)
      late  = Blatn(ii)
      lons  = Blonw(ii)
      lone  = Blone(ii)
      ;draw box for the region
      yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
      xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line
      resb                  = True ; polyline mods desired
      resb@gsLineColor      = "Magenta" ; color of lines
      resb@gsLineThicknessF = 2.5 ; thickness of lines
      do jj = 0,3
       regx(i,ii,jj)=gsn_add_polyline(wks1,plot(i),xptsa(jj:jj+1),yptsa(jj:jj+1),resb)
      end do
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
  gsn_panel(wks1,(/plot/),(/1,3/),pres1)

end
