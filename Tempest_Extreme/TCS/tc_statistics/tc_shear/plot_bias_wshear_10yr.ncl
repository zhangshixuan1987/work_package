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
  l_plot_latlon_box = True ; if True, add lat-lon box for region
  Basin  = (/"NA","NEP","NWP","NI","SI","SWP"/) ;,"SA1","SA2"/)
  Blats  = (/   0,    0,    0,   0, -90,  -90/) ;,  -90,   90/)
  Blatn  = (/  90,   90,   90,  90,   0,    0/) ;,    0,    0/)
  Blonw  = (/ 265,  180,  100,  30,  20,  135/) ;,  295,    0/)
  Blone  = (/ 360,  265,  180, 100, 135,  295/) ;,  360,   20/)
  nbsin  = dimsizes(Basin)
  odatdir = ""
  mdatdir = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups  = (/"ERA5", "CLIM",        "NDGUV",         "NDGUVTQ"/)
  LabList = (/"ERA5", "E3SM (CLIM)", "E3SM (NDG UV)", "E3SM (NDG UV+TQ)"/)
  timtag  = (/"1979-2017", "2007-2017","2007-2017","2007-2017","2007-2017"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  figkey  = "trkdens_10year"

  Figstr  = (/"(a)", "(b)", "(c)", "(d)", "(e)", "(f)","(g)","(h)"/)
  ystr    = 2008
  yend    = 2017
  nyr     = yend - ystr +1

  Varnam  = "trkdens"
  Varstr  = "Track density"
  Varunt  = "storm trainsts per month"
  vbtrk   = new ((/nbsin/),double)

  do i = 0,ngrps-1,1

   fna  = datdir+"/"+Groups(i)+"_"+sset+"_Track_statistics_"+timtag(i)+".nc"
   fla  = addfile(fna, "r")
   lat  = fla->lat
   lon  = fla->lon
   time = fla->time 

   ;for NH calculate May-November ;for SH calculate October-May.
   tim1 = new(nyr*7,integer)
   tim2 = new(nyr*8,integer)
   do ii = 0,nyr-1 
     kk1 = 0 
     kk2 = 0
     do jj = 1, 12
       if(jj.ge.5.and.jj.le.11)then 
         tim1(ii*7+kk1) = toint(sprinti("%04d",ystr + ii) + sprinti("%02d",jj))
         kk1 = kk1 + 1
       end if 
       if(.not.(jj.ge.6.and.jj.le.9))then
         tim2(ii*8+kk2) = toint(sprinti("%04d",ystr + ii) + sprinti("%02d",jj))
         kk2 = kk2 + 1
       end if
     end do
   end do 

   indx1 = get1Dindex(time,tim1) 
   indx2 = get1Dindex(time,tim2)
   var0  = dim_avg_n_Wrap(fla->$Varnam$(indx1,:,:),0)
   var1  = dim_avg_n_Wrap(fla->$Varnam$(indx2,:,:),0)
   var   = var0
   var({0:90},:)  = var0({0:90},:)
   var({-90:0},:) = var1({-90:0},:)
   delete([/indx1,indx2,var0,var1,tim1,tim2/])
   delete([/fna,fla,lat,lon,time/])

   do j = 0,nbsin-1,1
    if(Basin(j).ne."SA") then
      lmsk = var
      lmsk = 0
      lmsk({Blats(j):Blatn(j)},{Blonw(j):Blone(j)}) = 1
    else
      lmsk = var
      lmsk = 0
      lmsk({Blats(j):Blatn(j)},{Blonw(j):Blone(j)}) = 1
      lmsk({Blats(j+1):Blatn(j+1)},{Blonw(j+1):Blone(j+1)}) = 1
    end if
    var0     = where(lmsk.eq.0,var@_FillValue,var)
    vbtrk(j) = dim_sum_n_Wrap(dim_sum_n_Wrap(var0,1),0)
    delete([/var0,lmsk/])
   end do 
   print(" ")
   print("........................... ")
   print("igp : " + i) 
   print("........................... ")

   if(i.eq.0) then 
     vref = var 
   end if 
   
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
    text_id  = new ((/ngrps,nbsin/), graphic) 

    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")

   end if 

   ;;;set up the plot information;;;;
   lndcolor   = "Black" ;"Grey50"
   cnlevels  := (/0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25/) ;1.5,2.0/)
   cnlevel2  := (/-0.6,-0.4,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.4,0.6/)
   cnlevel3  := (/-1.2,-0.8,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.8,1.2/)

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False
   res@vpXF           = 0.08                  ; position and sizes
   res@vpYF           = -0.2 ;0.98                  ; for XY plot

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

   res@gsnStringFontHeightF  = FontHeightF*2.5
   res@lbLabelFontHeightF    = FontHeightF*0.9
   res@lbTitleFontHeightF    = FontHeightF*0.8

   res@cnFillOn              = True     ; color plot desired
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = False   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.50 ;61
   res@pmLabelBarOrthogonalPosF = 0.22
  ;res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.16 ;0.60
   res@pmLabelBarWidthF         = 0.60 ;0.06
   res@lbTitleString            = ""
   res@lbTitlePosition          = "Top" ;"Right"         ; title location
  ;res@lbTitleDirection         = "Across"        ; letter angle
  ;res@lbTitleAngleF            = 90.          ; title angle
   res@lbBoxLinesOn                = True
  ;res@lbBoxSeparatorLinesOn       = False
  ;res@lbLabelStride               = 15

   res@gsnCenterString    = ""
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
   mpres@mpGeophysicalLineColor = "Grey25"  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = 0.05          ; thickness of outlines

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
   mpres@mpMinLatF              = -65  
   mpres@mpMaxLatF              =  65 
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

   res1@gsnLeftString               = Figstr(i) + "  " +LabList(i)
   res1@gsnLeftStringOrthogonalPosF = 0.04
   res1@gsnRightString              = Varstr 
   res1@gsnRightStringOrthogonalPosF= 0.04
   res1@cnFillPalette               = "wh-bl-gr-ye-re";"precip2_17lev"
   res1@cnLevels                   := cnlevels
   res1@lbTitleString               = "("+Varunt+")";

   plot_org1  = gsn_csm_contour(wks1,var,res1)
   overlay(plot(i) ,plot_org1)

   ;;;;;add a box for selected region;;;;
   if(l_plot_latlon_box) then 
    do ii = 0,nbsin-1,1
      lats  = Blats(ii)
      late  = Blatn(ii)
      lons  = Blonw(ii)
      lone  = Blone(ii) ;-2.0
      ;draw box for the region
      yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
      xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line
      resb                  = True ; polyline mods desired
      resb@gsLineColor      = "Black";"Magenta" ; color of lines
      resb@gsLineThicknessF = 0.5 ; thickness of lines
      do jj = 0,3
       regx(i,ii,jj) = gsn_add_polyline(wks1,plot(i),xptsa(jj:jj+1),yptsa(jj:jj+1),resb)
      end do
      ;;add text 
      txres               = True
      txres@txFontHeightF = FontHeightF*0.8
      txres@txFontColor   = "Red"
      txres@txBackgroundFillColor = "White"
      txres@txJust        = "CenterLeft"
      text = Basin(ii)+ " (" +sprintf("%.1f", vbtrk(ii)) +")"
      latx = 0.63*(late+lats)
      lonx = lons+5.0
      if(Basin(ii).eq."NEP") then
         lonx = lons+35.0
      end if 
      if(Basin(ii).eq."SWP") then
         lonx = lons+90.0
         latx = -10.0
      end if
    ;  if(latx.lt.0) then 
    ;     latx = latx + 1.0
    ;  else
    ;     latx = latx - 1.0
    ;  end if 
      text_id(i,ii) = gsn_add_text(wks1,plot(i),text,lonx,latx,txres)
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
  pres1@gsnMaximize              = True 
  pres1@gsnPanelLabelBar         = True
  pres1@pmLabelBarParallelPosF   = 0.012
  pres1@pmLabelBarOrthogonalPosF = -0.02
 ;pres1@lbOrientation            = "vertical"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.03 ;0.21
  pres1@pmLabelBarWidthF         = 0.35 ;0.05
  pres1@lbLabelFontHeightF       = FontHeightF * 0.45
  pres1@lbTitleFontHeightF       = FontHeightF * 0.50
  pres1@lbTitleString            = "Track density (storms/month/4 degree box)"
  pres1@lbTitlePosition          = "Top" ;"Right"         ; title location
 ;pres1@lbTitleDirection         = "Across"        ; letter angle
 ;pres1@lbTitleAngleF            = 90.          ; title angle
  pres1@lbBoxLinesOn             = True
 ;pres1@lbBoxSeparatorLinesOn    = False
 ;pres1@lbLabelStride            = 15
  pres1@gsnPanelBottom           = 0.02	
  pres1@gsnPanelTop              = 0.95
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = 0.02
  gsn_panel(wks1,(/plot/),(/ngrps,1/),pres1)

end
