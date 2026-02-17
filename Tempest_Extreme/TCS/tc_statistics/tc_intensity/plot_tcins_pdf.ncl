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
load "./add_label_bar.ncl"

begin

 ;TC Basin definitions.(Knutson et. al.,2020,BAMS)
  l_plot_latlon_box = False ; if True, add lat-lon box for region
  Fignm   = (/"(a)", "(b)", "(c)","(d)","(e)","(f)","(g)","(h)"/)
  Basin   = (/"GLB","NA","NWP","NEP","NI","SI","SWP","SA"/)
  Basnm   = (/"Global","North Atlatic","Northwest Pacific","Northest Pacific",\
              "North Indian","South Indian","Southwest Pacific","South Atlatic"/)
  Blats   = (/  -90,   0,    0,    0,   0, -90,  -90,  -90,  -90/)
  Blatn   = (/   90,  90,   90,   90,  90,   0,    0,    0,    0/)
  Blonw   = (/    0, 265,  100,  180,  30,  20,  135,  295,    0/)
  Blone   = (/  360, 360,  180,  265, 100, 135,  295,  360,   20/)
  nbsin   = dimsizes(Basin)

  Ymax    = (/500, 180, 200, 200,  60, 200, 200, 100, 20/)
  Ymin    = (/  0,   0,   0,   0,   0,   0,   0,   0,  0/)
  month   = (/  1,   2,   3,   4,   5,   6,   7,   8,  9, 10, 11, 12/)
  nmons   = dimsizes(month)

  datdir  = "/pscratch/sd/z/zhan391/DARPA_project/TempestExtremes/TCS/track_info"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList = (/"ERA5",      "ERA5",      "E3SM (CLIM)", "E3SM (NDG UV)", "E3SM (NDG UV+TQ)"/)
  Groups  = (/"ERA5",      "ERA5",      "CLIM",        "NDGUV",         "NDGUVTQ"/)
  timtag  = (/"1979-2017", "1979-2017", "2007-2017",   "2007-2017",     "2007-2017"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  figkey  = "trkdens_anaycyc"

  ystr    = 2008
  yend    = 2017
  nyr     = yend - ystr +1

  Varnam  = "trkdens"
  Figstr  = (/"(d)", "(e)", "(f)"/)
  do j = 0,nbsin-1,1

   do i = 0,ngrps-1,1

    fna  = datdir+"/"+Groups(i)+"_"+sset+"_Track_info_"+ timtag(i)+ ".nc"
    fla  = addfile(fna, "r")
    latx = fla->lat(:,0)
    lonx = fla->lon(:,0)
    timx = cd_calendar(fla->stim(:,0),-1)
    tim1 = toint(ystr+"01")
    tim2 = toint(yend+"12") 
    lonx = where(lonx.lt.0,lonx+360.0,lonx)
    if(Basin(j).eq."SA") then 
      inds = ind ((latx.ge.Blats(j).and.latx.le.Blatn(j) .and. \
                   lonx.ge.Blonw(j).and.lonx.le.Blone(j) ).or. \
                  (latx.ge.Blats(j+1).and.latx.le.Blatn(j+1).and. \
                   lonx.ge.Blonw(j+1).and.lonx.le.Blone(j+1)).and. \
                   timx.ge.tim1.and.timx.le.tim2)
    else 
      inds = ind (latx.ge.Blats(j).and.latx.le.Blatn(j) .and. \
                  lonx.ge.Blonw(j).and.lonx.le.Blone(j) .and. \
                  timx.ge.tim1.and.timx.le.tim2)
    end if
    vsmc = fla->vsmc(inds,:)
    vpmc = fla->vpmc(inds,:)
    delete([/inds/])

    if(i.eq.0) then 
      if(Basin(j).eq."SA") then
        inds = ind ((latx.ge.Blats(j).and.latx.le.Blatn(j) .and. \
                     lonx.ge.Blonw(j).and.lonx.le.Blone(j) ).or. \
                    (latx.ge.Blats(j+1).and.latx.le.Blatn(j+1).and. \
                     lonx.ge.Blonw(j+1).and.lonx.le.Blone(j+1)).and. \
                    timx.ge.198801.and.timx.le.201712)
      else
        inds = ind (latx.ge.Blats(j).and.latx.le.Blatn(j) .and. \
                    lonx.ge.Blonw(j).and.lonx.le.Blone(j) .and. \
                    timx.ge.198801.and.timx.le.201712)
      end if
      vsmr = fla->vsmc(inds,:)
      vpmr = fla->vpmc(inds,:)
      delete([/inds/])
    end if 
   
    if(i.eq.0) then 
      nBin        = 51 
      bin         = fspan(0, 100, nBin)  ; 20 is arbitrary
      opt         = True
      opt@bin_min = 0
      opt@bin_max = 100.
      a1d     = ndtooned (vsmr)
      indx    = ind(a1d.ne.vsmr@_FillValue)
      apdf    = pdfx(a1d(indx), nBin, opt)
      xx      = new ((/ngrps+1, nBin/), typeof(apdf))
      yy      = new ((/ngrps+1, nBin/), typeof(apdf))
      xx(0,:) = apdf@bin_center     ; assign appropriate "x" axis values
      yy(0,:) = (/apdf/) /100.0
      delete([/a1d,indx,apdf,vsmr,vpmr/])
    end if 
    a1d       = ndtooned (vsmc)
    indx      = ind(a1d.ne.vsmc@_FillValue)
    apdf      = pdfx(a1d(indx), nBin, opt)
    xx(i+1,:) = apdf@bin_center     ; assign appropriate "x" axis values
    yy(i+1,:) = (/apdf/) /100.0      
    delete([/a1d,indx,apdf/])
    delete([/fna,fla,latx,lonx,timx,tim1,tim2,vsmc,vpmc/])
   end do 

  ;------------------
  ; Plotting section
  ;------------------
   wtype = "pdf"
   wtype@wkPaperSize  = "A4"
  ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
   plotname = "track_ins_pdf_" + Basin(j) 
   wks      = gsn_open_wks(wtype,plotname)
   cmap1    = read_colormap_file("srip_reanalysis")
   colors   = (/cmap1(18,:),cmap1(17,:),cmap1(15,:),cmap1(14,:),cmap1(13,:),cmap1(11,:)/)
   patcol   = (/cmap1(18,:),cmap1(17,:),cmap1(15,:),cmap1(14,:),cmap1(13,:),cmap1(11,:)/)

   BarChartFillScale = 0.4

   barWidth = (/0.8,0.8,0.8,0.8,0.8,0.8/)/4
   pattn    = (/0,0,0,0,0,0/)
   bopc     = (/0.8,0.8,0.8,0.8,0.8,0.8/)
   markers  = (/4,4,4,4,4,4/) ;11,8,7/)
   mkszf    = (/0.8,0.8,0.8,0.8,0.8,0.8/)*0.01
   dashind  = (/0,0,0,0,0,0,0/)
   lnthick  = (/5.0,5.0,5.0,5.0,5.0,5.0,5.0/)*0.8                                   

   ;---Resources to share between both plots
   res                 = True     ; Plot modes desired.
   res@gsnDraw         = False
   res@gsnFrame        = False
   res@gsnMaximize     = False     ; Maximize plot

   res@vpHeightF       = 0.50        ; Changes the aspect ratio
   res@vpWidthF        = 0.80
   res@vpXF            = 0.10        ; change start locations
   res@vpYF            = 0.75        ; the plot

   res@tmXTOn          = False
   res@tmYROn          = True
   res@tmYRBorderOn    = True         ; Habilita a borda do eixo y direito (YR). 
   res@tmXTBorderOn    = True         ; Habilita a borda do eixo x superior (XB). 

   FontHeight               = 0.024
   res@gsnStringFontHeightF = FontHeight*1.0
   res@tmXBLabelFontHeightF = FontHeight*1.0
   res@tmYLLabelFontHeightF = FontHeight*1.0
   res@tiMainFontHeightF    = FontHeight
   res@tiXAxisFontHeightF   = FontHeight
   res@tiYAxisFontHeightF   = FontHeight
 
   res@tmXMajorGrid         = True
   res@tmYMajorGrid         = True
   res@tmXMinorGrid         = False
   res@tmYMinorGrid         = False
   griddash = 2
   res@tmXMajorGridLineDashPattern = griddash
   res@tmYMajorGridLineDashPattern = griddash

   if(wtype.eq."png") then
    gridthickness = 8.0
   else
    gridthickness = 1.0
   end if
   res@tmBorderThicknessF       = gridthickness
   res@tmXBMajorThicknessF      = gridthickness
   res@tmXBMinorThicknessF      = gridthickness
   res@tmYLMajorThicknessF      = gridthickness
   res@tmYLMinorThicknessF      = gridthickness
   res@tmYRMajorThicknessF      = gridthickness
   res@tmYRMinorThicknessF      = gridthickness
   res@tmXMinorGridThicknessF   = gridthickness/2.0
   res@tmYMinorGridThicknessF   = gridthickness/2.0
   res@tmXMajorGridThicknessF   = gridthickness/2.0
   res@tmYMajorGridThicknessF   = gridthickness/2.0

   BorderThick = 1.5
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   res@gsnLeftString         = Fignm(j) + " " + Basnm(j)
   res@gsnRightString        = ""
   res@tiMainString          = ""
   res@tiYAxisString         = "Normalized pdf"
   res@tiXAxisString         = "Maximum wind speed (m s~S~-1~N~)"

   res@trXMinF               = 0.0 
   res@trXMaxF               = 120.0 
  ;res@tmXBMode              = "Explicit"
  ;res@tmXBValues            = month
  ;res@tmXBMinorValues       = month
  ;res@tmXBLabels            = month
   res@tmXBOn                = True            ; Turn off top tickmarks
   res@tmXBLabelsOn          = True            ; Turn off top tickmarks

   res@trYMinF               = 0 
   res@trYMaxF               = 0.3 ;30

   resL = res
   resL@tmYLOn                          = True
   resL@tmYROn                          = False
   resL@gsnXYBarChart                   = True ; create bar chart
   resL@gsnXYBarChartOutlineThicknessF  = 0.5
   resL@gsnXYBarChartFillLineThicknessF = 1.0
   resL@gsnXYBarChartOutlineOnly        = False
   resL@gsnXYBarChartBarWidth          := 0.7;(ixp)          ; change bar widths
   resL@gsnXYBarChartColors            := "grey60" ;(/colors(0,:)/)
   resL@gsnXYBarChartFillOpacityF      := (/bopc(0)/)
   resL@gsnXYBarChartPatterns          := (/pattn(0)/)
   resL@xyLineThicknesses              := (/lnthick(0)/)
   resL@xyLineColors                   := (/-1,-1/)
   plot = gsn_csm_xy(wks,xx(0,:),yy(0,:),resL)
   xlabel= LabList(0) +" (30yr)"  
   add_labelbar(wks,plot,xlabel,"grey60","grey60",bopc(0),pattn(0),pattn(0))

   resR                        = res
   resR@tmYLOn                 = False
   resR@tmYROn                 = True
   resR@tmXBOn                 = False
   resR@pmLegendDisplayMode    = "Always"            ; turn on legend
   resR@pmLegendSide           = "Top"               ; Change location of 
   resR@pmLegendParallelPosF   = .20                 ; move units right
   resR@pmLegendOrthogonalPosF = -0.52               ; move units down
   resR@pmLegendWidthF         = 0.10
   resR@pmLegendHeightF        = 0.16
   resR@lgPerimOn              = False               ; turn off box around
   resR@lgLabelFontHeightF     = FontHeight*.7       ; label font height
   resR@xyExplicitLegendLabels = " " + LabList
   resR@tmXMajorGrid           = False
   resR@tmYMajorGrid           = False
   resR@tmXMinorGrid           = False
   resR@tmYMinorGrid           = False
   resR@gsnLeftString          = ""
   resR@gsnRightString         = ""
   resR@gsnCenterString        = ""
   resR@tiXAxisString          = ""
   resR@gsnXYBarChart          = False
   resR@xyMarkLineMode        := "MarkLines" 
   resR@xyMarkers             := markers(1:ngrps)
   resR@xyMarkerSizes         := mkszf(1:ngrps)
   resR@xyMarkerThicknesses   := lnthick(1:ngrps)*0.5
   resR@xyMarkerColors        := colors(1:ngrps,:)
   resR@xyDashPatterns        := dashind(1:ngrps)
   resR@xyLineThicknesses     := lnthick(1:ngrps)
   resR@xyLineColors          := colors(1:ngrps,:)
   dumx = gsn_csm_xy(wks,xx(1:,:),yy(1:,:),resR)
   overlay(plot,dumx)

  ;---Draw both plots in a panel
   pnres                         = True
   pnres@gsnMaximize             = True
   pnres@gsnPanelMainString      = ""
   pnres@gsnPanelMainFontColor   = "Black"
   pnres@gsnPanelMainFontHeightF = 0.022
  ;pnres@gsnPanelBottom          = 0.05
   gsn_panel(wks,(/plot/),(/1,1/),pnres)

end do 

end
