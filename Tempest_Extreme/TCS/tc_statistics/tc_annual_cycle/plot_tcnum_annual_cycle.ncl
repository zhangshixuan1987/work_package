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
  datdir  = "/anvil/projects/x-phy220062/x-szhang3/tempest/track_info_ml/"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups  = (/"IBTrACS", "ERA5", "CLIM30",      "CLIM30ML" /)
  LabList = (/"IBTrACS", "ERA5", "CLIM-CTL",    "CLIM-LSTM"/)
  timtag  = (/"1979-2017","1979-2017", "1979-2014","1979-2014"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  figkey  = "trkdens_anaycyc_lstm"

  ystr    = 1979
  yend    = 2014
  nyr     = yend - ystr +1

  Varnam  = "trkdens"
  Figstr  = (/"(d)", "(e)", "(f)"/)
  do j = 0,nbsin-1,1

   do i = 0,ngrps-1,1

    fna   = datdir+"/"+Groups(i)+"_"+sset+"_Track_statistics_"+timtag(i)+".nc"
    fla   = addfile(fna, "r")
    lat   = fla->lat
    lon   = fla->lon
    time  = fla->time
    lmsk  = new((/dimsizes(lat),dimsizes(lon)/),float)
    lmsk!0="lat"
    lmsk!1="lon"
    lmsk&lat = lat 
    lmsk&lon = lon
    if(Basin(j).ne."SA") then 
      lmsk = 0
      lmsk({Blats(j):Blatn(j)},{Blonw(j):Blone(j)}) = 1
    else
      lmsk({Blats(j):Blatn(j)},{Blonw(j):Blone(j)}) = 1
      lmsk({Blats(j+1):Blatn(j+1)},{Blonw(j+1):Blone(j+1)}) = 1
    end if
    ;30-year climatology 
    if(i.eq.0) then 
      indx  = ind(time/100.ge.1988.and.time/100.le.2017)
      vref  = reshape(fla->$Varnam$(indx,:,:),(/30,nmons,dimsizes(lat),dimsizes(lon)/))
      vref  = vref*conform(vref,lmsk,(/2,3/))
      if(i.eq.0) then
        plotdat = new((/3,ngrps+1,nmons/),typeof(vref))
        vtc0    = dim_sum_n_Wrap(dim_sum_n_Wrap(vref,3),2)
        avgx    = dim_avg_n(vtc0,0)
        varx    = sqrt(dim_variance_n(vtc0,0))/sqrt(dimsizes(vtc0(:,0)))
        plotdat(0,0,:) = avgx
        plotdat(1,0,:) = avgx - cdft_t(0.95,dimsizes(vtc0(:,0)))*varx
        plotdat(2,0,:) = avgx + cdft_t(0.95,dimsizes(vtc0(:,0)))*varx
        delete([/avgx,varx,vref,vtc0/])
      end if
      delete(indx)
    end if

    indx = ind(time/100.ge.ystr.and.time/100.le.yend)
    var  = reshape(fla->$Varnam$(indx,:,:),(/nyr,nmons,dimsizes(lat),dimsizes(lon)/))
    var  = var*conform(var,lmsk,(/2,3/))
    vtc  = dim_sum_n_Wrap(dim_sum_n_Wrap(var,3),2)
    avgx = dim_avg_n(vtc,0)
    varx = sqrt(dim_variance_n(vtc,0))/dimsizes(vtc(:,0))
    plotdat(0,i+1,:) = avgx
    plotdat(1,i+1,:) = avgx - cdft_t(0.95,dimsizes(vtc(:,0)))*varx
    plotdat(2,i+1,:) = avgx + cdft_t(0.95,dimsizes(vtc(:,0)))*varx
    delete([/var,vtc,fna,fla,lat,lon,time,indx/])
   end do 

  ;----------------------------------------------------------------------
  ; Plotting section
  ;----------------------------------------------------------------------
   wtype = "pdf"
   wtype@wkPaperSize  = "A4"
  ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
   plotname = "track_dens_ts_" + Basin(j) 
   wks      = gsn_open_wks(wtype,plotname)
   cmap1    = read_colormap_file("srip_reanalysis")
   colors   = (/cmap1(18,:),cmap1(17,:),cmap1(15,:),cmap1(14,:),cmap1(13,:),cmap1(11,:)/)
   patcol   = (/cmap1(18,:),cmap1(17,:),cmap1(15,:),cmap1(14,:),cmap1(13,:),cmap1(11,:)/)

   BarChartFillScale = 0.4

   barWidth = (/0.8,0.8,0.8,0.8,0.8,0.8/)/1.2
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
   res@tiYAxisString         = "Storm transits per month"
   res@tiXAxisString         = "Month"

   res@trXMinF               = 0.0 
   res@trXMaxF               = max(month)+1
   res@tmXBMode              = "Explicit"
   res@tmXBValues            = month
   res@tmXBMinorValues       = month
   res@tmXBLabels            = month
   res@tmXBOn                = True            ; Turn off top tickmarks
   res@tmXBLabelsOn          = True            ; Turn off top tickmarks

   res@trYMinF                = Ymin(j)
   res@trYMaxF                = Ymax(j)

   resL = res
   resL@tmYLOn                          = True
   resL@tmYROn                          = False
   resL@gsnXYBarChart                   = True ; create bar chart
   resL@gsnXYBarChartOutlineThicknessF  = 0.5
   resL@gsnXYBarChartFillLineThicknessF = 1.0
   resL@gsnXYBarChartOutlineOnly        = False
   resL@gsnXYBarChartBarWidth          := 0.80 ;barWidth ;(ixp)          ; change bar widths
   resL@gsnXYBarChartColors            := "grey60" ;(/colors(0,:)/)
   resL@gsnXYBarChartFillOpacityF      := (/bopc(0)/)
   resL@gsnXYBarChartPatterns          := (/pattn(0)/)
   resL@xyLineThicknesses              := (/lnthick(0)/)
   resL@xyLineColors                   := (/-1,-1/)
   plot = gsn_csm_xy(wks,month,plotdat(0,0,:),resL)
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
   dumx = gsn_csm_xy(wks,month,plotdat(0,1:,:),resR)
   overlay(plot,dumx)

  ;************************************************
  ; add error bars
  ;************************************************
   polyres                   = True                       ; marker resources
   polyres@gsLineThicknessF  = lnthick(0)
   polyres@gsLineOpacityF    = bopc(0)
   polyres@gsLineColor       = "Gray25";colors(0,:)
   dum = new((/3,dimsizes(month)/),graphic)
   do t=0,dimsizes(month)-1
     xx  = (/month(t),month(t)/)
     yy  = (/plotdat(1,0,t),plotdat(2,0,t)/)
     xx2 = (/month(t)-0.1,month(t)+0.1/)
     yy1 = (/plotdat(1,0,t),plotdat(1,0,t)/)
     yy2 = (/plotdat(2,0,t),plotdat(2,0,t)/)
     dum(0,t) = gsn_add_polyline(wks,plot,xx,yy,polyres)
     dum(1,t) = gsn_add_polyline(wks,plot,xx2,yy1,polyres)
     dum(2,t) = gsn_add_polyline(wks,plot,xx2,yy2,polyres)
   end do

  ; reate a plot with the area between both curves filled in blue.
   dumy = new(ngrps+1,graphic)
   res2 = resR
   delete(res2@xyLineColors)
   res2@pmLegendDisplayMode  = "Never" 
   res2@xyMarkLineMode      := "Lines"
   res2@xyLineColor          = -1    
   do t = 1,ngrps,1
     res2@gsnXYFillColors    = colors(t,:)
     res2@gsnXYFillOpacities = bopc(t)*0.3
     dumy(t)  = gsn_csm_xy (wks,month,plotdat(1:2,t,:),res2)  ; Create filled XY plot.
     overlay(plot,dumy(t))
   end do 
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
