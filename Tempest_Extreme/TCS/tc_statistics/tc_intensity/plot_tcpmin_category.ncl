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
;Saffir-Simpson Hurricane Wind Scale
;Category 4 Hurricane: Winds 130-156 mph (113-136 kt or 209-251 km/hr)
;Category 5 Hurricane: Winds 157 mph or higher (137 kt or higher or 252 km/hr or higher).
;Category 3 Hurricane: Winds 111-129 mph (96-112 kt or 178-208 km/hr).
;Category 2 Hurricane: Winds 96-110 mph (83-95 kt or 154-177 km/hr).
;Category 1 Hurricane: Winds 74-95 mph (64-82 kt or 119-153 km/hr)
;Tropical storm: Winds 39-73 mph (34-63 kt or 63-118 km/hr)
;Tropical depression: Winds < 38 mph (33 kt or 62 km/hr)
;Atkinson and Holliday (1977): Vmax = 6.7(1010-pmin)^0.644
;Knaff and Zehr (2007): Vmax = 6.6(1010-MSLP)^0.65, (appendix A)
;Knaff and Zehr (2007) with updated coefficient: 
;vmax = 4.4*(1010.0 - pmin)^0.76

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

  Class   = (/"TS",   "Cat1",   "Cat2",   "Cat3",   "Cat4",   "Cat5"/)
  Thresh  = (/"994-994", "980-994","965-980","945-965","920-945","860-920"/)
  nclas   = dimsizes(Class)
  kt2ms   = 0.51444

  Ymax    = (/ 45,  12,  12,  12,  12, 12,   12,  12,   12/)
  Ymin    = (/  0,   0,   0,   0,   0,   0,   0,   0,    0/)
  datdir  = "/pscratch/sd/z/zhan391/DARPA_project/TempestExtremes/TCS/track_info"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList = (/"IBTrACS (30yr)", "IBTrACS", "E3SM (CLIM)", "E3SM (NDG UV)", "E3SM (NDG UV+TQ)"/)
  Groups  = (/"IBTrACS",   "IBTrACS",      "CLIM",        "NDGUV",         "NDGUVTQ"/)
  timtag  = (/"1979-2017", "1979-2017",    "2007-2017",   "2007-2017",     "2007-2017"/)
  period  = (/"1988-2017", "2008-2017",    "2008-2017",   "2008-2017",     "2008-2017"/)

  ngrps   = dimsizes(Groups)
  sset    = "set1"
  figkey  = "trkpmin_anaycyc"
  Figstr  = (/"(a)", "(b)", "(d)", "(e)", "(f)", "(g)","(h)", "(i)", "(j)"/)

  do j = 0,nbsin-1,1

   do i = 0,ngrps-1,1
    yrng = str_split(period(i),"-")
    nyr  = toint(yrng(1)) - toint(yrng(0)) + 1
    fna  = datdir+"/"+Groups(i)+"_"+sset+"_Track_info_"+ timtag(i)+ ".nc"
    fla  = addfile(fna, "r")
    latx = fla->lat(:,0)
    lonx = fla->lon(:,0)
    timx = cd_calendar(fla->stim(:,0),-1)
    tim1 = toint(yrng(0)+"01")
    tim2 = toint(yrng(1)+"12") 
    lonx = where(lonx.lt.0,lonx+360.0,lonx)
    inds = ind (latx.ge.Blats(j).and.latx.le.Blatn(j) .and. \
                lonx.ge.Blonw(j).and.lonx.le.Blone(j) .and. \
                timx.ge.tim1.and.timx.le.tim2)

    if(Groups(i).eq."IBTrACS") then 
      vsmc = short2flt(fla->vsmc(inds,:)) * kt2ms
      vpmc = short2flt(fla->vpmc(inds,:))
    else
      vsmc = fla->vsmc(inds,:)
      vpmc = fla->vpmc(inds,:)
      vpmc = vpmc/100.0
     ;vsmc = 4.4*(1010.0 - vpmc)^0.76
    end if 
    lfvmx    = dim_min_n(vpmc,1)

    if(i.eq.0) then 
      vcat = new((/nclas,ngrps/),typeof(lfvmx))
    end if 
    
    do ff = 0,nclas-1,1
     xx = toint(str_split(Thresh(ff),"-"))
     if(ff.ne.0) then
       vcat(ff,i) = num(lfvmx.gt.xx(0).and.lfvmx.le.xx(1))*1.0 /nyr
     else
       vcat(ff,i) = num(lfvmx.gt.xx(0)) *1.0/ nyr 
     end if
    end do
    delete([/fna,fla,latx,lonx,timx,tim1,tim2,vsmc,vpmc,lfvmx,inds,yrng/])
   end do

  ;----------------------------------------------------------------------
  ; Plotting section
  ;----------------------------------------------------------------------
   wtype = "pdf"
   wtype@wkPaperSize  = "A4"
  ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
   plotname  = "track_pmin_cat_" + Basin(j) 
   wks       = gsn_open_wks(wtype,plotname)

   cmap_name = "srip_reanalysis"
   gsn_define_colormap(wks,cmap_name)
   cmap1     = read_colormap_file(cmap_name)
   levels    = (/cmap1(18,:),cmap1(18,:),cmap1(17,:),cmap1(15,:),cmap1(14,:),cmap1(13,:),cmap1(11,:)/)
   colors    = (/19,20,17,15,14,13,11/) ;get_color_index(cmap_name,levels,0.0) 
   patcol    = colors 

   BarChartFillScale = 0.4

   barWidth = (/0.8,0.8,0.8,0.8,0.8,0.8/)/4.0
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
   res@tiYAxisString         = "TCs Per year"
   res@tiXAxisString         = "Month"

   res@trXMinF               = 0
   res@trXMaxF               = dimsizes(Class)*2+1
   res@tmXBMode              = "Explicit"
   res@tmXBValues            = ispan(1,dimsizes(Class)*2,2)
   res@tmXBMinorValues       = ispan(1,dimsizes(Class)*2,1)
   res@tmXBLabels            = Class
   res@tmXBOn                = True            ; Turn off top tickmarks
   res@tmXBLabelsOn          = True            ; Turn off top tickmarks

   res@trYMinF               = Ymin(j)
   res@trYMaxF               = Ymax(j)

   resL = res
   resL@tmYLOn                          = True
   resL@tmYROn                          = False
   resL@gsnXYBarChart                   = True ; create bar chart
   resL@gsnXYBarChartOutlineThicknessF  = 0.5
   resL@gsnXYBarChartFillLineThicknessF = 1.0
   resL@gsnXYBarChartOutlineOnly        = False
   resL@gsnXYBarChartBarWidth          := (/barWidth(0)/) ;(ixp)          ; change bar widths
   resL@gsnXYBarChartColors            := (/colors(0)/)
   resL@gsnXYBarChartFillOpacityF      := (/bopc(0)/)
   resL@gsnXYBarChartPatterns          := (/pattn(0)/)
   resL@xyLineThicknesses              := (/lnthick(0)/)
   resL@xyLineColors                   := (/-1/)
   dum  = new(ngrps,graphic)
   xx1  = ispan(1,dimsizes(Class)*2,2)-0.7
   plot = gsn_csm_xy(wks,xx1,vcat(:,0),resL)

   do ff = 1,ngrps-1,1
     resR = resL
     resR@tmYLOn                = False
     resR@tmYROn                = False
     resR@tmXBOn                = False
     resR@tmXTOn                = False
     resR@gsnLeftString         = "" 
     resR@gsnRightString        = ""
     resR@tiMainString          = ""
     resR@tiYAxisString         = ""
     resR@tiXAxisString         = ""
     resR@gsnXYBarChartBarWidth       := (/barWidth(ff)/) ;(ixp)          ; change bar widths
     resR@gsnXYBarChartColors         := (/colors(ff)/)
     resR@gsnXYBarChartFillOpacityF   := (/bopc(ff)/)
     resR@gsnXYBarChartPatterns       := (/pattn(ff)/)
     resR@xyLineThicknesses           := (/lnthick(ff)/)
     resR@xyLineColors                := (/-1/)
     xx1     = ispan(1,dimsizes(Class)*2,2)-0.7+0.25*ff
     print(vcat(:,ff))
     dum(ff) = gsn_csm_xy(wks,xx1,vcat(:,ff),resR)
     overlay(plot,dum(ff))
     delete(resR)
   end do 
   delete([/xx1,res,resL/])

   xlabel = LabList
   add_labelbar(wks,plot,xlabel,colors(0:ngrps-1),colors(0:ngrps-1),bopc(0:ngrps-1),pattn(0:ngrps-1),pattn(0:ngrps-1))

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
