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
load "./add_label_bar_midlat.ncl"

begin

  regions = (/"MidLat"/) ;(/"Polar","MidLat","Tropics","Global"/)
  nreg    = dimsizes(regions)
  freq    = "monthly"
  datdir  = "../fig_data/"
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups  = (/"ERA5","CLIM30", "CLIM30ML"/)
  LabList = (/"RD","CR",   "ML(CR)"/)
  ngrps   = dimsizes(Groups)
  fignam  = "fig_pdf_cr_vs_ml_monthly" 
  Varnam  = "trkdens"
  Figstr  = (/"(a)", "(b)","(c)","(d)"/)

  ilev    = 72 
  VarList = (/"U","V","T","Q"/) 
  VarStrs = (/"U","V","T", "Q"/) + " (" + ilev + ")"
  UntList = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~"/)
  StrList = (/"PDF of zonal wind at level " , \
              "PDF of meridional wind at level ", \
              "PDF of temperautre at level " , \
              "PDF of humidity at level " /) + ilev
  nvars   = dimsizes(VarList)
  yst     = 1979 
  yed     = 2014

  do k = 0,nreg-1,1
   reg = regions(k)

   do j = 0,nvars-1,1
    lftstr  = StrList(j) 

    do i = 0,ngrps-1,1
     if(.not.isvar("xx")) then 
       xx = new((/ngrps,50/),double)
     end if 
     if(.not.isvar("yy")) then
       yy = new((/ngrps,50/),double)
     end if
     model   = Groups(i)
     Varname = VarList(j)
     Varunit = UntList(j)
     fl1  = datdir+model+"_"+reg+"_"+Varname+"_lev"+ilev+"_pdf_"+freq+"_"+yst+"-"+yed+".nc"
     f1   = addfile(fl1, "r")
     vnam = Varname+"_"+reg
     var  = f1->$vnam$
     xx(i,:)   = var@bin_center     ; assign appropriate "x" axis values
     yy(i,:)   = (/var/) ;/100.0 
     delete([/fl1,f1,var/])
    end do 

  ;------------------
  ; Plotting section
  ;------------------
   if (j.eq.0) then 
     wtype = "pdf"
     wtype@wkPaperSize   = "A4"
    ;wtype@wkOrientation = "landscape" ;;pictures displayed along the length
     wks      = gsn_open_wks(wtype,fignam+"_"+reg)
     plot     = new(4,graphic)
     gsn_merge_colormaps(wks,"srip_reanalysis","cb_9step")
     gsn_merge_colormaps(wks,"srip_reanalysis","t2m_29lev")
    ;colors   = (/17,39,69,14,21/)
    ;colors   = (/17,39,12,14,21/)
     colors   = (/17,14,13,4,21/)
     colors   = (/17,14,13,29,46/)
    ;colors   = (/17,29,46,29,46/)
     colors   = (/17,14,13,14,13/)
     dashind  = (/0,0,0,0,0/)
     lnthick  = (/5.0,5.0,5.0,4.0,4.0/)*0.5
     markers  = (/4,4,4,4,4,4/) 
     mkszf    = (/0.8,0.8,0.8,0.8,0.8/)*0.01
     bopc     = (/0.8,0.8,0.8,0.8,0.8/)
     pattn    = (/0,0,0,0,0,0/)
   end if 

   ;---Resources to share between both plots
   res                 = True     ; Plot modes desired.
   res@gsnDraw         = False
   res@gsnFrame        = False
   res@gsnMaximize     = False     ; Maximize plot

   res@vpHeightF       = 0.45        ; Changes the aspect ratio
   res@vpWidthF        = 0.85
   res@vpXF            = 0.10        ; change start locations
   res@vpYF            = 0.75        ; the plot

  ;res@vpHeightF       = 0.65        ; Changes the aspect ratio
  ;res@vpWidthF        = 0.80
  ;res@vpXF            = 0.10        ; change start locations
  ;res@vpYF            = 0.75        ; the plot


   res@tmXTOn          = False
   res@tmYROn          = True
   res@tmYRBorderOn    = True         ; Habilita a borda do eixo y direito (YR). 
   res@tmXTBorderOn    = True         ; Habilita a borda do eixo x superior (XB). 

   FontHeight               = 0.030
   res@gsnStringFontHeightF = FontHeight*1.15
   res@tmXBLabelFontHeightF = FontHeight*1.1
   res@tmYLLabelFontHeightF = FontHeight*1.1
   res@tiMainFontHeightF    = FontHeight
   res@tiXAxisFontHeightF   = FontHeight*1.15
   res@tiYAxisFontHeightF   = FontHeight*1.15
 
   res@tmXMajorGrid         = False
   res@tmYMajorGrid         = False
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

   res@gsnLeftString         = Figstr(j) + " " + lftstr 
   res@gsnRightString        = ""
   res@tiMainString          = ""
   res@tiYAxisString         = "PDF (%) " ;+ Varname + " (%)" 

   if(Varname.eq."U") then 
     res@trXMinF                = -10. ; 239. ;200.0 
     res@trXMaxF                =  12. ;306. ;320.0 
     res@trYMinF                = 0.
     res@trYMaxF                = 8.
     res@tiXAxisString          = "Zonal Wind (" +Varunit +")" 
     res@gsnXYBarChartBarWidth := 0.39
   end if 

   if(Varname.eq."V") then
     res@trXMinF                = -12. ; 239. ;200.0
     res@trXMaxF                =  12. ;306. ;320.0
     res@trYMinF                = 0.
     res@trYMaxF                = 12.
     res@tiXAxisString          = "Meridonal wind (" +Varunit +")"
     res@gsnXYBarChartBarWidth := 0.39
   end if

   if(Varname.eq."T") then
     res@trXMinF                = 258. ; 239. ;200.0
     res@trXMaxF                = 310. ;306. ;320.0
     res@trYMinF                = 0.
     res@trYMaxF                = 6.
     res@tiXAxisString          = "Temperature (" +Varunit +")"
     res@gsnXYBarChartBarWidth := 1.0
   end if

   if(Varname.eq."Q") then
     res@trXMinF                = 0.0 ;0.0
     res@trXMaxF                = 16.0 ;20.0
     res@trYMinF                = 0. 
     res@trYMaxF                = 6.
     res@tiXAxisString          = "Specific humidity (" +Varunit +")"
     res@gsnXYBarChartBarWidth := 0.31
   end if 

  ;res@tmXBMode              = "Explicit"
   res@tmXBOn                = True            ; Turn off top tickmarks
   res@tmXBLabelsOn          = True            ; Turn off top tickmarks

   resL = res
   resL@tmYLOn                          = True
   resL@tmYROn                          = False
   resL@gsnXYBarChart                   = True ; create bar chart
   resL@gsnXYBarChartOutlineThicknessF  = 0.1
   resL@gsnXYBarChartFillLineThicknessF = 0.1
   resL@gsnXYBarChartOutlineOnly        = False
   resL@gsnXYBarChartColors            := "grey35" ;(/colors(0,:)/)
   resL@gsnXYBarChartFillOpacityF      := (/bopc(0)/)/10.0
   resL@gsnXYBarChartPatterns          := (/pattn(0)/)
   resL@xyLineThicknesses              := (/lnthick(0)/)/10.0
   resL@xyLineColors                   := (/-1,-1/)
   plot(j) = gsn_csm_xy(wks,xx(0,:),yy(0,:),resL)
   xlabel  = LabList(0) ;+" (Obs.)"  
   add_labelbar(wks,plot(j),xlabel,"grey60","grey60",bopc(0),pattn(0),pattn(0))

   resR                        = res
   resR@tmYLOn                 = False
   resR@tmYROn                 = True
   resR@tmXBOn                 = False
   resR@tmXMajorGrid           = False
   resR@tmYMajorGrid           = False
   resR@tmXMinorGrid           = False
   resR@tmYMinorGrid           = False
   resR@gsnLeftString          = ""
   resR@gsnRightString         = ""
   resR@gsnCenterString        = ""
   resR@tiXAxisString          = ""
   resR@gsnXYBarChart          = False
   resR@xyMarkLineMode        := "Lines";"MarkLines" 
   resR@xyMarkers             := markers(1:ngrps-1)
   resR@xyMarkerSizes         := mkszf(1:ngrps-1)
   resR@xyMarkerThicknesses   := lnthick(1:ngrps-1)*0.5
   resR@xyMarkerColors        := colors(1:ngrps-1)
   resR@xyDashPatterns        := dashind(1:ngrps-1)
   resR@xyLineThicknesses     := lnthick(1:ngrps-1)
   resR@xyLineColors          := colors(1:ngrps-1)
   resR@xyExplicitLegendLabels = "   " + LabList(1:ngrps-1)
   resR@pmLegendDisplayMode    = "Always"            ; turn on legend
   resR@pmLegendSide           = "Top"               ; Change location of
   resR@pmLegendWidthF         = 0.12
   resR@pmLegendHeightF        = 0.12   
   resR@pmLegendParallelPosF   = 0.835                ; move units right
   resR@pmLegendOrthogonalPosF = -0.65               ; move units down
   resR@lgPerimOn              = False               ; turn off box around
   resR@lgLabelFontHeightF     = FontHeight*0.9       ; label font height
   dumx = gsn_csm_xy(wks,xx(1:ngrps-1,:),yy(1:ngrps-1,:),resR)
   overlay(plot(j),dumx)

 end do 

;---Draw both plots in a panel
 pnres                         = True
 pnres@gsnMaximize             = False
 pnres@gsnPanelMainString      = ""
 pnres@gsnPanelMainFontColor   = "Black"
 pnres@gsnPanelMainFontHeightF = 0.022
;pnres@gsnPanelBottom          = 0.05
 pnres@gsnPanelYWhiteSpacePercent = 5
 pnres@gsnPanelXWhiteSpacePercent = 5

 gsn_panel(wks,(/plot/),(/2,2/),pnres)

end do 

end
