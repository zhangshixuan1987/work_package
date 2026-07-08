load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

  ldatalocal   = True
  ldatalocal   = False
  getminmaxind = True
  getminmaxind = False

  if (ldatalocal) then
     ModelOutRoot = "./data/"
  else
     ModelOutRoot = "/pic/projects/uq_climate/zhan391/cnvg_condens_master_revisit/run/"
  end if

;;;;;;;;;;select the regions to plot;;;;;;;;;;;;;
  Regions      = (/"Global"/) ;,"Land","Ocean"/)
  nregs        = dimsizes(Regions)

  Groups       = (/"RKZ_A1_B1_C2_ql19_lmt4_fmin2","RKZ_A1_B3_C2_ql19_lmt4_reint"/)
  ngrps        = dimsizes(Groups)
  
  Expnam       = (/"A+B1+C2","A+B3+C2"/) 
  nexps        = dimsizes(Expnam)
 
  Hash         = "E3SM_CONVG"
  compset      = "FC5AQUAP"
  configs      = (/"_dycore_mac","_dycore_mac"/)

  dtimeList    = (/1,1800/)
  ntimes       = dimsizes(dtimeList)

  nh = 13 ; nhours 
  stanum = new((/ntimes,ngrps,nh/),float)
  xx = ispan(0,nh-1,1)
;;;;;plot resources;;;;;;;
  colors    = (/"blue","red","green3","darkorange","Magenta","Magenta4"/)
  lthick    = (/5.0,5.0,5.0,5.0,5.0,5.0,5.0/)

;####################
do irg = 0,nregs-1

;ensemble members;
  irs = 1
  ire = 1
  nr = ire - irs + 1
 do iv = 0,0;nv-1

  do ir = 0,nr-1
    rid = "_"+sprinti("%02d",ir+irs)

   do igp = 0, ngrps-1

   do ih = 0,nh-1,1
    h2sec        = sprinti("%0.5i", ih*3600)
    ncase        = dimsizes(dtimeList)
    h0Suf        = ".cam.h0.0001-01-01-"+h2sec+".nc"
    itime        = 0   ;h0 file has only 1 time step in it

    do it        = 0,ntimes-1

    case1s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(it)) + rid + configs(igp) 
    FileName1s = systemfunc("ls -1 " + ModelOutRoot+case1s+"/"+case1s+h0Suf)
    File1s    = addfiles(FileName1s,"r")
    nfiles    = dimsizes(FileName1s)

    ;;read the information shared with all CAM files;;; 
    ;;you only need to do this once;;;;;;;;;;;;;;;;;;;;
     ifil = 0
     area     = File1s[ifil]->area
     lat       = File1s[ifil]->lat
     lon       = File1s[ifil]->lon
     hyam      = File1s[ifil]->hyam ; read from a file the mid-layer coef
     hybm      = File1s[ifil]->hybm ; read from a file
     ps        = File1s[ifil]->PS   ; surface pressure [Pa]
     pres      = File1s[ifil]->lev
;    landfrac  = File1s[ifil]->LANDFRAC(itime,:)     ; (hdim)
     p0        = 100000. ; since ps is in Pa or [ f->P0]
     pm        = pres_hybrid_ccm(ps,p0,hyam,hybm)

     ql        = File1s[:]->CLDLIQ
     f         = File1s[:]->RKZ_f
     al        = File1s[:]->RKZ_Al
     termB     = File1s[:]->RKZ_term_B
 
     mkflg     = ql
     mkflg     = where((f.eq.0).and.(al.lt.0.0),1.0,0.0)
     ;mkflg     = where((f.eq.0).and.(al.lt.0.0),1.0,0.0)
     stanum(it,igp,ih) = num(mkflg.eq.1)*100.0/num(mkflg.ge.0.0) 
    delete([/ql,f,al,termB/])

    end do ;;end it

   end do ;;end of ih 
   
 end do   ;; end of igrp
   
   print(stanum)
   
;;;;;define plot sorces;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    fignam="fig_time_1h_cond1_"+Regions(irg)+rid
    wks      = gsn_open_wks("png",fignam)
;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = True              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False
   ;res@tmYLPrecision        = 2        ; set the precision to 2 significant digits
   ;res@tmYLFormat           = "#0^se"   ; don't use "e" symbol, always show +/-

    res@trYLog               = False
    res@trXLog               = False
    res@trYMaxF              = max(stanum)+1.0 ;6.0
    res@trYMinF              = 0.0
    res@trXMaxF              = 12.0
    res@trXMinF              = 0.0

    FontHeight = 0.025
    numberScaleF = 0.9
    res@gsnStringFontHeightF = FontHeight*0.8
    res@tmXBLabelFontHeightF = FontHeight*0.9
    res@tmYLLabelFontHeightF = FontHeight*0.9
    res@tiMainFontHeightF    = FontHeight
    res@tiXAxisFontHeightF   = FontHeight
    res@tiYAxisFontHeightF   = FontHeight


    res@tmXMajorGrid         = True
    res@tmYMajorGrid         = True
 ;  res@tmXMinorGrid         = True
 ;  res@tmYMinorGrid         = True
    gridthickness = 0.5
    res@tmXMinorGridThicknessF = gridthickness
    res@tmYMinorGridThicknessF = gridthickness
    res@tmXMajorGridThicknessF = gridthickness
    res@tmYMajorGridThicknessF = gridthickness

    griddash = 2
    res@tmXMajorGridLineDashPattern = griddash
    res@tmYMajorGridLineDashPattern = griddash

    res@xyLineThicknesses        = lthick
    res@xyLineColors             = colors
    res@tiYAxisString            = "Ratio to total grid points (%)"
    res@tiXAxisString            = "Time (h)"
    res@gsnRightString           = "";"A+B1+C2" ;Expnam(igp)  ; "Total # of grid points = "+ngridpts
   ;res@gsnLeftString            = "Histogram of ql_incld (unit: kg/kg)"
    res@gsnCenterString          = "Condition: f = 0, Al < 0, term B ~F34~9~N~~F~ 0"

   ;res@gsnXYBarChart            = True              ; Create bar plot
   ;res@gsnXYBarChartOutlineOnly = True
    res@pmLegendDisplayMode    = "Always"            ; turn on legend
    res@pmLegendSide           = "Top"               ; Change location of 
    res@pmLegendWidthF         = 0.1               ; Change width and
    res@pmLegendHeightF        = 0.08               ; height of legend.
    res@lgPerimOn              = False                ; turn off/on box around
    res@lgLabelFontHeightF     = .018                ; label font height

    res0 = res
    res0@xyExplicitLegendLabels = " " + Expnam + " (      " + dtimeList(0) + "s)"
    res0@xyDashPattern  = 0                         ; Make curves all solid
    res0@xyMarkLineMode = "MarkLines"               ; Markers *and* lines
    res0@xyMarkers       = (/6,11,16/)               ; 3 different markers
    res0@xyMarkerColors := res@xyLineColors
    res0@pmLegendParallelPosF   = .25                 ; move units right
    res0@pmLegendOrthogonalPosF = -0.25               ; move units down
    plot0 = gsn_csm_xy (wks, xx, stanum(0,:,:), res0)

    res1 = res
    res1@xyExplicitLegendLabels = " " + Expnam + " (" + dtimeList(1) + "s)"
    res1@xyDashPattern   = 1                         ; Make curves all solid
    res1@xyMarkLineMode  = "MarkLines"               ; Markers *and* lines
    res1@xyMarkers       = (/6,11,16/)               ; 3 different markers
    res1@xyMarkerColors := res@xyLineColors
    res1@pmLegendParallelPosF   = .25                 ; move units right
    res1@pmLegendOrthogonalPosF = -0.37               ; move units down
    plot1 = gsn_csm_xy (wks, xx, stanum(1,:,:), res1)
    
    overlay(plot0,plot1)
;---Draw both plots in a panel
    pnres                    = True
    pnres@gsnMaximize        = True
    pnres@gsnPanelMainString =" " ;"A+B1+C2" ;+dtimeList(it) + "s - "+dtimeList(0)+ "s"  ;Expnam(igp) ;"Time integration: 1-h"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.028
;   pnres@gsnPanelBottom     = 0.05
    gsn_panel(wks,(/plot0/),(/1,1/),pnres)
   
  end do ; ensemble loop

 end do ; variable loop 

end do ; end of region loop
end
