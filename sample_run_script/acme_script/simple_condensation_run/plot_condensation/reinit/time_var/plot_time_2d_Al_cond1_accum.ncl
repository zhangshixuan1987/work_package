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

  Groups       = (/"RKZ_A1_B1_C2_ql19_lmt4_fmin2"/);,"RKZ_A1_B3_C2_ql19_lmt4"/)
  ngrps        = dimsizes(Groups)
  
  Expnam       = (/"A+B1+C2 (opt 19)","A+B3+C2 (opt 19)"/) 
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

   do it        = 0,ntimes-1

   do ih = 0,nh-1,1
    h2sec        = sprinti("%0.5i", ih*3600)
    ncase        = dimsizes(dtimeList)
    h0Suf        = ".cam.h0.0001-01-01-"+h2sec+".nc"
    itime        = 0   ;h0 file has only 1 time step in it

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

     ql        = (/File1s[:]->CLDLIQ(itime,:,:)/) 
     f         = (/File1s[:]->RKZ_f(itime,:,:)/)
     al        = (/File1s[:]->RKZ_Al(itime,:,:)/)
     termB     = (/File1s[:]->RKZ_term_B(itime,:,:)/)
 
     mkflg1     = ql
     mkflg1     = where((f.eq.0).and.(al.lt.0.0).and.(termB.ne.0.0),1.0,0.0)
   
     if (ih.eq.0)then
      mkflg     = mkflg1
     else
      mkflg     = mkflg + mkflg1
     end if 

    delete([/ql,f,al,termB/])

   
;;;;;define plot sorces;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    plot_org   = new(2,graphic)
    fignam="fig_cross_cond1_"+h2sec+"_"+Groups(igp)+"_"+sprinti("%04d",dtimeList(it))+"_"+Regions(irg)+rid
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
    res@vpHeightF= 0.6                    ; change aspect ratio of plot
    res@vpWidthF = 0.9


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


    res@tmYLMode             = "Explicit"                ; explicit labels
    res@tmYLValues           = ispan(30,0,5);(/ 1000, 700,500,400,300,200,100,50,30,10,5,1/)
    res@tmYLLabels           = ""+res@tmYLValues         ; make strings

    res@tiXAxisString        = "Columns"
    res@tiYAxisString        = "Model Vertical Levels"

    res@gsnLeftString        = ""                 ; no text at top of plots
    res@gsnRightString       = ""
    res@trYReverse           = True

    res@cnFillMode           = "RasterFill"      ; turn raster on
;   res@cnFillPalette        = "WhiteBlueGreenYellowRed";"WhViBlGrYeOrRe"; "BlueWhiteOrangeRed" ;"ViBlGrWhYeOrRe" ;"BlueWhiteOrangeRed"

    res@lbLabelBarOn          = True ;False      ; Will turn on in panel later
;   res@pmLabelBarDisplayMode = "Always"
;   res@pmLabelBarOrthogonalPosF=0.02
;   res@pmLabelBarParallelPosF =0.5
    res@pmLabelBarWidthF = 0.1
;   res@pmLabelBarHeightF= 0.10
    res@lbOrientation               = "Vertical"

    res@cnFillOn             = True              ; color plot desired
    res@cnLinesOn            = False             ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour labels
    res@cnMissingValFillColor = "White"

    res1 = res
    res1@gsnCenterString       = "Condition: f = 0, Al < 0, term B ~F34~9~N~~F~ 0"
    ;res1@cnFillColors          = (/"White","Blue"/);(/0,3,5,8,10,13,16/);(/0,8,7,5,3,2,11,12,14,15,17/)
    res1@cnFillPalette         = "WhBlGrYeRe";"precip3_16lev"
    ;res1@cnFillPalette        = "amwg";"precip3_16lev" ;"amwg_blueyellowred" ;"ViBlGrWhYeOrRe" ;"BlueWhiteOrangeRed"
    res1@tiMainString          = Expnam(igp)
    res1@cnLevelSelectionMode  = "ExplicitLevels"   ; set explicit contour levels
    res1@lbLabelAlignment      = "BoxCenters"
    res1@cnLevels              = fspan(0.5,11.5,12)
    res1@lbLabelStrings        = ispan(0,12,1)+""  ; need 1 more label than we have levels
    ;res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    ;res1@cnMinLevelValF       = 0               ; set min contour level
    ;res1@cnMaxLevelValF       = 1.0               ; set max contour level
    ;res1@cnLevelSpacingF      = 0.2               ; set contour spacing
    res1@cnInfoLabelOn         = False
    print(min(mkflg)+ " "+max(mkflg))
    plot_org(0)= gsn_csm_contour(wks,mkflg(:,:),res1)

;---Draw both plots in a panel
    pnres                    = True
    pnres@gsnMaximize        = True
    pnres@gsnPanelLabelBar   = True
    pnres@gsnPanelMainString = " " ;"Time integration: 1 hour"
    pnres@gsnPanelLabelBar   = False ;True                ; add common colorbar
    pnres@gsnPanelMainFontColor   = "Black"

    pnres@gsnPanelMainFontHeightF = 0.025
    pnres@lbLabelFontHeightF      = 0.02               ; make labels smaller

;   pnres@gsnPanelBottom     = 0.05

    gsn_panel(wks,(/plot_org(0)/),(/1,1/),pnres)

  end do ;;end of ih

  end do ;;end it

  end do   ;; end of igrp

  end do ; ensemble loop

 end do ; variable loop 

end do ; end of region loop
end
