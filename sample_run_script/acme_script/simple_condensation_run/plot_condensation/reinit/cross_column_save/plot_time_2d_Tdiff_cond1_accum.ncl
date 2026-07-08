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

  dtimeList    = (/1,2/)
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
   do it = 0,1 ;ntimes-1, 1

   Group        = Groups(igp)+"_ne30_ne30_FC5AQUAP_intel_constance"+"_"

   if (ldatalocal) then
     FileName1s = Group+"DT"+sprinti("%04d",dtimeList(0))+rid+configs(igp)+h0Suf
     FileName2s = Group+"DT"+sprinti("%04d",dtimeList(it))+rid+configs(igp)+h0Suf
   else
     FileName1s = ModelOutRoot+Group+"DT"+sprinti("%04d",dtimeList(0))+rid+configs(igp)+"/"+Group+"DT"+sprinti("%04d",dtimeList(0))+rid+configs(igp)+h0Suf
     FileName2s = ModelOutRoot+Group+"DT"+sprinti("%04d",dtimeList(it))+rid+configs(igp)+"/"+Group+"DT"+sprinti("%04d",dtimeList(it))+rid+configs(igp)+h0Suf
   end if
   File1s       = addfile(FileName1s,"r")
   File2s       = addfile(FileName2s,"r")
   ;;;pick the conditional points;;;;;;;;;;;;;;;;;;;;
   T1s    = File1s->T(itime,:,:)
   T2s    = File2s->T(itime,:,:)
   Tdf    = log10(abs(T2s - T1s))
   ql     = File2s->CLDLIQ(itime,:,:)
   f      = File2s->RKZ_f(itime,:,:)
   al     = File2s->RKZ_Al(itime,:,:)
   termB  = File2s->RKZ_term_B(itime,:,:)

   mkflg     = (/ql/)
   mkflg     = 0.0
  ;mkflg     = where((f.eq.0).and.(al.lt.0.0).and.(termB.ne.0.0),1.0,0.0)
  ;mkflg     = where((al.lt.0.0).and.(f.le.1e-3),1.0,0.0)
   mkflg     = where((Tdf.ge.-3),1.0,0.0)
   delete([/ql,f,al,termB/])

   end do ; end of it loop
   
;;;;;define plot sorces;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    plot_org   = new(2,graphic)
    fignam="fig_cross_Tdiff_"+h2sec+"_"+Groups(igp)+"_"+Regions(irg)+rid
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
    res1@gsnCenterString       = "Condition: ~F8~D~N~~F~T ~F34~3~N~~F~ 10~S~-3~N~ "
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

  end do   ;; end of igrp

  end do ; ensemble loop

 end do ; variable loop 

end do ; end of region loop
end
