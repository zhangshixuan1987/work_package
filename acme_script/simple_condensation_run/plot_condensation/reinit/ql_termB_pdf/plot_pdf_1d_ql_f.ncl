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

  Groups       = (/"RKZ_A1_B1_C2_ql19_lmt4_fmin1"/) ;,"RKZ_A1_B1_C2_ql19_lmt4_fmin1"/)
  ngrps        = dimsizes(Groups)
  
  Expnam       = (/"  A1+B1+C2 (opt 19 + fmin = 1E-12)"/) ;,"  A1+B1+C2 (opt 19 + fmin = 1E-12)"/) 
  nexps        = dimsizes(Expnam)
 
  Hash         = "E3SM_CONVG"
  compset      = "FC5AQUAP"
  config       = (/"_dycore_mac","_dycore_mac"/)
 
  dtimeList    = (/1,2,8,30,120,450,1800/)
  ntimes       = dimsizes(dtimeList)

  labels       = (/"  1s", "  2s","  8s","  30s","  120s","  450s","  1800s"/)
  nlabs        = dimsizes(labels)

  varNameList  = (/"RKZ_ql_incld","RKZ_f"/)
  

;;;;;plot resources;;;;;;;
  colors    = (/"black","blue","red","green3","darkorange3","magenta","Brown"/)
  lthick    = (/5.0,5.0,5.0,5.0,5.0,5.0,5.0/)

;;;;;;;;pdf control variables;;;;;;;
  maxbin    =   2.0
  minbin    = -12.0
  nbins     = 240
  xx        = new ( (/nlabs, nbins/), double)
  yy        = new ( (/nlabs, nbins/), double)
  ratio     = new (nlabs, double) ;;;save the ratio to total cells for each type
  labplt    = new (nlabs, string) ;;;labels for plot

do ih = 12,12

  h2sec        = sprinti("%0.5i", ih*3600)
  ncase        = dimsizes(dtimeList)
  h0Suf        = ".cam.h0.0001-01-01-"+h2sec+".nc"
  itime        = 0   ;h0 file has only 1 time step in it

;####################
do igp = 0, ngrps-1

do irg = 0,nregs-1

;ensemble members;
  irs = 1
  ire = 1       
  nr = ire - irs + 1

 do iv = 0,0;nv-1

  do ir = 0,nr-1

    rid = "_"+sprinti("%02d",ir+irs)
    casex1s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(0)) + rid + config(igp) ;1s simulation
    casex2s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(1)) + rid + config(igp) ;2s simulation
    FileNamex1s = systemfunc("ls -1 " + ModelOutRoot+casex1s+"/"+casex1s+h0Suf)
    FileNamex2s = systemfunc("ls -1 " + ModelOutRoot+casex2s+"/"+casex2s+h0Suf)
    
    ;print(FileNamex1s)
    ;print(FileNamex2s)
   
    Filex1s    = addfile(FileNamex1s,"r")
    Filex2s    = addfile(FileNamex2s,"r")
    nfiles    = dimsizes(FileNamex1s)

    tmp1s     = Filex1s->T(itime,:,:)
    tmp2s     = Filex2s->T(itime,:,:)
    tmpdf     = log10(abs(tmp2s - tmp1s))

    dfdt2s    = Filex2s->RKZ_dfdt(itime,:,:)
    f2s       = Filex2s->RKZ_f(itime,:,:)
    al2s      = Filex2s->RKZ_ql(itime,:,:)

    tmpdf     = where((dfdt2s.eq.0.0).and.(f2s.eq.0.0).and.(al2s.ne.0.0),tmpdf,-32766)
    a1D       = ndtooned(tmpdf)
    dsizes_a  = dimsizes(tmpdf)
    print(min(a1D)+ " "+max(a1D))
    indices   = ind_resolve(ind(a1D.ge.-1.0),dsizes_a)
    delete([/a1D/])
    dim_ida  = dimsizes(indices)     
    npts     = dim_ida(0)       ; number of elements > 5 (here 9)
    ndim     = dim_ida(1)       ; rank of "a"  (here 3)
  
    print(npts)
 
  end do ; ensemble loop

end do ; variable loop 

end do ; end of regions loop
end do ; end of times loop

end do ;end of ih loop


;####################
h0Suf1       = ".cam.h0.0001-01-01-*.nc";
itime        = 0   ;h0 file has only 1 time step in it

do igp = 0, ngrps-1

do irg = 0,nregs-1

;ensemble members;
  irs = 1
  ire = 1
  nr = ire - irs + 1
 do iv = 0,0;nv-1

  do ir = 0,nr-1

    rid = "_"+sprinti("%02d",ir+irs)

   do it        = 1,1 ;ntimes-1

    case1s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(0)) + rid + config(igp) ;1s simulation
    case2s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(it)) + rid + config(igp) ;2s simulation
    
    FileName1s = systemfunc("ls -1 " + ModelOutRoot+case1s+"/"+case1s+h0Suf1)
    FileName2s = systemfunc("ls -1 " + ModelOutRoot+case2s+"/"+case2s+h0Suf1)
    ;print(FileName1s)

    File1s    = addfiles(FileName1s(0:23),"r")
    File2s    = addfiles(FileName2s(0:23),"r")
    nfiles    = dimsizes(FileName1s)

    vql1s     = File1s[:]->RKZ_ql*1000.0 ;*3600.0*24.0 ;g/kg
    vql2s     = File2s[:]->RKZ_ql*1000.0 ;*3600.0*24.0 ;g/kg
    unit1     = "g/kg"

    vf1s      = File1s[:]->RKZ_f  ;RKZ_ql*1000.0*3600.0*24.0
    vf2s      = File2s[:]->RKZ_f  ;RKZ_ql*1000.0*3600.0*24.0
    unit2     = "fraction"

    print(min(vf1s)+ " " +max(vf1s))
    ntd      = dimsizes(vql1s(:,0,0))
    qlplt    = new((/2,ntd/),double)
    fplt     = new((/2,ntd/),double)
    hplt     = ispan(0,ntd-1,1)

    print(indices)
    do ipt = 0,npts-1 
  
     printVarSummary(indices)
     nx = indices(ipt,0)
     ny = indices(ipt,1)
     locstr="[lev,column]=["+nx+","+ny+"]"

     qlplt(0,:) = vql1s(:,nx,ny)     
     qlplt(1,:) = vql2s(:,nx,ny)
     fplt(0,:)  = vf1s(:,nx,ny)
     fplt(1,:)  = vf2s(:,nx,ny)
     printVarSummary(vql1s)
     print(nx+ " " +ny)
;;;;;define plot sorces;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    fignam="fig_pdf_lmt4_1h_ql_f_"+sprinti("%0.5i",ipt)+"_"+Groups(igp)+"_"+Regions(irg)+rid
    wks      = gsn_open_wks("png",fignam)
;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot
    res@vpWidthF             = 0.6            ; Change the aspect ratio, but 
    res@vpHeightF            = 0.3            ; make plot as large as possible.

;    res@tmXTOn               = False
;    ;res@tmYROn               = False
;    ;res@tmYLPrecision        = 2        ; set the precision to 2 significant digits
;    ;res@tmYLFormat           = "#0^se"   ; don't use "e" symbol, always show +/-
;
;    ;res@trYLog               = True
;    ;res@trXLog               = False
;    ;res@trYMinF              = 1e-1
;    ;res@trXMaxF              = max(hplt)
;    ;res@trXMinF              = 0
;
;    FontHeight = 0.025
;    numberScaleF = 0.9
;    res@gsnStringFontHeightF = FontHeight*0.8
;    res@tmXBLabelFontHeightF = FontHeight*0.9
;    res@tmYLLabelFontHeightF = FontHeight*0.9
;    res@tiMainFontHeightF    = FontHeight
;    res@tiXAxisFontHeightF   = FontHeight
;    res@tiYAxisFontHeightF   = FontHeight
;
;
;    res@tmXMajorGrid         = True
;    res@tmYMajorGrid         = True
; ;  res@tmXMinorGrid         = True
; ;  res@tmYMinorGrid         = True
;    gridthickness = 0.5
;    res@tmXMinorGridThicknessF = gridthickness
;    res@tmYMinorGridThicknessF = gridthickness
;    res@tmXMajorGridThicknessF = gridthickness
;    res@tmYMajorGridThicknessF = gridthickness
;
;    griddash = 2
;    res@tmXMajorGridLineDashPattern = griddash
;    res@tmYMajorGridLineDashPattern = griddash
;
;   ;res@xyLineThicknesses        = lthick
;   ;res@xyLineColors             = colors
;   ;res@xyMonoDashPattern        = True              ; all solid 
;   ;res@tiYAxisString            = "Number of grid points (#)"
;    res@tiXAxisString            = "Time (h)";"log~B~10~N~(|~F8~D~F~T|)"
;   ;res@gsnRightString           = "" ;"A+B1+C2" ;Expnam(igp)  ; "Total # of grid points = "+ngridpts
;   ;res@gsnLeftString            = "Histogram of |~F8~D~F~T| (unit: K)"
;   ;res@gsnCenterString            = "Histogram of "+varNameList(iv) ;+" (unit: kg/kg)"

    res@gsnDraw          = False
    res@gsnFrame         = False

    FontHeight = 0.025
    numberScaleF = 0.9
    res@gsnStringFontHeightF = FontHeight*0.8
    res@tmXBLabelFontHeightF = FontHeight*0.9
    res@tmYLLabelFontHeightF = FontHeight*0.9
    res@tiMainFontHeightF    = FontHeight
    res@tiXAxisFontHeightF   = FontHeight
    res@tiYAxisFontHeightF   = FontHeight

    res@tmXTOn               = False

    res@pmLegendDisplayMode    = "Always"            ; turn on legend
    res@pmLegendSide           = "Top"               ; Change location of
    res@pmLegendParallelPosF   = .25                 ; move units right
    res@pmLegendOrthogonalPosF = -0.55               ; move units down
    res@pmLegendWidthF         = 0.10               ; Change width and
    res@pmLegendHeightF        = 0.15              ; height of legend.
    res@lgPerimOn              = False                ; turn off/on box around
    res@lgLabelFontHeightF     = .018                ; label font height
   ;res@xyExplicitLegendLabels = (/"  RKZ_ql (2s)","  RKZ_f (2s)"/) ;labplt
    res@tiMainString           = Expnam(igp)
    res@tiXAxisString          = "Model integration time (h)"

; resources for "left" variable
    resL                  = res
    resL@trYMinF          = -max(abs(qlplt))          ; Set min/max of left Y axis
    resL@trYMaxF          =  max(abs(qlplt))
    resL@xyLineColor      = "Brown"
    resL@xyMarkLineMode   = "MarkLines"
    resL@xyMarker         = 7
    resL@xyMarkerSizeF    = 0.012
    resL@xyMarkerColor    = resL@xyLineColor
    resL@xyLineThicknessF = 5.0
    resL@xyDashPattern    = 0                   ; dashed line for 2nd
    resL@tiYAxisFontColor = resL@xyLineColor
   ;resL@tmYRFormat       = "f"            ; Remove trailing zeros from labels
    resL@tiYAxisString    = unit1 ;"Liquid water (g/kg)" 
    resL@pmLegendParallelPosF   = .21                 ; move units right
    resL@pmLegendOrthogonalPosF = -1.1               ; move units down
    resL@xyExplicitLegendLabels = (/"  RKZ_ql (2s)"/) ;labplt

; resources for "right" variable
    resR                  = res                    
    resR@trYMinF          = -1.0 ;min(fplt)           ; Set min/max of right Y axis
    resR@trYMaxF          = 1.0  ;max(fplt)
    resR@xyLineColor      = "NavyBlue"
    resR@xyMarkLineMode   = "MarkLines"
    resR@xyMarker         = 16
    resR@xyMarkerSizeF    = 0.012
    resR@xyMarkerColor    = resR@xyLineColor
    resR@xyDashPattern    = 0              ; Dashed line for 2nd curve
    resR@xyLineThicknessF = 5.0
    ;resR@tmYRFormat      = "f"            ; Remove trailing zeros from labels
    resR@tiYAxisFontColor = resR@xyLineColor
    resR@tiYAxisString    = unit2 ;"Cloud fraction"
    resR@pmLegendParallelPosF   = .20                 ; move units right
    resR@pmLegendOrthogonalPosF = -1.2               ; move units down
    resR@xyExplicitLegendLabels = (/"  RKZ_f (2s)"/) ;labplt
    resR@gsnCenterString        = locstr

    plot = gsn_csm_xy2(wks,hplt,qlplt(1,:),fplt(1,:),resL,resR)

;
; Add two new lines to the XY plot, one associated with the left
; Y axis, and one with the right Y axis. Use the special attribute
; "xy2" returned from gsn_csm_xy2 for the right Y axis.
;

; New line for left axis.

    lnres                   = True
    lnres@gsLineColor       = resL@xyLineColor
    lnres@gsLineThicknessF  = resL@xyLineThicknessF
    lnres@gsLineDashPattern = 2 ;resR@xyDashPattern
    ;dum0 = gsn_add_polyline(wks,plot,hplt,qlplt(1,:),lnres)


; New line for right axis, note "plot@xy2".

    lnres@gsLineColor       = resR@xyLineColor
    lnres@gsLineThicknessF  = resR@xyLineThicknessF
    lnres@gsLineDashPattern = 2 ;resL@xyDashPattern
    ;dum1 = gsn_add_polyline(wks,plot@xy2,hplt,fplt(1,:),lnres)

;
   ; res@pmLegendDisplayMode    = "Always"            ; turn on legend
   ; res@pmLegendSide           = "Top"               ; Change location of 
   ; res@pmLegendParallelPosF   = .25                 ; move units right
   ; res@pmLegendOrthogonalPosF = -0.55               ; move units down
   ; res@pmLegendWidthF         = 0.15               ; Change width and
   ; res@pmLegendHeightF        = 0.22              ; height of legend.
   ; res@lgPerimOn              = False                ; turn off/on box around
   ; res@lgLabelFontHeightF     = .018                ; label font height
 
;;;;;;;generate labels for the plot;;;;;;;;;
   ; do ilab = 0, nlabs-1,1
   ;  if(ilab.le.1)then
   ;   rstr = sprintf("8.0f",ratio(ilab))
   ;   labplt(ilab)= labels(ilab) ;+ " (" +str_squeeze(rstr)+ "%)"
   ;  else
   ;   rstr = sprintf("%8.4f",ratio(ilab))
   ;   labplt(ilab)= labels(ilab) ;+ " (" +str_squeeze(rstr)+ "%)"
   ;  end if

   ;end do
   ;res@xyExplicitLegendLabels = (/"1s","2s"/) ;labplt
   ;res@tiMainString           = Expnam(igp)
   ;plot = gsn_csm_xy (wks, xx, yy, res)
    
;---Draw both plots in a panel
    pnres                    = True
    pnres@gsnMaximize        = True
    pnres@gsnPanelMainString = "";"A+B1+C2: "+dtimeList(it) + "s - "+dtimeList(0)+ "s"  ;Expnam(igp) ;"Time integration: 1-h"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.028
;   pnres@gsnPanelBottom     = 0.05
    gsn_panel(wks,(/plot/),(/1,1/),pnres)


  end do ;; end of ipt loop

  end do ;;end of group loop
   
  end do ; ensemble loop

end do ; variable loop 

end do ; end of regions loop
end do ; end of times loop


end
