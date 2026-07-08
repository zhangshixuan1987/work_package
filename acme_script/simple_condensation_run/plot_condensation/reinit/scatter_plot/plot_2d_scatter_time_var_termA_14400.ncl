load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

  ldatalocal = True
  ldatalocal = False

  if (ldatalocal) then
     ModelOutRoot = "./data/"
  else
     ModelOutRoot = "/pic/projects/uq_climate/zhan391/cnvg_condens_master_revisit/run/"
  end if

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Global"/)
  nregs        = dimsizes(Regions)
  hash         = "E3SM_CONVG"

  Group0       = "RKZ_A1_B1_C2_ql19_lmt4_fmin2_ne30_ne30_FC5AQUAP_intel_constance"+"_"
  Group1       = "RKZ_A1_B1_C2_ql19_lmt4_fmin2_ne30_ne30_FC5AQUAP_intel_constance"+"_"

  config       = "_dycore_mac"
  config1      = "_dycore_mac"
 ;dtimeList    = (/1,2,8,30,120,450,1800/)
  dtimeList    = (/1,450,1800/)
  ntimes       = dimsizes(dtimeList)
  tailnc       = (/"14400.nc","14400.nc","14400.nc"/)
  Expnam       = "A+B1+C2 (opt 19)"
  Strnam       = "Grid points  with ~F8~D~N~~F~T ~F34~3~N~~F~ 10~S~-3~N~ "
  Grpnam       = "14400"
;------
  ncase        = dimsizes(dtimeList)
  h0Suf        = ".cam.h0.0001-01-01-"
  h0Suf0       = ".cam.h0.0001-01-01-14400.nc"

  itime        = 0 ;h0 file has only 1 time step in it

  varNameList  = (/"T","CLDLIQ","CLDICE","CLOUD"/)
  varNameList0 = (/"RKZ_Al","RKZ_term_A","RKZ_term_C"/)
  varNameList1 = (/"CLDLIQ","RKZ_fac"/)

  nvar         = dimsizes(varNameList)

  irs = 1
  ire = 1
  nr = ire - irs + 1
  diffnorm0 = new( (/3,ncase/), "double" )
  diffnormtable = new( (/nr,ncase-1/),"double" )
  slopes = new( nr,"double" )

;;;;;;;;;;;;;;
;##########
do it = 0,ntimes-1

 do irg = 0,nregs-1

  do ir = 0,nr-1

  rid = "_"+sprinti("%02d",ir+irs)
     case1s   = Group0 + "DT"+sprinti("%04d",1)+rid+config
     case2s   = Group0 + "DT"+sprinti("%04d",2)+rid+config
     case00   = Group0 + "DT"+sprinti("%04d",dtimeList(it))+rid+config
     case01   = Group1 + "DT"+sprinti("%04d",dtimeList(it))+rid+config1

     if (ldatalocal) then
        FileName00 = ModelOutRoot+case00+h0Suf
        FileName01 = ModelOutRoot+case01+h0Suf
        FileName1s = ModelOutRoot+case1s+h0Suf
        FileName2s = ModelOutRoot+case2s+h0Suf

     else
        FileName00 = systemfunc("ls -1 " + ModelOutRoot+case00+"/"+case00+h0Suf+tailnc(0))
        FileName01 = systemfunc("ls -1 " + ModelOutRoot+case01+"/"+case01+h0Suf+tailnc(0))
        FileName1s = systemfunc("ls -1 " + ModelOutRoot+case1s+"/"+case1s+h0Suf+tailnc(0))
        FileName2s = systemfunc("ls -1 " + ModelOutRoot+case2s+"/"+case2s+h0Suf+tailnc(0))
     end if

     File00   = addfile(FileName00,"r")
     File01   = addfile(FileName01,"r")
     File1s   = addfile(FileName1s,"r")
     File2s   = addfile(FileName2s,"r")
     T1s    = File1s->T(itime,:,:)
     T2s    = File2s->T(itime,:,:)
     Tdf    = log10(abs(T2s - T1s))
     mkflg     = (/T1s/)
     mkflg     = 0.0
     mkflg     = where((Tdf.ge.-3),1.0,0.0)

     tmp1s0   = File00->$varNameList0(0)$(itime,:,:)*1000.0*86400.0
     f1s0     = File00->$varNameList0(1)$(itime,:,:)*1000.0*86400.0
     x1s0     = File00->$varNameList0(2)$(itime,:,:)*1000.0*86400.0

     tmp1s1   = File01->$varNameList1(0)$(itime,:,:)*1000.0*86400.0
     f1s1     = File01->$varNameList1(1)$(itime,:,:)*1000.0*86400.0
     print(min(tmp1s0)+ " "+max(tmp1s0))
    
     printVarSummary(tmp1s1)

     tmp1s0   = where((mkflg.gt.0.0),tmp1s0,-32766)
     f1s0     = where((mkflg.gt.0.0),f1s0,-32766)
     x1s0     = where((mkflg.gt.0.0),x1s0,-32766)

     tmp1s0@_FillValue=-32766
     f1s0@_FillValue=-32766

     ratio1   = num((f1s0.eq.0.0).and.(x1s0.eq.0.0))*100.0/num(.not.ismissing(tmp1s0))
     ratio2   = num((f1s0.ne.0.0).or.(x1s0.ne.0.0))*100.0/num(.not.ismissing(tmp1s0))

     tmp1d1s0 = ndtooned(tmp1s0)  ;(:,int2d(:,0),int2d(:,1)))
     f1d1s0   = ndtooned(f1s0)  ;(:,int2d(:,0),int2d(:,1)))

     tmp1d1s1 = ndtooned(tmp1s1)  ;(:,int2d(:,0),int2d(:,1)))
     f1d1s1   = ndtooned(f1s1)  ;(:,int2d(:,0),int2d(:,1)))

     tmp1d1s0@_FillValue=-32766
     f1d1s0@_FillValue=-32766
     ind1s   = ind((.not.ismissing(tmp1d1s0))) ;.and.(.not.ismissing(f1d1s0)))
;************************************************
; plotting parameters
;************************************************
   fignam   = "fig_2D_A_B1_C2_scatter_RKZ_termA_"+Grpnam+"_"+sprinti("%04d",dtimeList(it))
   wks      = gsn_open_wks("png",fignam)
   plot     = new(3,graphic)
;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = True              ; Maximize plot
    res@vpHeightF            = 0.6                    ; change aspect ratio of plot
    res@vpWidthF             = 0.6

    res@tmXTOn               = False
    res@tmYROn               = False
;   res@tmYLPrecision        = 5                ; set the precision to 2 significant digits
;   res@tmYLFormat           = "#0^e"   ; don't use "e" symbol, always show +/-


    res@tiMainFontHeightF    = 0.025
    res@tiYAxisFontHeightF   = 0.025
    res@tiXAxisFontHeightF   = 0.025
    res@tmYLLabelFontHeightF = 0.022
    res@tmXBLabelFontHeightF = 0.022
    res@gsnStringFontHeightF = 0.02

    leftstr                     = Strnam 
;   res@gsnLeftString           = leftstr ;var_name0    ; long_name is too long!
;   res@gsnRightString          = leftstr ;var_name0    ; long_name is too long!
    res@gsnCenterString          = leftstr ;var_name0    ; long_name is too long!

    res@tiXAxisString          = "term A (g/kg/day) "              ;"Term A";"Global Temperature (~S~o~N~C)"
    res@tiYAxisString          = " Al (g/kg/day)"


    res@gsnMaximize         = True                   ; maximize plot in frame
    res@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
    res@xyMarkers           = 16                     ; choose type of marker 
    res@xyMarkerSizeF       = 0.01                   ; Marker size (default 0.01)
    res@xyDashPatterns      = 1                      ; solid line 
    res@xyLineThicknesses   = (/1,2/)                ; set second line to 2

    res@trYMinF              = -max(abs(tmp1d1s0)) ;-0.01;min(data(ncol-1,1:))                   ; min value on y-axis
    res@trYMaxF              = max(abs(tmp1d1s0)) ;max(data(ncol-1,1:))    
    ;res@tmYLMode            = "Explicit"
    ;res@tmYLValues          = fspan(0,0.5,6)
    ;res@tmYLLabels          = "" + res@tmYLValues
    ;res@tmYLMinorValues     = fspan(0,0.5,21)

    res@trXMinF              = -max(abs(f1d1s0)) ;-10.0 ;min(data(ncol-1,1:))                   ; min value on y-axis
    res@trXMaxF              =  max(abs(f1d1s0))  ;10.0;max(data(ncol-1,1:))
    ;res@tmXBMode            = "Explicit"	
    ;res@tmXBValues          = fspan(0,1,6) 
    ;res@tmXBLabels          = "" + res@tmXBValues 
    ;res@tmXBMinorValues     = fspan(0,1,21)

    ;res@tmXBFormat          = "#0^e"   ; don't use "e" symbol, always show +/-
    ;res@xyXStyle            = "Irregular"
    ;res@xyXIrregularPoints  = (/1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1/)
    ;res@trXAxisType         = "LogAxis "
    ;res@xyComputeXMin = True
    ;res@tmXBStyle           = "Irregular"
    ;res@tmXTIrregularPoints = (/1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1/)
    ;res@trXMinF              = 0 ;min(data(ncol-1,1:))                   ; min value on y-axis
    ;res@trXMaxF              = 1.0;max(data(ncol-1,1:))

    ;res@gsnLeftString       = dtimeList(it)+"s" ;var_name0    ; long_name is too long!
    res@tiMainString        = "~F8~D~F~t~F~ = " + dtimeList(it)+"s"

    res@xyMarkerColor       = "red"                  ; Marker color
    plot(0)  = gsn_csm_xy (wks,f1d1s0(ind1s),tmp1d1s0(ind1s),res)        ; create plot
    ;res@xyMarkerColor       = "Blue"                  ; Marker color
    ;plt0     = gsn_csm_xy (wks,f1d1s1(ind1s),tmp1d1s1(ind1s),res)        ; create plot
    ;overlay( plot(0),plt0)

    res_text               = True
    res_text@txFontHeightF = 0.02                       ; change font size
    txid1 = gsn_create_text(wks, "termA &  termC = 0 ("+sprintf("%3.1f",ratio1)+"%)", res_text)
    txid2 = gsn_create_text(wks, "termA || termC ~F34~9~N~~F~ 0 ("+sprintf("%3.1f",ratio2)+"%)", res_text)
    amres                  = True
    amres@amJust           = "BottomLeft"
    amres@amParallelPosF   = -0.45    ; This is the center of the plot.
    amres@amOrthogonalPosF =  0.45   ; This is above the top edge of the plot.
    annoid_top = gsn_add_annotation(plot(0), txid1, amres)
    amres@amOrthogonalPosF =  0.35  ; This is just below the previous title.
    annoid_mid = gsn_add_annotation(plot(0), txid2, amres)


;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = True
    pnres@gsnPanelMainString      = "";Expnam;"Time series of DT (2s-1s)"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.022
    pnres@gsnPanelYWhiteSpacePercent = 3
;   pnres@gsnPanelBottom          = 0.05
    gsn_panel(wks,(/plot(0)/),(/1,1/),pnres)
    
   delete([/ind1s/])
  end do ; ir

 end do ; region loop (different Regions)

end do ; end of time loop
end
