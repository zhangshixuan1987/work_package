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
  config       = "_dycore_mac"

  dtimeList    = (/2,8,450,1800/)
  ntimes       = dimsizes(dtimeList)

  Expnam       = "A+B1+C2 (opt 19)"
  Strnam       = "Grid points  with ~F8~D~N~~F~T ~F34~3~N~~F~ 10~S~-3~N~ "

  varNameList0 = (/"T"/) ;(/"RKZ_Al","RKZ_f"/)
  nvar         = dimsizes(varNameList0)

;;;;;plot resources;;;;;;;
  colors    = (/"blue","red","green3","darkorange","Magenta","Magenta4"/)
  lthick    = (/5.0,5.0,5.0,5.0,5.0,5.0,5.0/)

  irs = 1
  ire = 1
  nr = ire - irs + 1
do ih = 0,12,1

  h2sec        = sprinti("%0.5i", ih*3600)
  tailnc       = h2sec+".nc"
  Grpnam       = h2sec
;------
  ncase        = dimsizes(dtimeList)
  h0Suf        = ".cam.h0.0001-01-01-"

  itime        = 0 ;h0 file has only 1 time step in it


  do ir = 0,nr-1

  rid = "_"+sprinti("%02d",ir+irs)
     case1s   = Group0 + "DT"+sprinti("%04d",1)+rid+config
     case2s   = Group0 + "DT"+sprinti("%04d",2)+rid+config

     if (ldatalocal) then
        FileName1s = ModelOutRoot+case1s+h0Suf
        FileName2s = ModelOutRoot+case2s+h0Suf

     else
        FileName1s = systemfunc("ls -1 " + ModelOutRoot+case1s+"/"+case1s+h0Suf+tailnc(0))
        FileName2s = systemfunc("ls -1 " + ModelOutRoot+case2s+"/"+case2s+h0Suf+tailnc(0))
     end if

     File1s    = addfile(FileName1s,"r")
     File2s    = addfile(FileName2s,"r")
     ql2s      = File2s->RKZ_Al(itime,:,:)
     ql1s      = File1s->RKZ_Al(itime,:,:)
     f2s      = File2s->RKZ_f(itime,:,:)
     f1s      = File1s->RKZ_f(itime,:,:)
     mktmp     = (/ql1s/)
     mktmp     = 0.0
     mktmp     = where((ql2s.lt.0.0).or.(ql1s.lt.0.0),1.0,0.0)

     if (ih.eq.0)then
      mkfal = mktmp
     else
      mkfal = mkfal + mktmp
     end if 
    
 end do

end do 


do ih = 1,12,1

  h2sec        = sprinti("%0.5i", ih*3600)
  tailnc       = h2sec+".nc"
  Grpnam       = h2sec
;------
  ncase        = dimsizes(dtimeList)
  h0Suf        = ".cam.h0.0001-01-01-"

  itime        = 0 ;h0 file has only 1 time step in it

;;;;;;;;;;;;;;
;##########
do it = 0,0;ntimes-1

  do ir = 0,nr-1

  rid = "_"+sprinti("%02d",ir+irs)
     case1s   = Group0 + "DT"+sprinti("%04d",1)+rid+config
     case2s   = Group0 + "DT"+sprinti("%04d",2)+rid+config
     case00   = Group0 + "DT"+sprinti("%04d",dtimeList(it))+rid+config

     if (ldatalocal) then
        FileName00 = ModelOutRoot+case00+h0Suf
        FileName1s = ModelOutRoot+case1s+h0Suf
        FileName2s = ModelOutRoot+case2s+h0Suf

     else
        FileName00 = systemfunc("ls -1 " + ModelOutRoot+case00+"/"+case00+h0Suf+tailnc(0))
        FileName1s = systemfunc("ls -1 " + ModelOutRoot+case1s+"/"+case1s+h0Suf+tailnc(0))
        FileName2s = systemfunc("ls -1 " + ModelOutRoot+case2s+"/"+case2s+h0Suf+tailnc(0))
     end if

     File00    = addfile(FileName00,"r")
     File1s    = addfile(FileName1s,"r")
     File2s    = addfile(FileName2s,"r")
     ps1s      = File1s->PS(itime,:)
     T1s       = File1s->T(itime,:,:)
     T2s       = File2s->T(itime,:,:)
     Tdf       = log10(abs(T2s - T1s))
     mkflg     = (/T1s/)
     mkflg     = 0.0
     mkflg     = where((Tdf.ge.-3).and.(mkfal.gt.0.0),1.0,0.0)
     mkavg     = dim_sum_n(mkflg,0)

     lat       = File00->lat
     lon       = File00->lon

     var3d     = Tdf  ;File00->$varNameList0(0)$(itime,:,:) ;*1000.0*86400.0
     var3d     = where(mkflg.gt.0.0, var3d, 32766 )
     var3d@_FillValue = 32766
     print(num(ismissing(var3d)))  
     print(num(mkflg.eq.32766))
         
     ;indplt     = ind(mkavg.gt.0.0)
     varplt      = (/var3d(:,:)/)  ;*1000.0*24*3600.0 ;convert to g/kg/day
     varplt@_FillValue = 32766
     varplt@units = "K"
   
;;;;;define plot sorces;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    fignam="fig_cross_lev_colum_"+varNameList0(0)+"_"+h2sec+rid
    wks      = gsn_open_wks("png",fignam)
    gsn_define_colormap(wks,"GMT_wysiwyg")
    plot_org   = new(2,graphic)

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

    leftstr                  = "log~B~10~N~(~F8~D~N~~F~T) ("+ varplt@units+ ") ("+dtimeList(it)+"s)"
    print(leftstr)
    res@gsnLeftString        = leftstr ;var_name0    ; long_name is too long!
    res@gsnRightString       = "Condition: ~F8~D~N~~F~T ~F34~3~N~~F~ 10~S~-3~N~; Al < 0"

    res@trYReverse           = True

    res@cnFillMode           = "RasterFill"      ; turn raster on
;   res@cnFillPalette        = "WhiteBlueGreenYellowRed";"WhViBlGrYeOrRe"; "BlueWhiteOrangeRed" ;"ViBlGrWhYeOrRe" ;"BlueWhiteOrangeRed"

    res@lbLabelBarOn          = True ;False      ; Will turn on in panel later
;   res@pmLabelBarDisplayMode = "Always"
;   res@pmLabelBarOrthogonalPosF=0.02
;   res@pmLabelBarParallelPosF =0.5
;   res@pmLabelBarWidthF = 0.1
;   res@pmLabelBarHeightF= 0.10
    res@lbOrientation               = "Vertical"

    res@cnFillOn             = True              ; color plot desired
    res@cnLinesOn            = False             ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour labels
    res@cnMissingValFillColor = "White"


    res1 = res
    res1@cnFillPalette         = "WhBlGrYeRe";"precip3_16lev"
   ;res1@cnFillPalette        = "amwg";"precip3_16lev" ;"amwg_blueyellowred" ;"ViBlGrWhYeOrRe" ;"BlueWhiteOrangeRed"
    res1@tiMainString          = Expnam
    res1@cnLevelSelectionMode  = "ExplicitLevels"   ; set explicit contour levels
    ;res1@lbLabelAlignment      = "BoxCenters"
    res1@cnLevels              = fspan(-4,-1,9)
    ;res1@lbLabelStrings        = fspan(-4,-1,9)+""  ; need 1 more label than we have levels
    ;res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    ;res1@cnMinLevelValF       = 0               ; set min contour level
    ;res1@cnMaxLevelValF       = 1.0               ; set max contour level
    ;res1@cnLevelSpacingF      = 0.2               ; set contour spacing
    res1@cnInfoLabelOn         = False
    print(min(mkflg)+ " "+max(mkflg))
    plot_org(0)= gsn_csm_contour(wks, varplt(:,:),res1)

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
   
  end do ; ir

end do ; end of time loop

end do ; end of ih loop

end
