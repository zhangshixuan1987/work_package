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

  ;varNameList  = (/"T","CLDLIQ","CLDICE","CLOUD"/)
  varNameList0 = (/"T"/) ;(/"RKZ_Al","RKZ_f"/)
  ;varNameList1 = (/"CLDLIQ","RKZ_fac"/)
  nvar         = dimsizes(varNameList0)

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
     mktmp     = where(((ql2s.lt.0.0).and.(f2s.eq.0.0)).or.((ql1s.lt.0.0).and.(f1s.eq.0.0)),1.0,0.0)

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
do it = 0,ntimes-1

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
     ql2s      = File2s->RKZ_Al(itime,:,:)
     ql1s      = File1s->RKZ_Al(itime,:,:)
     Tdf       = log10(abs(T2s - T1s))
     mkflg     = (/T1s/)
     mkflg     = 0.0
     mkflg     = where((Tdf.ge.-3).and.(mkfal.gt.0.0),1.0,0.0)
     mkfla     = dim_sum_n(mkflg,0)

     lat       = File00->lat
     lon       = File00->lon
     p0        = File00->P0
     ps        = File00->PS(itime,:)
     hyai      = File00->hyai
     hybi      = File00->hybi
     nlev      = dimsizes(hyai) -1
     daAll     = hyai(1:nlev) - hyai(0:nlev-1)
     dbAll     = hybi(1:nlev) - hybi(0:nlev-1)

     da = daAll( : )
     db = dbAll( : )

     var3d   = File00->$varNameList0(0)$(itime,:,:) ;*1000.0*86400.0
     print(min(var3d)+ " "+max(var3d))
     printVarSummary(var3d)

     dp       = conform( var3d, da, 0 )* p0 \
              + conform( var3d, db, 0 )* conform( var3d, ps, 1 )
     dp1s     = conform( var3d, da, 0 )* p0 \
              + conform( var3d, db, 0 )* conform( var3d, ps1s, 1 )

     dpw = 0.5*( dp1s + dp )
     psw = 0.5*( ps1s + ps )

     diffvar3d=var3d
     diffvar3d=(var3d-T1s)^2*dpw
     vint = dim_sum_n(diffvar3d, 0 )
     pint = dim_sum_n(dpw, 0 )

     var1d   = log10(vint/pint)
     ind1s   = ind(mkfla.gt.0.0)
     varplt  = var1d
     varplt  = 32766
     varplt(ind1s)=var1d(ind1s)
     varplt@_FillValue = 32766

;************************************************
; plotting parameters
;************************************************
   fignam   = "fig_2D_A_B1_C2_stdv_f_Al_all_"+varNameList0(0)+"_"+Grpnam+"_"+sprinti("%04d",dtimeList(it))
   wks      = gsn_open_wks("png",fignam)
   plot_org   = new(1,graphic) 
   gsn_define_colormap(wks,"GMT_wysiwyg")
;----------------------------------------------------------------------
; Plotting section
    SUBREGION = False
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = True              ; Maximize plot

    res@cnFillOn             = True              ; color plot desired
    res@cnLinesOn            = False             ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour labels
    res@cnFillMode           = "RasterFill"      ; turn raster on
;   res@cnFillPalette        = "WhViBlGrYeOrRe"; "BlueWhiteOrangeRed" ;"ViBlGrWhYeOrRe" ;"BlueWhiteOrangeRed"
    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels             = fspan(-8,2,11)   ; set levels
    res@lbLabelBarOn         = False      ; Will turn on in panel later

    res@tmXTOn               = False
    res@tmYROn               = False
    res@tiMainFontHeightF    = 0.025
    res@tmYLLabelFontHeightF = 0.022
    res@tmXBLabelFontHeightF = 0.022
    res@gsnStringFontHeightF = 0.018
    res@lbLabelFontHeightF   = 0.020
    res@gsnMajorLatSpacing   = 30.0
    res@gsnMajorLonSpacing   = 60.0
   

    res@mpFillOn            = False

    if(SUBREGION) then
      res@mpMinLatF         = minlat
      res@mpMaxLatF         = maxlat
      res@mpMinLonF         = minlon
      res@mpMaxLonF         = maxlon
    else
      res@mpMinLonF         =   0
      res@mpMaxLonF         = 360
      res@mpCenterLonF      = 180     ; This is necessary to get the correct map
    end if

    leftstr                 = "log10(~F33~D~F~"+ varNameList0(0) +"|~B~RMSD~N~)("+dtimeList(it)+"s - 1s)"
    print(leftstr)

    res@gsnLeftString       = leftstr ;var_name0    ; long_name is too long!

    res@gsnAddCyclic        = False

;---Resources for plotting original data
    res@gsnAddCyclic   = False
    res@sfXArray       = lon
    res@sfYArray       = lat
;   res@tiMainString   = "Original CAM-SE grid"

    res@gsnRightString = Expnam 
    plot_org(0)= gsn_csm_contour_map(wks,varplt,res)

;---Draw markers on the plot in the lat/lon locations.
    mkres               = True
    mkres@gsMarkerIndex = 10     ; Filled circle
    mkres@gsMarkerSizeF = 0.01
    mkres@gsMarkerColor = "Black"
    ;dum0 = gsn_add_polymarker(wks,plot_org(0),lon(ind1s),lat(ind1s),mkres)
    ;mkres@gsMarkerColor = "Magenta"
    ;dum1 = gsn_add_polymarker(wks,plot_org(0),lon(ind1s),lat(ind1s),mkres)

    ;res_text               = True
    ;res_text@txFontHeightF = 0.02                       ; change font size
    ;txid1 = gsn_create_text(wks, "Al < 0 ("+sprintf("%3.1f",ratio1)+"%)", res_text)
    ;txid2 = gsn_create_text(wks, "Al > 0 ("+sprintf("%3.1f",ratio2)+"%)", res_text)
    ;amres                  = True
    ;amres@amJust           = "TopRight"
    ;amres@amParallelPosF   =  0.45    ; This is the center of the plot.
    ;amres@amOrthogonalPosF = -0.45   ; This is above the top edge of the plot.
    ;annoid_top = gsn_add_annotation(plot(0), txid1, amres)
    ;amres@amOrthogonalPosF = -0.35  ; This is just below the previous title.
    ;annoid_mid = gsn_add_annotation(plot(0), txid2, amres)


;---Draw both plots in a panel
    pnres                    = True
    pnres@gsnMaximize        = True
    pnres@gsnPanelLabelBar   = True
    pnres@gsnPanelMainString = "Time integration: "+ih+"-h"
    pnres@gsnPanelLabelBar   = True                ; add common colorbar
    pnres@gsnPanelMainFontColor   = "Black"

    pnres@gsnPanelMainFontHeightF = 0.025
    pnres@lbLabelFontHeightF      = 0.02               ; make labels smaller

;   pnres@gsnPanelBottom     = 0.05

    gsn_panel(wks,(/plot_org/),(/1,1/),pnres)
    
   delete([/ind1s/])

  end do ; ir

end do ; end of time loop

end do ; end of ih loop

end
