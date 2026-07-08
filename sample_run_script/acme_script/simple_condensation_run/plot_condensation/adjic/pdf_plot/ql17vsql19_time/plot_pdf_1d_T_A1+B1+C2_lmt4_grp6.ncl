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
     ModelOutRoot  = "/pic/projects/uq_climate/zhan391/cnvg_condens_master_adqldens/run/"
  end if

;;;;;;;;;;select the regions to plot;;;;;;;;;;;;;
  Regions      = (/"Global"/) ;,"Land","Ocean"/)
  nregs        = dimsizes(Regions)

  Groups       = (/"RKZ_A1_B1_C2_ql17_lmt4","RKZ_A1_B1_C2_ql19_lmt4"/)

  ngrps        = dimsizes(Groups)
  
  Expnam       = (/"A+B1+C2_lmt4","A+B1+C2_lmt5"/) 
  nexps        = dimsizes(Expnam)
 
  Hash         = "E3SM_CONVG"
  compset      = "FC5AQUAP"
  configs      = (/"_dycore_mac","_dycore_mac_nodens"/)

  dtimeList    = (/1,2,8,30,120,450,1800/)
  ntimes       = dimsizes(dtimeList)

  labels       = (/"  1s","  2s","  8s","  30s","  120s","  450s","  1800s"/)  ; create explicit labels
  nlabs        = dimsizes(labels)


;;;;;plot resources;;;;;;;
  colors    = (/"black","blue","red","green3","darkorange","Magenta","Magenta4"/)
  lthick    = (/5.0,5.0,5.0,5.0,5.0,5.0,5.0/)

;;;;;;;;;;variables for condensation model plot;;;;;;;;;;;;;;;
  varNameList  = (/"ql_incld"/) 
  h0Suf        = ".cam.h0.0001-01-01-03600.nc";
  h0Suf0       = ".cam.h2.0001-01-01-03600.nc";
  itime        = 0   ;h0 file has only 1 time step in it

;;;;;;;;pdf control variables;;;;;;;
  maxbin    = 2.0
  minbin    = -12.0
  nbins     = 200
  xx        = new ( (/nlabs, nbins/), double)
  yy        = new ( (/nlabs, nbins/), double)
  ratio     = new (nlabs, double) ;;;save the ratio to total cells for each type
  labplt    = new (nlabs, string) ;;;labels for plot

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
   do it        = 0,ntimes-1
    case1s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(0)) + rid + configs(igp) ;1s simulation
    case2s  = Groups(igp) + "_ne30_ne30_FC5AQUAP_intel_constance_DT" + sprinti("%04d",dtimeList(it)) + rid + configs(igp) ;2s simulation
    
    FileName1s = systemfunc("ls -1 " + ModelOutRoot+case1s+"/"+case1s+h0Suf)
    FileName2s = systemfunc("ls -1 " + ModelOutRoot+case2s+"/"+case2s+h0Suf)


    File1s    = addfiles(FileName1s,"r")
    File2s    = addfiles(FileName2s,"r")
    ;print(File1s)
    ;print(File2s)
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


     varName   = varNameList(iv)
     tmp1s     = File1s[:]->RKZ_ql_incld
     tmp2s     = File2s[:]->RKZ_ql_incld

     tmpdf     = log10(tmp2s)
     tmpdf1d   = ndtooned(tmpdf)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     opt         = True
     opt@bin_min = minbin
     opt@bin_max = maxbin
     tpdforg = pdfx(tmpdf1d, nbins, opt) 
     ;;;convert the pdf to grid point numbers;;;;;;;;;
     ;;;note that pdf function gives you percentage;;;
     ;;;the ratio=percentage/100;;;;;;;;;;;;;;;;;;;;;;
     tpdforg = tpdforg*dimsizes(tmpdf1d)/100.0
;;;;;;;check the results;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     print("Total grid cells "+dimsizes(tmpdf1d)+ " "+num(tmp2s.ge.0)+ " "+dim_sum_n(tpdforg,0))
     print("Total grid cells "+dimsizes(lat)*dimsizes(pres))
;;;;;;;;;refine the pdf for the plots;;;;;;;;
     tpdforg@_FillValue=32766
     tpdforg = where(tpdforg.gt.0,tpdforg,tpdforg@_FillValue)
     print(min(tpdforg)+" "+max(tpdforg)) 
;;;;;;;parpare the plot arrays;;;;;;;;;;;;;;;
    nBin    = nbins          
     xx(it,:) = tpdforg@bin_center
     yy(it,:) = (/ tpdforg /)
     ratio(it)= dim_sum_n(tpdforg,0)/dimsizes(tmpdf1d) *100.0
    delete([/tmpdf,tmpdf1d,tpdforg/])

    end do ;;end of group loop

;;;;;define plot sorces;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    fignam="fig_pdf_lmt4_1h_"+Groups(igp)+"_"+Regions(irg)+"_"+varNameList(iv)+rid
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
    res@tmYLPrecision        = 2        ; set the precision to 2 significant digits
    res@tmYLFormat           = "#0^se"   ; don't use "e" symbol, always show +/-

    res@trYLog               = True
    res@trXLog               = False
    res@trYMaxF              = 5e4
    res@trYMinF              = 1e-1
    res@trXMaxF              =  2
    res@trXMinF              = -14

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
    res@xyMonoDashPattern        = True              ; all solid 
    res@tiYAxisString            = "Number of grid points (#)"
    res@tiXAxisString            = "log~B~10~N~(ql_incld)"
    res@gsnRightString           = "";"A+B1+C2" ;Expnam(igp)  ; "Total # of grid points = "+ngridpts
   ;res@gsnLeftString            = "Histogram of ql_incld (unit: kg/kg)"
    res@gsnCenterString            = "Histogram of ql_incld (unit: kg/kg)"


    res@gsnXYBarChart            = True              ; Create bar plot
    res@gsnXYBarChartOutlineOnly = True

    res@pmLegendDisplayMode    = "Always"            ; turn on legend
    res@pmLegendSide           = "Top"               ; Change location of 
    res@pmLegendParallelPosF   = .55                 ; move units right
    res@pmLegendOrthogonalPosF = -0.8               ; move units down
    res@pmLegendWidthF         = 0.1               ; Change width and
    res@pmLegendHeightF        = 0.18               ; height of legend.
    res@lgPerimOn              = False                ; turn off/on box around
    res@lgLabelFontHeightF     = .018                ; label font height
 
;;;;;;;generate labels for the plot;;;;;;;;;
    do ilab = 0, nlabs-1,1
     if(ilab.le.1)then
      rstr = sprintf("%5.0f",ratio(ilab))
      labplt(ilab)= labels(ilab) ;+ " (" +str_squeeze(rstr)+ "%)"
     else
      rstr = sprintf("%8.4f",ratio(ilab))
      labplt(ilab)= labels(ilab) ;+ " (" +str_squeeze(rstr)+ "%)"
     end if

    end do
   res@xyExplicitLegendLabels = labplt
;   res@tiMainString           = Expnam(igp)
    plot = gsn_csm_xy (wks, xx, yy, res)
    
;---Draw both plots in a panel
    pnres                    = True
    pnres@gsnMaximize        = True
    pnres@gsnPanelMainString =" " ;"A+B1+C2" ;+dtimeList(it) + "s - "+dtimeList(0)+ "s"  ;Expnam(igp) ;"Time integration: 1-h"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@gsnPanelMainFontHeightF = 0.028
;   pnres@gsnPanelBottom     = 0.05
    gsn_panel(wks,(/plot/),(/1,1/),pnres)
   
  end do ; ensemble loop

end do ; variable loop 

end do ; end of regions loop
end do ; end of times loop

end
