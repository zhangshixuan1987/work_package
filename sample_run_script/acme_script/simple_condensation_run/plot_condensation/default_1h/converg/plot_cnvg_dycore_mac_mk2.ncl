load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

 plot_logdt_logrms = False
 plot_logdt_logrms = True


;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Global"/) ;,"Land","Ocean"/)
  nregs        = dimsizes(Regions)

  varNameList  = (/"T"/)
  nvars        = dimsizes(varNameList)

  hash         = "E3SM_CONVG"
  compset      = "FC5AQUAP"
  config       = "_dycore_mac"
  config1      = "_dycore_mac"
  Group        = "RKZ_A1_B3_C2_ql17_lmt4_fmin1_ne30_ne30_"+compset+"_intel_constance"+"_"

  leftstring = "Resolution: ne30 (1 degree) with 30 levels~C~"+\
               "compset: "+compset
  plotFileLabel = "_FC5AQUAP_intel_constance"

 ;---------------------
 ; time step lengths
 ;---------------------
 dtimeList = (/2,8,30,120,450,1800/)
 x         = log10(dtimeList*1.) 
 nx        = dimsizes(x)

;===============================================================

  ny = 5   ; # of configurations plotted by this script.
  nr = 1
  lglabel = new( ny, "string" )
  ldoregr = new( ny, "logical" )
  rawdata = new( (/ny,nr,nx/), "float" )
  colors  = (/"Black","Red", "Blue", "Green3", "Magenta"/) ;(/1, 8 , 26, 15/)
  markers = (/16,16,16,16,16/)
  mkszf   = (/0.8,0.8,0.8,0.8,0.8/)

do iv  = 0, nvars-1

 do irg = 0, nregs-1

  filnm1 = "rmse_reft_1s_03600_"+Regions(irg)+"_"+Group+varNameList(iv)+config+".txt"
  filnm2 = "rmse_reft_1s_07200_"+Regions(irg)+"_"+Group+varNameList(iv)+config1+".txt"
  filnm3 = "rmse_reft_1s_14400_"+Regions(irg)+"_"+Group+varNameList(iv)+config1+".txt"
  filnm4 = "rmse_reft_1s_21600_"+Regions(irg)+"_"+Group+varNameList(iv)+config1+".txt"
  filnm5 = "rmse_reft_1s_43200_"+Regions(irg)+"_"+Group+varNameList(iv)+config1+".txt"

  iy = 0 

  lglabel(iy) = "  1-h model run"

  ldoregr(iy) = True
  rawdata(iy,:,:) = readAsciiTable(filnm1, nx, "float", 0)

  iy = iy + 1 

  lglabel(iy) = "  2-h model run"
  ldoregr(iy) = True
  rawdata(iy,:,:) = readAsciiTable(filnm2, nx, "float", 0)

  iy = iy + 1

  lglabel(iy) = "  4-h model run"
  ldoregr(iy) = True
  rawdata(iy,:,:) = readAsciiTable(filnm3, nx, "float", 0)

  iy = iy + 1

  lglabel(iy) = "  6-h model run"
  ldoregr(iy) = True
  rawdata(iy,:,:) = readAsciiTable(filnm4, nx, "float", 0)

  iy = iy + 1

  lglabel(iy) = "  12-h model run"
  ldoregr(iy) = True
  rawdata(iy,:,:) = readAsciiTable(filnm5, nx, "float", 0)

;===============================================================
 rawdata = log10(rawdata)
 printVarSummary(rawdata)
 davg = dim_avg_n( rawdata, 1 )
 dstd = dim_stddev_n( rawdata, 1 )

 npanel = 1
 rmsd = new( (/3,ny,nx/),"float" )
 rmsd(0,:,:) = davg
 rmsd(1,:,:) = davg - dstd
 rmsd(2,:,:) = davg + dstd

 rmsd@_FillValue = -999.

 ;---------------------
 ; for plotting
 ;---------------------
  dum    = new((/ny,nx/),"graphic")
  dumtx  = new( ny,"graphic")

  if (plot_logdt_logrms) then
     plotFileSuf = ""
  else
     plotFileSuf = "_lin"
  end if

  wks  = gsn_open_wks("png","fig_cnvg_"+varNameList(iv)+"_"+Regions(irg)+plotFileLabel+plotFileSuf) ; specifies a ps plot
  gsn_define_colormap(wks,"StepSeq25")
 
  res          = True                   ; plot mods desired
 ;res@vpHeightF = 0.4
 ;res@vpWidthF =  0.6 
  res@gsnDraw  = False
  res@gsnFrame = False

  FontHeight = 0.022
  res@gsnStringFontHeightF = FontHeight*0.8
  res@tmXBLabelFontHeightF = FontHeight
  res@tmYLLabelFontHeightF = FontHeight
  res@tmYRLabelFontHeightF = FontHeight
  res@tiMainFontHeightF    = FontHeight
  res@tiXAxisFontHeightF   = FontHeight
  res@tiYAxisFontHeightF   = FontHeight

  res@xyMarkLineModes     = new( ny*2, "string" )
  res@xyMarkLineModes(:ny-1) =  "Lines"
  res@xyMarkLineModes(ny:)   =  "Markers"

  res@xyMarkers           = new( ny*2, "integer" )
  res@xyMarkers(:ny-1)    = markers
  res@xyMarkers(ny:)      = markers
  res@xyMarkerColors      = new( ny*2, "string" )
  res@xyMarkerColors(:ny-1) = colors 
  res@xyMarkerColors(ny:) = colors 
  res@xyMarkerSizes       = new( ny*2, "float" )
  res@xyMarkerSizes(:ny-1) = 0.015*mkszf
  res@xyMarkerSizes(ny:) = 0.015*mkszf
  res@xyDashPattern       = 5 
  res@xyLineThicknessF    = 5
  res@xyLineColors        = colors 

  res@tmXTOn = False

     res@trXMaxF = 3.5
     res@trXMinF = 0
     res@tmXBMode = "Manual"
     res@tmXBTickSpacingF = 1
     res@tmXBFormat = "f"

     res@tmYLLabelJust = "CenterRight"
     res@tmYRLabelJust = "CenterRight"
     res@tmYUseLeft = True
    ;res@tmYRLabelsOn = True

    ;res@tmYRMode = "Manual"
    ;res@tmYRTickSpacingF = 1
    ;res@tmYRFormat = "f"
    ;res@tmYRLabelJust = "CenterRight"
     xx = log10(dtimeList*1.)

     res@tiXAxisString = "log~B~10~N~(~F8~D~F10~t~F~)"

   if (plot_logdt_logrms) then

     res@tiYAxisString = "log~B~10~N~(RMSE)"
    ;res@trYMaxF = 0.
    ;res@trYMinF = -6.
    ;res@trYMinF = -1.5  ; for v1
    ;res@trYMinF = -2.5

    res@trYMaxF = -3.0 ;;for Q or CLW or CIE
    res@trYMinF = -10.0
    
    if(varNameList(iv).eq."U".or.varNameList(iv).eq."V")then
    res@trYMaxF = 0.0
    res@trYMinF = -6.0
    end if 
    
    if(varNameList(iv).eq."T")then
    res@trYMaxF = 0.0
    res@trYMinF = -6.0
    end if

     res@tmYLMode = "Automatic"
    ;res@tmYLMode = "Manual"
    ;res@tmYLTickSpacingF = 1
    ;res@tmYLFormat = "f"

   else
     res@tiYAxisString = "RMSE"
     res@trYMaxF =  0.16
     res@trYMinF = -0.04 
     res@tmYLTickSpacingF = 0.04
   end if

 ;************************************************
 ; do regression and plot results
 ;************************************************
 ; data is the tmp array for plotting 
 ; holding the original rmsd and 
 ; the fitted linear function

 data = new( (/2,ny,nx/),typeof(rmsd(0,:,:)) )
 data@_FillValue = -999.
 data = -999.

 slope = new( ny,typeof(rmsd(0,:,:)) )
 slope@_FillValue = -999.
 slope = -999.

 ; regression
 do iy = 0,ny-1

    if (ldoregr(iy)) then

       y  = rmsd(0,iy,:)
       rc = regline(x,y)

       data(1,iy,:) = y 
       data(0,iy,:) = rc*(x-rc@xave) + rc@yave

       slope(iy) = rc

    else
       data(1,iy,:) = rmsd(0,iy,:) 
    end if

 end do

 ; plotting 

 ;datatmp = reshape(data,(/2*ny,nx/)  )
  res@gsnLeftString = leftstring 
  datatmp = onedtond(ndtooned(data),(/2*ny,nx/)  )

  if (plot_logdt_logrms) then
     plot = gsn_csm_xy (wks,xx,datatmp,res)        ; create plot
  else
     plot = gsn_csm_xy (wks,xx,10^datatmp,res)        ; create plot
  end if

 ; add uncertainty bars

  plres = True
  plres@gsLineThicknessF = 3
  do iy = 0,ny-1

     plres@gsLineColor = colors(iy)

     do ix = 0,nx-1
        barxx = (/xx(ix),xx(ix)/)

        if (plot_logdt_logrms) then
        baryy = (/rmsd(1,iy,ix),rmsd(2,iy,ix)/)
        else
        baryy = 10^(/rmsd(1,iy,ix),rmsd(2,iy,ix)/)
        end if
        dum(iy,ix) =  gsn_add_polyline( wks, plot, barxx, baryy, plres )
     end do
     
  end do

;;manual legend

  txx = new( ny,"float" ) 

  if (plot_logdt_logrms) then
    if(varNameList(iv).eq."Q")then
     txx = (/ 3.4, 3.4, 3.4/)
     txy = (/-9.3,-9.0,-8.7/) 
    ;txy = (/-0.5/) 
    else
     txx = (/ 3.0, 3.0, 2.0, 2.0, 2.0/)
     txy = (/-4.1, -3.8, -1.5, -1.2, -0.9/)
    ;txy = (/-0.5/)
    end if 
  else
     txx = (/ 1.6/) + 0.22
     txy = (/0.14/) 
     txy = txy(::-1) 
  end if

  txres = True
  txres@txFontHeightF = FontHeight*0.93
  txres@txJust = "CenterRight"
  do iy=0,ny-1
     txres@txFontColor = res@xyMarkerColors(iy) 

     txt = "~F22~"+lglabel(iy)
     if (ldoregr(iy)) then
        txt = txt +" ("+sprintf("%3.1f",slope(iy))+")~F~"
     end if

     dumtx(iy) = gsn_add_text(wks,plot, txt, txx(iy),txy(iy),txres)
  end do


  draw(plot)
  frame(wks)

 end do ; end of regions loop

end do ; end of variables loop

end

