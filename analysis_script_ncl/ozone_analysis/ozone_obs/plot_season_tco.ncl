  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;Following  NASA Ozone Watch website and Tang (2020), two metrics are derived
;(i) the area (in millions of square kilometers)  with less than 220 DU in total column ozone
;(ii) the minimum total column ozone (in DU).
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/data/CVDP_RGD"
  data_di0   = "../data/ozone_vint"
  exp_name   = (/"C3S-L4_OZONE", "v2.LR.historical","v1.LR.historical"/)
  labels     = (/"OBS.",         "v2 Coupled"      , "v1 Coupled"/)
  nexps      = dimsizes(exp_name)
  ltrend     = new(nexps,string)

  regions    = (/"Antarctic"/)
  nregs      = dimsizes(regions)
  latS       = -90.0
  latN       = -40.0

 ;define regions for eof analysis
  varnam     = (/"TCO"/)
  varnamtr   = (/"Total column Ozone"/)
  facmod     = (/1.0/)
  varunt     = (/"DU"/)
  thresh     = 220.0 ;unit: Du 

  ystr       = 1979
  yend       = 2014
  time_tag   = ystr+"01-"+yend+"12" 
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  Min1       = (/120,240,250,220,250,150/)
  Max1       = (/360,320,300,320,300,360/)

  Min2       = (/0,   0,  0,  0,  0, 0/)
  Max2       = (/20, 20, 20, 20, 20, 20/)

  Min3       = (/260,280,280,280,270,240/)
  Max3       = (/400,350,360,370,320,400/)

  do k = 0,nexps-1

    if (k.eq.0) then 
       fl   = systemfunc("ls " + data_dir+"/"+exp_name(k)+"/"+exp_name(k)+"*.nc") 
       f    = addfile(fl,"r")
       lat  = f->lat({latS:latN})
       lon  = f->lon
       date = cd_calendar(f->time,0)
       indy = ind (date(:,0).ge.ystr.and.date(:,0).le.yend)
       time = yyyymm_to_yyyyfrac(cd_calendar(f->time(indy),-1),0.0)
       var  = f->TCO(indy,{latS:latN},:)  ; total column ozone
       delete([/fl,f,date,indy/])

       ozi      = new((/5,dimsizes(time),nexps/),double)
       ozi!0    = "ens"
       ozi!1    = "time" 
       ozi!2    = "exp"
       ozi&time = time
       oza      = ozi
       ozm      = ozi

       ;;calcualte the area
       jlat = dimsizes(lat)
       rad  = 4.0*atan(1.0)/180.0
       cost = lat
       cost = cos(lat*rad)       ; cosine weights
       re   = 6371220.0
       rr   = re*rad
       dlon = abs(lon(2)-lon(1))*rr
       dx   = dlon*cost
       dy   = new (jlat, typeof(dx))
       dy(0)         = abs(lat(2)-lat(1))*rr
       dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5
       dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr
       area          = dx*dy ; cell area function of latitude only
       area          = area / 1e6 ; convert to km^2
       area          = area * 1e-6
       area@units    = "10~S~6~N~ km~S~2~N~"
       delete([/jlat,rad,re,rr,dlon,dx,dy/])
       lmask = where(var.lt.thresh,1,0)

       nens = 1
       area1 = conform(var,area,(/1/))
       ozi(0:nens-1,:,k) = (/dim_min_n_Wrap(dim_min_n_Wrap(var,2),1)/)
       oza(0:nens-1,:,k) = (/dim_sum_n_Wrap(dim_sum_n_Wrap(area1*lmask,2),1)/)
       ozm(0:nens-1,:,k) = (/dim_avg_n_Wrap(dim_avg_wgt_n_Wrap(var(:,{-90:-60},:),cost({-90:-60}),1,1),1)/)
       delete(area1)
     else
      ;fl0  = systemfunc("ls " + data_dir+"/"+exp_name(k)+"/"+exp_name(k)+".en*."+varnam+"."+time_tag+".nc")  
      ;fl1  = fl0
      ; do ii = 0,dimsizes(fl0)-1,1
      ;   fl1(ii) = str_sub_str(fl0(ii),"TCO","SCO")
      ; end do 
      ; f0   = addfiles(fl0,"r") 
      ; f1   = addfiles(fl1,"r")
      ; ListSetType (f0, "join")
      ; ListSetType (f1, "join")
       time_tag = ystr+"-"+yend
       fl0  = systemfunc("ls " + data_di0+"/"+exp_name(k)+"_en*_"+varnam+"_MON_"+time_tag+".nc")
       f0   = addfiles(fl0,"r")
       ListSetType (f0, "join")

       lat  = f0[0]->lat({latS:latN})
       lon  = f0[0]->lon
       time = yyyymm_to_yyyyfrac(cd_calendar(f0[0]->time,-1),0.0)

      ;Total column ozone
       var   = f0[:]->$varnam$(:,:,{latS:latN},:)
       var!0 = "ens"
       printMinMax(var,False)
       printVarSummary(var)
       
       ;;calcualte the area
       jlat = dimsizes(lat)
       rad  = 4.0*atan(1.0)/180.0
       cost = lat
       cost = cos(lat*rad)       ; cosine weights
       re   = 6371220.0
       rr   = re*rad
       dlon = abs(lon(2)-lon(1))*rr
       dx   = dlon*cost
       dy   = new (jlat, typeof(dx))
       dy(0)         = abs(lat(2)-lat(1))*rr
       dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5
       dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr
       area          = dx*dy ; cell area function of latitude only
       area          = area / 1e6 ; convert to km^2
       area          = area * 1e-6
       area@units    = "10~S~6~N~ km~S~2~N~"
       delete([/jlat,rad,re,rr,dlon,dx,dy/])

       lmask = where(var.lt.thresh,1,0)
       area1 = conform(var,area,(/2/))
       nens  = dimsizes(var(:,0,0,0))
       printVarSummary(var)
       ozi(0:nens-1,:,k) = (/dim_min_n_Wrap(dim_min_n_Wrap(var,3),2)/)
       oza(0:nens-1,:,k) = (/dim_sum_n_Wrap(dim_sum_n_Wrap(area1*lmask,3),2)/)
       ozm(0:nens-1,:,k) = (/dim_avg_n_Wrap(dim_avg_wgt_n_Wrap(var(:,:,{-90:-60},:),cost({-90:-60}),1,2),2)/)
       delete(area1)
       delete([/fl0,f0/])
     end if 

     rc    = new((/3,nens/),double)
     do ie = 0,nens-1,1
       rc(0,ie) = regline(time(:),ozi(ie,:,k)) * 10.0
       rc(1,ie) = regline(time(:),oza(ie,:,k)) * 10.0
       rc(2,ie) = regline(time(:),ozm(ie,:,k)) * 10.0
     end do

     print(avg(rc)+ " " +stddev(rc))
     ltrend(k) = sprintf("%.2f",avg(rc)) + " ~F34~1~N~~F~ " + sprintf("%.3f",stddev(rc))
     delete([/var,rc,cost,lmask/])
   
   end do 

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;

   do i = 0,nseas-1,1

   wtype           = "png"
   wtype@wkWidth   = 1200
   wtype@wkHeight  = 1200
   fignam          = regions+"_"+varnam+"_"+seasons(i)+"_"+time_tag
   wks             = gsn_open_wks(wtype,fignam)
   gsn_define_colormap(wks,"t2m_29lev")
   colors = (/1,27,10,19,31/)
   lthick = (/3.0,3.0,3.0,5.0,5.0/)
   ldash  = (/0,0,0,0,0/)
   opcity = (/0.2,0.2,0.2,0.2,0.2/)
   dum1   = new(nexps,graphic)
   dum2   = new(nexps,graphic)
   dum3   = new(nexps,graphic)
   plot   = new(3,graphic)

   if (seasons(i).eq."MON") then 
     pltdat = new((/3,nexps,6,dimsizes(time)/),double)
     time_new = yyyymm_to_yyyyfrac(time,0.0)
   else 
     pltdat = new((/3,nexps,6,dimsizes(time)/12/),double)
     time_new = ispan(ystr,yend,1)*1.0 
   end if 

   do ff = 0,2
    if(ff.eq.0) then
      tmp = ozi(exp|:,ens|:,time|:)
    else if (ff.eq.0) then
      tmp = oza(exp|:,ens|:,time|:)
    else 
      tmp = ozm(exp|:,ens|:,time|:)
    end if 
    end if 
    tmp_sea        = runave_n_Wrap(tmp, 3, 0, 2)
    tmp_ann        = runave_n_Wrap(tmp, 12, 0, 2)
    tmp_sea(:,:,0) = dim_avg_n(tmp(:,:,0:1),2)
    printVarSummary(pltdat)
    printVarSummary(tmp_ann)
    if(seasons(i).eq."ANN") then
      pltdat(ff,:,1:,:) = tmp_ann(:,:,5::12)
    else if(seasons(i).eq."DJF") then
      pltdat(ff,:,1:,:) = tmp_sea(:,:,0::12)
    else if(seasons(i).eq."MAM") then
      pltdat(ff,:,1:,:) = tmp_sea(:,:,3::12)
    else if(seasons(i).eq."JJA") then
      pltdat(ff,:,1:,:) = tmp_sea(:,:,6::12)
    else if(seasons(i).eq."SON") then
      pltdat(ff,:,1:,:) = tmp_sea(:,:,9::12)
    else 
      pltdat(ff,:,1:,:) = tmp
    end if
    end if
    end if
    end if 
    end if 
    pltdat(ff,:,0,:)  = dim_avg_n(pltdat(ff,:,1:,:),1)
    delete([/tmp,tmp_sea,tmp_ann/])
   end do 

   pltdat1 = pltdat(0,:,:,:)
   pltdat2 = pltdat(1,:,:,:)
   pltdat3 = pltdat(2,:,:,:)

   delete([/pltdat/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then
     wtype           = "png"
     wtype@wkWidth   = 1200
     wtype@wkHeight  = 1200
     fignam          = regions+"_"+varnam+"_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/1,27,10,19,31/)
     lthick = (/3.0,3.0,3.0,5.0,5.0/)
     ldash  = (/0,0,0,0,0/)
     opcity = (/0.2,0.2,0.2,0.2,0.2/)
     dum1   = new(nexps,graphic)
     dum2   = new(nexps,graphic)
     dum3   = new(nexps,graphic)
     plot   = new(3,graphic)
   end if


   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True         
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling
 
    FontHeightF = 0.030
    res@tiMainFontHeightF     = FontHeightF*1.5
    res@tmYLLabelFontHeightF  = FontHeightF*1.5
    res@tmXBLabelFontHeightF  = FontHeightF*1.5
    res@gsnStringFontHeightF  = FontHeightF*1.5

    if(wtype.eq."png") then
     bthick = 3.0
    else
     bthick = 1.0
    end if

    res@tmBorderThicknessF    = bthick
    res@tmXBMajorThicknessF   = bthick
    res@tmXBMinorThicknessF   = bthick
    res@tmYLMajorThicknessF   = bthick
    res@tmYLMinorThicknessF   = bthick
    res@tmYRMajorThicknessF   = bthick
    res@tmYRMinorThicknessF   = bthick

    res@vpHeightF             = 0.80        ; Changes the aspect ratio
    res@vpWidthF              = 0.90
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@pmLegendDisplayMode    = "Never"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of 
    res@pmLegendParallelPosF   = 0.35                   ; move units right
    res@pmLegendOrthogonalPosF = -0.45                  ; more neg = down
    res@pmLegendWidthF         = 0.10                  ; Change width and
    res@pmLegendHeightF        = 0.12                  ; height of legend.
    res@lgLabelFontHeightF     = .024                   ; change font height
    res@lgPerimOn              = False                 ; no box around

    res@tmXTOn                 = False       ; have tick marks

    res@gsnLeftString           = ""
    res@gsnRightString          = ""
    res@tiMainString            = ""
    res@tiYAxisString           = "" 
    res@tiXAxisString           = ""

    res@gsnYRefLine             = (/-9999./)              ; reference line
    res@gsnYRefLineColors       = (/"Grey25"/)
    res@gsnYRefLineThicknesses  = (/0.5/)
    res@gsnYRefLineDashPatterns = (/0/)
    res@xyMonoLineColor         = False             ; want colored lines

    res@trXMinF                 = 1978 ;floor(min(tsm03))
    res@trXMaxF                 = 2015 ;ceil(max(tsm03))

    do k = 0,nexps-1,1

    ;*******************************************
    ; First plot
    ;*******************************************

    res1  = res
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = new(nens+1,integer)
    res1@xyLineThicknesses       = new(nens+1,float)
    res1@xyDashPatterns          = new(nens+1,integer)
    res1@xyLineOpacities         = new(nens+1,float)
    res1@xyLineColors(:)         = colors(k)
    if(seasons(i).eq."MON") then 
      res1@xyLineThicknesses(:)    = lthick(k)*0.1      ; line thicknesses
    else
      res1@xyLineThicknesses(:)    = lthick(k)*1.5
    end if 
    res1@xyDashPatterns(:)       = ldash(k)      ; make all lines solid
    res1@xyLineOpacities(:)      = 1.0
    res1@xyLineThicknesses(1:)   = lthick(k)/2.0
    res1@xyLineOpacities(1:)     = 0.5

   ;*******************************************
   ; First plot
   ;*******************************************
    if(k.eq.0) then
      res1@trYMinF                =  min( pltdat1) - abs(max( pltdat1))/20.0
      res1@trYMaxF                =  max( pltdat1) + abs(max( pltdat1))/20.0
    end if
    if(k.eq.0) then
      res1@tiYAxisString          = "Min. total column Ozone ("+ varunt+")"      ; y-axis label
      res1@tiXAxisString          = "Year"
      res1@gsnLeftString          = "Min. O3 in ozonal hole over Antarctic"
      res1@gsnRightString         = seasons(i)
      res1@tiMainString           = ""
    else
      res1@gsnLeftString          = ""
      res1@gsnRightString         = ""
      res1@tiMainString           = ""
      res1@tiYAxisString          = ""
      res1@tiXAxisString          = ""
    end if

    dum1(k)                      = gsn_csm_xy (wks, time_new, pltdat1(k,:,:), res1)

    ;res2 = res
    ;res2@gsnLeftString          = ""
    ;res2@gsnRightString         = ""
    ;res2@tiMainString           = ""
    ;res2@tiYAxisString          = ""
    ;res2@tiXAxisString          = ""
    ;res2@gsnXYFillColors        = colors(k)
    ;res2@gsnXYFillOpacities     = opcity(k)
    ;res2@xyLineColor            = -1          
    ;dum2(k) = gsn_csm_xy (wks, time_new, pltdat(k,1:,:), res2)
    ;overlay(dum1(k),dum2(k))


   ;*******************************************
   ; Second plot
   ;*******************************************
    res1  = res
    if(k.eq.0) then
      res1@trYMinF                =  min( pltdat2) - abs(max( pltdat2))/20.0
      res1@trYMaxF                =  max( pltdat2) + abs(max( pltdat2))/20.0
    end if
    if(k.eq.0) then
      res1@tiYAxisString          = "Area (km~S~2~N~)"      ; y-axis label
      res1@tiXAxisString          = "Year"
      res1@gsnLeftString          = "Ozone hole (< "+thresh+" DU) area over Antarctic"
      res1@gsnRightString         = seasons(i)
      res1@tiMainString           = ""
    else
      res1@gsnLeftString          = ""
      res1@gsnRightString         = ""
      res1@tiMainString           = ""
      res1@tiYAxisString          = ""
      res1@tiXAxisString          = ""
    end if

    dum2(k)                      = gsn_csm_xy (wks, time_new, pltdat2(k,:,:), res1)


   ;*******************************************
   ; Third plot
   ;*******************************************
    res1  = res
    if(k.eq.0) then 
      res1@trYMinF                 =  min( pltdat3) - abs(max( pltdat3))/20.0
      res1@trYMaxF                 =  max( pltdat3) + abs(max( pltdat3))/20.0
    end if 
    if(k.eq.0) then
      res1@tiYAxisString          = "TCO ("+ varunt+")"      ; y-axis label
      res1@tiXAxisString          = "Year"
      res1@gsnLeftString          = "Total Column Ozone"
      res1@gsnRightString         = seasons(i)
      res1@tiMainString           = ""
    else
      res1@gsnLeftString          = ""
      res1@gsnRightString         = ""
      res1@tiMainString           = ""
      res1@tiYAxisString          = ""
      res1@tiXAxisString          = ""
    end if

    dum3(k)                      = gsn_csm_xy (wks, time_new, pltdat3(k,:,:), res1)


    if(k.eq.0) then
      plot(0) = dum3(k)
      plot(1) = dum1(k)
      plot(2) = dum2(k)
    else
      overlay(plot(0),dum3(k))
      overlay(plot(1),dum1(k))
      overlay(plot(2),dum2(k))
    end if 
 
   delete([/res1/])

   end do 

   gres = True
   gres@YPosPercent          = 95.     
   gres@XPosPercent          = 65.    
   gres@Position             = "Top"        
   gres@ItemSpacePercent     = 8.          
   gres@LineLabelWhiteSpacePercent = 2 

   lineres = True
   lineres@lgLineColors       = colors(0:nexps-1)
  ;lineres@lgLineOpacity      = opcity(0:nexps-1)
   lineres@lgLineThicknesses  = lthick(0:nexps-1)
   lineres@lgDashIndexes      = ldash(0:nexps-1)
   lineres@LineLengthPercent  = 5      

   textres                    = True
   textres@lgLabels           = labels(0:nexps-1)
   textres@lgLabelFontHeights = FontHeightF*1.2
   textres@lgLabelColors      = colors(0:nexps-1)
  ;textres@lgLabelOpacity     = opcity(0:nexps-1)

   do ii = 0,2
    plot(ii) = simple_legend(wks,plot(ii),gres,lineres,textres)
   end do 

   ;;add text;;
   txid = new(3*(nexps+1),graphic)                  ; Nine text strings
   txres                  = True             
  ;txres@txFont           = "helvetica-bold" 
   txres@txFontHeightF    = FontHeightF*2.0
   txres@txFontThicknessF = 4.0
   txres@txJust           = "BottomLeft"

   do ii = 0,2
    do m = 0,nexps-1,1
      txres@txFontColor   = colors(m)
      txtstr = labels(m) + " ("+ltrend(m) +")"
      txid(m) = gsn_add_text (wks,plot(ii),txtstr, 1982.0, -15-4*m ,txres)
    end do 
   end do 
  
   ;panel plot only resources
   pres1                 = True                ; mods desired
   pres1@gsnFrame        = True               ; save panel until both ready
   gsn_panel(wks,plot(0),(/3,1/),pres1)          ; create first panel

   delete([/time_new,pltdat1,pltdat2,pltdat3,res,dum1,dum2,dum3,plot,txid,txres,pres1,lineres,gres/])

 end do 

end

