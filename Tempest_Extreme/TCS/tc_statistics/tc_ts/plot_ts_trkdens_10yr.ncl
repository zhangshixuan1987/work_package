;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;This script used to plot e3sm diagnostics for TCs;;;;;;
;TC Basin definitions.(Knutson et. al.,2020,BAMS)
;North Atlantic: 0-90N, 265-360E
;Northeast Pacific: 0-90N, 180-265E
;Northwest Pacific: 0-90N, 100-180E
;North Indian: 0-90N, 30-100E
;South Indian: 90S-0, 20-135E
;Southwest Pacific: 90S-0, 135-295E
;South Atlantic(SA1+SA2): 90S-0, South America (295-360E) to Africa (0-20E)
;The North Atlantic-northeast Pacific boundary is on a
;diagonal tracing a path through Mexico and Central America.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ;TC Basin definitions.(Knutson et. al.,2020,BAMS)
  l_plot_latlon_box = False ; if True, add lat-lon box for region
  Basin   = (/"GLB","NA","NEP","NWP","NI","SI","SWP","SA"/)
  Blats   = (/  -90,   0,    0,    0,   0, -90,  -90,  -90,  -90/)
  Blatn   = (/   90,  90,   90,   90,  90,   0,    0,    0,    0/)
  Blonw   = (/    0, 265,  180,  100,  30,  20,  135,  295,    0/)
  Blone   = (/  360, 360,  265,  180, 100, 135,  295,  360,   20/)
  nbsin   = dimsizes(Basin)

  datdir  = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS/track_info"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList = (/"OBS",       "ERA5",     "E3SM",     "NDGUV",      "NDGUVTQ"/)
  Groups  = (/"IBTrACS",   "ERA5",     "CLIM",     "NDGUV",      "NDGUVTQ"/)
  timtag  = (/"1979-2017", "1979-2017", "2007-2017","2007-2017","2007-2017"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  figkey  = "trkdens_10year"

  Figstr  = (/"(d)", "(e)", "(f)"/)
  ystr    = 2008
  yend    = 2017
  nyr     = yend - ystr +1

  vtrk   = new ((/nyr,ngrps,nbsin/),double)
  Varnam = "trkdens"

  do j = 0,nbsin-1,1

   do i = 0,ngrps-1,1

    fna  = datdir+"/"+Groups(i)+"_"+sset+"_Track_statistics_"+timtag(i)+".nc"
    fla  = addfile(fna, "r")

    lat   = fla->lat
    lon   = fla->lon
    time  = fla->time
    var   = fla->$Varnam$(0:nyr-1,:,:)
    do k = 0,nyr-1,1
     mon1  = (/5,6,7,8,9,10,11/)
     yy1   = mon1
     yy1(:)= ystr + k
     tim1  = toint(sprinti("%04d",yy1) + sprinti("%02d",mon1))
     indx1 = get1Dindex(time,tim1)
     mon2  = (/1,2,3,4,10,11,12/)
     yy2   = mon2
     yy2(:)= ystr + k
     tim2  = toint(sprinti("%04d",yy2) + sprinti("%02d",mon2))
     indx2 = get1Dindex(time,tim2)
     var(k,{5:90},:)  = dim_avg_n_Wrap(fla->$Varnam$(indx1,{5:90},:),0) 
     var(k,{-90:5},:) = dim_avg_n_Wrap(fla->$Varnam$(indx2,{-90:5},:),0) 
     delete([/yy1,mon1,tim1,indx1,yy2,mon2,tim2,indx2/])
    end do 

    if(Basin(j).ne."SA") then 
      lmsk = var
      lmsk = 0
      lmsk(:,{Blats(j):Blatn(j)},{Blonw(j):Blone(j)}) = 1
    else 
      lmsk = var
      lmsk = 0
      lmsk(:,{Blats(j):Blatn(j)},{Blonw(j):Blone(j)}) = 1
      lmsk(:,{Blats(j+1):Blatn(j+1)},{Blonw(j+1):Blone(j+1)}) = 1
    end if 
    
    var = where(lmsk.eq.0,var@_FillValue,var)

    vtrk(:,i,j) = dim_sum_n_Wrap(dim_sum_n_Wrap(var,2),1)

    delete([/var,fna,fla,lat,lon,time/])

   end do 

  end do 

  do ii = 0,nbsin-1,1

   ;print(vtrk(:,0,ii) + " " + vtrk(:,1,ii)+ " " + vtrk(:,2,ii) + " " +vtrk(:,3,ii))

    xx = ispan(ystr,yend,1)
    yy = new((/ngrps,dimsizes(xx)/),double)
    do jj = 0,ngrps-1,1
      yy(jj,:) = vtrk(:,jj,ii)
    end do 
   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "track_dens_ts_" + Basin(ii) 
    wks      = gsn_open_wks(wtype,plotname)

   ;---Resources to share between both plots
   res                 = True     ; Plot modes desired.
   res@vpHeightF       = 0.4      ; change aspect ratio of plot
   res@vpWidthF        = 0.6
  ;res@gsnDraw         = False
  ;res@gsnFrame        = False
  ;res@gsnMaximize     = False     ; Maximize plot
   res@tmXTOn          = False
   res@tmYROn          = True
   res@tmYRBorderOn    = True         ; Habilita a borda do eixo y direito (YR). 
   res@tmXTBorderOn    = True         ; Habilita a borda do eixo x superior (XB). 

   res@gsnLeftString   = "Track Density"
   res@gsnRightString  = ""
   res@tiMainString    = ""
   res@tiYAxisString   = "storm trainsts/month"
   res@tiXAxisString   = "Year"
    
   BorderThick = 1.5
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeightF = 0.022
   res@tiMainFontHeightF     = FontHeightF
   res@tmYLLabelFontHeightF  = FontHeightF
   res@tmXBLabelFontHeightF  = FontHeightF

  ;res@tmXBLabels            = modes 
   res@xyDashPattern         = 0                         ; Make curves all solid
   res@xyLineColors          = (/"Grey25","black","blue","red","green"/)
   res@xyMarkLineMode        = "MarkLines"               ; Markers *and* lines
   res@xyMarkers             = (/16,16,16,16,16/)               ; 3 different markers
   res@xyMarkerColors       := res@xyLineColors  

  ; add a legend
   res@pmLegendDisplayMode    = "Always"              ; turn on legend
   res@pmLegendSide           = "Top"                 ; Change location of 
   res@pmLegendParallelPosF   = .20                   ; move units right
   res@pmLegendOrthogonalPosF = -0.8                  ; more neg = down
   res@pmLegendWidthF         = 0.12                  ; Change width and
   res@pmLegendHeightF        = 0.25                  ; height of legend.
   res@lgLabelFontHeightF     = .02                   ; change font height
   res@lgPerimOn              = False                 ; no box around
 ; labels for the legend
   res@xyExplicitLegendLabels = " " + LabList 

   res@trXMinF = 2007.5
   res@trXMaxF = 2017.5

   plot  = gsn_csm_xy (wks,xx,yy,res) ; create plot

   draw(wks)                                     ; box plot does not call
   frame(wks)                             

end do 

end
