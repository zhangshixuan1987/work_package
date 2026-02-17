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
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  ;TC Basin definitions.(Knutson et. al.,2020,BAMS)
  l_plot_latlon_box = False ; if True, add lat-lon box for region
  Basin = (/"NA","NEP","NWP","NI","SI","SWP","SA1","SA2"/)
  Blats = (/   0,    0,    0,   0, -90,  -90,  -90,  -90/)
  Blatn = (/  90,   90,   90,  90,   0,    0,    0,    0/)
  Blonw = (/ 265,  180,  100,  30,  20,  135,  295,    0/)
  Blone = (/ 360,  265,  180, 100, 135,  295,  360,   20/)
  nbsin = dimsizes(Basin)

  datdir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS/track_info"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  LabList   = (/"ERA5", "E3SM", "E3SM (UV-NDG)", "E3SM (UVT-NDG)", "E3SM (UVTQ-NDG)"/)
  Groups    = (/"ERA5", "CLIM", "NDGUV", "NDGUVT", "NDGUVTQ"/)
  timtag    = (/"1979-2017","2007-2017","2007-2017","2007-2017","2007-2017"/)
  ngrps     = dimsizes(Groups)
  sset      = "set1"
  figkey    = "trkdens_pie"

  Figstr    = (/"(a)", "(b)", "(c)", "(d)","(e)"/)
  ystr      = 2008
  yend      = 2017
  nyr       = yend - ystr +1

  do i = 0,ngrps-1,1

   fna    = datdir+"/"+Groups(i)+"_"+sset+"_Track_statistics_"+timtag(i)+".nc"
   fla    = addfile(fna, "r")
   lat    = fla->lat
   lon    = fla->lon
   time   = fla->time 

   rad    = 4.0*atan(1.0)/180.0
   clat   = lat 
   clat   = cos(lat*rad)
   
   Varnam = "trkdens"
   Varstr = "Track density"
   Varunt = "storm trainsts per month"

   ;for NH calculate May-November; for SH calculate November-April.
   tim1 = new(nyr*7,integer)
   tim2 = new(nyr*6,integer)
   do ii = 0,nyr-1 
     kk1 = 0 
     kk2 = 0
     do jj = 1, 12
       if(jj.ge.5.and.jj.le.11)then 
         tim1(ii*7+kk1) = toint(sprinti("%04d",ystr + ii) + sprinti("%02d",jj))
         kk1 = kk1 + 1
       end if 
       if(.not.(jj.ge.5.and.jj.lt.11))then
         tim2(ii*6+kk2) = toint(sprinti("%04d",ystr + ii) + sprinti("%02d",jj))
         kk2 = kk2 + 1
       end if
     end do
   end do 

   indx1          = get1Dindex(time,tim1) 
   indx2          = get1Dindex(time,tim2)
   var0           = dim_avg_n_Wrap(fla->$Varnam$(indx1,:,:),0) 
   var1           = dim_avg_n_Wrap(fla->$Varnam$(indx2,:,:),0) 
   var0           = var0 * 7.0  ; convert to per year 
   var1           = var1 * 6.0  ; vonvert to per year 
   var            = var0
   var({5:90},:)  = var0({5:90},:)
   var({-90:5},:) = var1({-90:5},:) 
   delete([/indx1,indx2,var0,var1,tim1,tim2/])
   delete([/fna,fla,lat,lon,time/])

   print(" ")
   print("........................... ")
   print("igp : " + i) 
   print("........................... ")

   printMinMax(var,False)

   ;calculate statistics on each basin 
   bsvar = new((/nbsin+1/),typeof(var))
   bname = new((/nbsin+1/),typeof(Basin))

   do ib = 0,nbsin-1,1
     bsvar(ib) = avg(var({Blats(ib):Blatn(ib)},{Blonw(ib):Blone(ib)})) ;,clat({Blats(ib):Blatn(ib)}),1.0,0)
     bname(ib) = Basin(ib)
   end do 

  ;SA basin = SA1 + SA2
   bsvar(nbsin-2) = bsvar(nbsin-2) + bsvar(nbsin-1)
   bname(nbsin-2) = "SA"

  ;NH and SH average 
   bsvar(nbsin-1) = avg(var({0:90},:)) ;,clat({0:90}),1.0,1)
   bsvar(nbsin)   = avg(var({-90:0},:)) ;,clat({-90:0}),1.0,1) 
   bname(nbsin-1) = "NH"
   bname(nbsin)   = "SH"
   bsvar(0:nbsin-2) = bsvar(0:nbsin-2)*100.0 /sum(bsvar(0:nbsin-2)) ;(bsvar(nbsin-1)+ bsvar(nbsin))

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(i.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "track_" + figkey 
    wks1     = gsn_open_wks(wtype,plotname)
    plot     = new (ngrps, graphic)
    regx     = new ((/ngrps,nbsin,4/), graphic) 

    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")

   end if 

   ;;;set up the plot information;;;;
   lndcolor   = "Black" ;"Grey50"
   cnlevels  := (/1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 20/)
   cnlevel2  := (/-5,-4,-3,-2,-1,1,2,3,4,5/)

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 0.5
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeightF = 0.018
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF
   res@gsnStringFontHeightF = FontHeightF*2.0

   res@gsnCenterString      = ""
   res@gsnLeftString        = ""
   res@gsnRightString       = "" 
   res@tiMainString         = "" 


   res@tiMainString         = Figstr(i) + "  " +LabList(i)
  ;res@gsnLeftString               = Figstr(i) + "  " +LabList(i)
  ;res@gsnLeftStringOrthogonalPosF = 0.04
  ;res@gsnRightString              = Varstr
  ;res@gsnRightStringOrthogonalPosF= 0.04
 

  ;cmap   = read_colormap_file("cb_9step")
   xcolor = (/ "red", "green", "orange", "yellow", "brown","magenta","blue","gray" /)

   percent = bsvar(0:nbsin-2) 

  ;res@pcLabelType = "block"
   percent = bsvar(0:nbsin-2)
   name    = bname(0:nbsin-2) 
   color   = xcolor(0:nbsin-2) 
   plot(i) = pie_chart(wks1, percent, name, color, res)

   delete([/res/])

end do ; end of experiment group loop

;---Draw both plots in a panel
;***************************************
; panel first set of plots
;***************************************
  pres1                          = True
  pres1@gsnFrame                 = True      ; don't advance frame yet
  pres1@gsnMaximize              = False
  pres1@gsnPanelLabelBar         = False
  pres1@pmLabelBarParallelPosF   = -0.010
  pres1@pmLabelBarOrthogonalPosF = 0.0
  pres1@lbLabelFontHeightF       = FontHeightF * 0.45
  pres1@lbOrientation            = "vertical"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.21
  pres1@pmLabelBarWidthF         = 0.05
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = 0.028
  gsn_panel(wks1,(/plot/),(/1,ngrps/),pres1)

end
