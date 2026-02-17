;********************************************************
; Concepts illustrated:
;   - Plotting best tracks for a given season's storms
;   - Manually creating a legend using markers and text
;   - Attaching text strings to a map
;   - Attaching polylines to a map plot
;   - Drawing a custom legend inside a map plot
;   - Creating a color map using named colors
;   - Customizing the fill color of various map areas
;
;********************************************************
; The purpose of this script is to plot the best tracks
; for a given season storms including all data (subtropical
; storms, depressions, extratropical lows, etc.)
;********************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  datdir = "../data/"
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir    = (/"OBS","ERA5","CLIM","NDGUVTQ3_tau06"/)
  Groups    = (/"OBS","set2","set2","set2"/)
  Labels    = (/"IBTrACS", "ERA5", "E3SM(CLIM)","E3SM(NDGUVTQ)"/)
  ngrps     = dimsizes(Groups)
  LabList   = Groups
  plot      = new (ngrps, graphic)
  ystr      = 2007
  yend      = 2017
  timtag    = ystr+"-"+yend
  
  ;********************************
  ; select sub-regions
  ;********************************
   minlat = -70.0 ;0.0    ; deg N
   maxlat =  70.0 ;57.0    ; deg N
   minlon =   0.0 ;-100.0    ; - deg E = deg W
   maxlon = 360.0 ;-20.0    ; - deg E = deg W
   cenlon = 200.0 ;

  do i = 0,ngrps-1,1

   fna  = datdir+Casdir(i)+"_"+Groups(i)+"_Track_info_"+ timtag+ ".nc"
   fla  = addfile(fna, "r")
   lata = fla->lat
   lona = fla->lon
   tima = fla->stim
   vsmc = fla->vsmc
   vpmc = fla->vpmc

   nstorms  = dimsizes(lata(:,0))
   max_npts = dimsizes(lata(0,:))

  ;********************************
  ; create plot
  ;********************************

   if(i.eq.0) then 
     wks = gsn_open_wks("pdf","track_traj_e3sm_set2")     ; send graphics to PNG file
     dum  = new((/ngrps,300000/),graphic)
   end if 

   res                  = True 
   res@gsnMaximize      = False
   res@gsnDraw          = False   ; so we can add poly stuff
   res@gsnFrame         = False   ; do not advance frame

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

   FontHeightF = 0.020
   res@tiMainFontHeightF     = FontHeightF
   res@tmYLLabelFontHeightF  = FontHeightF*0.8
   res@tmXBLabelFontHeightF  = FontHeightF*0.8
   res@gsnStringFontHeightF  = FontHeightF*0.9
  ;res@tmXBMajorLengthF      = -0.001

   res@mpLimitMode  = "LatLon"
   res@mpMinLatF    = minlat 
   res@mpMaxLatF    = maxlat
   res@mpMinLonF    = minlon
   res@mpMaxLonF    = maxlon
   res@mpCenterLonF = cenlon

   res@mpDataBaseVersion           = "MediumRes"     ; Alias 'MediumRes'
   res@mpDataSetName               = "Earth..1"
   res@mpFillOn                    = True    ; False to turn off gray continents
   res@mpOutlineOn                 = True    ; turn on continental outline
   res@mpOutlineBoundarySets       = "Geophysical" ;"AllBoundaries" 
   res@mpLandFillColor             = "burlywood" ;"Yellow2enRod1"
   res@mpInlandWaterFillColor      = "PaleTurquoise3"
   res@mpOceanFillColor            = "PaleTurquoise3"
   res@mpGeophysicalLineColor      = "Grey37"
   res@mpGeophysicalLineThicknessF = 0.5
   res@mpUSStateLineColor          = "Grey37"
   res@mpUSStateLineThicknessF     = 0.5
   res@mpNationalLineColor	   = "Grey37" 
   res@mpNationalLineThicknessF    = 0.5
   res@mpGridAndLimbOn             = "False"
   res@mpGridAndLimbDrawOrder      = "Draw"
   res@mpGridMaskMode              = "MaskLand"
   res@mpGridSpacingF              = 5.0
   res@mpGridLineColor             = "Grey37"
   res@pmTickMarkDisplayMode       = "Always"

   left_str  = Labels(i)
   right_str = "Total # of TCs: " + nstorms

   res@gsnLeftString                = left_str
  ;res@gsnLeftStringOrthogonalPosF  = 0.04
  ;res@gsnLeftStringParallelPosF    = 0.2
   res@gsnRightString               = right_str
  ;res@gsnRightStringOrthogonalPosF = 0.04
  ;res@gsnRightStringParallelPosF   = 0.2

   plot(i) = gsn_csm_map(wks,res)            ; draw one of eight map projections

  ;***************************************
  ; Draw best track history as polylines *
  ;***************************************   

   res_poly                    = True			; polyline mods desired

   ic = 0
   do nc = 0,nstorms-1

    utc_date = cd_calendar(tima(nc,:),0)
    year     = tofloat(utc_date(:,0))
    month    = tofloat(utc_date(:,1))
    day      = tofloat(utc_date(:,2))
    hour     = tofloat(utc_date(:,3))

    lat      = tofloat(lata(nc,:))
    lon      = tofloat(lona(nc,:))
   ;vmax     = tofloat(vsmc(nc,:))
    pmin     = tofloat(vpmc(nc,:))
 
    pmin@_FillValue = vpmc@_FillValue 
    lat@_FillValue  = lata@_FillValue
    lon@_FillValue  = lona@_FillValue

   ;convert to hPa
    if(.not.all(ismissing(pmin)).and.min(pmin).ge.1050.0) then
      pmin = pmin / 100.0
    end if

   ;Atkinson and Holliday (1977): Vmax = 6.7(1010-pmin)^0.644
   ;Knaff and Zehr (2007): Vmax = 6.6(1010-MSLP)^0.65, (appendix A)                          
   if(Groups(i).eq."OBS") then 
     vmax     = tofloat(vsmc(nc,:))
   else 
    ;vmax = 6.7*(1010.0 - pmin)^0.644
    ;vmax = 6.6*(1010.0 - pmin)^0.65
     vmax = 4.4*(1010.0 - pmin)^0.76
   end if 

   do k=0,max_npts-2

      if(.not.(ismissing(vmax(k)).or.ismissing(lat(k)).or.ismissing(lon(k)))) then 
        ;print(lat(k)+ " " +lon(k)+" "+vmax(k))
         res_poly@gsLineDashPattern  = 0
         res_poly@gsLineThicknessF   = 1.2   
         res_poly@gsMarkerColor      = "transparent"
         res_poly@gsMarkerIndex      = 1
         res_poly@gsMarkerSizeF      = 0.001
         ;Saffir-Simpson Hurricane Wind Scale
         ;Category 4 Hurricane: Winds 130-156 mph (113-136 kt or 209-251 km/hr)
         ;Category 5 Hurricane: Winds 157 mph or higher (137 kt or higher or 252 km/hr or higher). 
         if((vmax(k).ge.58)) then 
            res_poly@gsLineColor        = "Red"                        
            dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         end if

         ;Category 3 Hurricane: Winds 111-129 mph (96-112 kt or 178-208 km/hr). 
         if((vmax(k).lt.58).and.(vmax(k).ge.50)) then 
            res_poly@gsLineColor        = "salmon3"                 
            dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         end if  

         ;Category 2 Hurricane: Winds 96-110 mph (83-95 kt or 154-177 km/hr). 
         if((vmax(k).lt.50).and.(vmax(k).ge.43)) then
            res_poly@gsLineColor        = "Orange"
            dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         end if

         ;Category 1 Hurricane: Winds 74-95 mph (64-82 kt or 119-153 km/hr)
         if((vmax(k).lt.43).and.(vmax(k).ge.33)) then
            res_poly@gsLineColor        = "Yellow2"                 
            dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         end if 

         ;Tropical storm: Winds 39-73 mph (34-63 kt or 63-118 km/hr)
         if((vmax(k).lt.33).and.(vmax(k).ge.18)) then 
            res_poly@gsLineColor        = "DarkGreen"
            dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         end if

         ;Tropical depression: Winds < 38 mph (33 kt or 62 km/hr)
         ;if((vmax(k).lt.18)) then
         ;   res_poly@gsLineColor        = "DarkGreen"
         ;   dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         ;end if

       ic = ic + 1

      end if

   end do

   delete([/utc_date,year,month,day,hour,lat,lon,vmax/])

 end do
  
;******************
; Plot a legend   *
;******************
 lgres                       = True
 lgres@lgAutoManage          = False
 lgres@vpWidthF              = 0.15     ; was 0.08        ; width of legend (NDC)
 lgres@vpHeightF             = 0.12     ; was 0.08        ; height of legend (NDC)
;lgres@lgBottomMarginF       = 0.17     ; was 0.25
 lgres@lgPerimFill           = 0                     ; Use solid fill (0) instead of the default hollow fill
 lgres@lgPerimFillColor      = "Background"
;lgres@lgBoxMajorExtentF     = 0.4
 lgres@lgBoxMinorExtentF     = 0.2 
;lgres@lgBoxBackground       = "PaleTurquoise3"
 
 lgres@lgMonoItemType        = False                 ; indicates that we wish to set the item types individually
 lgres@lgMonoMarkerIndex     = False
 lgres@lgMonoLineThickness   = False
 lgres@lgMonoMarkerThickness = False
 lgres@lgMonoMarkerSize      = False

 ; Position fine elements of legend relative to some positional values:
;lgres@lgLabelFont             = 0
 lgres@lgLabelFontHeightF      = FontHeightF*4.0
;lgres@lgLabelFontAspectF      = 0.9
;lgres@lgLabelConstantSpacingF = 0.05

 lgres@lgLineDashSegLenF  = 0.1
 lgres@lgItemTypes        = (/"Markers","Markers","Markers","Markers","Markers"/)
 lgres@lgMarkerIndexes    = (/      16,        16,       16,     16,     16/)
 lgres@lgLineThicknesses  = (/      0.1,      0.1,      0.1,    0.1,    0.1/)
 lgres@lgMarkerColors     = (/  "Red", "salmon3", "Orange", "Yellow2", "DarkGreen"/)
 lgres@lgMarkerSizes      = (/   0.004,   0.004,    0.004, 0.004,    0.004/) * 2.0
 lgres@lgLineColors       = (/  "Red", "salmon3", "Orange", "Yellow2", "DarkGreen"/)
 lgres@lgDashIndexes      = (/       0,      0,     0, 0, 0/)
 legend_labels            = (/"Category 4/5",  "Category 3", "Category 2", "Category 1","Tropical Storm"/)
 legend                   = gsn_create_legend(wks,dimsizes(legend_labels),legend_labels,lgres)

 amres = True
 amres@amParallelPosF   = -0.49
 amres@amOrthogonalPosF = -0.49
 amres@amJust           = "TopLeft"
 annoid1 = gsn_add_annotation(plot(i),legend,amres)

 delete([/lata,lona,tima,vsmc,vpmc/])

end do

;---Draw both plots in a panel
;***************************************
; panel first set of plots
;***************************************
 pres1                          = True
 pres1@gsnFrame                 = True      ; don't advance frame yet
 pres1@gsnMaximize              = False
 pres1@gsnPanelLabelBar         = False
 pres1@gsnPanelMainString       = " "
 pres1@gsnPanelMainFontColor    = "Black"
 pres1@gsnPanelMainFontHeightF  = 0.028
 pres1@gsnPanelYWhiteSpacePercent = 5.0
 pres1@gsnPanelXWhiteSpacePercent = 5.0
 gsn_panel(wks,(/plot/),(/ngrps,1/),pres1)
end
exit

