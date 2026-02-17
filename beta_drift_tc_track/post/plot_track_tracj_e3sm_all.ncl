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


  datadir = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/out"
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups    = (/"CLIM","NDGUV"/) ;,"NDGUVT","NDGUVTQ"/)
  ngrps     = dimsizes(Groups)
  Labels    = Groups
  plot      = new (ngrps, graphic)
  ystr      = 2007
  yend      = 2017
  nyrs      = yend - ystr +1
  timtag    = ystr+"-"+yend

  plot      = new (ngrps, graphic)
  markers   = new (ngrps, graphic)
  texts     = new (ngrps, graphic)
  ystr      = 2007
  yend      = 2017
  timtag    = ystr+"-"+yend
  
  ;********************************
  ; select sub-regions
  ;********************************
   minlat =  0.0 ;0.0    ; deg N
   maxlat =  70.0 ;57.0    ; deg N
   minlon = 200.0 ;-100.0    ; - deg E = deg W
   maxlon = 350.0 ;-20.0    ; - deg E = deg W
   cenlon = 200.0 ;

  do i = 0,ngrps-1,1

   fls    = systemfunc("ls "+ datadir + "/sythetic_track_forecast_" + Groups(i)+"_*.nc")
   f      = addfiles(fls(0:1),"r")
   ListSetType (f, "join")        ; concatenate (=default)

   nyrs    = dimsizes(fls)
   vartrk  = f[:]->strks
   varnam  = str_squeeze(f[0]->sname)
   vartim  = f[:]->stime

  ;mark fill values
   vartrk@_FillValue = -99999
   vartrk = where(vartrk.eq."nan".or.vartrk.eq."-nan",vartrk@_FillValue,vartrk)

   dimx    = dimsizes(vartrk)
   nyrs    = dimx(0)
   nvars   = dimx(1)
   ntrks   = dimx(2) 
   nsteps  = dimx(3) 

   indlat  = ind(varnam.eq."lat")
   indlon  = ind(varnam.eq."lon")
   indtim  = ind(varnam.eq."time")
   lat            = reshape(vartrk(:,indlat,:,:),(/nyrs*ntrks,nsteps/))
   lon            = reshape(vartrk(:,indlon,:,:),(/nyrs*ntrks,nsteps/))
   tmp            = reshape(vartrk(:,indtim,:,:),(/nyrs*ntrks,nsteps/)) 
   stim           = dim_cumsum_n(tmp,0,0)
   delete(tmp)

   lat!0          = "storm"
   lat!1          = "time"
   lat@units      = "degree_north"
   lat@long_name  = "latitude"

   lon!0          = "storm"
   lon!1          = "time"
   lon@units      = "degree_east"
   lon@long_name  = "longitude"

   stim!0         = "storm"
   stim!1         = "time"
   stim@units     = "days since "+(ystr)+"-12-31 00:00:00"
   stim@long_name = "date"

   latx = lat(:,0)
   lonx = lon(:,0)
   inds = ind (latx.ge.minlat.and.latx.le.maxlat.and.lonx.ge.minlon.and.lonx.le.maxlon)
 
   lata = lat(inds,:)
   lona = lon(inds,:)
   tima = stim(inds,:)

   nstorms  = dimsizes(lata(:,0))
   max_npts = dimsizes(lata(0,:))

   print("number of storms: "+ nstorms)

   delete([/inds,latx,lonx,lat,lon,stim,dimx,indlat,indlon,indtim,vartrk,varnam,vartim/])

  ;********************************
  ; create plot
  ;********************************

   if(i.eq.0) then 
     wks = gsn_open_wks("pdf","track_traj_e3sm_"+Groups(i))     ; send graphics to PNG file
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
   res@mpOutlineBoundarySets       = "AllBoundaries" ;"Geophysical" ;"AllBoundaries" 
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
   do nc = 0,5; nstorms-1

    utc_date = cd_calendar(tima(nc,:),0)
    year     = tofloat(utc_date(:,0))
    month    = tofloat(utc_date(:,1))
    day      = tofloat(utc_date(:,2))
    hour     = tofloat(utc_date(:,3))

    lat      = tofloat(lata(nc,:))
    lon      = tofloat(lona(nc,:))
   ;vmax     = tofloat(vsmc(nc,:))
    lat@_FillValue  = tofloat(lata@_FillValue)
    lon@_FillValue  = tofloat(lona@_FillValue)

   do k=0,max_npts-2

      if((.not.all(ismissing(lat(k:k+1)))).and.(.not.all(ismissing(lon(k:k+1))))) then 
        ;print(lat(k)+ " " +lon(k)+" "+vmax(k))
         res_poly@gsLineDashPattern  = 0
         res_poly@gsLineThicknessF   = 2.0   
         res_poly@gsMarkerColor      = "transparent"
         res_poly@gsMarkerIndex      = 1
         res_poly@gsMarkerSizeF      = 0.001
         res_poly@gsLineColor        = "Black"
         dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)

         ;Tropical depression: Winds < 38 mph (33 kt or 62 km/hr)
         ;if((vmax(k).lt.18)) then
         ;   res_poly@gsLineColor        = "DarkGreen"
         ;   dum(i,ic) = gsn_add_polyline(wks,plot(i),lon(k:k+1),lat(k:k+1),res_poly)
         ;end if

       ic = ic + 1

      end if

   end do

   delete([/utc_date,year,month,day,hour,lat,lon/])

 end do
 
;******************
; Plot a legend   *
;******************
 lgres                       = True
 lgres@lgAutoManage          = False
 lgres@vpWidthF              = 0.18     ; was 0.08        ; width of legend (NDC)
 lgres@vpHeightF             = 0.14     ; was 0.08        ; height of legend (NDC)
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
;annoid1 = gsn_add_annotation(plot(i),legend,amres)

 delete([/lata,lona,tima/])

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
 gsn_panel(wks,(/plot/),(/2,2/),pres1)
end
exit

