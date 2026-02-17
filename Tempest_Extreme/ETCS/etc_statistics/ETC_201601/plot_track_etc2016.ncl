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

  datdir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/ETCS/ETC_201601/"
  tcnam  = "ETC2016" 
  timtag = "2016"
  date   = yyyymmddhh_time(2007, 2017, 6, "integer")
  timplt = date({2016012000:2016012412})    ; coordinate subscripting

 ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir    = (/"ERA5", "ERA5", "NDGUV3_tau06", "NDGUVTQ3_tau06"/)
  Labels    = (/"ERA5", "ERA5", "E3SM (NDGUV)", "E3SM (NDGUVTQ)"/)
  ngrps     = dimsizes(Casdir)
  
  ;********************************
  ; select sub-regions
  ;********************************
   minlat =  10.0  ;0.0    ; deg N
   maxlat =  45.0  ;57.0    ; deg N
   minlon = 260.0  ;-100.0    ; - deg E = deg W
   maxlon = 310.0  ;-20.0    ; - deg E = deg W
   cenlon = 200.0  ;

  do i = 0,ngrps-1,1

   fna  = datdir+tcnam+"_"+Casdir(i)+"_Track_info_"+ timtag+ ".nc"
   fla  = addfile(fna, "r")

   if(i.eq.0) then 
     tmp    = fla->stim
     indxx  = ind(.not.ismissing(tmp))
     tima   = cd_calendar(fla->stim(indxx),-3)
     lon1   = fla->lon(indxx)
     lat1   = fla->lat(indxx)
     vmax1  = fla->vsmc(indxx)
     pmin1  = fla->vpmc(indxx)

     time   = new((/ngrps,dimsizes(timplt)/),integer)
     lon    = new((/ngrps,dimsizes(timplt)/),float)  
     lat    = new((/ngrps,dimsizes(timplt)/),float)
     vmax   = new((/ngrps,dimsizes(timplt)/),float)
     pmin   = new((/ngrps,dimsizes(timplt)/),float)
 
     idate     = get1Dindex(tima, timplt)
     lon(i,:)  = (/lon1(idate)/)
     lat(i,:)  = (/lat1(idate)/)
     vmax(i,:) = (/vmax1(idate)/)
     pmin(i,:) = (/pmin1(idate)/)
     time(i,:) = (/tima(idate)/)
     delete([/tmp,indxx,tima,lat1,lon1,vmax1,pmin1,idate/])
   else 
     tmp    = fla->stim
     indxx  = ind(.not.ismissing(tmp))
     if(.not.all(ismissing(indxx))) then 
       tima   = cd_calendar(fla->stim(indxx),-3)
       timx   = time(0,:)
       idate  = get1Dindex(timx, tima)
       time(i,idate) = (/cd_calendar(fla->stim(indxx),-3)/)
       lon(i,idate)  = (/fla->lon(indxx)/)
       lat(i,idate)  = (/fla->lat(indxx)/)
       vmax(i,idate) = (/fla->vsmc(indxx)/)
       pmin(i,idate) = (/fla->vpmc(indxx)/)
       delete([/idate,timx,tima/])
     end if 
     delete([/tmp,indxx/])
   end if 
  
  end do 

  ;********************************
  ; create plot
  ;********************************
   wks = gsn_open_wks("pdf","etc2016_track_traj")     ; send graphics to PNG file
   dum  = new((/ngrps,300000/),graphic)

   res                  = True 
   res@gsnMaximize      = False
   res@gsnDraw          = False   ; so we can add poly stuff
   res@gsnFrame         = False   ; do not advance frame

   BorderThick = 2.0
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

  ;res@mpProjection        = "LambertConformal" 
  ;res@mpLambertParallel1F =  20.0
  ;res@mpLambertParallel2F =  40.0
  ;res@mpLambertMeridianF  = -60.0

   res@mpLimitMode    = "LatLon"
   res@mpMinLatF       = minlat 
   res@mpMaxLatF       = maxlat
   res@mpMinLonF       = minlon
   res@mpMaxLonF       = maxlon
   res@mpCenterLonF    = cenlon

   res@mpDataBaseVersion           = "MediumRes"     ; Alias 'MediumRes'
   res@mpDataSetName               = "Earth..1"
   res@mpFillOn                    = True    ; False to turn off gray continents
   res@mpOutlineOn                 = True    ; turn on continental outline
  ;res@mpOutlineBoundarySets       = "Geophysical" ;"AllBoundaries" 
   res@mpOutlineBoundarySets       = "AllBoundaries" 
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

   left_str  = "Hurricane Sandy (2012)" 
   right_str = "" 

   res@gsnCenterString              = left_str
   res@gsnLeftString                = "" ;left_str
   res@gsnLeftStringOrthogonalPosF  = 0.04
  ;res@gsnLeftStringParallelPosF    = 0.2
   res@gsnRightString               = right_str
  ;res@gsnRightStringOrthogonalPosF = 0.04
  ;res@gsnRightStringParallelPosF   = 0.2

   map = gsn_csm_map(wks,res)            ; draw one of eight map projections

  ;***************************************
  ; Draw best track history as polylines *
  ;***************************************   

 ; Set up some legend resources.
   lgres                    = True
   lgres@lgLineColors       = (/"black","red","cyan","blue"/)
   lgres@lgLineThicknessF   = 2. 
   lgres@lgLabelFontHeightF = .08            ; set the legend label font thickness
   lgres@vpWidthF           = 0.15           ; width of legend (NDC)
   lgres@vpHeightF          = 0.1            ; height of legend (NDC)
   lgres@lgMonoDashIndex    = True   
   lgres@lgPerimColor       = "Black"       ; draw the box perimeter in orange
  ;lgres@lgPerimThicknessF = 3.0            ; thicken the box perimeter

   lbid   = gsn_create_legend(wks,ngrps,"  "+Labels,lgres)         ; create legend

; Set up resources to attach legend to map.
   amres = True
   amres@amParallelPosF   =  0.27 	         ; positive move legend to the right
   amres@amOrthogonalPosF = -0.35                 ; positive move the legend down
   annoid1 = gsn_add_annotation(map,lbid,amres)   ; attach legend to plot

; Add trajectory lines.
   pres                  = True               ; polyline resources
   pres@gsLineThicknessF = 3.0                ; line thickness
   pres@gsLineColor      = "Black"
   line1 = gsn_add_polyline(wks,map,lon(0,:),lat(0,:),pres)      ; draw the traj

   pres                  = True               ; polyline resources
   pres@gsLineColor      = "Red"
   indxx =  ind(.not.ismissing(lon(1,:)))
   line2 = gsn_add_polyline(wks,map,lon(1,indxx),lat(1,indxx),pres)      ; draw the traj
   delete(indxx)

   pres                  = True               ; polyline resources
   pres@gsLineColor      = "cyan"
   indxx =  ind(.not.ismissing(lon(2,:)))
   if(.not.all(ismissing(indxx)))
     line3 = gsn_add_polyline(wks,map,lon(2,indxx),lat(2,indxx),pres)      ; draw the traj
   end if 
   delete(indxx)

   pres                  = True               ; polyline resources
   pres@gsLineColor      = "blue"
   indxx =  ind(.not.ismissing(lon(3,:)))
   line4 = gsn_add_polyline(wks,map,lon(3,indxx),lat(3,indxx),pres)      ; draw the traj
   delete(indxx)

; Add markers to the trajectories.
   mres                = True         ; marker resources for best track
   mres@gsMarkerIndex  = 16           ; marker style (filled circle)
   mres@gsMarkerSizeF  = 8.0          ; marker size
   mres@gsMarkerColor  = "Black"      ; maker color
   markers = gsn_add_polymarker(wks,map,lon(0,:),lat(0,:),mres)

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
 gsn_panel(wks,(/map/),(/1,1/),pres1)
end
exit

