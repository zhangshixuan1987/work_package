;read ascii data (OMI & GOME2) and write to netcdf file
; 5th Nov. 2014
; by Jianfeng Li

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
; This file still has to be loaded manually
load "./shapefile_utils.ncl"

undef("draw_data_grid_lines")
procedure draw_data_grid_lines(wks,map,var)
local lats, lons, numlon, numlat, lon_step, lat_step, ilat, ilon, x, y
begin
 BOXES = 1   ;-- 0: lat/lon grid lines in cell center
 ;-- 1: lat/lon grid lines around cells (=boxes)
 lats = var&lat
 lons = var&lon
 lat_step = lats(1) - lats(0)
 lon_step = lons(1) - lons(0)
 numlat = dimsizes(lats)
 numlon = dimsizes(lons)

 lonadd = new(numlon+1,typeof(lons))   ;-- max. longitude less 360. so we have to add it
 lonadd(0:numlon-1) = lons
 lonadd(numlon) = lons(numlon-1)+lon_step

;-- grid line resources
 gres = True
 gres@gsLineColor = "grey50"   ;-- grid line color
 gres@gsLineThicknessF = 1.7   ;-- grid line thickness

;-- assign x and x arrays for the grid lines
 x = new(numlon+1,typeof(lons))
 y = new(numlon+1,typeof(lats))

;-- draw cell centers
 if(BOXES .eq. 0) then
    dx = 0.   ;-- cell center
    dy = 0.   ;-- cell center
 else
    dx = lon_step/2.   ;-- boxes
    dy = lat_step/2.   ;-- boxes
 end if

;-- draw horizontal lines
 x = lonadd + dx
 do ilat = 0,numlat-1
    y = lats(ilat) + dy
    str = unique_string("gridline") ;-- create a unique name
    map@$str$ = gsn_add_polyline(wks, map, x, y, gres) ;-- add polyline to map
 end do

;-- draw vertical lines
 y := fspan(min(lats),max(lats),numlon) + dy
 do ilon = 0,numlon
    x = lonadd(ilon) + dx
    str = unique_string("gridline")   ;-- create unique name
    map@$str$ = gsn_add_polyline(wks, map, x, y, gres)   ;-- add polyline to map
 end do

end ;-- end procedure draw_data_grid_lines


undef("read_base_map_data")
function read_base_map_data(filename)
local filename,Band1, Band2, Band3,latb,lonb,nlat,nlon
begin
  jpg_filename = filename + ".jpg"
  nc_filename  = filename + ".nc"
  system("rm -rvf " + nc_filename)
  system("gdal_translate -ot Int16 -of netCDF " + jpg_filename + " " + nc_filename)
  f     = addfile(nc_filename,"r")
  Band1 = where(f->Band1.gt.255, 255, f->Band1)  ; red channel
  Band2 = where(f->Band2.gt.255, 255, f->Band2)  ; green channel
  Band3 = where(f->Band3.gt.255, 255, f->Band3)  ; blue channel
  band_dims = dimsizes(Band3)
  nlat      = band_dims(0)
  nlon      = band_dims(1)
  latb      = fspan( -90, 90,nlat)
  lonb      = fspan(-180,180,nlon)
  latb@units = "degrees_north"
  lonb@units = "degrees_east"
  Band1!0   = "lat"
  Band1!1   = "lon"
  Band2!0   = "lat"
  Band2!1   = "lon"
  Band3!0   = "lat"
  Band3!1   = "lon"
  Band1&lat = latb
  Band1&lon = lonb
  Band2&lat = latb
  Band2&lon = lonb
  Band3&lat = latb
  Band3&lon = lonb
  delete([/band_dims,nlat,nlon,latb,lonb/])
  return[/Band1,Band2,Band3/]
end

undef("read_ocean_land_colormap")
function read_ocean_land_colormap(num_ocean)
local cmap_ocn, cmap_lnd
begin
  cmap_ocn = read_colormap_file("GMT_ocean")
  cmap_lnd = read_colormap_file("OceanLakeLandSnow")
  newcmap = array_append_record(cmap_ocn(0:num_ocean-1,:),cmap_lnd(2::2,:),0)
  return(newcmap)
end

;----------------------------------------------------------------------
; Given min/max elevation value and the number of colors you want
; for ocean, this function creates a custom color map AND calculates 
; levels for both land and ocean, based on an elevation pivot point.
;
; The assumption is that the end ocean index location is one
; less than the begining of the land index location.
;----------------------------------------------------------------------
undef("calc_levels_and_colors")
function calc_levels_and_colors(wks,emin,emax,split_elev,num_ocean_values)
local start_ocean, ocean_range, land_range, olevels, llevels, nol, nll, clen
begin
  cmap = read_ocean_land_colormap(num_ocean_values)
  clen = dimsizes(cmap(:,0))

  start_ocean = 0
  end_ocean   = num_ocean_values-1
  start_land  = end_ocean+1
  ocean_range = end_ocean-start_ocean+1
  land_range  = clen-start_land+1
  olevels     = fspan(emin,split_elev,ocean_range)
  llevels     = fspan(split_elev,emax,land_range)
  nol         = dimsizes(olevels)
  nll         = dimsizes(llevels)
  levels      = new((nol-1)+(nll-2),float)
  levels(0:nol-2) = olevels(1:)
  levels(nol-1:)  = llevels(1:nll-2)
  return([/levels,cmap/])
end

;----------------------------------------------------------------------
; This function reads a binary file containing elevation data and
; generates the necessary lat/lon coordinate arrays for plotting later.
; The information on the binary file is provided at the beginning of
; this script.
;
; The binary file was downloaded from:
;     http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO5/TOPO/ETOPO5/
;----------------------------------------------------------------------
undef("read_elev_data")
function read_elev_data(topo_file)
local nlat, nlon, topo_file, lat, lon
begin
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  f2        = addfile (topo_file, "r")
  elev      = tofloat(f2->PHIS)
  elev      = elev/9.80665
  return(elev)
end

begin

  month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

  datdir    = "/pscratch/sd/z/zhan391/DARPA_project/post_process/AR_201702/"
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Explist   = (/"E3SM_ARNEP_CTRL"/)
  Labels    = (/"EAMxx/SCREAM"/)
  ngrps     = dimsizes(Explist)
  timtag    = "2017-02-03_2017-02-19"

  fcord     = addfile("/pscratch/sd/z/zhan391/DARPA_project/post_process/AR_201702/terrian/E3SM.coord.nc","r")
  lat       = fcord->lat
  lon       = fcord->lon
 ;lon       = where(lon.gt.180.0,lon-360.0,lon)

  topo_file = "/pscratch/sd/z/zhan391/DARPA_project/post_process/AR_201702/terrian/E3SM.PHIS.nc"
  elev      = read_elev_data(topo_file)
  elev@lat  = lat 
  elev@lon  = lon 

  aa=addfile("../domain/ne32x32_ARNEP_scrip.nc","r")
  grid_corner_lat=aa->grid_corner_lat
  grid_corner_lon=aa->grid_corner_lon
  grid_center_lat=aa->grid_center_lat
  grid_center_lon=aa->grid_center_lon
  delete(aa)
  grid_corner_lon=mod(grid_corner_lon+180.,360.)-180.
  grid_center_lon=mod(grid_center_lon+180.,360.)-180.
  dim=dimsizes(grid_corner_lon)
  ngrid=dim(0)
  ncorner=dim(1)
  delete(dim)

  latmin = -90.
  latmax = 90.
  lonmin = 0   ;-180.
  lonmax = 360 ;180.
  edgell = 0.

 ;Feather River Basin region   
  lats      = 38
  latn      = 41
  lonw      = -122.5 + 360 
  lone      = -119.5 + 360 
  ;Oroville Dam 
  plot_oroville_dam = True 
  lat_od = 39 + (32.0+20.0/60.0)/60.0
  lon_od = 121 + (29.0+ 8.0/60.0)/60.0
  lon_od = 360.0 - lon_od

;---Construct RGBA colormaps...
  mapfile = "EarthMap_2500x1250"
  bnddat  = read_base_map_data(mapfile)
  Band1   = bnddat[0]
  Band2   = bnddat[1]
  Band3   = bnddat[2]
  band_dims = dimsizes(Band3)
  ramp   = fspan(0., 1., 255)
  reds   = new((/255, 4/), float)
  greens = new((/255, 4/), float)
  blues  = new((/255, 4/), float)
  reds   = 0
  greens = 0
  blues  = 0
  reds(:,0)   = ramp
  greens(:,1) = ramp
  blues(:,2)  = ramp

 ; The red contour map is plotted fully opaque; the green and blue
 ; are plotted completely transparent. When overlain, the colors 
 ; combine (rather magically).
  reds(:,3)   = 1.
  greens(:,3) = 0 
  blues(:,3)  = 0

 do i = 0,ngrps-1,1

  expnam     = Explist(i)
  flnm0      = datdir + expnam+"/"+expnam +".TUQ."+timtag+".nc"
  f0         = addfile (flnm0, "r")
  flnm1      = datdir + expnam+"/"+expnam +".TVQ."+timtag+".nc"
  f1         = addfile (flnm1, "r")
  flnmx0     = datdir + expnam+"/SE/"+expnam +".TUQ."+timtag+".nc"
  fx0        = addfile (flnmx0, "r")
  flnmx1     = datdir + expnam+"/SE/"+expnam +".TVQ."+timtag+".nc"
  fx1        = addfile (flnmx1, "r")
  flnmy0     = datdir + expnam+"/SE/"+expnam +".PRECT."+timtag+".nc"
  fy0        = addfile (flnmy0, "r")

  tim0       = cd_calendar(fx0->time,-3)
  indtim     = ind(tim0.ge.2017020300.and.tim0.le.2017021200)
 ;indtim     = ind(tim0.ge.2017020712.and.tim0.le.2017081200)
  delete(tim0)

  do k = 0,dimsizes(indtim)-1,1

  ivt = (/sqrt(fx0->TUQ(indtim(k),:)*fx0->TUQ(indtim(k),:)+ \
               fx1->TVQ(indtim(k),:)*fx1->TVQ(indtim(k),:))/)
  ivt@lat        = elev@lat
  ivt@lon        = elev@lon
  ivt@_FillValue = 1e20

  prc            = fy0->PRECT(indtim(k),:)
  prc@lat        = elev@lat
  prc@lon        = elev@lon
  prc@_FillValue = 1e20
  prc            = prc * 1000.0 * 3600.0
  prc@units      = "mm hr~S~-1~N~"

  var        = new((/3,dimsizes(f0->lat),dimsizes(f0->lon)/),typeof(lat))
  var(0,:,:) = f0->TUQ(indtim(k),:,:)
  var(1,:,:) = f1->TVQ(indtim(k),:,:)
  var(2,:,:) = (/sqrt(var(0,:,:)^2+var(1,:,:)^2)/)
  var!0      = "nvars"
  var!1      = "lat"
  var!2      = "lon"
  var&lat    = f0->lat
  var&lon    = f0->lon

  tdate     = cd_calendar(fx0->time(indtim(k)),0) 
  year      = toint(tdate(:,0))
  month     = toint(tdate(:,1))
  day       = toint(tdate(:,2))
  hour      = toint(tdate(:,3))
  mnt       = toint(tdate(:,4))
  
  timstr    = sprinti("%0.2i", hour)+sprinti("%0.2i", mnt)+"UTC " + sprinti("%0.2i ", day) + \
              month_abbr(month) + " "  + sprinti("%0.4i", year)
  delete([/tdate,year,month,day,hour/])

  timplt    = sprintf("%.12f",cd_calendar(fx0->time(indtim(k)),4))
  print("working on " + timstr)

 ;********************************
 ; RRM regions
 ;********************************
  minlat = 18.0  ;0.0    ; deg N
  maxlat = 60.0  ;57.0    ; deg N
  minlon = 180.0 ;-100.0    ; - deg E = deg W
  maxlon = 270.0 ;-20.0    ; - deg E = deg W
  
  indl   = ind(.not.(lat.le.maxlat.and.lat.ge.minlat.and.lon.ge.minlon.and.lon.le.maxlon))
  rivt   = ivt 
  rprc   = prc
  rlat   = lat
  rlon   = lon
  rivt(indl)   = rivt@_FillValue 
  rprc(indl)   = rprc@_FillValue 
  indr   = ind(lat.le.maxlat.and.lat.ge.minlat.and.lon.ge.minlon.and.lon.le.maxlon)
  latx   = var&lat
  lonx   = var&lon
  latx!0 = "lat"
  lonx!0 = "lon"
  latx&lat = latx
  lonx&lon = lonx
  rlat2d   = conform(var(0,{minlat:maxlat},{minlon:maxlon}),latx({minlat:maxlat}),0) 
  rlon2d   = conform(var(0,{minlat:maxlat},{minlon:maxlon}),lonx({minlon:maxlon}),1) 
  vdum2d   = var(0,{minlat:maxlat},{minlon:maxlon})
  vdum2d   = 1.0

  if (i.eq.0) then 
    figtype="png"
    figtype@wkWidth  = 2500
    figtype@wkHeight = 2500
    DIR = "ivt_sub_fig"
    system("if ! test -d " + DIR +" ; then mkdir " + DIR + " ; fi")
    fignam = DIR+"/ivt_ar_" + timplt 
    wks = gsn_open_wks(figtype,fignam)
    plot      = new (ngrps, graphic)
    vector    = new (ngrps, graphic)
    contour  = new (ngrps, graphic)
    contourx = new (ngrps, graphic)
   ;cmap      = read_colormap_file("precip4_11lev")
    cmap      = read_colormap_file("precip3_16lev")
   ;cmap      = read_colormap_file("WhBlGrYeRe")
    cnlevs    = (/250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000/)*1.0
   ;cnlevs    = fspan(200,800,101) 
  end if 
  
  res                  = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@gsnMaximize      = False
  res@tmXTOn           = False
  res@tmYROn           = False
  res@gsnLeftString    = ""
  res@gsnRightString   = ""
  res@tiMainString     = ""

  FontHeightF = 0.030
  res@tiMainFontHeightF    = FontHeightF
  res@tmYLLabelFontHeightF = FontHeightF
  res@tmXBLabelFontHeightF = FontHeightF
  res@gsnStringFontHeightF = FontHeightF
  res@tiXAxisFontHeightF   = FontHeightF
  res@tiYAxisFontHeightF   = FontHeightF

  res@tmYLLabelsOn      = False
  res@tmYRLabelsOn      = False
  res@tmXBLabelsOn      = False
  res@tmXTLabelsOn      = False

  res@tmYLOn            = False
  res@tmYROn            = False
  res@tmXBOn            = False
  res@tmXTOn            = False

  res@tmXTBorderOn      = False 
  res@tmXBBorderOn      = False
  res@tmYRBorderOn      = False
  res@tmYLBorderOn      = False

  resx = res
  resx@cnLinesOn       = False              ; Turn off contour lines      .
  resx@cnLineLabelsOn  = False              ; Turn off contour labels
  resx@cnInfoLabelOn   = False              ; Turn off info label
  resx@gsnTickMarksOn  = False              ; Turn off tickmarks
  resx@lbLabelBarOn    = False              ; Turn off labelbar
  resx@gsnRightString  = ""                 ; Turn off subtitles
  resx@gsnLeftString   = ""
  resx@cnFillOn        = True         ; Turn on filled rasters
  resx@cnFillMode      = "RasterFill"
  resx@cnLevelSelectionMode  = "EqualSpacedLevels"
  resx@cnMaxLevelCount       = 254
  resx@cnFillBackgroundColor = (/ 1., 1., 1., 1./)
  resx@cnFillColors = greens
  resx@gsnAddCyclic             = False
  greenMap = gsn_csm_contour(wks, Band2, resx)
  resx@cnFillColors = blues
  blueMap = gsn_csm_contour(wks, Band3, resx)

  mresP = res
  mresP@mpPerimOn                   = False             ; turn off box around plot(i)
  mresP@mpFillOn                    = False
  mresP@mpLimitMode                 = "angles" ;-- angles represent angular
  mresP@mpLeftAngleF                = 25.     ;-- deviation from the line
  mresP@mpRightAngleF               = 25.    ;-- of sight from the satellite
  mresP@mpTopAngleF                 = 25.      ;-- to the projection center
  mresP@mpBottomAngleF              = 25.
  mresP@mpDataBaseVersion           = "MediumRes"    ; better map outlines
  mresP@mpProjection                = "Satellite"    ;-- map projection
  mresP@mpSatelliteAngle1F          = 7.*57.2957795130823*asin(1./30.)/8. ;-- 1.671437 to 3.58404
  mresP@mpSatelliteAngle2F          = 85.   ;-- orthogonal view
  mresP@mpSatelliteDistF            = 2.32    ;-- n-times Earth distance
  mresP@mpEllipticalBoundary        = True
  mresP@mpCenterLonF                = 235.0       ; choose center lon
  mresP@mpCenterLatF                = 30.         ; choose center lat
  mresP@mpGridAndLimbOn             = True           ; turn on limb and grid lines
  mresP@mpGridLineColor             = -1             ; ...but don't draw grid lines...
  mresP@mpGridAndLimbDrawOrder      = "PreDraw"   ; draw lines first, so
  mresP@mpOutlineOn                 = True
  mresP@mpOutlineDrawOrder          = "PostDraw"      ; draw continental outline last
  mresP@mpOutlineBoundarySets       = "GeophysicalAndUSStates" ; state boundaries
  mresP@pmTickMarkDisplayMode       = "Always"; use NCL default lat/lon labels
  mresP@mpGeophysicalLineThicknessF = 5.0 
  mresP@mpUSStateLineThicknessF     = 5.0 
  mresP@mpNationalLineThicknessF    = 5.0 
  mresP@mpPerimLineColor            = -1
  mresP@mpPerimOn                   = False          ; turn off map perimetera
  mresP@mpLimbLineThicknessF        = 5.0
  mresP@mpOceanFillColor            = "skyblue1"
  mresP@mpLandFillColor             = "transparent"
  mresP@mpInlandWaterFillColor      = "skyblue1" 
  mresP@mpFillDrawOrder             = "PreDraw"

  mpres  = mresP
  mpres@cnLinesOn       = False              ; Turn off contour lines      .
  mpres@cnLineLabelsOn  = False              ; Turn off contour labels
  mpres@cnInfoLabelOn   = False              ; Turn off info label
  mpres@gsnTickMarksOn  = False              ; Turn off tickmarks
  mpres@lbLabelBarOn    = False              ; Turn off labelbar
  mpres@gsnRightString  = ""                 ; Turn off subtitles
  mpres@gsnLeftString   = ""
  mpres@cnFillOn        = True         ; Turn on filled rasters
  mpres@cnFillMode      = "RasterFill"
  mpres@cnLevelSelectionMode  = "EqualSpacedLevels"
  mpres@cnMaxLevelCount       = 254
  mpres@cnFillBackgroundColor = (/ 1., 1., 1., 1./)
  mpres@cnFillColors          = reds
  mpres@gsnAddCyclic          = False
  mpres@mpFillOn              = False
  mpres@gsnAddCyclic          = False
  plot(i) = gsn_csm_contour_map(wks, Band1, mpres) 
 ;---Overlay everything and draw
  overlay( plot(i) , greenMap)
  overlay( plot(i) , blueMap)

  ref_res                          = res
  ref_res@cnFillOn                 = True
  ref_res@cnLevelSelectionMode     = "ExplicitLevels"
  ref_res@cnFillMode               = "RasterFill"
  ref_res@cnRasterSmoothingOn      = True
  ref_res@cnOutOfRangeFillColor    = "white"
  ref_res@cnRasterSampleFactorF    = 5.0
  ref_res@cnRasterMinCellSizeF     = 1e-5
  ref_res@cnOutOfRangeFillColor    = "white"
  ref_res@cnInfoLabelOn            = False
  ref_res@cnLineLabelsOn           = False
  ref_res@cnLinesOn                = False
  ref_res@cnLineLabelBackgroundColor = -1    ; transparent
  ref_res@cnFillDrawOrder          = "PostDraw"
  ref_res@lbLabelBarOn             = True   ; Will turn on in panel later
  ref_res@pmLabelBarParallelPosF   = 0.50 ;61
  ref_res@pmLabelBarOrthogonalPosF = 0.05
  ref_res@lbOrientation            = "vertical"  ; vertical label bar
  ref_res@pmLabelBarHeightF        = 0.55
  ref_res@pmLabelBarWidthF         = 0.06
  ref_res@lbLabelFontHeightF       = FontHeightF*0.45
  ref_res@lbTitleFontHeightF       = FontHeightF*0.50
  ref_res@lbTitleString            = "Integrated Vapor Transoport (kg m~S~-1~N~ s~S~-1~N~)"
  ref_res@lbTitlePosition          = "Right"         ; title location
  ref_res@lbTitleDirection         = "Across"        ; letter angle
  ref_res@lbTitleAngleF            = 90.          ; title angle
  ref_res@lbBoxLinesOn             = True 
  ref_res@lbBoxSeparatorLinesOn    = False
  ref_res@cnLevelSelectionMode    := "ExplicitLevels" ;"ManualLevels"
  ref_res@cnLineThicknessF        := 1.0
  ref_res@gsnRightString          := "" 
  ref_res@gsnLeftString           := ""
  ref_res@sfXArray                      = ivt@lon
  ref_res@sfYArray                      = ivt@lat
  ref_res@gsnCenterString               = timstr
  ref_res@gsnCenterStringFontColor      = "Black"
  ref_res@gsnCenterStringFontHeightF    = FontHeightF*0.43
  ref_res@gsnCenterStringOrthogonalPosF = -0.07
  ref_res@gsnCenterStringParallelPosF   = 0.88
  ncols      = 101
  xnlevsy    = fspan(200,1200,ncols)
  xmapy      = read_colormap_file("MPL_BuPu") ;"MPL_BuPu";"MPL_YlGn"; "MPL_GnBu" ;"MPL_PuBuGn" ;cmap
  xmapz      = span_color_rgba(xmapy(0:,:),ncols)
  ref_res@cnFillPalette           := xmapz 
  ref_res@cnLevels                := xnlevsy
  ref_res@lbLabelStride            = 10
  ivt = where(ivt.le.min(cnlevs),ivt@_FillValue,ivt)
  contour(i) = gsn_csm_contour(wks,ivt,ref_res)

  res3    = ref_res
  ;---Calculate "nice" contour levels, and create a color map to match
  ncols      = 101
  cnlevsy    = fspan(0,20,101)
  cnlevsy(0) = 0.1
  levs       = ispan(1,11,1)
  cmapy      = read_colormap_file("BlAqGrYeOrReVi200") ;percent_11lev") ;("WhBlGrYeRe") ;BlAqGrYeOrReVi200")
  cmapz      = span_color_rgba(cmapy(8:,:),ncols)
  res3@cnLevelSelectionMode     = "ExplicitLevels"
  res3@cnLevels                := cnlevsy
  res3@cnFillPalette           := cmapz
  res3@cnFillMode               = "RasterFill"
  res3@cnRasterSmoothingOn      = True
  res3@cnOutOfRangeFillColor    = "white"
  res3@cnRasterSampleFactorF    = 5.0
  res3@cnRasterMinCellSizeF     = 1e-4
  res3@cnFillOpacityF          := 1.0     
  res3@lbLabelStride            = 10                 ; skip every other label
  res3@gsnRightString          := ""
  res3@gsnLeftString           := ""        
  res3@sfXArray                 = prc@lon   
  res3@sfYArray                 = prc@lat
  res3@lbLabelBarOn             = True   ; Will turn on in panel later
  res3@pmLabelBarParallelPosF   = 0.5
  res3@pmLabelBarOrthogonalPosF = 0.05
  res3@lbOrientation            = "horizontal"  
  res3@pmLabelBarHeightF        = 0.06
  res3@pmLabelBarWidthF         = 0.65
  res3@lbTitleString            = "Total Precipitation Rate (mm hr~S~-1~N~)"
  res3@lbTitlePosition          = "Top"         ; title location
  res3@lbTitleDirection         = "Across"        ; letter angle
  res3@lbTitleAngleF            = 0.0          ; title angle
  if(.not.isvar("contou2")) then 
    contou2 = new(100,graphic)
  end if 
  prc = where(prc.le.min(res3@cnLevels),prc@_FillValue,prc)
  contou2(i) = gsn_csm_contour(wks,prc,res3)
  
 ;set vector;;
  res_vc                            = res
  res_vc@vcGlyphStyle               = "CurlyVector"
  res_vc@vcLineArrowThicknessF      = 4
  res_vc@vcLineArrowHeadMaxSizeF    = 0.003
  res_vc@vcMinDistanceF             = 0.015
  res_vc@vcRefLengthF               = 0.01
  res_vc@vcRefAnnoOn                = False ;True
  res_vc@vcRefAnnoSide              = "Top"
  res_vc@vcRefAnnoString2On         = True
  res_vc@vcRefAnnoPerimOn           = True
  res_vc@vcRefAnnoOrthogonalPosF    = -0.40
  res_vc@vcRefAnnoParallelPosF      = 0.15
  res_vc@vcRefAnnoBackgroundColor   = "White";"Transparent"; "Purple"
  res_vc@vcVectorDrawOrder          = "PostDraw"
  res_vc@gsnRightString             = ""
  res_vc@vcRefAnnoFontHeightF       = FontHeightF*0.65
  res_vc@vcRefMagnitudeF            = 200
  res_vc@vcRefAnnoString1           = res_vc@vcRefMagnitudeF
  res_vc@vcRefAnnoString2           = "(kg/m/s)"
  tuq = var(0,:,:)
  tvq = var(1,:,:)
  wrf_smooth_2d(tuq,3)
  wrf_smooth_2d(tvq,3)
  iv1 = sqrt(tuq*tuq+tvq*tvq)
 ;tuq = where(iv1.lt.100.0,0,tuq)
 ;tvq = where(iv1.lt.100.0,0,tvq)
  if(.not.isvar("dum1")) then
    dum1 = new (100,graphic)
  end if
  res_vc@vcGlyphOpacityF     = 0.15
  res_vc@vcLineArrowColor    = "Grey40"           ; change vector color
  dum1(i) = gsn_csm_vector(wks,tuq,tvq,res_vc)
 ;vector(i) = gsn_csm_vector(wks,tuq,tvq,res_vc)
  res_vcx = res_vc
  tuq({minlat:maxlat},{minlon:maxlon}) = tuq@_FillValue
  tvq({minlat:maxlat},{minlon:maxlon}) = tvq@_FillValue
  res_vcx@vcGlyphOpacityF     = 0.5
  res_vcx@vcLineArrowColor    = "Grey60"           ; change vector color
  vector(i) = gsn_csm_vector(wks,tuq,tvq,res_vcx)
  overlay(vector(i),dum1(i))
  delete([/tuq,tvq,iv1/])

  ;-- draw edges
  plres = True
  plres@tmYLOn            = False
  plres@tmYROn            = False
  plres@tmXBOn            = False
  plres@tmXTOn            = False
  plres@gsLineColor       = "Grey50"
  plres@gsLineThicknessF  = 1e-6
  plres@gsSegments        = ispan(0,ngrid * ncorner,ncorner)
  plres@tfPolyDrawOrder   = "PostDraw"
  if(.not.isvar("dum3")) then
    dum3 = new (100,graphic)
  end if
  dum3(i) = gsn_add_polyline(wks, contour(i), ndtooned(grid_corner_lon), ndtooned(grid_corner_lat), plres)

 ;overlay(plot(i),contour(i))
 ;overlay(plot(i),contou2(i))
  overlay(plot(i),vector(i))

; plot regional data in high resolution 

;-- regional contour resources
 cnres2 = ref_res
 delete(cnres2@sfXArray)
 delete(cnres2@sfYArray)
 cnres2@trYReverse            = True   ;-- reverse y-axis
 cnres2@tfDoNDCOverlay        = False  ;-- transform to standard lat/lon
 cnres2@sfXArray              = rlon   ;-- longitude array
 cnres2@sfYArray              = rlat   ;-- latitude array
 cnres2@lbLabelBarOn          = False   ;-- don't draw a labelbar
 cnres2@cnFillDrawOrder       = "PostDraw"   ;-- draw last; default: "Draw"
 cnres2@cnMissingValFillColor = "transparent" ;res_vc@vcLineArrowColor  
 cnres2@cnFillOpacityF        = 1.0
 ncols      = 101
 fnlevsy    = fspan(50,1150,ncols) ;min(rivt),max(rivt),ncols)
 fmapy      = read_colormap_file("cmocean_gray") ;"MPL_BuPu";"MPL_YlGn"; "MPL_GnBu" ;"MPL_PuBuGn" ;cmap
 fmapz      = span_color_rgba(fmapy,ncols)
 do kk = 0,ncols-1,1
   fmapz(kk,3) =  fmapz(kk,3) * (kk/ncols)
 end do 
 cnres2@cnFillPalette           := fmapz
 cnres2@cnLevels                := fnlevsy
 cnres2@lbLabelStride            = 10
 if(.not.isvar("contou5")) then
    contou5 = new (100,graphic)
 end if
;rivt       = where(rivt.le.100,rivt@_FillValue,rivt)
;rivt       = where(rprc.le.0.01.and.rivt.le.50,rivt@_FillValue,rivt)
 contou5(i) = gsn_csm_contour(wks,rivt,cnres2)   ;-- create Europe plot

 cnres3 = cnres2
 ;---Calculate "nice" contour levels, and create a color map to match
 ncols      = 101
 cnlevsy    = fspan(0,20,101)
 cnlevsy(0) = 0.1
 levs       = ispan(1,11,1)
 cmapy      = read_colormap_file("BlAqGrYeOrReVi200") ;percent_11lev") ;("WhBlGrYeRe") ;BlAqGrYeOrReVi200")
 cmapz      = span_color_rgba(cmapy(8:,:),ncols)
; do kk = 0,ncols-1,1
;   cmapz(kk,3) =  cmapz(kk,3) * (kk/ncols)
; end do
 cnres3@cnLevelSelectionMode     = "ExplicitLevels"
 cnres3@cnLevels                := cnlevsy
 cnres3@cnFillPalette           := cmapz
 cnres3@cnFillMode               = "RasterFill"
 cnres3@cnRasterSmoothingOn      = True
;cnres3@cnOutOfRangeFillColor    = "white"
 cnres3@cnRasterSampleFactorF    = 5.0
 cnres3@cnRasterMinCellSizeF     = 1e-4
 cnres3@cnFillOpacityF          := 1.0
 cnres3@lbLabelStride           := 10                 ; skip every other label
 cnres3@gsnRightString          := ""
 cnres3@gsnLeftString           := ""
 cnres3@sfXArray                 = prc@lon
 cnres3@sfYArray                 = prc@lat
 cnres3@lbLabelBarOn             = True   ; Will turn on in panel later
;cnres3@pmLabelBarParallelPosF   = 0.5
;cnres3@pmLabelBarOrthogonalPosF = 0.05
;cnres3@lbOrientation            = "horizontal"
;cnres3@pmLabelBarHeightF        = 0.06
;cnres3@pmLabelBarWidthF         = 0.65
 cnres3@lbTitleString            = "Total Precipitation Rate (mm hr~S~-1~N~)"
;cnres3@lbTitlePosition          = "Top"         ; title location
;cnres3@lbTitleDirection         = "Across"        ; letter angle
;cnres3@lbTitleAngleF            = 0.0          ; title angle
 if(.not.isvar("contou6")) then
   contou6 = new(100,graphic)
 end if
 rprc = where(rprc.le.0.05,rprc@_FillValue,rprc)
 contou6(i) = gsn_csm_contour(wks,rprc,cnres3)

;-- edges for the shadow box and the borderline of plot2a
 nlat = dimsizes(rlat2d(:,0))  ;-- number of latitudes
 nlon = dimsizes(rlon2d(0,:))  ;-- number of longitudes
 lon_val_upper = rlon2d(nlat-1,:)
 lat_val_upper = rlat2d(nlat-1,:)
 lon_val_lower = rlon2d(0,:)
 lat_val_lower = rlat2d(0,:)
 lon_val_left  = rlon2d(:,0)
 lat_val_left  = rlat2d(:,0)
 lon_val_right = rlon2d(:,nlon-1)
 lat_val_right = rlat2d(:,nlon-1)
;-- draw a "shadow" rectangle around the regional map for better visibility
 dx = 1.5   ;-- x-offset of shadow rectangle
 dy = 1.0   ;-- y-offset of shadow rectangle
;-- x-array for shadow rectangle
 x = array_append_record((/lon_val_lower+dx/), (/lon_val_right+dx/), 0)
 x := array_append_record((/x/), (/lon_val_upper(::-1)+dx/), 0)
 x := array_append_record((/x/), (/lon_val_left(::-1)+dx/), 0)
 x := array_append_record((/x/), (/rlon2d(0,0)+dx/), 0)

;-- y-array for shadow rectangle
 y = array_append_record((/lat_val_lower-dy/), (/lat_val_right-dy/), 0)
 y := array_append_record((/y/), (/lat_val_upper(::-1)-dy/), 0)
 y := array_append_record((/y/), (/lat_val_left(::-1)-dy/), 0)
 y := array_append_record((/y/), (/rlat2d(0,0)-dy/), 0)

;-- draw the transparent shadow rectangle (filled polygon)
 pgres = True
 pgres@tfPolyDrawOrder = "PostDraw"  ;-- draw the rectangle last
 pgres@gsLineColor     = "grey40" ; grey30"   ;-- rectangle line color
 pgres@gsFillColor     = "grey40" ; "grey30"   ;-- rectangle fill color
 pgres@gsFillOpacityF  = 0.3   ;-- light transparency

 shadow = gsn_add_polygon(wks, plot(i), x, y, pgres)   ;-- attach shadow rectangle to map
 print("---> regional plot2 shadow done")
 delete([/x,y/])

 overlay(plot(i),contou5(i))
 overlay(plot(i),contou6(i))

;-- create a simple global map with geophysical outlines only
 mpres2 = mresP
 mpres2@mpGridAndLimbOn = False   ;-- don't draw map grid lines
 mpres2@mpFillOn        = False   ;-- don't fill map
 map2 = gsn_csm_map(wks,mpres2) ;-- create simple map
 print("---> simple map2 outlines done")

;-- attach map2 to map to draw the geophysical outlines ontop of the regional map, too.
 anres = True
 anid = gsn_add_annotation(plot(i),map2,anres)   ;-- attach map2 to map
 print("---> add map2 annotation done")

;-- draw the edges of the regional map in black
 plres = True
 plres@tfPolyDrawOrder  = "PostDraw"  ;-- draw line last
 plres@gsLineThicknessF = 3.0 ;1.5   ;-- line thickness
 plres@gsLineColor     := "grey50"   ;-- line color

 upper  = gsn_add_polyline(wks, map2, lon_val_upper, lat_val_upper-0.5, plres)   ;-- attach upper line
 lower  = gsn_add_polyline(wks, map2, lon_val_lower, lat_val_lower, plres)   ;-- attach lower line
 left   = gsn_add_polyline(wks, map2, lon_val_left, lat_val_left, plres)      ;-- attach left line
 right  = gsn_add_polyline(wks, map2, lon_val_right, lat_val_right, plres)   ;-- attach right line
 print("---> black rectangle/edge done")

 ;---Draw a box using gsn_add_polygon. See how the edges are curved.
 if (plot_oroville_dam) then
   pmres                        = True
   pmres@gsMarkerColor          = "Red"       ;-- marker color
   pmres@gsMarkerSizeF          = 0.012        ;-- set size of marker
   pmres@gsMarkerThicknessF     = 8.          ;-- marker line thickness
   pmres@gsMarkerIndex          = 12
   pmres@tfPolyDrawOrder = "PostDraw"
   if(.not.isvar("dum2")) then
     dum2 = new (100,graphic)
   end if
   dum2(i) = gsn_add_polymarker(wks,contou6(i),lon_od,lat_od,pmres) ; Draw the polymarkers.
 end if

 resP                            = True
 resP@gsnMaximize                = True
 resP@gsnPaperOrientation        = "portrait"
 resP@gsnPanelLabelBar           = False 
 resP@gsnPanelBottom             = 0.1
 resP@gsnPanelYWhiteSpacePercent = 0.
 resP@gsnPanelXWhiteSpacePercent = 3.
 resP@gsnFrame                   = True
 gsn_panel(wks,plot,(/1,1/),resP)

 delete([/latx,lonx,rlat2d,rlon2d,vdum2d/])

end do 
end do

end
