;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  datadir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/OBS/tc_track"
  outdir  = "./"
  grid    = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output/landmsk_coord/landmsk_coord.nc"

 ;;; flags to control the plot style
  plot_contour      = False ; if True, plot the contour in the figure
  plot_leftstr      = True  ; if True, the lefadd igure captions as left string
  l_plot_latlon_box = False ; if True, add lat-lon box for SEP region

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups  = (/"OBS"/)
  Casdir  = (/"IBTrACS"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  ystr    = 1979
  yend    = 2017
  timtag  = ystr+"-"+yend 
  nyr     = yend - ystr + 1
  nmon    = 12  

  ;construct 4x4 global grid to bin storms 
  nlon    = toint(360./4)
  nlat    = toint(180./4)
  latm    = latGlobeFo(nlat,"lat","latitude","degrees_north")
  lonm    = lonGlobeFo(nlon,"lon","longitude","degrees_east")

  ;;create array to save the track counts 
  trkdens      = new((/nyr,nmon,nlat,nlon/),float)
  trkdens!0    = "year"
  trkdens!1    = "month"
  trkdens!2    = "lat"
  trkdens!3    = "lon"
  trkdens&year = ispan(ystr,yend,1)
  trkdens&month= ispan(1,nmon,1)
  trkdens&lat  = latm
  trkdens&lon  = lonm
  tgedens      = trkdens
  tlydens      = trkdens
  tcltime      = trkdens
  tcints       = trkdens
  tcspd        = trkdens

  trkdens@long_name  = "track density (number density of storm tracks passing through a region)" 
  trkdens@short_name = "tdens" 
  trkdens@units      = "transits per unit area equivalent to a 4-degree spherical cap"

  tgedens@long_name  = "genesis density (number density of starting points)"
  tgedens@short_name = "gdens"
  tgedens@units      = "transits per unit area equivalent to a 4-degree spherical cap"

  tlydens@long_name  = "lysis density (number density of ending points)"
  tlydens@short_name = "ldens"
  tlydens@units      = "transits per unit area equivalent to a 4-degree spherical cap"

  tcints@long_name   = "mean intensity"
  tcints@short_name  = "int"
  tcint@units        = "m/s"

  tcltime@long_name  = "mean lifetime"
  tcltime@short_name = "lft"
  tcltime@units      = "days"

  tcspd@long_name    = "mean track speed"
  tcspd@short_name   = "mspd"
  tcspd@units        = "m/s"

  do i = 0,ngrps-1,1

   fna  = systemfunc("ls -1 " + datadir+"/"+Casdir(i)+"/"+Groups(i)+"*_"+sset+"_TCS_Track_"+ timtag+ ".txt") 
   fla  = asciiread(fna,-1,"string")
   inds = str_match_ind_ic_regex(fla,"start")
   nst  = num(.not.ismissing(inds)) 

   do j = 0,nst-1,1
    
    ih      = inds(j)
    head    = fla(ih)
    delim   = "	" ;note: tab space 
    ntrk    = stringtointeger(str_get_field(head,2, delim))
    print("working on storm: " + head)

    is      = ih+1
    ie      = ih+ntrk
    hsub    = fla(is:ie)
    ilon    = stringtointeger(str_get_field(hsub,1, delim)) 
    ilat    = stringtointeger(str_get_field(hsub,2, delim))  
    slon    = stringtointeger(str_get_field(hsub,3, delim))
    slat    = stringtointeger(str_get_field(hsub,4, delim))

    is    = is + 1   
    ie    = ntlin + ntrk -1
    hsub  = fla(ns1:ne1)
    print(hsub)
  
    ;trkdens() = trkdens({})+ 
    print(head)
    exit

    ntrk  = stringtointeger(str_get_field(head, 2, "	"))
    ns1   = ntlin
    ne1   = ntlin + ntrk -1
    hsub  = fla(ns1:ne1)
    ntlin = ntlin + ntrk
    lon(ix,0:ntrk-1) = stringtofloat(str_get_field(hsub, 2,"	"))
    lat(ix,0:ntrk-1) = stringtofloat(str_get_field(hsub, 3,"	"))
    vsmc(ix,0:ntrk-1)= stringtofloat(str_get_field(hsub, 5,"	"))
    vpmc(ix,0:ntrk-1)= stringtofloat(str_get_field(hsub, 4,"	"))
    print(hsub)
    yyyy = stringtointeger(str_get_field(hsub, 7,"	"))
    mm   = stringtointeger(str_get_field(hsub, 8,"	"))
    dd   = stringtointeger(str_get_field(hsub, 9,"	"))
    hh   = stringtointeger(str_get_field(hsub, 10,"	"))
    delete([/ns1,ne1,ntrk,head,hsub/])
    printVarSummary(yyyy)
    print(yyyy)
   end do 

   delete([/fna,fla,inds/])
   printVarSummary(yyyy)
   exit
  ;loop each storm and calculate the AEI 
   
  ;loop each year and mont to obtain the binned storm density 
   lat1  = lat(ix,0:ntrk-1)
   lon1  = lon(ix,0:ntrk-1)
   lon1  = where(lon1.lt.0.0,lon1+360.0,lon1)
   nm = getind_latlon2d (lat2d,lon2d, lat1, lon1)
   tdens(nm(:,0),nm(:,1)) = tdens(nm(:,0),nm(:,1)) + 1

   print ("total number of TCs in " + Groups(i) + " is " + nst)
   tdens!0 = "lat"
   tdens!1 = "lon"
   tdens&lat = latgrd
   tdens&lon = longrd
   tdens@long_name = "track density (number of sttorms passed the 1x1 grid box)"

  ;bining storms into 4x4 box
   nlat    = toint(360./4)
   nlon    = toint(180./4)
   bin_lat = latGlobeFo(nlat,"lat","latitude","degrees_north")
   bin_lon = lonGlobeFo(nlon,"lon","longitude","degrees_east")

   ;;;;open netcdf file and save data for the pfigures;;;;;;;
   setfileoption("nc", "Format",  "NetCDF4")
   out_file_nam1 = outdir+Casdir(i)+"_"+Groups(i)+"_Track_info_"+ timtag+ ".nc"
   out_file_nam2 = outdir+Casdir(i)+"_"+Groups(i)+"_Track_dens_"+ timtag+ ".nc"
   system("rm " + out_file_nam1)
   system("rm " + out_file_nam2)
   fout1 = addfile(out_file_nam1,"cw")
   fout2 = addfile(out_file_nam2,"cw")
   fout1->lat  = lat
   fout1->lon  = lon
   fout1->stim = stim 
   fout1->vsmc = vsmc 
   fout1->vpmc = vpmc

   fout2->lat      = latgrd
   fout2->lon      = longrd
   fout2->density  = tdens  

   delete([/lat,lon,stim,vsmc,vpmc,tdens/])

 end do 


end
