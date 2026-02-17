;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  datadir  = "/anvil/scratch/x-szhang3/data_240x120"
  outdir   = "../fig_data/"
  fcoord   = "/anvil/projects/x-phy220062/x-szhang3/map_file/coord.nc"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5"/)
  Groups  = Casdir
  ngrps   = dimsizes(Groups)

  freq    = "monthly"
  yst     = 1979
  yed     = 2014
  yymm    = yyyymm_time(yst,yed, "integer")

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U","V","T","Q","RELHUM"/)
  VnmList  = (/"Zonal wind","Meridional wind", "Temperature", \
               "Specific humidity", "Relative humidity"/)
  UnitList = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~","unitless"/)
  nvars    = dimsizes(VarList)

  do i = 0,ngrps-1,1
   model    = Casdir(i)
   modeldir = datadir+ "/" + model + "/" + freq+ "/"
   flm      = systemfunc("ls "+modeldir+model+"*.monthly.*.nc")
   fm       = addfiles(flm, "r")
   tim0     = cd_calendar(fm[:]->time,-1)
   indt     = get1Dindex(tim0,yymm) ; index for each desired data
   if(dimsizes(dimsizes(fm[0]->lat)).gt.1) then
     lat    = fm[0]->lat(:,0)
     lon    = fm[0]->lon(0,:)
   else
     lat    = fm[0]->lat
     lon    = fm[0]->lon
   end if
   lat@units  = "degree_north"
   lon@units  = "degree_east"
   lat!0      = "lat"
   lon!0      = "lon"
   lat&lat    = lat
   lon&lon    = lon
   
   fc         = addfile(fcoord,"r")
   hya        = fc->hyam                                 ; levels to which the variable
   hyb        = fc->hybm     
   lev        = fc->lev 

   p0         = 100000.
   plev       = fm[0]->lev
   plev       = plev*100.                                ; match units of p0 and psfc
   plev@units = "Pa"
   ps         = fm[:]->PS(indt,:,:)
   ps@_FillValue = -99999 

   ;calculate the global/regional mean values
   latrg    = lat
   rad      = 4.*atan(1.)/180.    ;
   cost     = cos(latrg*rad)       ; cosine weights

   do jj = 0,nvars-2, 1

    Varname  = VarList(jj)
    Varunt   = UnitList(jj)
    Vnmstr   = VnmList(jj)
    vtmp     = fm[:]->$Varname$(indt,:,:,:)
    vtmp@_FillValue = -99999. 

    x11      = pres2hybrid(plev,ps,p0,vtmp,hya,hyb,0)
    x11!0    = "time"
    x11!1    = "lev"
    x11!2    = "lat"
    x11!3    = "lon"
    x11&lat  = lat
    x11&lon  = lon
    x11&time = tim0(indt)
    x11&lev  = lev
    x11&lev@units = "hPa"

    if(Varname.eq."RELHUM") then
     x11 = x11/100.0
    end if

    if(Varname.eq."Q") then
     x11 = x11*1000.0
    end if

    y11     = month_to_annual(x11,1)
    y11!0   = "time"

    z11           = dim_avg_n_Wrap(y11,3)
    z11@long_name = VnmList(jj)
    z11@units     = UnitList(jj)

    ;calculate time average
    vmmod = dim_avg_n_Wrap(z11,0)
    vsmod = dim_variance_n_Wrap(z11,0)
    vsdmd = vmmod
    vsdmd = (/sqrt(vsmod)/)
    
   ;open netcdf file to save the data for plots;;;;;
    setfileoption("nc", "Format",  "NetCDF4")
    out_file_name = outdir+model+"_"+Varname+"_"+yst+"-"+yed+".nc"
    system("rm " + out_file_name)
    fout = addfile(out_file_name,"c")
    vo1  = Varname
    vo2  = "variance"
    vo3  = "stddev"
    fout->$vo1$ = vmmod
    fout->$vo2$ = vsmod
    fout->$vo3$ = vsdmd
    fout->lev   = lev
    fout->lat   = lat
    fout->hyam  = hya
    fout->hybm  = hyb 
    delete([/vmmod,vsmod,vsdmd/])
    delete([/x11,y11,z11,vtmp/])

   end do

  end do 

end
