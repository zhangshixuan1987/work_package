;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  datadir  = "/anvil/scratch/x-szhang3/data_240x120"
  outdir   = "../fig_data/"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5_MLEV"/)
  Groups  = Casdir
  ngrps   = dimsizes(Groups)

  freq    = "monthly"
  yst     = 2007
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
   print(flm)
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
   lat@units = "degree_north"
   lon@units = "degree_east"
   lat!0    = "lat"
   lon!0    = "lon"
   lat&lat  = lat
   lon&lon  = lon
   lev      = fm[0]->lev

   do jj = 0,nvars-2, 1

    Varname  = VarList(jj)
    Varunt   = UnitList(jj)
    Vnmstr   = VnmList(jj)

    ;calculate the global/regional mean values
    latrg    = lat
    rad      = 4.*atan(1.)/180.    ;
    cost     = cos(latrg*rad)       ; cosine weights
 
    x11      = fm[:]->$Varname$(indt,:,:,:)
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

    printVarSummary(x11)
    vtmp           = dim_avg_n_Wrap(x11,3)
    vtmp@long_name = VnmList(jj)
    vtmp@units     = UnitList(jj)

    ;calculate time average
    vmmod = dim_avg_n_Wrap(vtmp,0)
    vsmod = dim_variance_n_Wrap(vtmp,0)

   ;open netcdf file to save the data for plots;;;;;
    setfileoption("nc", "Format",  "NetCDF4")
    out_file_name = outdir+model+"_"+Varname+"_"+yst+"-"+yed+".nc"
    system("rm " + out_file_name)
    fout = addfile(out_file_name,"c")
    vo1  = Varname 
    vo2  = "variance"
    fout->lev= lev
    fout->lat= lat
    fout->$vo1$=vmmod
    fout->$vo2$=vsmod

    vsmo1 = dim_avg_n_Wrap(dim_stddev_n_Wrap(x11,0),2)
    fout->stddev = vsmo1
    delete([/vmmod,vsmod/])
    delete([/x11,vtmp/])

   end do

  end do 

end
