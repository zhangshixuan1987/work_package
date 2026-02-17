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
  Casdir  = (/"ERA5","CLIM30"/)
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

  do i = 1,ngrps-1,1
   obstr    = Casdir(0)
   obsdir   = datadir+ "/" + obstr + "/" + freq+ "/"
   flo      = systemfunc("ls "+obsdir+obstr+"*.monthly.*.nc")
   fo       = addfiles(flo, "r")
   tim0     = cd_calendar(fo[:]->time,-1)
   ind0     = get1Dindex(tim0,yymm) 
   lat      = fo[0]->lat
   lon      = fo[0]->lon
   p0       = 100000.
   plev     = fo[0]->lev
   plev     = plev*100.                                ; match units of p0 and psfc
   plev@units = "Pa"
   
   model    = Casdir(i)
   modeldir = datadir+ "/" + model + "/" + freq+ "/"
   flm      = systemfunc("ls "+modeldir+model+"*.monthly.*.nc")
   fm       = addfiles(flm, "r")
   tim1     = cd_calendar(fm[:]->time,-1)
   ind1     = get1Dindex(tim1,yymm) ; index for each desired data

   fc       = addfile(fcoord,"r")
    
   do jj = 0,nvars-2, 1

    Varname  = VarList(jj)
    Varunt   = UnitList(jj)
    Vnmstr   = VnmList(jj)

   ;interpolate obs to model level
    x10     = pres2hybrid_Wrap(plev,fo[:]->PS(ind0,:,:),p0, \
                               fo[:]->$Varname$(ind0,:,:,:),fc->hyam,fc->hybm,0)
    y10     = fm[:]->$Varname$(ind1,:,:,:)

    x11     = month_to_annual(x10,1)
    y11     = month_to_annual(y10,1)
    x11!0   = "time"
    y11!0   = "time"

    x11@_FillValue = -99999
    y11@_FillValue = -99999

    x11 = where(ismissing(y11),x11@_FillValue,x11)
    y11 = where(ismissing(x11),y11@_FillValue,y11)

    ;calculate the global/regional mean values
    latrg    = lat
    rad      = 4.*atan(1.)/180.    ;
    cost     = cos(latrg*rad)       ; cosine weights
 
    if(Varname.eq."RELHUM") then
     x11 = x11/100.0
     y11 = y11/100.0
    end if

    if(Varname.eq."Q") then
     x11 = x11*1000.0
     y11 = y11*1000.0
    end if

    votmp          = dim_avg_n_Wrap(x11,3)
    vmtmp          = dim_avg_n_Wrap(y11,3)

    vmtmp@long_name = VnmList(jj)
    vmtmp@units     = UnitList(jj)
    votmp@long_name = VnmList(jj)
    votmp@units     = UnitList(jj)

    ;calculate time average
    vmmod = dim_avg_n_Wrap(vmtmp,0)
    vsmod = dim_variance_n_Wrap(vmtmp,0)
    vmobs = dim_avg_n_Wrap(votmp,0)
    vsobs = dim_variance_n_Wrap(votmp,0)
    vsdmd = vmmod 
    vsdmd = (/sqrt(vsmod)/)
    
    sigr  = 0.05                        ; critical sig lvl for r
    xEqv  = equiv_sample_size (votmp(lev|:,lat|:,time|:), sigr,0)
    yEqv  = equiv_sample_size (vmtmp(lev|:,lat|:,time|:), sigr,0)
    xN    = wgt_areaave (xEqv, 1.0, cost,  0)    ; wgty could be gaussian weights 
    yN    = wgt_areaave (yEqv, 1.0, cost,  0) 

    iflag = False                        ; population variance similar
    prob  = ttest(vmobs,vsobs,xN,vmmod,vsmod,yN,iflag,False)
    prob  = where(prob.gt.sigr,prob@_FillValue,prob)
    copy_VarCoords(vmmod,prob)

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
    fout->prob  = prob

    fout->lev   = fc->lev
    fout->hybn  = fc->hybm
    fout->hyam  = fc->hyam
    fout->lat   = lat

    delete([/vmmod,vsmod,vsobs,vmobs,vsdmd/])
    delete([/x11,vmtmp,y11,votmp,x10,y10/])

   end do

  end do 

end
