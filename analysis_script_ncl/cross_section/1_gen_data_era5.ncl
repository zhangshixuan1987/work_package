;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
  outdir   = "../fig_data/"
  datadir  = "/anvil/scratch/x-szhang3/data"
  freq     = "monthly"
  ystr     = 1979
  yend     = 2014

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5"/)
  Groups  = Casdir
  ngrps   = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U","V","T","Q","RELHUM"/)
  VnmList  = (/"Zonal wind","Meridional wind", "Temperature", \
               "Specific humidity", "Relative humidity"/)
  UnitList = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~","unitless"/)
  nvars    = dimsizes(VarList)

  ;open the ERA5 data file;;;
  fl0      =  systemfunc("ls " + datadir + "/"+Casdir+"." + freq + "/" + Casdir+"*analysis*.nc")
  f0       = addfile(fl0,"r")
  utcdate  = cd_calendar(f0->time,0)

  do ii = 0,ngrps-1,1

  do year=ystr,yend,1

   indtim   = ind (utcdate(:,0).eq.year)

   do jj = 0,nvars-1, 1

    Varname  = VarList(jj)
    Varunt   = UnitList(jj)
    Vnmstr   = VnmList(jj)

    lat    = f0->lat
    lon    = f0->lon
    lev    = f0->lev

    ;calculate the global/regional mean values
    latrg  = lat
    rad    = 4.*atan(1.)/180.    ;
    cost   = cos(latrg*rad)       ; cosine weights

    xtmp   = f0->$Varname$
    x11    = xtmp(indtim,:,:,:)

    x11&lev@units = "hPa"

    if(Varname.eq."RELHUM") then
     x11 = x11/100.0
    end if

    if(Varname.eq."Q") then
     x11 = x11*1000.0
    end if

    printVarSummary(x11)
    pncre0300 = dim_avg_n_Wrap(x11,3)
    pncre0300@long_name = VnmList(jj)
    pncre0300@units     = UnitList(jj)

    system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
    fno  = outdir+ VarList(jj) + "_" + Groups(ii) + "_"+year+".nc"
    system("rm -f " + fno)
    print(" ")
    print("output : " + fno)
    print(" ")

    flo = addfile(fno,"c")

    vna = VarList(jj)
    flo->$vna$=pncre0300

    delete(flo)
    delete([/x11,pncre0300/])

   end do

  end do 

  end do

end
