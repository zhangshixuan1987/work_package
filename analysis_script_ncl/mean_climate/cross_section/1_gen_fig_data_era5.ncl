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
  outdir  = "../fig_data/"
  datadir = "/global/cfs/cdirs/e3sm/zhan391/data/"
  freq    = "monthly"
  ystr    = 200801
  yend    = 201712

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5"/)
  Groups  = Casdir
  ngrps   = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U","V","T","Q","RELHUM","CLOUD"/)
  VnmList  = (/"Zonal wind","Meridional wind", "Temperature", \
               "Specific humidity", "Relative humidity", "Cloud fraction"/)
  UnitList = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~","unitless","fraction"/)
  nvars    = dimsizes(VarList)

  ;open the ERA5 data file;;;
  fl0      =  systemfunc("ls " + datadir + "/"+Casdir+"/" + freq + "/" + Casdir+"*analysis*.nc")
  f0       = addfile(fl0,"r")
  utcdate  = cd_calendar(f0->time,-1)
  indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)

  do ii = 0,ngrps-1,1 

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
    fno    = outdir+ VarList(jj) + "_" + Groups(ii) + ".nc"
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
 
end
