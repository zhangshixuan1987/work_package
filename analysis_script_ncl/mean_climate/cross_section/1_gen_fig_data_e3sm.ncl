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
  outdir   = "../fig_data/"
  datadir  = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"
  freq     = "monthly"
  ystr     = 200801
  yend     = 201712

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = Casdir
  ngrps    = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U","V","T","Q","RELHUM","CLOUD"/)
  VnmList  = (/"Zonal wind","Meridional wind", "Temperature", \
               "Specific humidity", "Relative humidity", "Cloud fraction"/)
  UnitList = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~","unitless","fraction"/)
  nvars    = dimsizes(VarList)

  ;open the ERA5 data file;;;
  obsdir   = "/global/cfs/cdirs/e3sm/zhan391/data/"
  flo      =  systemfunc("ls " + obsdir + "/ERA5/" + freq + "/ERA5*analysis*.nc")
  fo       = addfile(flo,"r")
  lev      = fo->lev
  utcdate  = cd_calendar(fo->time,-1)
  indobs   = ind(utcdate.ge.ystr.and.utcdate.le.yend)
  delete([/utcdate/])

  do ii = 0,ngrps-1,1 

   ;open the ERA5 data file;;;
   fl0      = systemfunc("ls " + datadir + "/"+Casdir(ii)+"/" + freq + "/" +"*.nc")
   f0       = addfiles(fl0,"r")
   time     = f0[:]->time
   time     = time - 1
   utcdate  = cd_calendar(time,-1)
   indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)

   do jj = 0,nvars-1, 1

    Varname  = VarList(jj)
    Varunt   = UnitList(jj)
    Vnmstr   = VnmList(jj)

    lat    = f0[0]->lat
    lon    = f0[0]->lon

    ;calculate the global/regional mean values
    latrg  = lat
    rad    = 4.*atan(1.)/180.    ;
    cost   = cos(latrg*rad)       ; cosine weights

    P0mb   = 1000.0
    intyp  = 1
    kxtrp  = False 

    if(Varname.eq."T") then 
      varflg = 1  
    else if (Varname.eq."Z3") then
      varflg = -1
    else 
      varflg = 0
    end if 
    end if 
    xtmp = vinth2p_ecmwf(f0[:]->$Varname$,f0[0]->hyam,f0[0]->hybm,lev,f0[:]->PS,intyp,P0mb, \
                           1,kxtrp,varflg,f0[:]->TS,f0[:]->PHIS)
 
    x11  = xtmp(indtim,:,:,:)
    x11&lev@units = "hPa"
    y11  = fo->$Varname$(indobs,:,:,:)

    if(Varname.eq."RELHUM") then
     x11 = x11/100.0
     y11 = y11/100.0
    end if

    if(Varname.eq."Q") then
     x11 = x11*1000.0
     y11 = y11*1000.0
    end if

    y11  = where(ismissing(x11),y11@_FillValue,y11)
    x11  = where(ismissing(y11),x11@_FillValue,x11)
    xdif = x11
    xdif = (/x11 - y11/)
 
    vmod = dim_avg_n_Wrap(x11,3)
    vmod@long_name = VnmList(jj)
    vmod@units     = UnitList(jj)

    verr = dim_avg_n_Wrap(xdif,3)
    verr@long_name = VnmList(jj) + " bias"
    verr@units     = UnitList(jj)

    system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
    fnx    = outdir+ VarList(jj) + "_" + Groups(ii) + ".nc"
    system("rm -f " + fnx)

    print(" ")
    print("output : " + fnx)
    print(" ")

    flx = addfile(fnx,"c")
    vna = VarList(jj) 
    flx->$vna$=vmod
    vnb = "d"+VarList(jj)
    flx->$vnb$=verr
    
    delete([/flx,x11,y11,vmod,verr,xdif/])

   end do 
  end do 
 
end
