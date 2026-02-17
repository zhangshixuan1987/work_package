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
  Casdir   = (/"CLIM30","CLIM30ML"/)
  Groups   = Casdir
  ngrps    = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U","V","T","Q","RELHUM"/)
  VnmList  = (/"Zonal wind","Meridional wind", "Temperature", \
               "Specific humidity", "Relative humidity"/)
  UnitList = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~","unitless"/)
  nvars    = dimsizes(VarList)

  ;open the ERA5 data file;;;
  flo      = systemfunc("ls " + datadir + "/ERA5.monthly/ERA5*analysis*.nc")
  fo       = addfile(flo,"r")
  udate    = cd_calendar(fo->time,0)
  yyyy     = udate(:,0)
  lev      = fo->lev 

  do ii = 0,ngrps-1,1

   do year  = ystr,yend,1
    
    ;open the E3SM data file;;;
    fl0     = systemfunc("ls " + datadir + "/"+Casdir(ii)+"/" + freq + "/*" +year+"*.nc")
    f0      = addfile(fl0,"r")
    indo    = ind(yyyy.eq.year)
   
    do jj = 0,nvars-1, 1

     Varname  = VarList(jj)
     Varunt   = UnitList(jj)
     Vnmstr   = VnmList(jj)

     lat      = f0->lat
     lon      = f0->lon
     tbot     = f0->TBOT
     phis     = conform(tbot,f0->PHIS,(/1,2/))

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

     if(Varname.eq."RELHUM") then
      pres = f0->T 
      do n=0,dimsizes(f0->hyam)-1
        pres(:,n,:,:) = tofloat(f0->P0*f0->hyam(n) + f0->hybm(n)*f0->PS)
      end do
      vtmp = f0->T
      vtmp = (/relhum (f0->T, f0->Q, pres)/)
      vtmp@long_name="Relative humidity"
      vtmp@units="%"
      delete([/pres/])
     else
      vtmp = f0->$Varname$   
     end if 
     x11  = vinth2p_ecmwf(vtmp,f0->hyam,f0->hybm,lev,f0->PS,intyp,P0mb, \
                           1,kxtrp,varflg,tbot,phis)
     x11&lev@units = "hPa"
     y11  = fo->$Varname$(indo,:,:,:)

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
     fnx    = outdir+ VarList(jj) + "_" + Groups(ii) + "_"+year+".nc"
     system("rm -f " + fnx)

     print(" ")
     print("output : " + fnx)
     print(" ")

     flx = addfile(fnx,"c")
     vna = VarList(jj)
     flx->$vna$=vmod
     vnb = "d"+VarList(jj)
     flx->$vnb$=verr

     delete([/flx,x11,y11,vmod,verr/])

    end do
   end do 
  end do

end

