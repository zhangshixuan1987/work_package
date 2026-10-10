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
  obsdir   = "/global/cfs/cdirs/e3sm/zhan391/data"
  datadir  = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"
  freq     = "monthly"
  ystr     = 200801
  yend     = 201712

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"ERA5","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = Casdir
  ngrps    = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList   = (/"CLDTOT", "PSL",   "U10",   "TS",   "TREFHT", "LHFLX", "SHFLX", \
                "STRESS_MAG", "TGCLDLWP", "TGCLDIWP","PRECL","PRECC", \
                "U",           "U",       "U",        "Z3"/) 
  VstList   = (/"CLDTOT", "PSL",   "U10",   "TS",   "TREFHT", "LHFLX", "SHFLX", \
                "STRESS_MAG", "LWP", "IWP","PRECL","PRECC", \
                "U200", "U500", "U850", "Z500"/) 
  nvars     = dimsizes(VarList)

  do ii = 0, ngrps-1,1

   do jj = 0,0 ;nvars-1, 1

   Varname = VarList(jj)
   Varstr  = VstList(jj)

   if(ii.eq.0) then 
     fl0      =  systemfunc("ls " + obsdir + "/"+Casdir(ii)+"/" + freq + "/" + Casdir(ii)+"*analysis*.nc")
     f0       = addfile(fl0,"r")
     time     = f0->time
     utcdate  = cd_calendar(time,-1)
     indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)
     
     if (Varstr.eq."STRESS_MAG") then
      var     = f0->TAUX(indtim,:,:)
      var     = sqrt(f0->TAUX(indtim,:,:)*f0->TAUX(indtim,:,:) + \
                     f0->TAUY(indtim,:,:)*f0->TAUY(indtim,:,:))
     else if (Varstr.eq."U200".or.Varstr.eq."U500".or.Varstr.eq."U850".or.Varstr.eq."Z500") then 
      if(Varname.eq."Z3") then 
        vnm = "Z"
      else 
        vnm = Varname
      end if 
      plev    = tofloat(str_sub_str(Varstr,vnm,""))
      indlev  = ind(f0->lev.eq.plev)
      var     = f0->$Varname$(indtim,indlev,:,:) 
      delete([/plev,indlev,vnm/])
     else 
      var      = f0->$Varname$(indtim,:,:)
     end if 
     end if 
   else 
     fl0      = systemfunc("ls " + datadir + "/"+Casdir(ii)+"/" + freq + "/" +"*.nc")
     f0       = addfiles(fl0,"r")
     time     = f0[:]->time
     time     = time - 1
     utcdate  = cd_calendar(time,-1)
     indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)
     if (Varstr.eq."STRESS_MAG") then 
      var     = f0[:]->TAUX(indtim,:,:) 
      var     = sqrt(f0[:]->TAUX(indtim,:,:)*f0[:]->TAUX(indtim,:,:) + \
                     f0[:]->TAUY(indtim,:,:)*f0[:]->TAUY(indtim,:,:))
     else if (Varstr.eq."U200".or.Varstr.eq."U500".or.Varstr.eq."U850".or.Varstr.eq."Z500") then
      if(Varname.eq."Z3") then
        vnm = "Z"
      else
        vnm = Varname
      end if
      plev    = tofloat(str_sub_str(Varstr,vnm,""))
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
      xtmp = vinth2p_ecmwf(f0[:]->$Varname$,f0[0]->hyam,f0[0]->hybm,plev,f0[:]->PS,intyp,P0mb, \
                             1,kxtrp,varflg,f0[:]->TS,f0[:]->PHIS)
      var  = xtmp(indtim,0,:,:)
      delete([/plev,P0mb,intyp,kxtrp,varflg,xtmp,vnm/])
     else 
       var    = f0[:]->$Varname$(indtim,:,:)
     end if 
     end if 

   end if   
    
   delete([/indtim,utcdate,fl0,f0,time/])

   system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
   fno = outdir+ Varstr + "_" + Groups(ii) + ".nc"
   system("rm -f " + fno)

   print(" ")
   print("output : " + fno)
   print(" ")

   flo = addfile(fno,"c")
   flo->$Varstr$=var
   delete([/var,flo/])

 end do 

end do ; end of experiment group loop

end
