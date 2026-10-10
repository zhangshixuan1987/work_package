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
  obsdir   = "/global/cfs/cdirs/e3sm/zhan391/data/CVDP_RGD"
  datadir  = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"
  freq     = "monthly"
  ystr     = 200801
  yend     = 201712

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"GPCP","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = Casdir
  ngrps    = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"PRECT"/)
  VstList  = (/"PRECT"/)
  UnitList = (/"mm day~S~-1~N~"/)
  nvars    = dimsizes(VarList)

  do ii = 0, ngrps-1,1

   do jj = 0,nvars-1, 1

   Varname = VarList(jj)
   Varunt  = UnitList(jj)
   Varstr  = VstList(jj)

   if(ii.eq.0) then 
     fl0      = systemfunc("ls " + obsdir + "/"+Casdir(ii)+"/" + Casdir(ii)+"*.nc")
     print(fl0)
     f0       = addfile(fl0,"r")
     time     = f0->time
     utcdate  = cd_calendar(time,-1)
     indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)
     var      = f0->$Varname$(indtim,:,:)
   else 
     fl0      = systemfunc("ls " + datadir + "/"+Casdir(ii)+"/" + freq + "/" +"*.nc")
     f0       = addfiles(fl0,"r")
     time     = f0[:]->time
     time     = time - 1
     utcdate  = cd_calendar(time,-1)
     indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)
     if(Varstr.eq."PRECT") then 
      var     = f0[:]->PRECC(indtim,:,:)
      var     = f0[:]->PRECC(indtim,:,:) + f0[:]->PRECL(indtim,:,:)
      var@long_name = "Total Precipitation"
     else 
       var    = f0[:]->$Varname$(indtim,:,:)
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
