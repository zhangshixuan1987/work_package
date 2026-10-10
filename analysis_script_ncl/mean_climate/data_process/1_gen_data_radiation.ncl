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
  obsdir   = "/global/cfs/cdirs/e3sm/zhan391/data/Observation"
  datadir  = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"
  freq     = "monthly"
  ystr     = 200801
  yend     = 201712

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"CERES_EBAF","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = Casdir
  ngrps    = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"CLDTOT",   "SWCF",        "FLUT",        "LWCF",        "NETCF", \
               "FLNT",     "FSNT",        "FNET",        "FLDS",        "FLNS" , \
               "FLUTC",    "FLNTC",       "FSNS",        "FSNTOAC",     "FSNTOA"/)
  VstList  = (/"CLDTOT",   "SWCRE",       "FLUT",        "LWCRE",        "CRE", \
               "FLNTOA",   "FSNTOA",      "FNET",        "FLDS",        "FLNS" , \
               "FLUTC",    "FLNTC",       "FSNS",        "FSNTOAC",     "FSNTOA" /)
  nvars    = dimsizes(VarList)

  do ii = 0, ngrps-1,1

   do jj = 0,nvars-1, 1

   Varname = VarList(jj)
   Varstr  = VstList(jj)

   if(ii.eq.0) then 
     fl0      = systemfunc("ls " + obsdir + "/"+Casdir(ii)+"/" + Casdir(ii)+"*1x1.nc")
     print(fl0)
     f0       = addfile(fl0,"r")
     time     = f0->time
     utcdate  = cd_calendar(time,-1)
     indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)
     if(Varname.eq."FLNT") then 
       var    = f0->FLNTOA(indtim,:,:)
     else if (Varname.eq."FLNTC") then
       var    = f0->FLNTOAC(indtim,:,:)  
     else if (Varname.eq."FSNT") then
       var    = f0->FSNTOA(indtim,:,:)
     else if (Varname.eq."FSNTC") then
       var    = f0->FSNTOAC(indtim,:,:)
     else if (Varname.eq."FNET") then
       var    = f0->FSNTOA(indtim,:,:)
       var    = f0->FSNTOA(indtim,:,:) + f0->FLNTOA(indtim,:,:) 
     else if(Varname.eq."CLDTOT") then
       var    = f0->$Varname$(indtim,:,:)
       var    = var /100.0
       var@units = "fraction"
     else 
       var    = f0->$Varname$(indtim,:,:)
     end if 
     end if 
     end if
     end if 
     end if 
     end if 
   else 
     fl0      = systemfunc("ls " + datadir + "/"+Casdir(ii)+"/" + freq + "/" +"*.nc")
     f0       = addfiles(fl0,"r")
     time     = f0[:]->time
     time     = time - 1
     utcdate  = cd_calendar(time,-1)
     indtim   = ind (utcdate.ge.ystr.and.utcdate.le.yend)
     if(Varstr.eq."CRE") then 
      var     = f0[:]->SWCF(indtim,:,:)
      var     = f0[:]->SWCF(indtim,:,:) + f0[:]->LWCF(indtim,:,:)
      var@long_name = "Cloud radiative effects"
     else if (Varstr.eq."FNET") then
      var     = f0[:]->FSNT(indtim,:,:)
      var     = f0[:]->FLNT(indtim,:,:) + f0[:]->FSNT(indtim,:,:)
      var@long_name = "All-sky net radiative flux at top of model"
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
