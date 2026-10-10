load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

undef("nice_string")
  
function nice_string(var)

begin

  string_len = sizeof(stringtochar(var)) - 1

  if (string_len .lt. 25) then
      ilen = 16
      do it = 0, ilen - string_len
         var = var + " "
      end do
  end if

  return(var)

end

;;;;;;;;;;;;;;;;
; main program ;
;;;;;;;;;;;;;;;;

begin

  model     = getenv("MODEL_NAME")
  obs       = getenv("OBS_NAME")
  region    = getenv("REG_NAME")
  Input     = getenv("INPUT")
  var       = getenv("VAR")
  vname     = var

  yst       = toint(getenv("yst"))
  yed       = toint(getenv("yed"))

 ; setup regions 
  
  if(region.eq."Global") then 
   minlat = -90.
   maxlat =  90.
   minlon =   0.
   maxlon = 360.
  else if(region.eq."Polar")
   minlat = -90.
   maxlat = -55.
   minlon =   0.
   maxlon = 360.
  else
   print("please specify the region for calculation....")
   exit 
  end if 
  end if 

  modeldir  = Input+"/E3SMv2/"+model+"/"
  flm       = systemfunc("ls "+modeldir+"*.nc")
  fm        = addfiles(flm, "r")
  modtime   = cd_calendar(fm[0]->time, -1)
  modind    = ind(modtime.ge.yst.and.modtime.le.yed)

  obsdir    = Input+"/"+obs+"/monthly/"
  flo       = systemfunc("ls "+obsdir+"*_1x1.nc")
  fo        = addfiles(flo, "r")
  obstime   = cd_calendar(fo[0]->time, -1)
  obsind    = ind(obstime.ge.yst.and.obstime.le.yed)
 
  ListSetType (fm, "join")
  ListSetType (fo, "join")

 ;print (modtime(modind) + " " + obstime(obsind))
  if (var .eq. "NETCF") then
      vmod  = fm[:]->LWCF(:,modind,:,:)
      vmod  = vmod + fm[:]->SWCF(:,modind,:,:)
      if(isfilevar(fo[0],"LWCF")) then
        vobs  = fo[:]->LWCF(:,obsind,:,:)
        vobs  = vobs + fo[:]->SWCF(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
  else if (var .eq. "FNET") then
      vmod  = fm[:]->FSNT(:,modind,:,:)
      vmod  = vmod - fm[:]->FLNT(:,modind,:,:)
      if(isfilevar(fo[0],"FLNTOA")) then 
        vobs  = fo[:]->FSNTOA(:,obsind,:,:)
        vobs  = vobs - fo[:]->FLNTOA(:,obsind,:,:)
      else 
        vobs  = vmod
        vobs  = -99999
      end if 
  else if (var .eq. "FLNTC") then
      vmod  = fm[:]->FLNTC(:,modind,:,:)
      if(isfilevar(fo[0],"FLNTOAC")) then
        vobs  = fo[:]->FLNTOAC(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
  else if (var .eq. "FLNT") then
      vmod  = fm[:]->FLNT(:,modind,:,:)
      if(isfilevar(fo[0],"FLNTOA")) then
        vobs  = fo[:]->FLNTOA(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
  else if (var .eq. "FSNT") then
      vmod  = fm[:]->FSNT(:,modind,:,:)
      if(isfilevar(fo[0],"FSNTOA")) then
        vobs  = fo[:]->FSNTOA(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
  else if (var .eq. "PRECT") then
      vmod  = fm[:]->PRECL(:,modind,:,:) 
      vmod  = vmod + fm[:]->PRECC(:,modind,:,:) 
      vmod  = vmod * 24.0 * 3600.0 * 1000.0 
      if(isfilevar(fo[0],var)) then
        vobs  = fo[:]->$var$(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
      vmod@units = "mm/day"
      vobs@units = "mm/day"
  else if (var .eq. "STRESS_MAG") then
      vmod  = fm[:]->TAUX(:,modind,:,:)
      vmod  = sqrt(vmod^2 + (fm[:]->TAUY(:,modind,:,:))^2)
      if(isfilevar(fo[0],"TAUX")) then
        vobs  = fo[:]->TAUX(:,obsind,:,:)
        vobs  = sqrt(vobs^2 + (fo[:]->TAUY(:,obsind,:,:))^2)
      else
        vobs  = vmod
        vobs  = -99999
      end if
  else if (var .eq. "U200")  then 
      xtmp  = vinth2p(fm[:]->U,fm[0]->hyam,fm[0]->hybm,200.0,fm[:]->PS,1,fm[0]->P0/100.0,1,False)
      vmod  = xtmp(:,modind,0,:,:)  
      if(isfilevar(fo[0],"U")) then
        vobs  = fo[:]->U(:,obsind,{200},:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
      vmod@units = "m/s"
      vobs@units = "m/s"
     ;printMinMax(vmod,False)
     ;printMinMax(vobs,False)
      delete([/xtmp/])
 
  else if (var .eq. "U500")  then
      xtmp  = vinth2p(fm[:]->U,fm[0]->hyam,fm[0]->hybm,500.0,fm[:]->PS,1,fm[0]->P0/100.0,1,False)
      vmod  = xtmp(:,modind,0,:,:)
      if(isfilevar(fo[0],"U")) then
        vobs  = fo[:]->U(:,obsind,{500},:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
      vmod@units = "m/s"
      vobs@units = "m/s"
     ;printMinMax(vmod,False)
     ;printMinMax(vobs,False)
      delete([/xtmp/])

  else if (var .eq. "U850")  then
      xtmp  = vinth2p(fm[:]->U,fm[0]->hyam,fm[0]->hybm,850.0,fm[:]->PS,1,fm[0]->P0/100.0,1,False)
      vmod  = xtmp(:,modind,0,:,:)
      if(isfilevar(fo[0],"U")) then
        vobs  = fo[:]->U(:,obsind,{850},:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
     ;printMinMax(vmod,False)
     ;printMinMax(vobs,False)
      delete([/xtmp/])
      vmod@units = "m/s"
      vobs@units = "m/s"
  else if (var .eq. "Z500")  then
      xtmp  = vinth2p(fm[:]->Z3,fm[0]->hyam,fm[0]->hybm,500.0,fm[:]->PS,1,fm[0]->P0/100.0,1,False)
      vmod  = xtmp(:,modind,0,:,:)
      if(isfilevar(fo[0],"Z3")) then
        vobs  = fo[:]->Z3(:,obsind,{500},:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
      vobs = vobs / 9.8 
      vobs@units = "gpm"
      vmod@units = "gpm"
     ;printMinMax(vmod,False)
     ;printMinMax(vobs,False)
      delete([/xtmp/])

  else if (var .eq. "LWP")  then
      vmod  = fm[:]->TGCLDLWP(:,modind,:,:)
      if(isfilevar(fo[0],"TGCLDLWP")) then
        vobs  = fo[:]->TGCLDLWP(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
      vmod = vmod * 1000.0
      vobs = vobs * 1000.0
      vmod@units = "g/m2"
      vobs@units = "g/m2"
     ;printMinMax(vmod,False)
     ;printMinMax(vobs,False)

  else if (var .eq. "IWP")  then
      vmod  = fm[:]->TGCLDIWP(:,modind,:,:)
      if(isfilevar(fo[0],"TGCLDIWP")) then
        vobs  = fo[:]->TGCLDIWP(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
      vmod = vmod * 1000.0
      vobs = vobs * 1000.0
      vmod@units = "g/m2"
      vobs@units = "g/m2"
     ;printMinMax(vmod,False)
     ;printMinMax(vobs,False)

  else
      vmod  = fm[:]->$var$(:,modind,:,:)
      if(isfilevar(fo[0],var)) then
        vobs  = fo[:]->$var$(:,obsind,:,:)
      else
        vobs  = vmod
        vobs  = -99999
      end if
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 
  end if 

  printVarSummary(vobs)
  printVarSummary(vmod)
  exit

  vmod@_FillValue = -99999
  vobs@_FillValue = -99999

  if (var .eq. "PRECC" .or. var .eq. "PRECL") then
    vmod = vmod * 24.0 * 3600.0 * 1000.0
    vobs = vobs * 24.0 * 3600.0 * 1000.0
    vmod@units = "mm/day"
    vobs@units = "mm/day"
  end if 

  if (var .eq. "SHFLX" .or. var .eq. "LHFLX") then
    if (obs.eq."ERA5") then 
      vobs = vobs * -1.0 
    end if 
  end if 

  if(var.eq."AODVIS") then 
    vmod@units = "1"
    vobs@units = "1"
  end if 

  if (var .eq. "PSL") then
    vmod = vmod / 100.0
    vobs = vobs / 100.0
    vmod@units = "hPa"
    vobs@units = "hPa"
  end if

 ;printMinMax(vmod,False)
 ;printMinMax(vobs,False)

 ;mask fill values 
  if (.not.all(ismissing(vobs))) then 
    vmod = where (ismissing(vobs), vmod@_FillValue, vmod)
    vobs = where (ismissing(vmod), vobs@_FillValue, vobs)
  end if 

 ;select regional data 
  
 ;calculate time average
  vmmod = dim_avg_n_Wrap(vmod(:,:,{minlat:maxlat},{minlon:maxlon}),0)
  vmobs = dim_avg_n_Wrap(vobs(:,:,{minlat:maxlat},{minlon:maxlon}),0)
  vmdif = vmmod
  vmdif = vmmod - vmobs
  vm1d  = ndtooned(vmmod)
  vo1d  = ndtooned(vmobs)
  vm1d@_FillValue = vmod@_FillValue
  vo1d@_FillValue = vobs@_FillValue
  varunt = vmod@units
  printVarSummary(vmmod)
  exit

 ;process metrics
  metrics    = (/"Variable", "Observation", "Mean(obs)", "RMS(obs)", "Mean(model)", "RMS(model)", \
                 "Bias(m-o)", "RMSE(m-o)", "PCORR(m-o)"/)

  nmetric    = dimsizes(metrics)
  vmetric    = new(nmetric-2,double) 

  lat        = fm[0]->lat
  rad        = 4.0*atan(1.0)/180.0
  clat       = cos(lat*rad)

  vmetric(0) = wgt_areaave(vmobs, clat, 1.0, 0)
  vmetric(1) = dim_stddev(vo1d)
  vmetric(2) = wgt_areaave(vmmod, clat, 1.0, 0)
  vmetric(3) = dim_stddev(vm1d)
  vmetric(4) = wgt_areaave(vmdif, clat, 1.0, 0)
  vmetric(5) = wgt_arearmse(vmmod, vmobs, clat, 1.0, 0)
  vmetric(6) = pattern_cor(vmmod, vmobs,clat, 1)

  if (all(ismissing(vmobs)))then 
    vmetric(0) = -999.0
    vmetric(1) = -999.0
    vmetric(4) = -999.0
    vmetric(5) = -999.0
    vmetric(6) = -999.0
  end if 

  delete([/modtime,modind,flm,fm/])

  delete([/vmod,vobs,vmmod,vmobs,vmdif,vm1d,vo1d/])

;;;;;;;;;;
; Output ;
;;;;;;;;;;
  header = ""
  str    = ""
  do i = 0, nmetric-1
   header = header + str_insert(metrics(i)," ", 12)
   if(i.eq.0)then 
     tmp    = var + "("+varunt+")"
     str    = str + str_insert(tmp," ",18)
   else if (i.eq.1)then
     str    = str + str_insert(obs," ",10)
   else
     str    = str + sprintf("%10.3f",vmetric(i-2))
   end if 
   end if 
  end do 

  fname   = "./output/"+model+"_clim_output.txt"
  if (fileexists(fname)) then
    write_table(fname, "a", [/str/], "%s")
  else
    asciiwrite(fname, header)
    write_table(fname, "a", [/str/], "%s")
  end if 

end
