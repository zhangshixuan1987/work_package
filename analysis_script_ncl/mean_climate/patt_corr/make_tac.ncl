;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; This script is used to compare the correlation between climatology and nudged data ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "./area_rectlinear_grid.ncl"

begin

  model     = getenv("MODEL_NAME")
  obsnm     = getenv("OBS_NAME")
  Input     = getenv("INPUT")
  Output    = getenv("OUTPUT")
  Obsin     = getenv("OBSDIR")
  var       = getenv("VAR")
  plev      = toint(getenv("PLEV"))
  year      = toint(getenv("year"))
  yms       = toint(year+"01")
  yme       = toint(year+"12")
 
  if(var.ne."OLR".and.var.ne."PRECT") then
    vnam1   = var+sprinti("%03d",plev)
    vnam2   = var+sprinti("%03d",plev)
    vname   = var+sprinti("%03d",plev)
    freq    = "6hourly"
  else
    if(var.eq."OLR") then
      vnam1 = "FLUT"
      if(obsnm.eq."ERA5")
        vnam2 = "OLR"
      else
        vnam2 = "OLR"
      end if
    else
      vnam1 = var
      if(obsnm.eq."ERA5")
        vnam2 = "PRECT"
      else
        vnam2 = vnam1
      end if
    end if
    vname   = var
    freq    = "daily"
  end if

  print("work on variable: " + vname)

  ;specify regions 
  region_lobnd = (/-90, -60, -20, 30, 60, -90/)
  region_hibnd = (/-60, -30, 20,  60, 90, 90/)
  region_name  = (/1, 2, 3, 4, 5, 6/)
  region_name@long_name = "1 - SH Polar, 90S-60S; " + \
                          "2 - SH Mid-lat, 60S-30S; " + \
                          "3 - Tropical zone, 20S-20N; " + \
                          "4 - NH Mid-lat, 30N-60N; " + \
                          "5 - NH Polar, 60N-90N; " + \
                          "6 - Global"
  dim0     = dimsizes(region_name)
  mon_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  ;;;process model simulations
  freqm     = "6hourly"
  modeldir  = Input+"/"+model+"/"+freqm+"/"
  flm       = systemfunc("ls "+modeldir+vnam1+"_"+year+"_*.nc")
  fm        = addfiles(flm, "r")
  modtime   = cd_calendar(fm[:]->time, -1)
  modind    = ind(modtime.ge.yms.and.modtime.le.yme)

  if(obsnm.eq."TRMM") then
    freqo   = "daily"
    obsdir  = Obsin+"/Observation/"+obsnm+"/"+freqo+"/"
  else if(obsnm.eq."ERA5") then
    freqo   = "6hourly"
    obsdir  = Obsin+"/CVDP_RGD/"+obsnm+"."+freqo+"/"
  else
    freqo   = "daily"
    obsdir  = Obsin+"/Observation/"+obsnm+"/"+freqo+"/"
  end if
  end if

  if(obsnm.eq."ERA5") then
    flo     = systemfunc("ls "+obsdir+"*"+vnam2+"*"+yms+"-"+yme+".nc")
  else
    flo     = systemfunc("ls "+obsdir+"*"+year+"*.nc")
  end if

  fo        = addfiles(flo, "r")
  if(isleapyear(year)) then
    obstime   = cd_calendar(fo[:]->time, -1)
    ymdh      = cd_calendar(fo[:]->time, 0)
    obsind    = ind(obstime.ge.yms.and.obstime.le.yme.and.\
                    (.not.(ymdh(:,1).eq.2 .and. ymdh(:,2).eq.29)))
    delete(ymdh)
  else
    obstime   = cd_calendar(fo[:]->time, -1)
    obsind    = ind(obstime.ge.yms.and.obstime.le.yme)
  end if
 ;print(modtime(modind))
 ;print(obstime(obsind))

 ;print(flo)
 ;print(flm)

  ;output data directory
  OUTDir    = Output+"/"+model+"-"+obsnm
  ret       = systemfunc("test -d "+OUTDir+"; echo $?")
  if (ret .ne. 0) then
      ret   = systemfunc("mkdir -p "+OUTDir)
  end if

  lat  = fm[0]->lat
  lon  = fm[0]->lon

  if(dimsizes(dimsizes(fm[:]->$vnam1$)).eq.4) then
    vm1  = fm[:]->$vnam1$(modind,0,:,:)
  else
    vm1  = fm[:]->$vnam1$(modind,:,:)
  end if
  if(dimsizes(dimsizes(fo[:]->$vnam2$)).eq.4) then
    vo1  = fo[:]->$vnam2$(obsind,0,:,:)
  else
    vo1  = fo[:]->$vnam2$(obsind,:,:)
  end if
  vm1@_FillValue = 9.e36
  vo1@_FillValue = 9.e36
  delete([/obstime,obsind,flo,flm,fm,fo/])

  if (freq.eq."daily") then 
 
   if(freqm.eq."daily")
     v1 = vm1
   else
     opt     = True
     opt@nval_crit = 4
     v1 = calculate_daily_values (vm1, "avg", 0, opt)
     v1!0 = "time"
   end if 
   if(freqo.eq."daily")
     v2 = vo1
   else 
     opt     = True
     opt@nval_crit = 4
     v2 = calculate_daily_values (vo1, "avg", 0, opt)
     v2!0 = "time"
   end if 

   delete([/vm1,vo1/])

  else

   tim1  = cd_calendar(vm1&time, -3)
   tim2  = cd_calendar(vo1&time, -3)
   ntim1 = dimsizes(tim1)
   ntim2 = dimsizes(tim2)
  ;print(ntim1 + " " +ntim2)

   if (ntim1 .lt. ntim2) then
     idate = get1Dindex(tim2, tim1) ; index for each desired data
     v1    = vm1
     v2    = vo1(idate,:,:)
   else
     idate = get1Dindex(tim1, tim2) ; index for each desired data
     v2    = vo1
     v1    = vm1(idate,:,:)
   end if
   delete([/tim1,tim2,vm1,vo1,ntim1,ntim2,idate/])

  end if 

  dim1 = dimsizes(v1(:,0,0))
  dim2 = dimsizes(lat)
  dim3 = dimsizes(lon)

  if (var.eq."OMEGA") then 
   ncase = 3   
  else
   ncase = 1
  end if 

  do i = 1 , ncase

    if  (i .eq. 1) then
      fname = OUTDir+"/tac_"+vname+"_"+year+".nc"
      system("rm -rvf "+fname)
      fout  = addfile(fname,"cw")
      xv1   = v1 
      xv2   = v2 
    else if (i .eq. 2) then
      fname = OUTDir+"/tac_"+vname+"_NEGATIVE"+"_"+year+".nc"
      system("rm -rvf "+fname)
      fout  = addfile(fname,"cw")
      xv1   = mask(v1, v1 .lt. 0., True)
      xv2   = mask(v2, v2 .lt. 0., True)
    else
      fname = OUTDir+"/tac_"+vname+"_POSITIVE"+"_"+year+".nc"
      system("rm -rvf "+fname)
      fout = addfile(fname,"cw")
      xv1  = mask(v1, v1 .gt. 0., True)
      xv2  = mask(v2, v2 .gt. 0., True)
    end if
    end if

   ;printMinMax(xv1,False)
   ;printMinMax(xv2,False)

    if (var .eq. "PRECT") then
      xv1 = mask(xv1, xv1 .gt. -1.e-20, True)
      xv2 = mask(xv2, xv2 .gt. -1.e-20, True)
    end if
    xv1 = where(ismissing(xv2),xv1@_FillValue,xv1)
    xv2 = where(ismissing(xv1),xv2@_FillValue,xv2)

   ;printMinMax(xv1,False)
   ;printMinMax(xv2,False)

   ;;;;;;;;;;;;;;;;;;;;;;;;
   ;; create the anomaly ;;
   ;;;;;;;;;;;;;;;;;;;;;;;;
    opt = True
    opt@nval_crit = 12  ; require at least 12 values for the "avg" is calculated.
    v1Mon = calculate_monthly_values (xv1, "avg", 0,opt)
    v2Mon = calculate_monthly_values (xv2, "avg", 0,opt)
    v1Mon!0 = "time"
    v2Mon!0 = "time"
   ;utc_date1 = cd_calendar(v2Mon&time, -2)
   ;utc_date2 = cd_calendar(v1Mon&time, -2)
   ;print(utc_date1+ " " +utc_date2)
   ;exit
    utc_date  = cd_calendar(v2Mon&time, 0)
    month     = tointeger(utc_date(:,1))
    ids = 0
    do i = 1, dimsizes(month)
      if  (freq .eq. "daily") then
        ide = ids + mon_days(month(i-1)-1)
      else if (freq .eq. "6hourly") then
        ide = ids + mon_days(month(i-1)-1)*4
      else
        print("Unsupported FREQ for SAC = "+freq)
        exit
      end if
      end if
     ;print(ids +  " " +ide)
      xv1(ids:ide-1,:,:) = xv1(ids:ide-1,:,:) - conform(xv1(ids:ide-1,:,:),v1Mon(i-1,:,:),(/1,2/))
      xv2(ids:ide-1,:,:) = xv2(ids:ide-1,:,:) - conform(xv2(ids:ide-1,:,:),v2Mon(i-1,:,:),(/1,2/))
      ids = ide
    end do
    delete([/utc_date,month,v1Mon,v2Mon,opt/])
   ;printMinMax(xv1,False)
   ;printMinMax(xv2,False)

    ;;;calculate the correlation;;;;;;;; 
    v11 = conform_dims((/1,dim1,dim2,dim3/), xv1, (/1,2,3/))
    v21 = conform_dims((/1,dim1,dim2,dim3/), xv2, (/1,2,3/))
    v11!0 = "ens"
    v11!1 = "time"
    v11!2 = "lat"
    v11!3 = "lon"
    v21!0 = "ens"
    v21!1 = "time"
    v21!2 = "lat"
    v21!3 = "lon"
    v11@_FillValue = v1@_FillValue  
    v21@_FillValue = v2@_FillValue

   ;printVarSummary(v21)
   ;printVarSummary(v11)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; https://www.ncl.ucar.edu/Document/Functions/Contributed/pattern_cor.shtml ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    corr = pattern_cor(v11(lat|:,lon|:,ens|:,time|:), v21(lat|:,lon|:,ens|:,time|:), 1.0, 1)
   ;corr = pattern_cor(vt1,vt2,1.0,1)
   ;printVarSummary(corr)
   ;print(corr)
   ;exit
    corr!0 = "lat"
    corr!1 = "lon"
    corr&lat = lat
    corr&lon = lon
    corr@long_name = "Temporal anomaly uncentered correlation"
    corr@units = "unitless"
    corr@_FillValue = 9.e+36
    fout->corr = corr
    delete([/corr,xv1,xv2,v11,v21/])
 
  end do ; end of case 
  delete([/fout,fname,v1,v2/])

end
exit

