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
  region_lobnd = (/-90, -60, -20, -90 /)
  region_hibnd = (/ 90,  60,  20,  90 /)
  region_name = (/1, 2, 3, 4/)
  region_name@long_name = "1 - Polar region, 90S-60S, 60N-90N; " + \
                          "2 - Mid-latitude, 60S-30S, 30N-60N; " + \
                          "3 - Tropical zone, 20S-20N; " + \
                          "4 - Global"
  dim0     = dimsizes(region_name)
  mon_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  ;;;process model simulations
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

  ;output data directory 
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

  if (var.eq."OMEGA") then 
   ncase = 3   
  else
   ncase = 1
  end if 

  do i = 1 , ncase

    if  (i .eq. 1) then
      fname = OUTDir+"/sac_3regions_"+vname+"_"+year+".nc"
      system("rm -rvf "+fname)
      fout = addfile(fname,"cw")
      xv1  = v1
      xv2  = v2
    else if (i .eq. 2) then
      fname = OUTDir+"/sac_3regions_"+vname+"_NEGATIVE"+"_"+year+".nc"
      system("rm -rvf "+fname)
      fout = addfile(fname,"cw")
      xv1 = mask(xv1, xv1 .lt. 0., True)
      xv2 = mask(xv2, xv2 .lt. 0., True)
    else
      fname = OUTDir+"/sac_3regions_"+vname+"_POSITIVE"+"_"+year+".nc"
      system("rm -rvf "+fname)
      fout = addfile(fname,"cw")
      xv1 = mask(xv1, xv1 .gt. 0., True)
      xv2 = mask(xv2, xv2 .gt. 0., True)
    end if
    end if

    if (var .eq. "PRECT") then
       xv1 = mask(xv1, xv1 .gt. -1.e-20, True)
       xv2 = mask(xv2, xv2 .gt. -1.e-20, True)
    end if

    xv1 = where(ismissing(xv2),xv1@_FillValue,xv1)
    xv2 = where(ismissing(xv1),xv2@_FillValue,xv2)

    dim1 = dimsizes(xv1(:,0,0))
    dim2 = dimsizes(lat)
    dim3 = dimsizes(lon)
    corr = new((/dim1,dim0/), double)

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
    do j = 0, dim0-1 

      ;;;calculate the correlation;;;;;;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; https://www.ncl.ucar.edu/Document/Functions/Contributed/pattern_cor.shtml ;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      gridarea = tofloat(area_global_rectilinear_grid(lat, lon, False))
      gridarea@_FillValue = 9.e36

      if  (j .eq. 0) then
        a    = closest_val(-90.,lat)
        b    = closest_val(-60.,lat)
        c    = closest_val(60.,lat)
        d    = closest_val(90.,lat)
        nx1  = dimsizes(gridarea(a:b,0))
        nx2  = dimsizes(gridarea(c:d,0))
        dimx = nx1 + nx2
        wgt  = new ((/dimx,dim3/),typeof(gridarea))
        v11  = new ((/dim1,dimx,dim3/),typeof(v1))
        v21  = new ((/dim1,dimx,dim3/),typeof(v2))
        wgt(0:nx1-1,:)      = gridarea(a:b,:)
        wgt(nx1:dimx-1,:)   = gridarea(c:d,:)
        v11(:,0:nx1-1,:)    = xv1(:,a:b,:)
        v11(:,nx1:dimx-1,:) = xv1(:,c:d,:)
        v21(:,0:nx1-1,:)    = xv2(:,a:b,:)
        v21(:,nx1:dimx-1,:) = xv2(:,c:d,:)
        wgt@_FillValue = gridarea@_FillValue
        v11@_FillValue = xv1@_FillValue
        v21@_FillValue = xv2@_FillValue

        avgarea  = sum(wgt) / (num(.not. ismissing(wgt)) * 1.0)
        v11  = v11 * conform(v11,wgt,(/1,2/)) / avgarea
        v21  = v21 * conform(v21,wgt,(/1,2/)) / avgarea
        corr(:,j) = pattern_cor(v11, v21, 1.0, 1)

        delete([/nx1,nx2,dimx,wgt,v11,v21,avgarea,gridarea/])

      else if  (j .eq. 1) then
        a = closest_val(-60.,lat)
        b = closest_val(-30.,lat)
        c = closest_val(30.,lat)
        d = closest_val(60.,lat)
        nx1  = dimsizes(gridarea(a:b,0))
        nx2  = dimsizes(gridarea(c:d,0))
        dimx = nx1 + nx2
        wgt  = new ((/dimx,dim3/),typeof(gridarea))
        v11  = new ((/dim1,dimx,dim3/),typeof(v1))
        v21  = new ((/dim1,dimx,dim3/),typeof(v2))
        wgt(0:nx1-1,:)      = gridarea(a:b,:)
        wgt(nx1:dimx-1,:)   = gridarea(c:d,:)
        v11(:,0:nx1-1,:)    = xv1(:,a:b,:)
        v11(:,nx1:dimx-1,:) = xv1(:,c:d,:)
        v21(:,0:nx1-1,:)    = xv2(:,a:b,:)
        v21(:,nx1:dimx-1,:) = xv2(:,c:d,:)
        wgt@_FillValue = gridarea@_FillValue
        v11@_FillValue = v1@_FillValue
        v21@_FillValue = v2@_FillValue

        avgarea  = sum(wgt) / (num(.not. ismissing(wgt)) * 1.0)
        v11  = v11 * conform(v11,wgt,(/1,2/)) / avgarea
        v21  = v21 * conform(v21,wgt,(/1,2/)) / avgarea
        corr(:,j) = pattern_cor(v11, v21, 1.0, 1)
  
        delete([/nx1,nx2,dimx,wgt,v11,v21,avgarea,gridarea/])

      else 

        a = closest_val(region_lobnd(j),lat)
        b = closest_val(region_hibnd(j),lat)
        avgarea = sum(gridarea(a:b,:)) / (num(.not. ismissing(gridarea(a:b,:))) * 1.0)
        v11  = xv1(:,a:b,:)
        v21  = xv2(:,a:b,:)
        wgt  = gridarea(a:b,:)
        v11  = v11 * conform(v11,wgt,(/1,2/)) / avgarea
        v21  = v21 * conform(v21,wgt,(/1,2/)) / avgarea
        corr(:,j) = pattern_cor(v11, v21, 1.0, 1)
        delete([/v11,v21,wgt,avgarea,gridarea/])
      end if 
      end if 

    end do 

    corr!0 = "time"
    corr!1 = "num_region"
    corr&num_region = region_name
    corr@long_name = "Spatial anomaly uncentered correlation"
    corr@units = "unitless"
    corr@_FillValue = 9.e+36 
    fout->corr = corr
    delete([/corr,xv1,xv2/])

 end do 
 delete([/fout,fname,v1,v2/])

end
exit

