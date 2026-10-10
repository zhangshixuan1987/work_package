load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;;;;;;;;;;;;;;;;
; main program ;
;;;;;;;;;;;;;;;;

begin

  data_dir = "/global/cfs/cdirs/e3sm/zhan391/data/CVDP_RGD/"
  out_dir  = "./fig_data"
  latS     = -90
  latN     = -60
  lonW     = 0 
  lonE     = 360
  
  var      = (/"CLDTOT","TGCLDCWP","TGCLDIWP","PRECST",\
               "PSL","STRESS_MAG","TREFHT","TS",\
               "U200","T850","U850","Z500",\
               "LHFLX","SHFLX","PRECT","TMQ",\
               "OLR","FSNTOA","FNET","CRE",\
               "LWCRE","SWCRE","FLDS","FLNS",\
               "FLUTC","FLDSC","FSNS","FSNTOAC","SICONC"/)
  obdat    = (/"CERES-OAFlux","ERA5","ERA5","ERA5",\
               "ERA5","ERA5","ERA5","ERA5",\
               "ERA5","ERA5","ERA5","ERA5",\
               "ERA5","ERA5","GPCP","ERA5",\
               "CERES-OAFlux","CERES-OAFlux","CERES-OAFlux","CERES-OAFlux",\
               "CERES-OAFlux","CERES-OAFlux","CERES-OAFlux","CERES-OAFlux",\
               "CERES-OAFlux","CERES-OAFlux","CERES-OAFlux","CERES-OAFlux", "HadSST"/)
  nvar     = dimsizes(var)

  models   = (/"CMIP6_MME","v1.LR.historical","v2.LR.historical","v2.LR.amip","v1.LR.amip"/)
  nmod     = dimsizes(models)
 
  do im = 0,nmod-1,1

   if (models(im).eq."CMIP6_MME") then 
     Input = data_dir + models(im)+"/historical/*/"
   else
     Input = data_dir + models(im)+"/"
   end if 

   do iv   = 0,nvar-1,1
 
    if(obdat(iv).eq."CERES-OAFlux") then
      yst  = 200101
      yed  = 201412
    else
      yst  = 197901
      yed  = 201412
    end if 

    obsdir = data_dir + obdat(iv) +"/"

    if (models(im).eq."CMIP6_MME") then
      flm     = systemfunc("ls "+Input+"*en00*."+var(iv)+".*.nc")
      fm      = addfiles(flm, "r")
      ListSetType (fm, "join")
    else
      flm     = systemfunc("ls "+Input+"*en*."+var(iv)+".*.nc")
      fm      = addfiles(flm, "r")
      ListSetType (fm, "join")
    end if 

    if (var(iv).eq."SICONC".and.(models(im).eq."v2.LR.amip".or.models(im).eq."v1.LR.amip")) then
      delete([/flm,fm/])
      flm     = systemfunc("ls "+Input+"*en*.TS.*.nc")
      fm      = addfiles(flm, "r")
      ListSetType (fm, "join")
    end if 

    if(var(iv).eq."FLNTC") then 
      flo     = systemfunc("ls "+obsdir+"*en00*.FLNTOAC.*.nc")
    else 
      flo     = systemfunc("ls "+obsdir+"*en00*."+var(iv)+".*.nc")
    end if 

    fo        = addfile(flo, "r")
    timm      = cd_calendar(fm[0]->time,1)
    timo      = cd_calendar(fo->time,1)
    indm      = ind(timm(:).ge.yst.and.timm(:).le.yed)
    indo      = ind(timo(:).ge.yst.and.timo(:).le.yed)

    if (var(iv).eq."SICONC".and.(models(im).eq."v2.LR.amip".or.models(im).eq."v1.LR.amip")) then
      vtmp0   = fm[:]->TS
      vtmp0   = -9999.0
      vtmp0@_FillValue = -9999.0
    else 
      vtmp0   = fm[:]->$var(iv)$
    end if 
    if(max(vtmp0&lon).le.180.0) then
      vtmp0 = lonFlip(vtmp0)
    end if
    vmdat  = vtmp0(:,indm,{latS:latN},{lonW:lonE})
   
    vtmp1     = fo->$var(iv)$
    if(max(vtmp1&lon).le.180.0) then
      vtmp1 = lonFlip(vtmp1)
    end if
    vodat     = vtmp1(indo,{latS:latN},{lonW:lonE})
    delete([/vtmp0,vtmp1,timm,timo,indm,indo/])

    if(var(iv).eq."TGCLDCWP".or.var(iv).eq."TGCLDLWP".or.var(iv).eq."TGCLDIWP".or.var(iv).eq."TMQ") then
      vodat  = where(vodat.lt.0.0,0.0,vodat)
      vmdat  = where(vmdat.lt.0.0,0.0,vmdat)
    end if 
 
    nens      = dimsizes(vmdat(:,0,0,0))

    do ie = 0,nens-1,1

      vmod = vmdat(ie,:,:,:)
      vobs = vodat

      if(var(iv).eq."CLDTOT") then
        if(max(vobs).le.1) then 
          vobs = vobs * 100.0
        end if 
        if(max(vmod).le.1) then
          vmod = vmod * 100.0
        end if
        vmod@units = "%"
        vobs@units = "%"
      end if

      if(var(iv).eq."SICONC") then 
        vmod@units = "%"
        vobs@units = "%"
      end if 

      if (var(iv).eq."PSL") then
        vmod = vmod / 100.0
        vobs = vobs / 100.0
        vmod@units = "hPa"
        vobs@units = "hPa"
      end if
 
     ;mask fill values 
      if (.not.all(ismissing(vobs))) then 
        vmod = where (ismissing(vobs), vmod@_FillValue, vmod)
        vobs = where (ismissing(vmod), vobs@_FillValue, vobs)
      end if 

     ;calculate time average
      vmmod = dim_avg_n_Wrap(vmod,0)
      vmobs = dim_avg_n_Wrap(vobs,0)
      vmdif = vmmod
      vmdif = vmmod - vmobs
      vm1d  = ndtooned(vmmod)
      vo1d  = ndtooned(vmobs)
      vm1d@_FillValue = vmod@_FillValue
      vo1d@_FillValue = vobs@_FillValue
      varunt = vmod@units

     ;process metrics
      metrics    = (/"Variable", "Observation", "Mean(obs)", "RMS(obs)", "Mean(model)", "RMS(model)", \
                     "Bias(m-o)", "RMSE(m-o)", "PCORR(m-o)"/)

      nmetric    = dimsizes(metrics)
      vmetric    = new(nmetric-2,double) 

      lat        = vmod&lat
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

      ;;;;;;;;;;
      ; Output ;
      ;;;;;;;;;;
      header = ""
      str    = ""
      do i = 0, nmetric-1
       header = header + str_insert(metrics(i)," ", 12)
       if(var(iv).eq."TGCLDLWP") then 
         varo = "LWP"
       else if (var(iv).eq."TGCLDIWP") then  
         varo = "IWP"
       else if (var(iv).eq."TGCLDCWP") then
         varo = "CWP"
       else
         varo = var(iv)
       end if 
       end if 
       end if 
       if(obdat(iv).eq."CERES-OAFlux") then 
         obso = "CERES-EBAF"
       else
         obso = obdat(iv)
       end if 
       if(i.eq.0)then 
         tmp    = varo + "("+varunt+")"
         str    = str + str_insert(tmp," ",18)
       else if (i.eq.1)then
         str    = str + str_insert(obso," ",10)
       else
         str    = str + sprintf("%10.3f",vmetric(i-2))
       end if 
       end if 
      end do 

      fname   = out_dir+"/"+models(im)+"_clim_output_"+sprinti("%02d",ie)+".txt"
      if (fileexists(fname)) then
        write_table(fname, "a", [/str/], "%s")
      else 
        asciiwrite(fname, header)
        write_table(fname, "a", [/str/], "%s")
      end if 
      delete([/vmod,vobs,vmmod,vmobs,vmdif,vm1d,vo1d/])
    end do 
    delete([/flm,flo,fm,fo,vmdat,vodat/])
  end do
  end do 
end
