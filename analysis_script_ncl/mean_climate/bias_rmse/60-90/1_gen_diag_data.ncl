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

  models   = (/"CMIP6_MME","v1.LR.historical","v2.LR.historical","v1.LR.amip","v2.LR.amip"/)
  nmod     = dimsizes(models)
 
  do im = 2,nmod-1,1

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

    if(var(iv).eq."STRESS_MAG") then
      dirs  = systemfunc("ls -d "+Input)
      do kk = 0,dimsizes(dirs)-1,1
       flxs = systemfunc("ls "+dirs(kk)+"*en*TAUX*.nc")
       flys = systemfunc("ls "+dirs(kk)+"*en*TAUY*.nc")
       do jj = 0,dimsizes(flxs)-1
        flx = flxs(jj)
        fly = flys(jj)
        flz = str_sub_str(flx, "TAUX", "STRESS_MAG")
        system("rm -rvf " + flz)
        if(isfilepresent(flx).and.isfilepresent(fly).and.(.not.isfilepresent(flz))) then
          fout = addfile(flz,"wc")
          fx   = addfile(flx,"r")
          fy   = addfile(fly,"r")
          vout = fx->TAUX
          vout = sqrt(fx->TAUX*fx->TAUX + fy->TAUY* fy->TAUY)
          vout@original_name = "STRESS_MAG"
          fout->$var(iv)$ = vout
          fout->time     = fx->time
          fout->lat      = fx->lat
          fout->lon      = fx->lon
          delete([/vout,fx,fy,fout/])
        end if
        delete([/flx,fly,flz/])
       end do 
       delete([/flxs,flys/])
      end do
      delete(dirs)
    end if
    
    if(models(im).eq."CMIP6_MME".and.var(iv).eq."TGCLDLWP") then 
      dirs  = systemfunc("ls -d "+Input)
      do kk = 0,dimsizes(dirs)-1,1
       flxs = systemfunc("ls "+dirs(kk)+"*en*TGCLDCWP*.nc")
       flys = systemfunc("ls "+dirs(kk)+"*en*TGCLDIWP*.nc")
       do jj = 0,dimsizes(flxs)-1
        flx = flxs(jj)
        fly = flys(jj)
        flz = str_sub_str(flx, "TAUX", "STRESS_MAG")
        system("rm -rvf " + flz)
        if(isfilepresent(flx).and.isfilepresent(fly).and.(.not.isfilepresent(flz))) then 
          fout = addfile(flz,"wc")
          fx   = addfile(flx,"r")
          fy   = addfile(fly,"r")
          vout = fx->TGCLDCWP
          vout = fx->TGCLDCWP - fy->TGCLDIWP
          vout@original_name = "TGCLDLWP"
          fout->$var(iv)$ = vout
          fout->time     = fx->time
          fout->lat      = fx->lat
          fout->lon      = fx->lon
          delete([/vout,fx,fy,fout/])
        end if 
        delete([/flx,fly,flz/])
       end do 
       delete([/flxs,flys/])
      end do 
      delete(dirs)
    end if 

    if(models(im).ne."CMIP6_MME".and.var(iv).eq."TGCLDCWP") then
      dirs  = systemfunc("ls -d "+Input)
      do kk = 0,dimsizes(dirs)-1,1
       flxs = systemfunc("ls "+dirs(kk)+"*en*TGCLDLWP*.nc")
       flys = systemfunc("ls "+dirs(kk)+"*en*TGCLDIWP*.nc")
       do jj = 0,dimsizes(flxs)-1
        flx = flxs(jj)
        fly = flys(jj)
        flz = str_sub_str(flx, "TGCLDLWP", "TGCLDCWP")
        system("rm -rvf " + flz)
        if(isfilepresent(flx).and.isfilepresent(fly).and.(.not.isfilepresent(flz))) then
          fout = addfile(flz,"wc")
          fx   = addfile(flx,"r")
          fy   = addfile(fly,"r")
          vout = fx->TGCLDLWP
          vout = fx->TGCLDLWP + fy->TGCLDIWP
          vout@original_name = "TGCLDCWP"
          fout->$var(iv)$ = vout
          fout->time     = fx->time
          fout->lat      = fx->lat
          fout->lon      = fx->lon
          delete([/vout,fx,fy,fout/])
        end if
        delete([/flx,fly,flz/])
       end do
       delete([/flxs,flys/])
      end do
      delete(dirs)
    end if

    if(var(iv).eq."STRESS_MAG") then
      dirs  = systemfunc("ls -d "+obsdir)
      do kk = 0,dimsizes(dirs)-1,1
        flx = systemfunc("ls "+dirs(kk)+"*en00*TAUX*.nc")
        fly = systemfunc("ls "+dirs(kk)+"*en00*TAUY*.nc")
        flz = str_sub_str(flx, "TAUX", "STRESS_MAG")
        system("rm -rvf " + flz)
        if(isfilepresent(flx).and.isfilepresent(fly).and.(.not.isfilepresent(flz))) then
          fout = addfile(flz,"wc")
          fx   = addfile(flx,"r")
          fy   = addfile(fly,"r")
          vout = fx->TAUX
          vout = sqrt(fx->TAUX*fx->TAUX + fy->TAUY* fy->TAUY)
          vout@original_name = "STRESS_MAG"
          fout->$var(iv)$ = vout
          fout->time     = fx->time
          fout->lat      = fx->lat
          fout->lon      = fx->lon
          delete([/vout,fx,fy,fout/])
        end if
        delete([/flx,fly,flz/])
      end do
      delete(dirs)
    end if

    if(var(iv).eq."TGCLDCWP") then
      dirs  = systemfunc("ls -d "+obsdir)
      do kk = 0,dimsizes(dirs)-1,1
        flx = systemfunc("ls "+dirs(kk)+"*en00*TGCLDLWP*.nc")
        fly = systemfunc("ls "+dirs(kk)+"*en00*TGCLDIWP*.nc")
        flz = str_sub_str(flx, "TGCLDLWP", "TGCLDCWP")
        system("rm -rvf " + flz)
        if(isfilepresent(flx).and.isfilepresent(fly).and.(.not.isfilepresent(flz))) then
          fout = addfile(flz,"wc")
          fx   = addfile(flx,"r")
          fy   = addfile(fly,"r")
          vout = fx->TGCLDLWP
          vout = fx->TGCLDLWP + fy->TGCLDIWP
          vout@original_name = "TGCLDCWP"
          fout->$var(iv)$ = vout
          fout->time     = fx->time
          fout->lat      = fx->lat
          fout->lon      = fx->lon
          delete([/vout,fx,fy,fout/])
        end if
        delete([/flx,fly,flz/])
      end do
      delete(dirs)
    end if

  end do
  end do 
end
