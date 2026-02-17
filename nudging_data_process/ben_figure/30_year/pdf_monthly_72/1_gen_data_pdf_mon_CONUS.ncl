load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;;;;;;;;;;;;;;;;
; main program ;
;;;;;;;;;;;;;;;;

begin

  region_name  = (/"1 - Polar region, 90S-60S, 60N-90N", \
                   "2 - Mid-latitude, 60S-30S, 30N-60N", \
                   "3 - Tropics, 20S-20N", \ 
                   "4 - Global", \
                   "5 - CONUS, 25-55N, 90-120W", \
                   "6 - NEUS,  25-55N, 60-90W", \
                   "7 - NEURO,  40-70N, 10-40E", \
                   "8 - NWP,   30-60N, 150-180E"/)

  regions      = (/"Polar","MidLat","Tropics","Global", \
                   "CONUS","NEUS","NEURO","NWP"/)
  nreg         = dimsizes(region_name)
  ireg         = 4

  BASE_DIR = "/anvil/scratch/x-szhang3/data_240x120"
  fcoord   = "/anvil/projects/x-phy220062/x-szhang3/map_file/coord.nc"
  models   = (/"ERA5","CLIM30","CLIM30ML"/)
  nmods    = dimsizes(models)
  outdir   = "../fig_data/"
  syear    = 1979
  eyear    = 2014
  freq     = "monthly"
  yymm     = yyyymm_time(syear,eyear, "integer")

  ilev     = 72 
  var_List = (/"U","V","T","Q"/)
  nvars    = dimsizes(var_List)
 
  do  i = 1,nmods-1,1

   model  = models(i) 
   Output = outdir 
   ym1    = toint(syear+"01")
   ym2    = toint(eyear+"12")

   do j = 0,nvars-1,1

    var      = var_List(j)
    fl1      = systemfunc("ls "+BASE_DIR +"/" + model+"/" +freq+"/"+model+".monthly.*.nc")
    f1       = addfiles(fl1, "r")
    tim1     = cd_calendar(f1[:]->time,-1)
    if(dimsizes(dimsizes(f1[0]->lat)).ge.2) then
      lat = f1[0]->lat(:,0)
      lon = f1[0]->lon(0,:)
      lon       = where(lon.ge.180.0,lon-360.0,lon)
      lat!0     = "lat"
      lat@units = "degree_north"
      lat&lat   = lat
      lon!0     = "lon"
      lon@units = "degree_east"
      lon&lon   = lon
    else
      lat = f1[0]->lat
      lon = f1[0]->lon
    end if
    ind1     = get1Dindex(tim1,yymm) ; index for each desired data
    if(model.eq."ERA5") then
      plev       = f1[0]->lev
      plev       = plev*100.                                ; match units of p0 and psfc
      plev@units = "Pa"
      fc         = addfile(fcoord,"r")
      hya        = fc->hyam                                 ; levels to which the variable
      hyb        = fc->hybm
      lev        = fc->lev
      p0         = 100000.
      x11        = pres2hybrid(plev,f1[:]->PS(ind1,:,:),p0,f1[:]->$var$(ind1,:,:,:),hya,hyb,2)
      x11!0      = "time"
      x11!1      = "lev"
      x11!2      = "lat"
      x11!3      = "lon"
      vmod     = x11(:,ilev-1,:,:)
      delete([/x11,plev,fc,hya,hyb,lev,p0/])
    else
      vmod     = f1[:]->$var$(ind1,ilev-1,:,:)
    end if

    vmod@_FillValue = -99999
    vmod!0   = "time"
    vmod!1   = "lat"
    vmod!2   = "lon"
    vmod&lat = lat 
    vmod&lon = lon 

    if(isStrSubset(var,"Q")) then 
      vmod = vmod * 1000.0
      vmod@units = "g/kg"
    end if 
    if ( isStrSubset(var,"U") .or. isStrSubset(var,"V")) then
       vmod@units = "m/s"
    end if
    if ( isStrSubset(var,"T")) then 
      vmod@units = "K"
    end if
    if (var .eq. "PS".or.var .eq. "PSL") then
      vmod = vmod / 100.0
      vmod@units = "hPa"
    end if

   ;loop each region
    do ii = ireg,ireg
      regnm = regions(ii)
      vtmp0 = vmod
      ;open netcdf file to save the data for plots;;;;;
      setfileoption("nc", "Format",  "NetCDF4")
      out_file_name = Output+model+"_"+regnm+"_"+var+"_lev"+ilev+"_pdf_"+freq+"_"+syear+"-"+eyear+".nc"
      system("rm " + out_file_name)
      fout = addfile(out_file_name,"c")
      if(regnm.eq."Polar") then
        vtmp0(:,{-60:60},:) = vmod@_FillValue
      end if
      if(regnm.eq."MidLat") then
        vtmp0(:,{-90:-60},:) = vmod@_FillValue
        vtmp0(:,{60:90},:)   = vmod@_FillValue
        vtmp0(:,{-30:30},:)  = vmod@_FillValue
      end if
      if(regnm.eq."Tropics") then
        vtmp0(:,{-90:-20},:) = vmod@_FillValue
        vtmp0(:,{20:90},:)   = vmod@_FillValue
      end if
      if(regnm.eq."Global") then
        vtmp0(:,{-90:-75},:) = vmod@_FillValue
        vtmp0(:,{75:90},:)   = vmod@_FillValue
      end if
      if(regnm.eq."CONUS") then
        lmsk  = vmod
        vtmp0 = vmod
        lmsk(:,{25:55},{-120:-90}) = vmod@_FillValue
        lmsk  = where(ismissing(lmsk),1.0,0.0) 
        vtmp0 = where(lmsk.eq.0,vmod@_FillValue,vmod) 
        delete([/lmsk/])
      end if
      if(regnm.eq."NEUS") then
        lmsk  = vmod
        vtmp0 = vmod
        lmsk(:,{25:55},{-90:-60}) = vmod@_FillValue
        lmsk  = where(ismissing(lmsk),1.0,0.0) 
        vtmp0 = where(lmsk.eq.0,vmod@_FillValue,vmod)  
        delete([/lmsk/])
      end if
      if(regnm.eq."NEURO") then
        lmsk  = vmod
        vtmp0 = vmod
        lmsk(:,{40:70},{10:40}) = vmod@_FillValue
        lmsk  = where(ismissing(lmsk),1.0,0.0)
        vtmp0 = where(lmsk.eq.0,vmod@_FillValue,vmod)
        delete([/lmsk/])
      end if
      if(regnm.eq."NWP") then
        lmsk  = vmod
        vtmp0 = vmod
        lmsk(:,{30:60},{150:180}) = vmod@_FillValue
        lmsk  = where(ismissing(lmsk),1.0,0.0)
        vtmp0 = where(lmsk.eq.0,vmod@_FillValue,vmod)
        delete([/lmsk/])
      end if
     ;calculate pdfs 
     ;opt          = True
     ;opt@bin_nice = True
     ;vpdf = pdfx(vmod, 50, opt)  

     if (isStrSubset(var,"T")) then 
       opt         = True
       opt@bin_min = 270
       opt@bin_max = 320
     end if 
     if (isStrSubset(var,"Q")) then
       opt         = True
       opt@bin_min = 0
       opt@bin_max = 20
     end if

     if (isStrSubset(var,"U")) then
       opt         = True
       opt@bin_min = -10
       opt@bin_max =  10
     end if

     if (isStrSubset(var,"V")) then
       opt         = True
       opt@bin_min = -10
       opt@bin_max =  10
     end if

     vtmp1 = ndtooned(vtmp0) 
     vtmp1@_FillValue = vtmp0@_FillValue
     inmis = ind(.not.ismissing(vtmp1)) 
     vpdf = pdfx(vtmp1(inmis), 50, opt)
     vo1  = var+"_"+regnm
     vpdf@region = region_name(ii)
     fout->$vo1$=vpdf
     delete([/vtmp0,vtmp1,vpdf,inmis/])
    end do 
    delete([/tim1,fl1,f1,ind1,vmod/])
  end do 
end do 

end
