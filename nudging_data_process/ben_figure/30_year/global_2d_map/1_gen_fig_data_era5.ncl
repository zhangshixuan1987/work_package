load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;;;;;;;;;;;;;;;;
; main program ;
;;;;;;;;;;;;;;;;

begin

  BASE_DIR = "/anvil/scratch/x-szhang3/data_240x120"
  models   = ("ERA5")
  nmods    = dimsizes(models)
  fcoord   = "/anvil/projects/x-phy220062/x-szhang3/map_file/coord.nc"
  frefx    = "/anvil/scratch/x-szhang3/data_240x120/ERA5/monthly/ERA5.monthly.2011.nc"
  outdir   = "../fig_data/"
  
  freq     = "monthly"
  yst      = 1979
  yed      = 2014
  yymm     = yyyymm_time(yst,yed, "integer")

  var_List = (/"PSL","UBOT","VBOT","QBOT", "TBOT", \
               "TUQ","TVQ"/)
  nvars   = dimsizes(var_List)
  
  do  i = 0,nmods-1,1

   model  = models(i) 
   Input  = BASE_DIR +"/" + models(i)
   Output = outdir 

   do j = 0, nvars-1,1

    var = var_List(j)
    modeldir = Input+ "/" + freq+ "/"
    if(var.ne."PSL".and.var.ne."TUQ".and.var.ne."TVQ") then 
      flm    = systemfunc("ls "+modeldir+model+"*.monthly.*.nc")
    else 
      flm    = systemfunc("ls "+modeldir+model+"*.sfc.*.nc")
    end if 
    print(flm)

    fm       = addfiles(flm, "r")
    tim0     = cd_calendar(fm[:]->time,-1)
    indt     = get1Dindex(tim0,yymm) ; index for each desired data

    if(model.eq."ERA5".and.(var.eq."TVQ".or.var.eq."TUQ"))then 
      if(var.eq."TUQ") then 
        vmod = fm[:]->VIWVE(indt,:,:)
      end if 
      if (var.eq."TVQ") then 
        vmod = fm[:]->VIWVN(indt,:,:)
      end if 
    else
      if(var.eq."PSL") then 
        vmod = fm[:]->$var$(indt,:,:)
      else 
        plev       = fm[0]->lev
        plev       = plev*100.                                ; match units of p0 and psfc
        plev@units = "Pa"
        fc         = addfile(fcoord,"r")
        fx         = addfile(frefx,"r")
        hya        = fc->hyam                                 ; levels to which the variable
        hyb        = fc->hybm
        lev        = fc->lev
        p0         = 100000.
        if (isStrSubset(var,"UBOT")) then 
          vnm = "U" 
        end if 
        if (isStrSubset(var,"VBOT")) then
          vnm = "V"
        end if
        if (isStrSubset(var,"TBOT")) then
          vnm = "T"
        end if
        if (isStrSubset(var,"QBOT")) then
          vnm = "Q"
        end if
        x11 = pres2hybrid(plev,fm[:]->PS(indt,:,:),p0,fm[:]->$vnm$(indt,:,:,:),hya(71),hyb(71),3)
        vmod = x11(:,0,:,:)
        vmod!0 = "time"
        vmod!1 = "lat"
        vmod!2 = "lon"
        vmod&lat = fx->lat
        vmod&lon = fx->lon
        delete([/plev,fc,hya,hyb,lev,p0,x11,fx/])
      end if 
    end if 

    vmod@_FillValue = -99999

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

    ;calculate time average
    vmmod = dim_avg_n_Wrap(vmod,0)
    vsmod = dim_variance_n_Wrap(vmod,0)
    
   ;open netcdf file to save the data for plots;;;;;
    setfileoption("nc", "Format",  "NetCDF4")
    out_file_name = Output+model+"_"+var+"_"+yst+"-"+yed+".nc"
    system("rm " + out_file_name)
    fout = addfile(out_file_name,"c")
    vo1=var
    vo2="variance"
    fout->$vo1$=vmmod
    fout->$vo2$=vsmod
    delete([/tim0,flm,fm,indt,vmmod,vsmod,vmod/])
    
  end do 
end do 

end
