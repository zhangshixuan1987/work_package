load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;;;;;;;;;;;;;;;;
; main program ;
;;;;;;;;;;;;;;;;

begin

  BASE_DIR = "/anvil/scratch/x-szhang3/data_240x120"
  models   = ("CLIM30ML")
  nmods    = dimsizes(models)
  outdir   = "../fig_data/"

  syear    = 1979
  eyear    = 2014
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
    flm      = systemfunc("ls "+modeldir+model+"*.sfc.*.nc")
    print(flm)

    ;do kk = 0,dimsizes(flm)-1 
    ; ff = addfile(flm(kk), "r")
    ; if(.not.isfilevar(ff,var))then 
    ;  print(flm(kk) + " " +var)
    ; end if 
    ;  delete(ff)
    ;end do 
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
      vmod = fm[:]->$var$(indt,:,:)
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
