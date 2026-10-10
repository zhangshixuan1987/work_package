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

  ModelOutRoot = "/compyfs/zhan391/darpa_scratch/post_process/generate_bias_rmse/output/"
  
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Global"/)
  nregs        = dimsizes(Regions)
  season       = (/"DJF"/)
  nseas        = dimsizes(season)

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir       = (/"v2.LR.amip", "v2.LR.historical"/)
  Groups       = (/"v2.LR.amip", "v2.LR.historical"/)
  ngrps        = dimsizes(Groups)
 
 ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  xlabnam = (/"FLDS","FLNS","FLUT","FLUTC","LWCF","FSNS","FSNTOA","FSNTOAC","SWCF","CLDTOT","U (200 hPa) ","Z3 (500 hPa)","PRECT","TMQ","PSL","LHFLX","SHFLX","STRESS MAG","TREFHT","TS"/)
  varnam  = (/"FLDS","FLNS","FLUT","FLUTC","LWCF","FSNS","FSNTOA","FSNTOAC","SWCF","CLDTOT","U200","Z500","PRECT","TMQ","PSL","LHFLX","SHFLX","STRESS_MAG","TREFHT","TS"/)
  nvars   = dimsizes(varnam)

 ;;;processe the zonal mean data;;;;
  plotmean  = new((/ngrps,nvars/),float)
  plotrmse  = new((/ngrps,nvars/),float)

 ;;open netcdf file to save the data for plots;;;;;
 setfileoption("nc", "Format",  "NetCDF4")
 out_file_name = "Bias_and_RMSE_1year_mean_with_respect_to_obs.nc"
 system("rm " + out_file_name)
 fout = addfile(out_file_name,"c")

 ovarnam                 = (/xlabnam/)
 ovarnam!0               = "nvars"
 ovarnam&nvars           = ispan(1,nvars,1)
 ovarnam&nvars@long_name = "number of variables"
 fout->Evaluated_variables=ovarnam

 ;;read and calculate the mean and rmsd;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 do igp    = 0,ngrps-1

  fname  = ModelOutRoot+ Groups(igp)+"_clim_output.txt"
  data   = asciiread(fname,-1,"string")
  vstr   = str_get_field(data(1::), 1," ")
  mobs   = stringtofloat(str_get_field(data(1::), 3," "))
  rmso   = stringtofloat(str_get_field(data(1::), 4," "))
  mmod   = stringtofloat(str_get_field(data(1::), 5," "))
  rmsd   = stringtofloat(str_get_field(data(1::), 8," "))
  
  do i = 0,nvars-1,1
    do j = 0,dimsizes(vstr)-1,1
      if (isStrSubset(vstr(j),varnam(i))) then 
        plotmean(igp,i) = (mmod(j)-mobs(j))/mobs(j)
        plotrmse(igp,i) = rmsd(j)/rmso(j)
      end if 
    end do 
  end do 

 ;;;save the data to netcdf file;;;;
 vnmbias = "Bias_"+Groups(igp)
 vnmrmse = "RMSE_"+Groups(igp)

;;;;open a netcdf file and save the data for plot;;;;;
 vbias   = plotmean(igp,:)
 vrmse   = plotrmse(igp,:)
 
 vbias!0     = "nvars"
 vbias&nvars = ispan(1,nvars,1)
 vbias&nvars@long_name = "number of variables"

 vrmse!0     = "nvars"
 vrmse&nvars = ispan(1,nvars,1)
 vrmse&nvars@long_name = "number of variables"
 
 fout->$vnmbias$=vbias
 fout->$vnmrmse$=vrmse
 delete([/vbias,vrmse/])
end do ; end igp loop

end

