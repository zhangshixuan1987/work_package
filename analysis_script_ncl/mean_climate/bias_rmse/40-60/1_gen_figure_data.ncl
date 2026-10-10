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

  ModelOutRoot = "/pscratch/sd/z/zhan391/DARPA_project/evaluation/paper_material/nudging_evaluation/fig_data/"
  Casdir       = (/"CLIM","NDGUV", "NDGUVT", "NDGUVTQ", "NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups       = Casdir
  ngrps        = dimsizes(Groups)
  ystr         = 2008
  yend         = 2017
  nyr          = yend - ystr + 1

 ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  xlabnam = (/"FLDS","FLNS","FLUT","FLUTC","LWCRE","FSNS","FSNTOA","FSNTOAC","SWCRE","CLDTOT",\
              "U (200 hPa) ","Z3 (500 hPa)","PRECT","TMQ","PSL","LHFLX","SHFLX","STRESS MAG","TREFHT","TS"/)
  varnam  = (/"FLDS","FLNS","FLUT","FLUTC","LWCRE","FSNS","FSNTOA","FSNTOAC","SWCRE","CLDTOT",\
              "U200","Z500","PRECT","TMQ","PSL","LHFLX","SHFLX","STRESS_MAG","TREFHT","TS"/)
  nvars   = dimsizes(varnam)

 ;;;processe the zonal mean data;;;;
  plotmean  = new((/ngrps,nvars/),float)
  plotrmse  = new((/ngrps,nvars/),float)
  merr      = new((/nyr,nvars/),float)
  mrms      = new((/nyr,nvars/),float)

 ;;open netcdf file to save the data for plots;;;;;
 setfileoption("nc", "Format",  "NetCDF4")
 out_file_name = "Bias_and_RMS_with_respect_to_obs_2008-2017.nc"
 system("rm " + out_file_name)
 fout = addfile(out_file_name,"c")

 ovarnam                 = (/xlabnam/)
 ovarnam!0               = "nvars"
 ovarnam&nvars           = ispan(1,nvars,1)
 ovarnam&nvars@long_name = "number of variables"
 fout->Evaluated_variables=ovarnam

 ;;read and calculate the mean and rmsd;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 do igp = 0,ngrps-1

  do iy = 0, nyr-1
    year   = ystr + iy 
    fname  = ModelOutRoot+ Groups(igp)+"_clim_output_"+year+".txt"
    data   = asciiread(fname,-1,"string")
    vstr   = str_get_field(data(1::), 1," ")
    mobs   = stringtofloat(str_get_field(data(1::), 3," "))
    rmso   = stringtofloat(str_get_field(data(1::), 4," "))
    mmod   = stringtofloat(str_get_field(data(1::), 5," "))
    rmsd   = stringtofloat(str_get_field(data(1::), 8," "))
    do i = 0,nvars-1,1
      do j = 0,dimsizes(vstr)-1,1
        if (isStrSubset(vstr(j),varnam(i))) then
          merr(iy,i) = (mmod(j)-mobs(j))/mobs(j)
          mrms(iy,i) = rmsd(j)/rmso(j)
        end if
      end do
    end do
  end do 

 ;;;;open a netcdf file and save the data for plot;;;;;
 vbias = dim_avg_n(merr,0)
 vrmse = dim_avg_n(mrms,0)
 sbias = dim_stddev_n(merr,0)
 srmse = dim_stddev_n(mrms,0)

 vbias!0     = "nvars"
 vbias&nvars = ispan(1,nvars,1)
 vbias&nvars@long_name = "number of variables"

 sbias!0     = "nvars"
 sbias&nvars = ispan(1,nvars,1)
 sbias&nvars@long_name = "number of variables"

 vrmse!0     = "nvars"
 vrmse&nvars = ispan(1,nvars,1)
 vrmse&nvars@long_name = "number of variables"
 
 srmse!0     = "nvars"
 srmse&nvars = ispan(1,nvars,1)
 srmse&nvars@long_name = "number of variables"

 ;;;save the data to netcdf file;;;;
 vnmbias = "Bias_"+Groups(igp)
 vnmrmse = "RMSE_"+Groups(igp)
 fout->$vnmbias$=vbias
 fout->$vnmrmse$=vrmse

 vnmbias = "Bias_SD_"+Groups(igp)
 vnmrmse = "RMSE_SD_"+Groups(igp)
 fout->$vnmbias$=sbias
 fout->$vnmrmse$=srmse

 delete([/vbias,vrmse,srmse,sbias/])
end do ; end igp loop

end

