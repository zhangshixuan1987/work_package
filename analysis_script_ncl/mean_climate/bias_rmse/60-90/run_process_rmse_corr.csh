#!/bin/csh

  set BASE_DIR = "/pscratch/sd/z/zhan391/DARPA_project/evaluation/paper_material/nudging_evaluation/fig_data"
  set models   = ("CLIM" "NDGUV" "NDGUVT" "NDGUVTQ" "NDGUVTQ_SRF1" "NDGUVTQ_SRF2")
  set lenm     = $#models

  set var      = ("CLDTOT" "LWP" "IWP" "PRECL" "PRECC" "PSL" "STRESS_MAG" \
                  "TREFHT" "TS" "U200" "U500" "U850" "Z500" \
                  "LHFLX" "SHFLX" "PRECT" "TMQ" "AODVIS" \
                  "FLNT" "FSNT" "FNET" "CRE" "LWCRE" "SWCRE" \
                  "FLUT" "FLDS" "FLNS" "FLUTC" "FLNTC" "FSNS" "FSNTOAC")
  set lenv     = $#var

  set obdat    = ("CERES_EBAF" "ERA5" "ERA5" "ERA5" "ERA5" "ERA5" "ERA5" \
                  "ERA5" "ERA5" "ERA5" "ERA5" "ERA5" "ERA5" \
                  "ERA5" "ERA5" "GPCP" "REMSS"  "ERA5" \
                  "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" \
                  "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" "CERES_EBAF" )
  set leno     = $#obdat

  set outdir = "../fig_data/" 
  if( ! -d $outdir ) then
    mkdir -p  $outdir
  endif

  foreach i (`seq 1 1 $lenm`)

    setenv MODEL_NAME $models[$i]
    setenv INPUT      $BASE_DIR
    setenv OUTPUT     $outdir
    setenv model      $models[$i]
    setenv yst        "2008"
    setenv yed        "2017"
    
    rm -rvf $outdir/${model}_clim_output_*.txt

    foreach j (`seq 1 1 $lenv`)
      setenv OBS_NAME $obdat[$j]
      setenv VAR      $var[$j]
      ncl 1_gen_clim_rmse_corr.ncl 
    end

    echo "Calculate climo mean done..."

  end 

