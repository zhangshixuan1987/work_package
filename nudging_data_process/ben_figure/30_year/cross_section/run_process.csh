#!/bin/csh

  set BASE_DIR = "/anvil/scratch/x-szhang3/data_240x120"
  set models   = ("CLIM30" "CLIM30ML")

  set lenm     = $#models
  set var      = ("PSL" \
                  "U200" "U500" "U850" \
                  "V200" "V500" "V850" \
                  "T200" "T500" "T850" \
                  "Q200" "Q500" "Q850" ) 
  set lenv     = $#var

  set obdat    = ("ERA5" \
                  "ERA5" "ERA5" "ERA5" \
                  "ERA5" "ERA5" "ERA5" \
                  "ERA5" "ERA5" "ERA5" \
                  "ERA5" "ERA5" "ERA5")
  set leno     = $#obdat

  set vlev     = ("sfc" \
                  "plev" "plev" "plev" \
                  "plev" "plev" "plev" \
                  "plev" "plev" "plev" \
                  "plev" "plev" "plev")
  set lenl     = $#vlev

  set outdir = "./output" 
  if( ! -d $outdir ) then
    mkdir -p  $outdir
  endif

  foreach i (`seq 1 1 $lenm`)
    setenv MODEL_NAME $models[$i]
    setenv INPUT      $BASE_DIR
    setenv model      $models[$i]
    setenv yst        "197901"
    setenv yed        "201412"
    set outfile = "${outdir}/${MODEL_NAME}*_clim_output_$yst-$yed.txt" 
    rm -rvf ${outfile}
    foreach j (`seq 1 1 $lenv`)
      setenv OBS_NAME $obdat[$j]
      setenv VAR      $var[$j]
      setenv VLEV     $vlev[$j]
      ncl 2D_CLIM_RMSE_CORR_CALC.ncl
    end
    echo "Calculate climo mean done..."
  end 
