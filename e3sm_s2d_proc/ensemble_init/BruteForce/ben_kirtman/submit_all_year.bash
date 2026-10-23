#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=m4849
#SBATCH  --constraint=cpu
#SBATCH  --qos=regular
#SBATCH  --nodes=1
#SBATCH  --time=04:00:00
#SBATCH  --job-name=eam_pert 

work_dir=${PWD}

group="BruteForce"
fosi="JRA55"
atm="ERA5-Tang"
#ocn="BenKirtman"
ocn="LukeVanRoekel"

# default years & months
years=($(seq 1980 2014))
months=(5 11)

if [[ ${ocn} == "BenKirtman" ]]; then
  syear=1982
  eyear=1992
  years=($(seq ${syear} ${eyear}))
  months=(11)
elif [[ ${ocn} == "LukeVanRoekel" ]]; then
  years=(1987 1990 1992 1993 1994 2001 2003 2006 2013 2015 2018)
  months=(11)
fi

years=(1987)
months=(11)

cfg="${work_dir}/create_and_setup_case.sh"

for year in "${years[@]}"; do
  for month in ${months[@]}; do

    cd "${work_dir}" || exit 1

    ymd=$(printf "%04d-%02d-01" "$year" "$month")
    echo "Running ${ymd}"

    sed -i -E \
      -e "s|^export[[:space:]]+my_casedate=.*$|export my_casedate=\"${ymd}\"|" \
      -e "s|^export[[:space:]]+my_runcase=.*$|export my_runcase=\"${group}\"|" \
      -e "s|^export[[:space:]]+my_fosopt=.*$|export my_fosopt=\"${fosi}\"|" \
      -e "s|^export[[:space:]]+my_eamopt=.*$|export my_eamopt=\"${atm}\"|" \
      -e "s|^export[[:space:]]+my_ocnopt=.*$|export my_ocnopt=\"${ocn}\"|" \
      "${cfg}"

    ./gen_ensemble_s2d_init.sh

  done
done

exit 0 
