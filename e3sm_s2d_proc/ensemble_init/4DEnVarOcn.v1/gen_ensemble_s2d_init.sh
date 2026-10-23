#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=e3sm
#SBATCH  --constraint=cpu
#SBATCH  --qos=debug
#SBATCH  --nodes=4
#SBATCH  --time=00:20:00
#SBATCH  --job-name=eam_pert 

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
if [[ ${GENERATOR_PREFLIGHT_ONLY:-false} != true ]]; then
  source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
fi

source ./create_and_setup_case.sh
if [[ ${GENERATOR_PREFLIGHT_ONLY:-false} != true ]]; then
  source "${my_dart_env}"
fi

die() {
  echo "ERROR: $*" >&2
  exit 1
}

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

# ---------------------
# Purpose
# ---------------------
#
# This script is used to generate initial ensembles with DART 
# The perturbations were added to Temperature fields and then the model 
# was run for 15-days to spin up the ensemble 
#
#*******************************************************************************
# Set paths
WORK_DIR=`pwd`
E3SM_ROOT=${my_e3sm_code}
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo} 
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}

BASE_ATMICFIL="${my_eami_template}"
DART_ATMICFIL="${my_refeam_ic}"
BASE_RUNCASE="${my_runcase}"
BASE_RUNCASE_LC="${BASE_RUNCASE,,}"

is_da_runcase() {
  [[ ${1} == "4denvarocn" ]]
}

# Run options
START_DATE=${my_casedate}
START_TOD=${my_casetod}
DART_CASE=${my_dart_case}
DART_ENSNUM=${my_ensnum}
DART_RUNDIR=${my_dart_runpath}
DART_NTASKS=${my_dart_ntask}
DART_ON_PGRID=${my_dart_pgrid}
DATA_ASSIMILATION_CYCLES=0
DATA_ASSIMILATION_WINDOW=${my_dart_window}
ATM_DATE_EXT=${START_DATE}-${START_TOD}

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

if is_da_runcase "${BASE_RUNCASE_LC}"; then
  my_foscase=""
  my_fosdir=""
  my_fosrefy=""
  my_fossimy=""
elif [[ ${my_fosopt} == "JRA55" ]]; then
  my_foscase=${my_foscase1}
  my_fosdir=${my_fosdir1}
  my_fosrefy=${my_fosrefy1}
  my_fossimy=${my_fossimy1} 
elif [[ ${my_fosopt} == "ERA5" ]]; then
  my_foscase=${my_foscase2}
  my_fosdir=${my_fosdir2}
  my_fosrefy=${my_fosrefy2}
  my_fossimy=${my_fossimy2} 
else
  die "option ${my_fosopt} for FOSI restart source is not supported"
fi

if is_da_runcase "${BASE_RUNCASE_LC}"; then
  my_eam_ic=""
elif [[ ${my_eamopt} == "ERA5-Tang1" ]]; then 
  my_eam_ic=${my_eam_ic1}
elif [[ ${my_eamopt} == "ERA5-Tang2" ]]; then
  my_eam_ic=${my_eam_ic2}
elif [[ ${my_eamopt} == "CUSTOMIZE" ]]; then
  my_eam_ic=${my_eam_icx:-}
else
  die "option ${my_eamopt} for atmosphere initial condition is not supported"
fi

if ! is_da_runcase "${BASE_RUNCASE_LC}"; then
  [[ -n ${my_eam_ic} ]] || die "atmosphere initial condition path is empty for my_eamopt=${my_eamopt}"
fi
[[ ${DART_ENSNUM} =~ ^[0-9]+$ ]] || die "DART_ENSNUM must be an integer, got '${DART_ENSNUM}'"
[[ ${START_TOD} =~ ^[0-9]+$ ]] || die "START_TOD must be seconds-of-day, got '${START_TOD}'"

# Ocean-analysis replacement is specific to BruteForce. FOSIRL cases retain
# the ocean restart selected above and do not require my_ocnopt.
if [[ ${BASE_RUNCASE} == "BruteForce" ]]; then
if [[ ${my_ocnopt:-} == "BenKirtman" ]];then

  #############################################
  # from Dr. Ben Kirtman
  # Available calendar year -> FOSIRL model year
  #############################################
  declare -A YEAR_MAP=(
    [1982]=366
    [1983]=360
    [1984]=379
    [1985]=391
    [1986]=386
    [1987]=367
    [1988]=365
    [1989]=374
    [1990]=388
    [1991]=376
    [1992]=363
  )
  START_YEAR=${START_DATE:0:4}
  if [[ -n "${YEAR_MAP[$START_YEAR]}" ]]; then
      MODEL_YEAR=${YEAR_MAP[$START_YEAR]}
  else
      echo "BruteForce: Year ${START_YEAR} not in supported list"
      exit 1
  fi
  BASE_OCNNAM="${my_ocncase1}"
  BASE_OCNDIR="${my_ocndir1}"
  BASE_OCNDATE=${START_DATE}

elif [[ ${my_ocnopt:-} == "LukeVanRoekel1" ]];then

  declare -A YEAR_MAP=(
    [1987]=360
    [1990]=363
    [1992]=365
    [1993]=366
    [1994]=367
    [2001]=374
    [2003]=376
    [2006]=379
    [2013]=386
    [2015]=388
    [2018]=391
  )
  START_YEAR=${START_DATE:0:4}
  if [[ -n "${YEAR_MAP[$START_YEAR]}" ]]; then
      MODEL_YEAR=${YEAR_MAP[$START_YEAR]}
  else
      echo "BruteForce: Year ${START_YEAR} not in supported list"
      exit 1
  fi
  BASE_OCNNAM="${my_ocncase2}"
  BASE_OCNDIR="${my_ocndir2}"
  BASE_OCNDATE=$(printf "%04d" "$MODEL_YEAR")${START_DATE:4}

elif [[ ${my_ocnopt:-} == "LukeVanRoekel2" ]];then

  declare -A YEAR_MAP=(
    [1980]=353 [1981]=354 [1982]=355 [1983]=356 [1984]=357
    [1985]=358 [1986]=359 [1987]=360 [1988]=361 [1989]=362
    [1990]=363 [1991]=364 [1992]=365 [1993]=366 [1994]=367
    [1995]=368 [1996]=369 [1997]=370 [1998]=371 [1999]=372
    [2000]=373 [2001]=374 [2002]=375 [2003]=376 [2004]=377
    [2005]=378 [2006]=379 [2007]=380 [2008]=381 [2009]=382
    [2010]=383 [2011]=384 [2012]=385 [2013]=386 [2014]=387
    [2015]=388 [2016]=389 [2017]=390 [2018]=391 [2019]=392
    [2020]=393 [2021]=394 [2022]=395 [2023]=396 [2024]=397
  )
  START_YEAR=${START_DATE:0:4}
  if [[ -n "${YEAR_MAP[$START_YEAR]}" ]]; then
      MODEL_YEAR=${YEAR_MAP[$START_YEAR]}
  else
      echo "BruteForce: Year ${START_YEAR} not in supported list"
      exit 1
  fi
  BASE_OCNNAM="${my_ocncase3}"
  BASE_OCNDIR="${my_ocndir3}"
  BASE_OCNDATE=$(printf "%04d" "$MODEL_YEAR")${START_DATE:4}

else

  die "BruteForce: option ${my_ocnopt:-<unset>} for ocean analysis initial condition not in supported list"

fi
fi

DA_PROVIDES_RESTARTS=false

if [[ ${BASE_RUNCASE_LC} == "jra55-fosirl" || ${BASE_RUNCASE_LC} == "era5-fosirl" ]]; then

  #match ic file for other component 
  BASE_BCNAM="${my_foscase}"
  BASE_BCDIR="${my_fosdir}"
  BASE_BCYREF1="${my_fosrefy}"
  BASE_BCYREF2="${my_fossimy}"
  BASE_REPLACE_OCN=false
  BASE_ATMICFIL="${my_eami_template}"
  DART_ATMICFIL="${my_eam_ic}"

elif [[ ${BASE_RUNCASE_LC} == "bruteforce" ]]; then

  #match ic file for other component 
  BASE_BCNAM="${my_foscase}"
  BASE_BCDIR="${my_fosdir}"
  BASE_BCYREF1="${START_YEAR}"
  BASE_BCYREF2="${MODEL_YEAR}"
  BASE_REPLACE_OCN=true
  BASE_ATMICFIL="${my_eami_template}"
  DART_ATMICFIL="${my_eam_ic}"

elif is_da_runcase "${BASE_RUNCASE_LC}"; then

  ##############################
  # BASE_RUNCASE = '4DEnVarOcn'
  ##############################
  case ${START_DATE:0:4} in
    198[0-9]) da_subdir="1980_1989" ;;
    199[0-9]) da_subdir="1990_1999" ;;
    200[0-9]) da_subdir="2000-2009" ;;
    201[0-9]) da_subdir="2010-2019" ;;
    *) die "no 4DEnVar source directory mapping for year ${START_DATE:0:4}" ;;
  esac
  BASE_BCNAM="${my_dacase}"
  BASE_BCDIR="${my_daroot}/${da_subdir}"
  BASE_BCYREF1="${START_DATE:0:4}"
  BASE_BCYREF2="${START_DATE:0:4}"
  BASE_REPLACE_OCN=false
  BASE_ATMICFIL="${my_eami_template}"
  DART_ATMICFIL="${BASE_BCDIR}/$(printf "${my_daeam}" "${START_DATE}")"
  DA_PROVIDES_RESTARTS=true
  DART_CASE="${BASE_BCNAM}"

else
  die "BASE_RUNCASE=${BASE_RUNCASE} is not supported"
fi

if [[ ${DA_PROVIDES_RESTARTS} == true ]]; then
  [[ -d "${BASE_BCDIR}" ]] || die "4DEnVar source directory not readable: ${BASE_BCDIR}"
  [[ -r "${DART_ATMICFIL}" ]] || die "4DEnVar atmosphere initial condition not readable: ${DART_ATMICFIL}"
else
  [[ -r "${BASE_ATMICFIL}" ]] || die "template atmosphere initial condition not readable: ${BASE_ATMICFIL}"
  [[ -r "${DART_ATMICFIL}" ]] || die "DART atmosphere initial condition not readable: ${DART_ATMICFIL}"
fi

# ==============================================================================
# standard commands:
# Make sure that this script is using standard system commands
# instead of aliases defined by the user.
# If the standard commands are not in the location listed below,
# The 'force' (-f) options listed are added to commands where they are used.
# The verbose (-v) argument has been separated from these command definitions
# because these commands may not accept it on some systems.  On those systems
# set VERBOSE = ''
# ==============================================================================
# machine-specific dereferencing
# suppress "rm" warnings if wildcard does not match anything
nonomatch=1
case ${my_machine} in
        "compy")
                VERBOSE='-v'
                MOVE='/usr/bin/mv'
                COPY='/usr/bin/cp --preserve=timestamps'
                LINK='/usr/bin/ln -fs'
                LINKV=TRUE
                LIST='/usr/bin/ls'
                REMOVE='/usr/bin/rm -fr'
                LAUNCHCMD="srun --mpi=pmi2 --ntasks=${DART_NTASKS} --kill-on-bad-exit -l --cpu_bind=cores -c 1 -m plane=40"
                ;;
        "pm-cpu")
                VERBOSE='-v'
                MOVE='/usr/bin/mv'
                COPY='/usr/bin/cp --preserve=timestamps'
                LINK='/usr/bin/ln -fs'
                LINKV=TRUE
                LIST='/usr/bin/ls'
                REMOVE='/usr/bin/rm -fr'
                LAUNCHCMD=mpirun.lsf 
                ;;
         *)
                VERBOSE='-v'
                MOVE='/usr/bin/mv'
                COPY='/usr/bin/cp --preserve=timestamps'
                LINK='/usr/bin/ln -fs'
                LINKV=TRUE
                LIST='/usr/bin/ls'
                REMOVE='/usr/bin/rm -fr'
                LAUNCHCMD="srun --mpi=pmi2 --ntasks=${DART_NTASKS} "
                ;;

esac

# ==============================================================================
# Match land,river, ocean and seaice initial conditions
# ==============================================================================
declare -A comp_map=(
  ["lnd"]="elm"
  ["rof"]="mosart"
  ["ocn"]="mpaso"
  ["ice"]="mpassi"
  ["drv"]="cpl"        # or whatever driver maps to in your case
)
# Reference mapping: model_year 353 == calendar_year 1980
# Match other component ic with atmospehre 
START_YEAR=$((10#${START_DATE:0:4}))
START_HOUR=$((10#${START_TOD} / 3600))
START_HOUR=$(printf "%02d" "${START_HOUR}")
REF_YEAR=$(( START_YEAR - 10#$BASE_BCYREF1 + 10#$BASE_BCYREF2 ))
REF_YEAR_PAD=$(printf "%04d" "$REF_YEAR")
REF_DATE="${REF_YEAR_PAD}${START_DATE:4}"
REF_TOD="${START_TOD}"

component_source() {
  local scomp=$1
  local smod=$2
  local rest_tag=$3
  local ref_date_ext=$4

  if [[ ${DA_PROVIDES_RESTARTS} == true && ${scomp} == "ocn" ]]; then
    printf '%s/%s\n' "${BASE_BCDIR}" "$(printf "${my_daocn}" "${START_DATE}")"
  elif [[ ${DA_PROVIDES_RESTARTS} == true && ${scomp} == "lnd" ]]; then
    printf '%s/%s\n' "${BASE_BCDIR}" "$(printf "${my_dalnd}" "${START_DATE}")"
  elif [[ ${DA_PROVIDES_RESTARTS} == true ]]; then
    printf '%s/%s.%s.%s.nc\n' "${BASE_BCDIR}" "${BASE_BCNAM}" "${rest_tag}" "${ref_date_ext}"
  else
    printf '%s/%s-%s/%s.%s.%s.nc\n' \
      "${BASE_BCDIR}" "${REF_DATE}" "${START_TOD}" \
      "${BASE_BCNAM}" "${rest_tag}" "${ref_date_ext}"
  fi
}

# Validate every source before creating or replacing an output directory.
for scomp in lnd rof ocn ice drv; do
  smod="${comp_map[$scomp]}"
  if [[ ${scomp} == "ocn" || ${scomp} == "ice" ]]; then
    REF_DATE_EXT=${REF_DATE}_${REF_TOD}
    REST_TAG="${smod}.rst"
  else
    REF_DATE_EXT=${REF_DATE}-${REF_TOD}
    REST_TAG="${smod}.r"
  fi
  REF_FILENAME=$(component_source "${scomp}" "${smod}" "${REST_TAG}" "${REF_DATE_EXT}")
  [[ -r "${REF_FILENAME}" ]] || die "initial condition file not readable: ${REF_FILENAME}"
done

if [[ ${DA_PROVIDES_RESTARTS} == true ]]; then
  for atm_rest_tag in eam.r eam.rs; do
    ATM_RESTART_SOURCE="${BASE_BCDIR}/${BASE_BCNAM}.${atm_rest_tag}.${ATM_DATE_EXT}.nc"
    [[ -r "${ATM_RESTART_SOURCE}" ]] || die "4DEnVar atmosphere restart not readable: ${ATM_RESTART_SOURCE}"
  done
fi

if [[ ${BASE_REPLACE_OCN} == true ]]; then
  ocnfile=$(printf "${BASE_OCNNAM}" "${BASE_OCNDATE}")
  [[ -r "${BASE_OCNDIR}/${ocnfile}" ]] || die "ocean analysis file not readable: ${BASE_OCNDIR}/${ocnfile}"
fi

if [[ ${GENERATOR_PREFLIGHT_ONLY:-false} == true ]]; then
  echo "Generator preflight passed: case=${BASE_RUNCASE}, date=${START_DATE}"
  exit 0
fi

if [[ -d "${DART_RUNDIR}" ]]; then
  shopt -s nullglob dotglob
  existing_outputs=("${DART_RUNDIR}"/*)
  shopt -u nullglob dotglob
  if (( ${#existing_outputs[@]} > 0 )); then
    [[ ${my_overwrite} == true ]] || die "output directory is not empty: ${DART_RUNDIR}; set my_overwrite=true to replace it"
    ${REMOVE} "${existing_outputs[@]}"
  fi
else
  mkdir -p "${DART_RUNDIR}"
fi

for scomp in lnd rof ocn ice drv; do
  cd "${DART_RUNDIR}"
  smod="${comp_map[$scomp]}"
  echo === E3SM component ${scomp} ${smod}===
  if [[ ${scomp} == "ocn" || ${scomp} == "ice" ]]; then 
     REF_DATE_EXT=${REF_DATE}_${REF_TOD}
     TGT_DATE_EXT=${START_DATE}_${START_TOD}
     REST_TAG="${smod}.rst"
  else
     REF_DATE_EXT=${REF_DATE}-${REF_TOD}
     TGT_DATE_EXT=${START_DATE}-${START_TOD}
     REST_TAG="${smod}.r"
  fi

  INITIAL_FILENAME="${BASE_BCNAM}.${REST_TAG}.${TGT_DATE_EXT}.nc"
  echo "process ${smod} ic: ${INITIAL_FILENAME}"
  REF_FILENAME=$(component_source "${scomp}" "${smod}" "${REST_TAG}" "${REF_DATE_EXT}")
  if [ ! -f "${REF_FILENAME}" ];then
    die "initial condition file not found: ${REF_FILENAME}"
  fi

  if [[ ${scomp} == "ocn" || ${scomp} == "ice" ]]; then
    ${COPY} "${REF_FILENAME}" "${INITIAL_FILENAME}.tmp" || exit 16
    ncks -O --hdr_pad=10000 "${INITIAL_FILENAME}.tmp" "${INITIAL_FILENAME}"
    ncrename -v xtime,xtime.orig "${INITIAL_FILENAME}"
    rm -rvf "${INITIAL_FILENAME}.tmp"
    echo "${START_DATE}_${START_HOUR}:00:00" > "rpointer.${scomp}"
  else
    pwd
    ${COPY} "${REF_FILENAME}" "${INITIAL_FILENAME}" || exit 17
    echo "./${INITIAL_FILENAME}" > "rpointer.${scomp}"
  fi

  if [[ ${scomp} == "ocn" && ${BASE_REPLACE_OCN} == true ]];then
    # update the ocean initial condition with analysis-based file
    echo "ocean data path: ${BASE_OCNDIR}"
    echo "ocean data file: ${BASE_OCNNAM}"
    echo "ocean data time: ${BASE_OCNDATE}"
    ocnfile=$(printf "${BASE_OCNNAM}" "${BASE_OCNDATE}")
    ncks -A -C -x -v xtime ${BASE_OCNDIR}/${ocnfile} ${INITIAL_FILENAME}
    ncatted -O -a history_of_appended_files,global,d,, ${INITIAL_FILENAME}
  fi

done

if [[ ${DA_PROVIDES_RESTARTS} == true ]]; then
  cd "${DART_RUNDIR}"
  for atm_rest_tag in eam.r eam.rs; do
    ATM_RESTART_FILENAME="${BASE_BCNAM}.${atm_rest_tag}.${ATM_DATE_EXT}.nc"
    ATM_RESTART_SOURCE="${BASE_BCDIR}/${ATM_RESTART_FILENAME}"
    echo "process atmosphere restart: ${ATM_RESTART_FILENAME}"
    ${COPY} "${ATM_RESTART_SOURCE}" "${ATM_RESTART_FILENAME}" || exit 20
  done
  echo "./${BASE_BCNAM}.eam.r.${ATM_DATE_EXT}.nc" > rpointer.atm
  echo "./${BASE_BCNAM}.eam.rs.${ATM_DATE_EXT}.nc" >> rpointer.atm
  echo "$(date) -- END 4DEnVar restart staging"
fi

echo "$(date) -- START EAM_DART_PERTURBATION"

# ==============================================================================
# Block 0: Set command environment
# Make sure the DART executables exist or build them if we can't find them.
# The DART input.nml in the model directory IS IMPORTANT during this part
# because it defines what observation types are supported.
# ==============================================================================
targetdir=${DART_ROOT}/models/${DART_MODEL}/work
if [ ! -x ${targetdir}/filter ]; then
   echo ""
   echo "WARNING: executable file 'filter' not found."
   echo "         Looking for: $targetdir/filter "
   echo "         Trying to rebuild all executables for $MODEL now ..."
   echo "         This will be incorrect, if input.nml:preprocess_nml is not correct."
   cd $targetdir 
   ./quickbuild.sh mpi
   if [ ! -x ${targetdir}/filter ]; then
      echo "ERROR: executable file 'filter' not found."
      echo "       Unsuccessfully tried to rebuild: $targetdir/filter "
      echo "       Required DART assimilation executables are not found."
      echo "       Stopping prematurely."
      exit 01
   fi
fi

# ==============================================================================
# Block 1: generate copies for eam.i file 
# ==============================================================================
cd ${DART_RUNDIR}

if [[ ${my_eamopt} == ERA5-Tang* ]]; then
  #Reconstruct IC file 
  VARS="CO,C2H6,C3H8,CH3COCH3,E90,N2OLNZ,NOYLNZ,CH4LNZ,H2OLNZ,DMS,SO2,H2SO4,so4_a1,so4_a2,so4_a3,so4_a5,pom_a1,pom_a3,pom_a4,bc_a1,bc_a3,bc_a4,dst_a1,dst_a3,ncl_a1,ncl_a2,ncl_a3,mom_a1,mom_a2,mom_a3,mom_a4,num_a1,num_a2,num_a3,num_a4,num_a5,O3,OH,HO2,H2O2,CH2O,CH3O2,CH3OOH,NO,NO2,NO3,N2O5,HNO3,HO2NO2,PAN,C2H5O2,C2H5OOH,CH3CHO,CH3CO3,C2H4,ROHO2,ISOP,ISOP_VBS,ISOPO2,C10H16,MVKMACR,MVKO2,soa_a1,soa_a2,soa_a3,SOAG0,SOAG15,SOAG24,SOAG31,SOAG32,SOAG33,SOAG34,SOAG35"
  ncks -O -x -v "$VARS" "${DART_ATMICFIL}" tmp.base.nc
  ncks -O -v "$VARS" "${DART_ATMICFIL}" tmp.conv.nc

  # Rename ncol -> ncol_d 
  if ! ncdump -h tmp.base.nc | grep -Eq '^[[:space:]]*ncol[[:space:]]*='; then
    echo "ERROR: tmp.base.nc does not define dimension 'ncol'. Nothing to rename."
    echo "Dims found:"
    ncdump -h tmp.base.nc | sed -n '/dimensions:/,/variables:/p'
    exit 02
  fi
  ncrename -O -d .ncol,ncol_d tmp.base.nc
  ncrename -O -v .area,area_d -v .lat,lat_d -v .lon,lon_d tmp.base.nc
  ncap2 -O \
    -s 'time=double(time)' \
    -s 'time_bnds=double(time_bnds)' \
    -s 'date=int(date)' \
    -s 'datesec=int(datesec)' \
    -s 'nbdate=int(nbdate)' \
    -s 'nbsec=int(nbsec)' \
    -s 'ndbase=int(ndbase)' \
    -s 'ndcur=int(ndcur)' \
    -s 'nsbase=int(nsbase)' \
    -s 'nscur=int(nscur)' \
    -s 'nsteph=int(nsteph)' \
    -s 'number=int(number)' \
    -s 'CLDICE=double(CLDICE)' \
    -s 'CLDLIQ=double(CLDLIQ)' \
    -s 'PS=double(PS)' \
    -s 'Q=double(Q)' \
    -s 'T=double(T)' \
    -s 'U=double(U)' \
    -s 'V=double(V)' \
    tmp.base.nc tmp.base.nc

  ncap2 -O \
    -s 'ndbase0=ndbase(0); nsbase0=nsbase(0); nbdate0=nbdate(0); nbsec0=nbsec(0);' \
    tmp.base.nc tmp.base.nc
  ncks -O -x -v ndbase,nsbase,nbdate,nbsec tmp.base.nc tmp.base.tmp.nc
  ncrename -O -v ndbase0,ndbase -v nsbase0,nsbase -v nbdate0,nbdate -v nbsec0,nbsec tmp.base.tmp.nc
  mv tmp.base.tmp.nc tmp.base.nc

  for i in `seq 0 ${DART_ENSNUM}`; do
    ENSTR=EN`printf "%02d" ${i}`
    ATM_INITIAL_FILENAME="${DART_RUNDIR}/${DART_CASE}.${ENSTR}.eam.i.${ATM_DATE_EXT}.nc"
    ${COPY} "${BASE_ATMICFIL}" ${ATM_INITIAL_FILENAME}
    # Append ncol_d vars first, then renamed base (usually want base to override)
    ncks -A -C -x -v date,datesec,nbdate,nbsec,ndbase,ndcur,nsbase,nscur,nsteph,number,time,time_bnds,lat,lon,area,lev,ilev,PS,hyai,hybi,hyam,hybm,P0,ICEFRAC,SICTHK,SNOWHICE tmp.conv.nc ${ATM_INITIAL_FILENAME}
    ncks -A -C -x -v ICEFRAC,SICTHK,SNOWHICE,TS,TS1,TS2,TS3,TS4 tmp.base.nc ${ATM_INITIAL_FILENAME}
    ncatted -O -a history_of_appended_files,global,d,, ${ATM_INITIAL_FILENAME}
  done
  rm -rvf tmp.*.nc 

elif [[ ${DA_PROVIDES_RESTARTS} == true ]]; then
  # The 4DEnVar atmosphere IC is already on the target EAM grid.
  for i in $(seq 0 "${DART_ENSNUM}"); do
    ENSTR=EN$(printf "%02d" "${i}")
    ATM_INITIAL_FILENAME="${DART_RUNDIR}/${DART_CASE}.${ENSTR}.eam.i.${ATM_DATE_EXT}.nc"
    ${COPY} "${DART_ATMICFIL}" "${ATM_INITIAL_FILENAME}"
  done

else 
  # simply copy the eam.i initial condition file to template 	
  for i in `seq 0 ${DART_ENSNUM}`; do
    ENSTR=EN`printf "%02d" ${i}`
    ATM_INITIAL_FILENAME="${DART_RUNDIR}/${DART_CASE}.${ENSTR}.eam.i.${ATM_DATE_EXT}.nc"
    ${COPY} "${BASE_ATMICFIL}" ${ATM_INITIAL_FILENAME}
    # Append ncol_d vars first, then renamed base (usually want base to override)
    ncks -A "${DART_ATMICFIL}"  ${ATM_INITIAL_FILENAME}
    ncatted -O -a history_of_appended_files,global,d,, ${ATM_INITIAL_FILENAME}
  done
fi 

# ==============================================================================
# Block 2: Populate a run-time directory with the input needed to run DART.
# ==============================================================================
CURRENT_DADIR="${DART_RUNDIR}/pert"
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
else
  ${REMOVE} ${CURRENT_DADIR}/*
fi
cd ${CURRENT_DADIR}

# members with perturbations 
for i in `seq 1 ${DART_ENSNUM}`;do
  ENSTR=EN`printf "%02d" ${i}`
  ATM_INITIAL_FILENAME="${DART_RUNDIR}/${DART_CASE}.${ENSTR}.eam.i.${ATM_DATE_EXT}.nc"
  if [ ! -f "${ATM_INITIAL_FILENAME}" ];then 
    echo "ERROR: required file missing ${ATM_INITIAL_FILENAME}"
    exit 03
  fi 
  inst_string=`printf _%04d ${i}`
  ATM_DART_FILENAME=${DART_CASE}.eam${inst_string}.i.${ATM_DATE_EXT}.nc
  if [ $LINKV == TRUE ]; then
    echo "Linking $ATM_INITIAL_FILENAME  $ATM_DART_FILENAME"
  fi
  $LINK ${ATM_INITIAL_FILENAME} ${ATM_DART_FILENAME} || exit 04
done

# ==============================================================================
# Block 3: Actually run the assimilation.
#
# DART namelist settings required:
# &perturb_single_instance_nml
#  ens_size               = 50
#  input_files            = ''
#  output_files           = ''
#  output_file_list       = ''
#  perturbation_amplitude = 0.0
#  single_restart_file_in = .false.
# /
# ==============================================================================
cd ${CURRENT_DADIR}
# Put a pared down copy (no comments) of input.nml in this assimilate_eam directory.
# The contents may change from one cycle to the next, so always start from
# the known configuration in the CASE_ROOT directory.
if [ -e "${DART_WORKDIR}/eam_pert_input.nml" ]; then
  ${COPY} ${DART_WORKDIR}/eam_pert_input.nml input.nml  || exit 05
  sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
  sed -i '1,1i\WARNING: Changes to this file will be ignored. \n Edit \$DART_WORKDIR/eam_pert_input.nml instead.\n\n\n'  input.nml
else
  echo "ERROR ... DART required file ${DART_WORKDIR}/eam_pert_input.nml not found ... ERROR"
  exit 06
fi

#This file is needed for DART code after 2025/02/01
if [ -e "${DART_WORKDIR}/qceff_table.csv" ]; then
  ${COPY} ${DART_WORKDIR}/qceff_table.csv qceff_table.csv  || exit 07
else
  echo "ERROR ... DART required file ${DART_WORKDIR}/qceff_table.csv not found ... ERROR"
  exit 08
fi

# Ensure that the input.nml ensemble size matches the number of instances.
# WARNING: the output files contain ALL ensemble members ==> BIG
ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${DART_ENSNUM};
g;num_output_state_members ;s;= .*;= ${DART_ENSNUM};
g;num_output_obs_members ;s;= .*;= ${DART_ENSNUM};
g;eam_use_pgrid ;s;= .*;= ${DART_ON_PGRID};
wq
ex_end

list=`grep '^[ ]*vertical_localization_coord' input.nml`
list=( `echo $list | sed -e "s#[=,']# #g"` )
if [ "${list[1]}" == "SCALEHEIGHT" ]; then
   list1=`grep '^[ ]*vert_normalization_scale_height' input.nml `
   list1=( `echo $list1 | sed -e "s#[=,]##g"` )
   if [ "${list1[1]}" != "1.5" ]; then
      echo "WARNING!  input.nml is not using 1.5 for vert_normalization_scale_height."
      echo "          Use a different value only if you definitely want to. "
   fi
else
   echo "WARNING!  input.nml is not using SCALEHEIGHT for vertical_localization_coord."
   echo "          SCALEHEIGHT is highly recommended for EAM"
fi

# If possible, use the round-robin approach to deal out the tasks.
# This facilitates using multiple nodes for the simultaneous I/O operations.
if [[ -n ${my_task_per_node:-} ]]; then
  sed -i "s#layout.*#layout = 2#" input.nml
  sed -i "s#tasks_per_node.*#tasks_per_node = ${my_task_per_node}#" input.nml
fi

#stage_dart_files 
${COPY} -f ${DART_WORKDIR}/perturb_single_instance ${CURRENT_DADIR} || exit 09

# always needs a eaminput.nc and a eam_phis.nc for geometry information, etc.
MYSTRING=`grep eam_phis_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAM_PHIS=${MYSTRING[1]}
${REMOVE} ${EAM_PHIS}
${LINK} ${BASE_PHIS} ${EAM_PHIS} || exit 10

MYSTRING=`grep eam_template_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAMINPUT=${MYSTRING[1]}
${REMOVE} ${EAMINPUT}
${LINK}  "${DART_RUNDIR}/${DART_CASE}.EN00.eam.i.${ATM_DATE_EXT}.nc" ${EAMINPUT} || exit 11

#Now, Link the grid information files 
if [ ! -f ${BASE_SEMAPS} ] && [ ! -f ${BASE_CSGRID} ]; then
  echo "ERROR ... no mapping file ${homme_map_file}"
  echo "ERROR ... no gridinfo file ${cs_grid_file}"
  echo "ERROR ... must provide either of them"
  exit 12
else
  if [ -f ${BASE_SEMAPS} ]; then
    ${REMOVE} ${homme_map_file}
    ${COPY} -r ${BASE_SEMAPS} ${homme_map_file} || exit 13
  fi
  if [ -f ${BASE_CSGRID} ]; then
    ${REMOVE} ${cs_grid_file}
    ${COPY} -r ${BASE_CSGRID} ${cs_grid_file} || exit 14
  fi
fi

line=( `grep output_file_list input.nml | sed -e "s#[=,'\.]# #g"` )
output_file_list_name=${line[1]}

${LIST} -1 ${DART_CASE}.eam_[0-9][0-9][0-9][0-9].i.${ATM_DATE_EXT}.nc > $output_file_list_name

${LAUNCHCMD} ${CURRENT_DADIR}/perturb_single_instance || exit 15

# clean the temporary run directory
${REMOVE} "${CURRENT_DADIR}"

echo "$(date) -- END EAM_DART_PERTURBATION"

exit 0
