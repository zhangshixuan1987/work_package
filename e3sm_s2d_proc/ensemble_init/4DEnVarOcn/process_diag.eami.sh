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
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

source ./create_and_setup_case.sh
source ${my_dart_env}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

is_da_runcase() {
  [[ ${1} =~ ^4denvar(ocn[12]?|lnd)$ ]]
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

map_file="${my_eami_rgdmap}"
homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

BASE_RUNCASE_LC="${my_runcase,,}"

if [[ ${my_eamopt} == "ERA5-Tang1" ]]; then
  my_eam_ic=${my_eam_ic1}
elif [[ ${my_eamopt} == "ERA5-Tang2" ]]; then
  my_eam_ic=${my_eam_ic2}
elif [[ ${my_eamopt} == "CUSTOMIZE" ]]; then
  my_eam_ic=${my_eam_icx:-}
elif [[ ${BASE_RUNCASE_LC} == "4denvarocn" || ${BASE_RUNCASE_LC} == "4denvarocn1" ]]; then
  if [[ ${my_ocnopt,,} == "4denvarocn2" || ( ${BASE_RUNCASE_LC} == "4denvarocn" && ${START_DATE:0:4} -ge 1990 ) ]]; then
    my_eam_ic=$(printf "${my_daeam2}" "${START_DATE}")
    DART_CASE="${my_dacase2}"
  else
    my_eam_ic=$(printf "${my_daeam1}" "${START_DATE}")
    DART_CASE="${my_dacase1}"
  fi
elif [[ ${BASE_RUNCASE_LC} == "4denvarocn2" ]]; then
  my_eam_ic=$(printf "${my_daeam2}" "${START_DATE}")
  DART_CASE="${my_dacase2}"
elif [[ ${BASE_RUNCASE_LC} == "4denvarlnd" ]]; then
  my_eam_ic=$(printf "${my_daeam3}" "${START_DATE}")
  DART_CASE="${my_dacase3}"
elif is_da_runcase "${BASE_RUNCASE_LC}"; then
  die "diagnostic atmosphere IC is not configured for my_runcase=${my_runcase}"
else
  die "option ${my_eamopt} for atmosphere initial condition is not supported"
fi

[[ -r "${my_eam_ic}" ]] || die "diagnostic atmosphere initial condition not readable: ${my_eam_ic}"

BASE_ICFIL="${my_eami_template}"
DART_ICFIL="${my_eam_ic}"

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
                REMOVE='/usr/bin/rm'
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

# target time tag
ATM_DATE_EXT=${START_DATE}-${START_TOD}

echo "`date` -- START EAM.I DIAG"

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

CURRENT_DADIR="${DART_RUNDIR}/diag"
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
fi
cd ${CURRENT_DADIR}

#Reconstruct IC file 
VAR_ALL="U,V,T,Q,PS,CLDLIQ,CLDICE,NUMLIQ,NUMICE,CO,C2H6,C3H8,CH3COCH3,E90,N2OLNZ,NOYLNZ,CH4LNZ,H2OLNZ,DMS,SO2,H2SO4,so4_a1,so4_a2,so4_a3,so4_a5,pom_a1,pom_a3,pom_a4,bc_a1,bc_a3,bc_a4,dst_a1,dst_a3,ncl_a1,ncl_a2,ncl_a3,mom_a1,mom_a2,mom_a3,mom_a4,num_a1,num_a2,num_a3,num_a4,num_a5,O3,OH,HO2,H2O2,CH2O,CH3O2,CH3OOH,NO,NO2,NO3,N2O5,HNO3,HO2NO2,PAN,C2H5O2,C2H5OOH,CH3CHO,CH3CO3,C2H4,ROHO2,ISOP,ISOP_VBS,ISOPO2,C10H16,MVKMACR,MVKO2,soa_a1,soa_a2,soa_a3,SOAG0,SOAG15,SOAG24,SOAG31,SOAG32,SOAG33,SOAG34,SOAG35"

VAR_DIAG="T"

for VAR in ${VAR_DIAG};do 
  tgt_file="${VAR}.base.eam.i.${ATM_DATE_EXT}_SE.nc"
  rgd_file="${VAR}.base.eam.i.${ATM_DATE_EXT}_360x180.nc"
  ncks -O -v "$VAR" "${DART_ICFIL}" ${tgt_file}
  ncap2 -O \
    -s "PS=double(PS)" \
    -s "${VAR}=double(${VAR})" \
    ${tgt_file} ${tgt_file}
  ncrename -O -d .ncol_d,ncol  ${tgt_file}
  ncrename -O -v .area_d,area -v .lat_d,lat -v .lon_d,lon ${tgt_file}
  ncremap -m ${map_file} -i ${tgt_file} -o ${rgd_file}
  for i in `seq 0 ${DART_ENSNUM}`; do
    ENSTR=EN`printf "%02d" ${i}`
    ATM_INITIAL_FILENAME="${DART_RUNDIR}/${DART_CASE}.${ENSTR}.eam.i.${ATM_DATE_EXT}.nc"
    echo $ATM_INITIAL_FILENAME
    tgt_file="${VAR}.${ENSTR}.eam.i.${ATM_DATE_EXT}_SE.nc"
    rgd_file="${VAR}.${ENSTR}.eam.i.${ATM_DATE_EXT}_360x180.nc"
    ncks -O -v "area_d,lat_d,lon_d,PS,$VAR" "${ATM_INITIAL_FILENAME}" ${tgt_file}
    ncks -O -C -x -v "lat,lon,area" ${tgt_file} ${tgt_file}
    ncrename -O -d .ncol_d,ncol  ${tgt_file}
    ncrename -O -v .area_d,area -v .lat_d,lat -v .lon_d,lon ${tgt_file}
    ncremap -m ${map_file} -i ${tgt_file} -o ${rgd_file}
  done
done 

echo "`date` -- END EAM.I DIAG"

exit 0
