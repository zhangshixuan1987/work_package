#!/bin/bash -l
#################################################################################
# This script attempts to prepare initial condition files for startup of DART DA 
#################################################################################
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=e3sm_dart_ensda_init 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_dart_ensda_init.%j 
#SBATCH  --exclusive 
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --qos=short

environment_ext=/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
source ${environment_ext}

source create_and_setup_case.sh

my_wkdir=${PWD}

# Run options
REF_DATE=${my_refdate}
REF_TOD=${my_reftod}
REF_CASE=${my_refcase}
REF_DIR=${my_refdir}
REF_HOUR=`echo "${REF_TOD}" \* 3600 | bc`
REF_HOUR=`printf "%02d" $REF_HOUR`
echo "valid time is $REF_DATE $REF_TOD (seconds) $REF_HOUR (hours)"

ARCHIVE_DIR=${my_modeldir}/archive/rest

# ==============================================================================
# machine-specific dereferencing
# suppress "rm" warnings if wildcard does not match anything
VERBOSE='-v'
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
LINKV=TRUE
LIST='/usr/bin/ls'
REMOVE='/usr/bin/rm -fr'
LAUNCHCMD="srun --mpi=pmi2 --ntasks=${DART_NTASKS} --kill-on-bad-exit -l --cpu_bind=cores -c 1 -m plane=40"
# ==============================================================================

# Loop over members
for i in `seq 1 ${my_ensnum}`;do
  echo === Member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  DART_CASE=${my_casename}.${ENSTR}
  CASE_ARCHIVE_DIR=${ARCHIVE_DIR}/${REF_DATE}-${REF_TOD}
  echo "Run Case: ${DART_CASE}"
  echo "Run Directory: ${CASE_ARCHIVE_DIR}"
  if [ ! -d "${CASE_ARCHIVE_DIR}" ];then
    mkdir -p "${CASE_ARCHIVE_DIR}"
  fi
  for scomp in "atm" "lnd" "rof" "ocn" "ice" "drv"; do 
     echo === E3SM component ${scomp} ===
     cd ${CASE_ARCHIVE_DIR}
     if [[ ${scomp} == "atm" ]]; then 
        smod="eam"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        ATM_INITIAL_FILENAME=${DART_CASE}.${smod}.i.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${ATM_INITIAL_FILENAME}"
        if [ -f "${my_refeam_in}" ];then
          ${COPY} ${my_refeam_in} ${ATM_INITIAL_FILENAME} || exit 01
          echo "reference file: ${my_refeam_ic}"
          ncks -A -v PS,U,V,T,Q,CLDLIQ,CLDICE ${my_refeam_ic} ${ATM_INITIAL_FILENAME} || exit 02
          ATM_DATE=( `echo $REF_DATE_EXT | sed -e "s#-# #g"` )
          mdate=`echo ${REF_DATE} | sed "s/-//g"`
          mtods=${REF_TOD}
          ncap2 -O -s "date=${mdate};datesec=${mtods}" ${ATM_INITIAL_FILENAME} ${ATM_INITIAL_FILENAME} || exit 03
        else
          echo "ERROR: initial condition file not found: ${my_refeam_in}"
          exit
        fi
        ATM_REST_FILENAME="${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc"
        echo "${ATM_REST_FILENAME}"   >  rpointer.atm
        if [ ! -f "${ATM_REST_FILENAME}" ];then 
          touch ${ATM_REST_FILENAME}
        fi 
     elif [[ ${scomp} == "lnd" ]]; then
        smod="elm"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        LND_INITIAL_FILENAME=${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${LND_INITIAL_FILENAME}"
        if [ -f "${my_refelm_in}" ];then
          ${COPY} ${my_refelm_in} ${LND_INITIAL_FILENAME} || exit 04
        else
          echo "ERROR: initial condition file not found: ${my_refelm_in}"
          exit
        fi 
        echo "./${LND_INITIAL_FILENAME}"   >  rpointer.lnd
     elif [[ ${scomp} == "rof" ]]; then
        smod="mosart"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        ROF_INITIAL_FILENAME=${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${ROF_INITIAL_FILENAME}"
        if [ -f "${my_refrof_in}" ];then
          ${COPY} ${my_refrof_in} ${ROF_INITIAL_FILENAME} || exit 05
        else
          echo "ERROR: initial condition file not found: ${my_refrof_in}"
          exit
        fi
        echo "./${ROF_INITIAL_FILENAME}"   >  rpointer.rof
     elif [[ ${scomp} == "ocn"  &&  ${my_runtype} == "Full-CPL" ]]; then
        smod="mpaso"
        REF_DATE_EXT=${REF_DATE}_${REF_TOD}
        OCN_INITIAL_FILENAME=${DART_CASE}.${smod}.rst.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME}"
        if [ -f "${my_refocn_in}" ];then
          ${COPY} ${my_refocn_in} ${OCN_INITIAL_FILENAME} || exit 06
          mv ${OCN_INITIAL_FILENAME} tmp.nc 
          ncks -O --hdr_pad=10000 tmp.nc ${OCN_INITIAL_FILENAME} || exit 07
          ncrename -v xtime,xtime.orig ${OCN_INITIAL_FILENAME} || exit 08
          rm -f tmp.nc
        else
          echo "ERROR: initial condition file not found: ${my_refocn_in}"
          exit
        fi
        echo "${REF_DATE}_`printf "%02d" ${REF_HOUR}`:00:00"  > rpointer.ocn
     elif [[ ${scomp} == "ocn"  && ${my_runtype} == "AMIP" ]]; then
        smod="docn"
        REF_DATE_EXT=${REF_DATE}_${REF_TOD}
        OCN_INITIAL_FILENAME1="${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc"
        OCN_INITIAL_FILENAME2="${DART_CASE}.${smod}.rs1.${REF_DATE_EXT}.bin"
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME1}"
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME2}"
        if [ ! -f "${OCN_INITIAL_FILENAME1}" ];then
          touch ${OCN_INITIAL_FILENAME1}
        fi
        if [ ! -f "${OCN_INITIAL_FILENAME2}" ];then
          touch ${OCN_INITIAL_FILENAME2}
        fi
        echo "${OCN_INITIAL_FILENAME1}"  >  rpointer.ocn
        echo "${OCN_INITIAL_FILENAME2}"  >> rpointer.ocn
     elif [[ ${scomp} == "ice" ]]; then
        smod="mpassi"
        REF_DATE_EXT=${REF_DATE}_${REF_TOD}
        ICE_INITIAL_FILENAME=${DART_CASE}.${smod}.rst.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${ICE_INITIAL_FILENAME}"
        if [ -f "${my_refice_in}" ];then
          ${COPY} ${my_refice_in} ${ICE_INITIAL_FILENAME} || exit 09
          # Modify MPAS restart files
          mv ${ICE_INITIAL_FILENAME} tmp.nc
          ncks -O --hdr_pad=10000 tmp.nc ${ICE_INITIAL_FILENAME} || exit 10
          ncrename -v xtime,xtime.orig ${ICE_INITIAL_FILENAME} || exit 11  
          rm -f tmp.nc
        else
          echo "ERROR: initial condition file not found: ${my_refice_in}"
          exit
        fi
        echo "${REF_DATE}_`printf "%02d" ${REF_HOUR}`:00:00"  > rpointer.ice
     else 
        smod="cpl"
        REF_DATE_EXT=${REF_DATE}-${REF_TOD}
        CPL_INITIAL_FILENAME=${DART_CASE}.${smod}.r.${REF_DATE_EXT}.nc
        echo "process ${smod} ic: ${CPL_INITIAL_FILENAME}"
        if [ -f "${my_refcpl_in}" ];then
          ${COPY} ${my_refcpl_in} ${CPL_INITIAL_FILENAME} || exit 12
        else
          echo "ERROR: initial condition file not found: ${my_refcpl_in}"
          exit
        fi
        echo "${CPL_INITIAL_FILENAME}"  >  rpointer.drv
     fi
  done
done

wait
 
exit
