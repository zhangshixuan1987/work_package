#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=24:00:00
#SBATCH  --partition=slurm
#SBATCH  --job-name=e3sm_compress_diag 
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_compress_diag.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

echo == Start of output file compressing ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

VERBOSE='-v'
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
LINKV=TRUE
LIST='/usr/bin/ls'
REMOVE='/usr/bin/rm'
LAUNCHCMD=mpirun.lsf

my_wkdir=${PWD}
scomp="eam"
cd ${my_wkdir}

source ./create_and_setup_case.sh

# Set paths
E3SM_ROOT=${my_e3sm_code}
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo}
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}

CASE_ROOT=${my_modelcase}
CASE_ROOT=${my_modeldir}/case_scripts
RUN_ROOT=${my_modeldir}/run
ARCHIVE_DIR="${my_modeldir}/archive"

#START_YMD="2011-12-09-00000"
#START_YMD="2011-12-10-64800"
#START_YMD="2011-12-12-43200"
#END_YMD="2011-12-13-00000"
START_YMD="2011-12-13-00000"
END_YMD="2011-12-29-00000"

if [ ${my_runtype} == 'AMIP' ];then
  comps=( rest atm cpl lnd rof )
else
  comps=( rest atm cpl ice lnd ocn rof )
fi
START_DATE=( `echo ${START_YMD} | sed -e "s#-# #g"` )
END_DATE=( `echo ${END_YMD} | sed -e "s#-# #g"` )
for year in `seq ${START_DATE[0]} ${END_DATE[0]}`;do 
  for month in `seq ${START_DATE[1]} ${END_DATE[1]}`;do
    for day in `seq ${START_DATE[2]} ${END_DATE[2]}`;do
       timstr=`printf "%04d" $year`-`printf "%02d" $month`-`printf "%02d" $day`
       echo $timstr
       for rundir in ${comps[@]};do
         if [[ ${rundir} != 'rest' ]];then 
           workdir=${ARCHIVE_DIR}/${rundir}/hist  
           if [ -d "${workdir}" ];then 
             cd ${workdir}
             echo $workdir
             xfiles=(${workdir}/*${timstr}*.nc)
             nfiles=${xfiles[0]}
             if [ ! -d "${workdir}/tmp_dir" ] && [ -f "${xfiles[0]}" ];then
                mkdir -p ${workdir}/tmp_dir
                mv ${workdir}/*${timstr}*.nc ${workdir}/tmp_dir/
             fi
             i=0
             for file in ${workdir}/tmp_dir/*.nc ; do
               if [ -f ${file} ] && [ -s ${file} ] ;then
                 fname=`basename ${file}`
                 nccopy -7 -d 5 ${workdir}/tmp_dir/${fname} ${workdir}/${fname} &
                 i=$((i+1))
               fi
               if [[ $i == 40 ]];then
                 wait
                 i=0
               fi
             done 
             wait
             all_complete=0
             for file in ${workdir}/tmp_dir/*.nc ; do
                fname=`basename ${file}`
                if [ ! -f "${workdir}/${fname}" ];then 
                  all_complete=1
                fi 
             done  
             if [[ ${all_complete} == 0 ]];then 
               #now clean up the temporary directory 
               if [ -d "${workdir}/tmp_dir" ];then
                  rm -rvf ${workdir}/tmp_dir
               fi
             fi 
           fi 
         else
           for tod in `seq 0 21600 64800`;do 
             workdir=${ARCHIVE_DIR}/${rundir}/${timstr}-`printf "%05d" ${tod}`
             if [ -d "${workdir}" ];then
               cd ${workdir}
               echo $workdir
               if [ ! -d "${workdir}/tmp_dir" ];then
                  mkdir -p ${workdir}/tmp_dir
                  mv ${workdir}/*${timstr}*.nc ${workdir}/tmp_dir/
               fi
               i=0
               for file in ${workdir}/tmp_dir/*.nc ; do
                 if [ -f ${file} ];then
                   fname=`basename ${file}`
                   nccopy -7 -d 5 ${workdir}/tmp_dir/${fname} ${workdir}/${fname} &
                   i=$((i+1))
                 fi
                 if [[ $i == 40 ]];then
                   wait
                   i=0
                 fi
               done
               wait
               all_complete=0
               for file in ${workdir}/tmp_dir/*.nc ; do
                  fname=`basename ${file}`
                  if [ ! -f "${workdir}/${fname}" ];then
                    all_complete=1
                  fi
               done
               if [[ ${all_complete} == 0 ]];then
                 #now clean up the temporary directory 
                 if [ -d "${workdir}/tmp_dir" ];then
                    rm -rvf ${workdir}/tmp_dir
                 fi
               fi
             fi
           done 
         fi 
       done
    done 
  done 
done 

wait

echo ===== End of output file compressing ==
date
echo =====================================
