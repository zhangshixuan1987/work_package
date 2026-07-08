#!/bin/bash -el

echo == Start of eam_boundle_cycling.sh ==
date

echo ============================================
NEXT_DATE=`printf "%04d" ${DART_YEAR}`-`printf "%02d" ${DART_MONTH}`-`printf "%02d" ${DART_DAY}`
NEXT_TOD=`printf "%05d" ${DART_SECONDS}`
# Loop over members
for i in `seq 1 ${my_ensnum}`;do
  echo === Starting member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  REF_DIR="${ARCHIVE_DIR}/rest/${NEXT_DATE}-${NEXT_TOD}"
  cd ${CASE_DIR}
  ./xmlchange RUN_STARTDATE="${NEXT_DATE}"
  ./xmlchange START_TOD="${NEXT_TOD}"
  ./xmlchange STOP_OPTION="nhours"
  ./xmlchange STOP_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange RUN_REFDATE="${NEXT_DATE}"
  ./xmlchange RUN_REFTOD="${NEXT_TOD}"
  cp -rp ${REF_DIR}/${CASE_NAME}.eam.i.${NEXT_DATE}-${NEXT_TOD}.nc  ${RUN_DIR}/
done

cd ${my_wkdir}
ex create_and_setup_case.sh <<ex_end
g;export my_dart_cycle=;s;.*;export my_dart_cycle=${DATA_ASSIMILATION_CYCLES};
wq
ex_end

if [[ ${NEXT_DATE} == ${my_dartymde} && ${NEXT_DATE} == ${my_darttode} ]] ; then
  echo "valid time for DA is $DART_YEAR $DART_MONTH $DART_DAY $DART_HOUR (end time )"
else
  sbatch 4_run_dart_eam_cycleda.sh
fi 

#===compress data to reduce size ===
sbatch eam_compress_data.sh ${NEXT_DATE} ${NEXT_TOD}

echo ===================================
