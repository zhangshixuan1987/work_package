#!/bin/bash
# Submit this script as : sbatch ./[script-name]
#SBATCH  --job-name=regrid0300
#SBATCH  --nodes=1
#SBATCH  --time=12:00:00
#SBATCH  --exclusive
#SBATCH -A condo
#SBATCH -p acme-small

source /lcrc/soft/climate/e3sm-unified/load_e3sm_unified_1.9.3_chrysalis.sh

jobdir=`pwd`
cd $jobdir

region="Arctic"
#location of work directory 
WORK_DIR=/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/data/mpas_ts

#Mapping file
MAP_FILE=/lcrc/group/acme/ac.szhang/acme_scratch/data/regrid_maps/map_IcoswISC30E3r4_to_1.0x1.0degree_conserve.nc

start_year=1850
end_year=2050
period="${start_year}01-${end_year}12"
time_tag=`printf "%04d" $start_year`01-`printf "%04d" $end_year`12
time_unt="days since `printf "%04d" $start_year`-01-01 00:00:0.0"

rundir="/lcrc/group/e3sm2/ac.wlin/E3SMv3"
if [[ ${end_year} -eq 2024 ]]; then
   subdir="post/analysis/mpas_analysis/ts_${start_year}-${end_year}_climo_1995-2024/timeseries/mpasTimeSeriesOcean.nc"
else
   subdir="post/analysis/mpas_analysis/ts_1850-2049_climo_2021-2050/timeseries/mpasTimeSeriesOcean.nc"
fi 

fvari="timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature,timeMonthly_avg_volumeCellGlobal,timeMonthly_avg_areaCellGlobal"
fvaro="sst,volumeCell,areaCell"

expstr=( "0051" "0101" "0151" "0201" "0251" "0301" \
	 "0111" "0121" "0131" "0141" "0161" "0171" \
         "0181" "0191" "0211" "0221" "0231" "0241" \
         "0261" "0271" "0281" "0291" "0311" "0321" \
         "0091" )

ensnum=( "en00" "en01" "en02" "en03" "en04" "en05" \
         "en06" "en07" "en08" "en09" "en10" "en11" \
         "en12" "en13" "en14" "en15" "en16" "en17" \
         "en18" "en19" "en20" "en21" "en22" "en23" \
         "en24" )

# Each line: key name
runs=$(cat <<'EOF'
v3.LR.historical.en00 v3.LR.historical_0051
v3.LR.historical.en01 v3.LR.historical_0101
v3.LR.historical.en02 v3.LR.historical_0151
v3.LR.historical.en03 v3.LR.historical_0201
v3.LR.historical.en04 v3.LR.historical_0251
v3.LR.historical.en05 v3.LR.historical_0301
v3.LR.historical.en06 v3.LR.historical_0111
v3.LR.historical.en07 v3.LR.historical_0121
v3.LR.historical.en08 v3.LR.historical_0131
v3.LR.historical.en09 v3.LR.historical_0141
v3.LR.historical.en10 v3.LR.historical_0161
v3.LR.historical.en11 v3.LR.historical_0171
v3.LR.historical.en12 v3.LR.historical_0181
v3.LR.historical.en13 v3.LR.historical_0191
v3.LR.historical.en14 v3.LR.historical_0211
v3.LR.historical.en15 v3.LR.historical_0221
v3.LR.historical.en16 v3.LR.historical_0231
v3.LR.historical.en17 v3.LR.historical_0241
v3.LR.historical.en18 v3.LR.historical_0261
v3.LR.historical.en19 v3.LR.historical_0271
v3.LR.historical.en20 v3.LR.historical_0281
v3.LR.historical.en21 v3.LR.historical_0291
v3.LR.historical.en22 v3.LR.historical_0311
v3.LR.historical.en23 v3.LR.historical_0321
v3.LR.historical.en24 v3.LR.historical_0091
EOF
)

if [ ! -d $WORK_DIR ];then
  mkdir -p $WORK_DIR
fi

while read -r key name; do
 #echo "$key => name=$name, period=$period"
 CASE_NAME=$key 
 exp_name=$name 
 #echo $CASE_NAME $exp_name
 RUN_FILE_DIR=${rundir}/${exp_name}/${subdir}
 if [ ! -f "${RUN_FILE_DIR}" ];then 
    alternate=`echo ${RUN_FILE_DIR} | sed "s/1995-2024/1985-2024/g"` 
    if [ -f "${alternate}" ];then 
      echo $alternate
      RUN_FILE_DIR=${alternate}
    fi 
 fi 
 echo $RUN_FILE_DIR

 # Output file path
 OUT_FILE_DIR="${WORK_DIR}/${CASE_NAME}.mpaso_ts.${period}.nc"
 
 # Remove old files
 rm -vf "${OUT_FILE_DIR}" "${OUT_FILE_DIR}.tmp"

 if [ -f ${RUN_FILE_DIR} ];then
   # Extract selected variables from input NetCDF file
   ncks -7 -O -v "${fvari}" "${RUN_FILE_DIR}" "${OUT_FILE_DIR}"
   ncatted -O -h -a ,global,d,, "${OUT_FILE_DIR}"

   # Rename Time dimension to time
   ncrename -v Time,time -d Time,time "${OUT_FILE_DIR}"
   # --- Split into arrays ---
   IFS=',' read -r -a vi <<< "$fvari"
   IFS=',' read -r -a vo <<< "$fvaro"
   # --- Rename variables one by one ---
   for i in "${!vi[@]}"; do
       echo "Renaming ${vi[$i]} → ${vo[$i]}"
       ncrename -v "${vi[$i]},${vo[$i]}" "${OUT_FILE_DIR}"
   done

   # Replace time with mid-month values (days since start, 365-day year)
   ncap2 -O -s 'time=array(0,1,$time)*365.0/12+15' "${OUT_FILE_DIR}" "${OUT_FILE_DIR}.tmp"

   # Set time units
   ncatted -O -a units,time,m,c,"${time_unt}" "${OUT_FILE_DIR}.tmp"

   # Set calendar attribute
   ncatted -a calendar,time,m,c,"365_day" "${OUT_FILE_DIR}.tmp"
 else 
    echo $RUN_FILE_DIR >> "/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3LE_paper/data/logs/missing.txt"
    echo "missing, try two segements...."
    alternate1=`echo ${RUN_FILE_DIR} | sed "s/ts_1850-2049_climo_2021-2050/ts_1850-2024_climo_1985-2024/g"`
    if [ ! -f "${alternate1}" ];then
      alternate1=`echo ${RUN_FILE_DIR} | sed "s/ts_1850-2049_climo_2021-2050/ts_1850-2024_climo_1995-2024/g"`
    fi
    alternate2=`echo ${RUN_FILE_DIR} | sed "s/ts_1850-2049_climo_2021-2050/ts_2025-2050_climo_2031-2050/g"`
    echo $alternate1
    echo $alternate2
    rm -rvf ${OUT_FILE_DIR}.tmp1 ${OUT_FILE_DIR}.tmp2 
    ncks -7 -O -v "${fvari}" "${alternate1}" "${OUT_FILE_DIR}.tmp1"
    ncks -7 -O -v "${fvari}" "${alternate2}" "${OUT_FILE_DIR}.tmp2"
    ncatted -O -h -a ,global,d,, "${OUT_FILE_DIR}.tmp1"
    ncatted -O -h -a ,global,d,, "${OUT_FILE_DIR}.tmp2"

    ncrename  -v Time,time -d Time,time "${OUT_FILE_DIR}.tmp1"
    ncrename  -v Time,time -d Time,time "${OUT_FILE_DIR}.tmp2"
    # --- Split into arrays ---
    IFS=',' read -r -a vi <<< "$fvari"
    IFS=',' read -r -a vo <<< "$fvaro"
    # --- Rename variables one by one ---
    for i in "${!vi[@]}"; do
        echo "Renaming ${vi[$i]} → ${vo[$i]}"
        ncrename -v "${vi[$i]},${vo[$i]}" "${OUT_FILE_DIR}.tmp1"
        ncrename -v "${vi[$i]},${vo[$i]}" "${OUT_FILE_DIR}.tmp2"
    done

    # Replace time with mid-month values (days since start, 365-day year)
    # Assign synthetic mid-month time values (0-based indexing)
    # Adjust lengths if your files differ; this assumes 2412 months total
    ncap2 -O -s 'time=array(0,1,$time)*365.0/12.0+15'    "${OUT_FILE_DIR}.tmp1" "${OUT_FILE_DIR}.tmp1"
    ncap2 -O -s 'time=array(2100,1,$time)*365.0/12.0+15' "${OUT_FILE_DIR}.tmp2" "${OUT_FILE_DIR}.tmp2"
    # Set time units
    ncatted -O -a units,time,m,c,"${time_unt}" "${OUT_FILE_DIR}.tmp1"
    ncatted -O -a units,time,m,c,"${time_unt}" "${OUT_FILE_DIR}.tmp2"

    # Set calendar attribute
    ncatted -a calendar,time,m,c,"365_day" "${OUT_FILE_DIR}.tmp1"
    ncatted -a calendar,time,m,c,"365_day" "${OUT_FILE_DIR}.tmp2"
    ncks -O --mk_rec_dmn time "${OUT_FILE_DIR}.tmp1" "${OUT_FILE_DIR}.tmp1"
    ncks -O --mk_rec_dmn time "${OUT_FILE_DIR}.tmp2" "${OUT_FILE_DIR}.tmp2"
    ncrcat -O "${OUT_FILE_DIR}.tmp1" "${OUT_FILE_DIR}.tmp2" "${OUT_FILE_DIR}.tmp"
 fi

 # Add time bounds (2-element bounds dimension)
 ncap2 -O -s 'defdim("bnds",2); time_bnds=make_bounds(time,$bnds,"time_bnds")' \
        "${OUT_FILE_DIR}.tmp" "${OUT_FILE_DIR}"

 ncatted -a units,sst,m,c,"degC"        "${OUT_FILE_DIR}"

 # Remove temporary file
 rm -rvf  "${OUT_FILE_DIR}.tmp" "${OUT_FILE_DIR}.tmp1" "${OUT_FILE_DIR}.tmp2"

echo "done"

done <<< "$runs"
