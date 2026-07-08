#!bin/csh

# Submit this script as : sbatch ./[script-name]

#SBATCH -A m3089
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -C knl,quad,cache
#SBATCH -o regrid0300
#SBATCH -t 03:00:00

#work dir
set workdir  = `pwd`
echo $workdir

#Script name and path
set script_name = ncrcat
set script_path = /global/common/sw/cray/cnl7/haswell/nco/4.7.9/gcc/8.2.0/unbt25h/bin

#Modify path so that ncclimo can find latest ncra and other nco utilities
set path = ( $script_path  $path )

set CASE_NAME  = AMIPNUDG_F20TRC5-CMIP6_NODRB_DT1800
set start_year = 2001
set end_year   = 2010

#location of regridded annual mean file
#set RUN_FILE_DIR  = /global/cscratch1/sd/zhan391/climate_analysis/data_process/regrid_climo/$CASE_NAME/climo
set RUN_FILE_DIR = /global/cscratch1/sd/zhan391/F20TRC5-CMIP6_201806_10YR_FRE/climo
#/10yr_${ystr1}to${ystr2}yrs/$CASE_NAME

#Output combined files in FV grid
#set OUT_DIR       = /global/cscratch1/sd/zhan391/climate_analysis/data_process/regrid_climo/$CASE_NAME/means

set OUT_DIR       = /global/cscratch1/sd/zhan391/F20TRC5-CMIP6_201806_10YR_FRE/climo/$CASE_NAME


#if ( -d ${OUT_DIR} ) then
#rm -rvf ${OUT_DIR}
#endif 
#mkdir -p ${OUT_DIR} 
#echo ${OUT_DIR}

foreach n (01 02 03 04 05 06 07 08 09 10 11 12 ANN DJF JJA MAM SON)
echo $n
set REGRIDDED_FILES = ` ls ${RUN_FILE_DIR}/10yr*/$CASE_NAME/*_${n}_climo.nc`
echo  $REGRIDDED_FILES

set OUT_FILES       = ${OUT_DIR}/${CASE_NAME}_${n}_means.nc
echo $REGRIDDED_FILES
echo ${OUT_FILES}

set CAT_SCRIPT = $script_path/$script_name
$CAT_SCRIPT -O $REGRIDDED_FILES  ${OUT_FILES}
#mv ${OUT_FILES}  ${workdir}/climo/${CASE_NAME}/

end
exit
