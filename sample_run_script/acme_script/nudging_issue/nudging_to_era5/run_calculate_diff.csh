#!/bin/csh
# Submit this script as : sbatch ./[script-name]

#SBATCH -A m3089
#SBATCH -N 1
##SBATCH -q regular
#SBATCH -q scavenger
#SBATCH -C knl,quad,cache
#SBATCH -o DT1800.cam.h1
#SBATCH -t 04:00:00

#Set the minimal environment
#module purge
#unsetenv LD_LIBRARY_PATH
#module load intel/16.0.2
#module load mvapich2/2.2b

setenv OMP_NUM_THREADS 1

set dtimeList = (1800)
set irS = 1
set irE = 1
set casestr = "AMIP_NDRB_NUDG_F20TRC5-CMIP6"
set rgdstr  = "AMIP_NDRB_NUDG_F20TRC5-CMIP6"

#Script name and path
set script_name1 = ncremap
set script_name2 = ncea #ncrcat
set script_path = /share/apps/nco/4.7.9/bin

set path = ( $script_path  $path )

set CLIMO_SCRIPT1 = $script_path/$script_name1
set CLIMO_SCRIPT2 = $script_path/$script_name2

#Mapping file
set MAP_FILE      = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/run/map_ne30np4_to_0.9x1.25_aave_110121.nc
#location of model history file
set RUN_FILE_DIR  = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/run
#Output climo files in FV grid
set REGRIDDED_DIR = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/mean_diff

if ( ! -d $REGRIDDED_DIR ) then
  mkdir -p $REGRIDDED_DIR
endif

set NDG1_DIR = /compyfs/zhan391/acme_init/NUDGING_DATA/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800
set NDG2_DIR = /compyfs/zhan391/F20TRC5-CMIP6_NUDG_ISSUE/run/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800

#${CLIMO_SCRIPT2} ${NDG1_DIR}/*cam.h1*2010* ${REGRIDDED_DIR}/NDRB_NUDG1_DATA_F20TRC5-CMIP6_DT1800_2010.nc 
${CLIMO_SCRIPT2} ${NDG2_DIR}/*cam.h1*2010* ${REGRIDDED_DIR}/NDRB_NUDG2_DATA_F20TRC5-CMIP6_DT1800_2010.nc

#set file = ${REGRIDDED_DIR}/NDRB_NUDG1_DATA_F20TRC5-CMIP6_DT1800_2010.nc
#$CLIMO_SCRIPT1 -i ${file} -m ${MAP_FILE} -o ${REGRIDDED_DIR}/RGD_NDRB_NUDG1_DATA_F20TRC5-CMIP6_DT1800_2010.nc

set file = ${REGRIDDED_DIR}/NDRB_NUDG2_DATA_F20TRC5-CMIP6_DT1800_2010.nc
$CLIMO_SCRIPT1 -i ${file} -m ${MAP_FILE} -o ${REGRIDDED_DIR}/RGD_NDRB_NUDG2_DATA_F20TRC5-CMIP6_DT1800_2010.nc


exit

set ndtime = $#dtimeList
set idtime = 1
while ( $idtime <= $ndtime )

  set dtime  = `printf "%04d" $dtimeList[$idtime]`
  set id = $irS
  while ( ${id} <= ${irE} )

   set idstr=`printf "%02d" ${id}`
   set CASE_NAME  = ${casestr}_DT${dtime}

   set ifcst = 1
   cd ${RUN_FILE_DIR}/${CASE_NAME}/

   foreach file (*cam.h1*)
    set nchar = `echo $file | awk '{print length($0)}'` 
    @ stri = $nchar - 25
    echo $stri
    set filestr = `echo $file | awk '{print substr($0,'$stri')}'` 
    set rgdfil = ${rgdstr}.${filestr}
    echo $rgdfil
    rm -rvf ${REGRIDDED_DIR}/${rgdfil}
    $CLIMO_SCRIPT1 -i ${file} -m ${MAP_FILE} -o ${REGRIDDED_DIR}/${rgdfil}
    @ ifcst++
   end  #----------------------------------------------

  @ id++
  end  #----------------------------------------------

@ idtime++
end  #----------------------------------------------

#Other mapping files:
#/lustre/atlas/scratch/bsingh/cli112/map_ne30np4_to_0.9x1.25_aave_110121.nc (this was used prior to 03/06/2016)
#map_ne30_to_0.9x1.25_bilinear.nc
#map_ne30np4_to_fv257x512_bilin.20150901.nc
