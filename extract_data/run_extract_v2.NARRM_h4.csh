#!/bin/csh
#SBATCH --account=e3sm 
#SBATCH -J era5-process
#SBATCH -q flex
#SBATCH -C knl
#SBATCH -N 1
#SBATCH --time=48:00:00      #the max walltime allowed for flex QOS jobs
#SBATCH --time-min=2:00:00   #the minimum amount of time the job should run
#SBATCH --error=%x%j.err 
#SBATCH --output=%x%j.out 

#module load python
#conda activate e3sm_analysis

set jobdir = `pwd`
cd $jobdir

set exp_name   = v2.NARRM.historical 
#set model_case = (v2.NARRM.historical_0101 v2.NARRM.historical_0151 v2.NARRM.historical_0201 v2.NARRM.historical_0251 v2.NARRM.historical_0301)
set model_case = (v2.NARRM.historical_0151 v2.NARRM.historical_0201 v2.NARRM.historical_0251)
set ncase = $#model_case

set var2d_list = ("FLUT" "U200" "PRECT" "OMEGA500")
set nvars      = $#var2d_list

set start_year = 1979
set end_year   = 2014
set time_tag   = `printf "%04d" $start_year`-`printf "%04d" $end_year`

#location of work directory 
set WORK_DIR   = "/global/cfs/cdirs/e3sm/zhan391/extract_data/${exp_name}"

if( ! -d $WORK_DIR )then
 mkdir -p $WORK_DIR
endif

if( ! -d ${WORK_DIR}/SE )then
  mkdir -p ${WORK_DIR}/SE
endif

set i = 1
while ( $i <= $ncase )

 set CASE_NAME  = $model_case[$i]
 set iy = $start_year
 set RUN_FILE_DIR = /global/cfs/projectdirs/e3sm/tang30/E3SMv2/${CASE_NAME}/archive/atm/hist

 while ($iy <= $end_year)

    @ ym1  = $iy - 1
    @ yp1  = $iy + 1

    set eam_files = ($RUN_FILE_DIR/*eam.h5*`printf "%04d" $ym1`-12*)
    set eam_files = ($eam_files $RUN_FILE_DIR/*eam.h5*`printf "%04d" $iy`*)
    if ( $iy < $end_year ) then
      set eam_files = ($eam_files $RUN_FILE_DIR/*eam.h5*`printf "%04d" $yp1`-01*)
    endif

    set j = 1
    while ( $j <= $nvars )
      set var  = $var2d_list[$j]
      @ ens    = $i - 1
      set ensr = en`printf "%02d" $ens`
      set SE_FILE = ${WORK_DIR}/SE/${exp_name}.${ensr}.${var}.${iy}.nc
      rm -rvf $SE_FILE
      set timstr1  = "${iy}-01-01 00:00:0.0"
      set timstr2  = "${iy}-12-31 21:00:0.0"
      ncrcat -d time,"${timstr1}","${timstr2}",1 -v $var $eam_files ${SE_FILE} & 
      exit

     @ j++
    end

    wait

  @ iy ++
 end

 @ i++
end

##finish the process and remove the temp data
rm -rvf ${WORK_DIR}/SE
echo "done"
exit
