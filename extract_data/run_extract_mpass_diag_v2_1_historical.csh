#!/bin/csh
# Submit this script as : sbatch ./[script-name]
#SBATCH -A m3525
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH -N 2
#SBATCH  --job-name=ncclimo_ctrl
#SBATCH  --output=job%j 

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.csh

set workdir = `pwd`
set syear = 1850
set eyear = 2014

set hpssdir = "/home/g/golaz/E3SMv2_1"
set case    = "v2_1.LR.historical_0101"
set ddir  = "${hpssdir}/${case}"
#zstash ls --hpss=$ddir

set subdir  = "post/analysis/mpas_analysis/ts_1850-2014_climo_1985-2014/timeseries"
set rundir  = $workdir/$case 

if ( ! -d $rundir ) then 
  mkdir -p $rundir 
endif
cd $rundir 

set ddir  = "${hpssdir}/${case}"
set fil1   = "$subdir/moc/mocTimeSeries_????-????.nc"
set fil2   = "$subdir/transport/transport_????-????.nc"
zstash extract --hpss=$ddir -v "$fil1" "$fil2" >& zstash.log & 

exit
