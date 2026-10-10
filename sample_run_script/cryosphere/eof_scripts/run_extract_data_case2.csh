#/bin/csh
set Cases       = "20180215.DECKv1b_H1.ne30_oEC.edison"
set ncas        = $#Cases

set SE_data_dir    = "/global/cfs/cdirs/e3sm/wlin/E3SM_simulations/data/${Cases}/post/atm/180x360/"
set outdir      = "/global/cfs/cdirs/m3089/shixuan/run_e3sm_crysphere/eof_scripts_for_wuyin/data"

set Vars        = ("T","TREFHT","U","V","PSL","PRECC","PRECL","Z3","PS","PHIS")
set nvars       = $#Vars

###duplicate the data######
mkdir -p $outdir/$Cases/cam_h0

ln -sf $SE_data_dir/*cam.h0* $outdir/$Cases/cam_h0/


#################
set script_path = /global/u1/z/zender/bin_cori #/share/apps/nco/4.7.9/bin 
set MAP_FILE    = /qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc

set l_extr_dat  = 'TRUE'
set l_regrid    = 'FALSE' # flag to process horizontal remapping 


set Syr   = 1978
set Eyr   = 2014

set icas = 1

while ( $icas <= $ncas )
 
 set case  = $Cases[$icas]
 echo $case

 if ( ! -d $outdir/$case )then
  mkdir -p $outdir/$case
 endif

 cd $outdir/$case 

 if ($l_extr_dat == 'TRUE') then

  rm -rvf tmp1_????.nc  tmp2_????.nc

  set ys = $Syr 
  while ( $ys <= $Eyr ) 
   set ystr  = `printf "%04d" $ys`
   echo $ystr
   set files = `ls ${SE_data_dir}/*.cam.h0*${ystr}*`
 
   echo $files
   echo $Vars
   
   $script_path/ncrcat -v lat,lon,$Vars  $files tmp1_${ystr}.nc

   @ ys ++
  end
 
  rm -rvf $outdir/$case/SE_${Syr}-${Eyr}.nc
  $script_path/ncrcat -O tmp1_????.nc $outdir/$case/SE_${Syr}-${Eyr}.nc
  rm -rvf tmp1_????.nc  tmp2_????.nc

 endif # end of extracting data on SE grid

 if ($l_regrid == 'TRUE') then

  rm -rvf $outdir/$case/FV_${Syr}-${Eyr}.nc
 
  $script_path/ncremap -i $outdir/$case/SE_${Syr}-${Eyr}.nc -m ${MAP_FILE} -o  $outdir/$case/FV_${Syr}-${Eyr}.nc
 
 else
  
  mv $outdir/$case/SE_${Syr}-${Eyr}.nc $outdir/$case/FV_${Syr}-${Eyr}.nc 
 endif

 @ icas ++ 
end 
