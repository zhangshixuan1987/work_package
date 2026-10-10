#/bin/csh
set Cases       = "20200910.alpha5_0.F2010.ne30pg2_r05_oECv3.compy"
set ncas        = $#Cases

set SE_data_dir = "/compyfs/zhen797/E3SM_simulations/${Cases}/archive/atm/hist/"
set outdir      = "./data"

set Vars        = ("T","TREFHT","U","V","PSL","PRECC","PRECL","Z3","PS","PHIS")
set nvars       = $#Vars

###duplicate the data######
mkdir -p $outdir/$Cases/cam_h0

ln -sf $SE_data_dir/*cam.h0* $outdir/$Cases/cam_h0/


#################
set script_path = /share/apps/nco/4.7.9/bin #/share/apps/netcdf/4.3.2/intel/15.0.1/bin # nco directory
set MAP_FILE    = /qfs/people/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc

set l_extr_dat  = 'TRUE'
set l_regrid    = 'TRUE' # flag to process horizontal remapping 


set Syr   = 2
set Eyr   = 6

set icas = 1

while ( $icas <= $ncas )
 
 set case  = $Cases[$icas]
 echo $case

 if ( ! -d $outdir/$case )then
  mkdir -p $outdir/$case
 endif


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
 
 endif

 @ icas ++ 
end 
