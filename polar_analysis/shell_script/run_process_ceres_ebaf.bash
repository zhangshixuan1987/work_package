#!/bin/bash

# Parameters: override from the environment for another machine/experiment,
# e.g.  DATA_DIR=/other/data bash run_process_ceres_ebaf.bash   (or: sbatch --export=ALL,DATA_DIR=... run_process_ceres_ebaf.bash)
DATA_DIR="${DATA_DIR:-/lcrc/group/e3sm/public_html/diagnostic_output/ac.szhang/v3_polar_analysis/data}"

outdir="${DATA_DIR}/climo"
if [ ! -d ${outdir} ];then 
  mkdir -p ${outdir}
fi 

data_dir="/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology"
syear=2001
eyear=2018
exp="CERES_EBAF"
period="${syear}-${eyear}"

sfiles=
pfiles=
for mon in `seq 1 12`;do 
  mstr=`printf "%02d" $mon`
  sfiles=( ${sfiles[@]} `echo ${data_dir}/ceres_ebaf_surface_v4.1/*_${mstr}_*.nc` )
  pfiles=( ${pfiles[@]} `echo ${data_dir}/ceres_ebaf_toa_v4.1/*_${mstr}_*.nc` )
done 

#for mon in `seq 1 12`;do
#  mstr=`printf "%02d" $mon`
#  sfiles=( ${sfiles[@]} `echo ${data_dir}/ceres_ebaf_surface_v4.0/*_${mstr}_*.nc` )
#  pfiles=( ${pfiles[@]} `echo ${data_dir}/ceres_ebaf_toa_v4.0/*_${mstr}_*.nc` )
#done

echo ${sfiles[@]}
echo ${pfiles[@]}

ens="observation"
outfile="${outdir}/${exp}.${ens}.climo.${period}.nc"
rm -vf ${outfile}.sfc.tmp ${outfile}.toa.tmp ${outfile}
ncrcat -d time,0, ${sfiles[@]} ${outfile}.sfc.tmp
ncrcat -d time,0, ${pfiles[@]} ${outfile}.toa.tmp
mv ${outfile}.sfc.tmp ${outfile}
ncks -h -A ${outfile}.toa.tmp  ${outfile}
rm -rvf ${outfile}.toa.tmp 

ncrename -v "rlds,FLDS"       ${outfile}
ncrename -v "rlus,FLUS"       ${outfile}
ncap2 -s "FLNS=FLUS-FLDS"     ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncrename -v "rldscs,FLDSC"    ${outfile}
ncrename -v "rluscs,FLUSC"    ${outfile}
ncap2 -s "FLNSC=FLUSC-FLDSC"  ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncrename -v "rsds,FSDS"       ${outfile}
ncrename -v "rsus,FSUS"       ${outfile}
ncap2 -s "FSNS=FSDS-FSUS"     ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncrename -v "rsdscs,FSDSC"    ${outfile}
ncrename -v "rsuscs,FSUSC"    ${outfile}
ncap2 -s "FSNSC=FSDSC-FSUSC"  ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncrename -v "rsdt,SOLIN"      ${outfile}
ncrename -v "rsut,FSUTOA"     ${outfile}
ncrename -v "rsutcs,FSUTOAC"  ${outfile}
ncrename -v "rlut,FLUT"       ${outfile}
ncrename -v "rlutcs,FLUTC"    ${outfile}

ncrename -v "toa_net_all_mon,RESTOA"        ${outfile}
ncrename -v "toa_net_clr_t_mon,RESTOAC"     ${outfile}
ncap2 -s "RESCOM=RESTOA"  ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "RESCOMC=RESTOAC"  ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncrename -v "toa_cre_sw_mon,SWCF"           ${outfile}
ncrename -v "toa_cre_lw_mon,LWCF"           ${outfile}
ncap2 -s "NETCF=SWCF+LWCF"  ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "SWCF_SRF=sfc_net_sw_all_mon-sfc_net_sw_clr_t_mon"   ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "LWCF_SRF=sfc_net_lw_all_mon-sfc_net_lw_clr_t_mon"   ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "NETCF_SRF=SWCF_SRF+LWCF_SRF"  ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "FSNTOA=SOLIN-FSUTOA"   ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FSNTOAC=SOLIN-FSUTOAC" ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "FLNT=FLUT"   ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "FLNTC=FLUTC" ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -s "ALBEDO_SRF=FSUS/FSDS"        ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "ALBEDOC_SRF=FSUSC/FSDSC"     ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "ALBEDO=FSUTOA/SOLIN"         ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}
ncap2 -s "ALBEDOC=FSUTOAC/SOLIN"       ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

ncap2 -Oh -s 'defdim("bnds",2);time[time]={15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5,319, 349.5} ; time_bnds[time,bnds]={0, 31, 31, 59, 59, 90, 90, 120, 120, 151, 151, 181, 181, 212, 212, 243, 243, 273, 273, 304, 304, 334, 334, 365.} ; time@units="days since 2001-01-01 00:00:00"; time@calendar="standard"' ${outfile} ${outfile}.tmp
mv ${outfile}.tmp ${outfile}

#ncrename -v "sfc_cre_net_sw_mon,SWCF_SRF"   ${outfile}
#ncrename -v "sfc_cre_net_lw_mon,LWCF_SRF"   ${outfile}
#ncrename -v "sfc_cre_net_tot_mon,NETCF_SRF" ${outfile}
#ncrename -v "sfc_net_lw_all_mon,FLNS"       ${outfile}
#ncrename -v "sfc_net_lw_clr_t_mon,FLNSC"    ${outfile}
#ncrename -v "sfc_net_sw_all_mon,FSNS"       ${outfile}
#ncrename -v "sfc_net_sw_clr_t_mon,FSNSC"    ${outfile}
#ncrename -v "toa_cre_sw_mon,SWCF"           ${outfile}
#ncrename -v "toa_cre_lw_mon,LWCF"           ${outfile}
#ncrename -v "toa_cre_net_mon,NETCF"         ${outfile}
#ncrename -v "toa_net_all_mon,RESTOA"        ${outfile}
#ncrename -v "toa_net_clr_t_mon,RESTOAC"     ${outfile}
#ncrename -v "sfc_lw_down_clr_c_mon,FLDSC"   ${outfile}
#ncrename -v "sfc_lw_up_clr_c_mon,FLUSC"     ${outfile}
#ncrename -v "toa_net_clr_c_mon,RESTOAC"     ${outfile}
#ncrename -v "toa_sw_clr_c_mon,FSDT"        ${outfile}
#ncrename -v "toa_lw_clr_c_mon,FLUCT"       ${outfile}
#ncrename -v "sfc_net_tot_all_mon,FSNET"     ${outfile}
#ncrename -v "sfc_net_tot_clr_t_mon,FSNETC"  ${outfile}
#ncrename -v "sfc_net_lw_clr_c_mon,FLNSC"    ${outfile}
#ncrename -v "sfc_net_sw_clr_c_mon,FSNSC"    ${outfile}
#ncrename -v "sfc_net_tot_clr_c_mon,FLNTC"   ${outfile}
#ncrename -v "sfc_sw_down_clr_c_mon,FSDSC"   ${outfile}     
#ncrename -v "sfc_sw_up_clr_c_mon,FSUSC"     ${outfile} 
