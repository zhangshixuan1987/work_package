# EAM ensemble members from one HICCUP ERA5 file (sourced by make_atm_ic.bash).
# Adapted from ../BruteForce/luke_cfsr/gen_ensemble_s2d_init.sh (ERA5-Tang
# reconstruction, then DART perturb_single_instance on T); EN00 is unperturbed.

# Settings, as my_* in ../BruteForce/luke_cfsr/create_and_setup_case.sh.
eami_template="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.amip_0101/1980-01-01-00000/v3.LR.amip_0101.eam.i.1980-01-01-00000.nc"
dart_code="/global/homes/z/zhan391/code/DART"
dart_env="${dart_code}/models/eam-se/mach_env/env_pm-cpu_specific.sh"
dart_workdir="${dart_code}/models/eam-se/work"
e3sm_topo="/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
e3sm_semap="/global/homes/z/zhan391/code/HOMME/SEMapping.nc"
e3sm_csgrid="${dart_workdir}/SEMapping_cs_grid_NE30.nc"
dart_pgrid=".true."
dart_tasks_per_node=64

# Inputs eam_ens_members needs besides the HICCUP file (one per line).
eam_ens_inputs() {
  printf '%s\n' "${eami_template}" "${dart_env}" "${dart_workdir}/eam_pert_input.nml" \
    "${dart_workdir}/qceff_table.csv" "${dart_workdir}/perturb_single_instance" "${e3sm_topo}" \
    "${e3sm_semap}" "${e3sm_csgrid}"
}

# Tracer fields kept from the HICCUP file on the template's ncol_d grid.
TRACER_VARS="CO,C2H6,C3H8,CH3COCH3,E90,N2OLNZ,NOYLNZ,CH4LNZ,H2OLNZ,DMS,SO2,H2SO4,so4_a1,so4_a2,so4_a3,so4_a5,pom_a1,pom_a3,pom_a4,bc_a1,bc_a3,bc_a4,dst_a1,dst_a3,ncl_a1,ncl_a2,ncl_a3,mom_a1,mom_a2,mom_a3,mom_a4,num_a1,num_a2,num_a3,num_a4,num_a5,O3,OH,HO2,H2O2,CH2O,CH3O2,CH3OOH,NO,NO2,NO3,N2O5,HNO3,HO2NO2,PAN,C2H5O2,C2H5OOH,CH3CHO,CH3CO3,C2H4,ROHO2,ISOP,ISOP_VBS,ISOPO2,C10H16,MVKMACR,MVKO2,soa_a1,soa_a2,soa_a3,SOAG0,SOAG15,SOAG24,SOAG31,SOAG32,SOAG33,SOAG34,SOAG35"

# Reconstruct the HICCUP file $1 on the eam.i template into member file $2.
eam_reconstruct() {
  local hiccup=$1 out=$2 work
  work=$(mktemp -d "$(dirname "${out}")/.eam_rebuild.XXXX") || return 1
  (
    set -e
    cd "${work}"
    ncks -O -x -v "${TRACER_VARS}" "${hiccup}" base.nc
    ncks -O -v "${TRACER_VARS}" "${hiccup}" conv.nc
    ncrename -O -d .ncol,ncol_d base.nc
    ncrename -O -v .area,area_d -v .lat,lat_d -v .lon,lon_d base.nc
    ncap2 -O \
      -s 'time=double(time)' -s 'time_bnds=double(time_bnds)' \
      -s 'date=int(date)' -s 'datesec=int(datesec)' \
      -s 'nbdate=int(nbdate)' -s 'nbsec=int(nbsec)' -s 'ndbase=int(ndbase)' \
      -s 'ndcur=int(ndcur)' -s 'nsbase=int(nsbase)' -s 'nscur=int(nscur)' \
      -s 'nsteph=int(nsteph)' -s 'number=int(number)' \
      -s 'CLDICE=double(CLDICE)' -s 'CLDLIQ=double(CLDLIQ)' -s 'PS=double(PS)' \
      -s 'Q=double(Q)' -s 'T=double(T)' -s 'U=double(U)' -s 'V=double(V)' \
      base.nc base.nc
    ncap2 -O -s 'ndbase0=ndbase(0); nsbase0=nsbase(0); nbdate0=nbdate(0); nbsec0=nbsec(0);' base.nc base.nc
    ncks -O -x -v ndbase,nsbase,nbdate,nbsec base.nc base.tmp.nc
    ncrename -O -v ndbase0,ndbase -v nsbase0,nsbase -v nbdate0,nbdate -v nbsec0,nbsec base.tmp.nc
    cp --preserve=timestamps "${eami_template}" member.nc
    ncks -A -C -x -v date,datesec,nbdate,nbsec,ndbase,ndcur,nsbase,nscur,nsteph,number,time,time_bnds,lat,lon,area,lev,ilev,PS,hyai,hybi,hyam,hybm,P0,ICEFRAC,SICTHK,SNOWHICE conv.nc member.nc
    ncks -A -C -x -v ICEFRAC,SICTHK,SNOWHICE,TS,TS1,TS2,TS3,TS4 base.tmp.nc member.nc
    ncatted -O -a history_of_appended_files,global,d,, member.nc
    mv -f member.nc "${out}"
  )
  local rc=$?
  rm -rf "${work}"
  return ${rc}
}

# Members EN00-EN<nmem> of HICCUP file $1 in directory $2, named
# <case>.ENxx.eam.i.<stamp>.nc ($3 case, $4 stamp YYYY-MM-DD-SSSSS, $5 nmem).
# EN00 is the reconstruction; DART perturbs T of copies EN01-EN<nmem>. Needs a
# Slurm allocation (srun).
eam_ens_members() {
  local hiccup=$1 dir=$2 case=$3 stamp=$4 nmem=$5 i en pert ntasks
  local en00="${dir}/${case}.EN00.eam.i.${stamp}.nc"
  eam_reconstruct "${hiccup}" "${en00}" || return 1
  for i in $(seq 1 "${nmem}"); do
    en=$(printf "EN%02d" "${i}")
    cp "${en00}" "${dir}/${case}.${en}.eam.i.${stamp}.nc" || return 1
  done
  (( nmem > 0 )) || return 0

  pert="${dir}/.pert"
  rm -rf "${pert}"; mkdir -p "${pert}"
  (
    set -e
    source "${dart_env}"
    cd "${pert}"
    for i in $(seq 1 "${nmem}"); do
      ln -fs "${dir}/${case}.$(printf "EN%02d" "${i}").eam.i.${stamp}.nc" \
        "${case}.eam$(printf "_%04d" "${i}").i.${stamp}.nc"
    done
    cp "${dart_workdir}/eam_pert_input.nml" input.nml
    sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
    sed -i -e "s/^\([ ]*ens_size \).*/\1= ${nmem},/" \
           -e "s/^\([ ]*num_output_state_members \).*/\1= ${nmem},/" \
           -e "s/^\([ ]*num_output_obs_members \).*/\1= ${nmem},/" \
           -e "s/^\([ ]*eam_use_pgrid \).*/\1= ${dart_pgrid},/" \
           -e "s#layout.*#layout = 2#" \
           -e "s#tasks_per_node.*#tasks_per_node = ${dart_tasks_per_node}#" input.nml
    cp "${dart_workdir}/qceff_table.csv" .
    cp -f "${dart_workdir}/perturb_single_instance" .
    nml_value() { grep "^[ ]*$1" input.nml | sed -e "s#[=,']# #g" | awk '{print $2}'; }
    ln -fs "${e3sm_topo}" "$(nml_value eam_phis_filename)"
    ln -fs "${en00}" "$(nml_value eam_template_filename)"
    cp "${e3sm_semap}" SEMapping.nc
    cp "${e3sm_csgrid}" SEMapping_cs_grid.nc
    ls -1 "${case}".eam_[0-9][0-9][0-9][0-9].i."${stamp}".nc > "$(nml_value output_file_list)"
    ntasks=$(( ${SLURM_JOB_NUM_NODES:-1} * dart_tasks_per_node ))
    srun --mpi=pmi2 --ntasks="${ntasks}" ./perturb_single_instance
  )
  local rc=$?
  rm -rf "${pert}"
  return ${rc}
}
