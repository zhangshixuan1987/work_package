#!/bin/bash

out_dir="/global/cfs/cdirs/e3smdata/simulations/S2S2D/post_process"
climo_dir="/global/cfs/cdirs/e3smdata/simulations/S2S2D/s2d_diag/climo"

run_dir="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_JRA55_FOSIRL"
comp="atm"
grid="180x360_aave"
freq="monthly"

years=$(seq 1980 2018)
members=$(seq 0 9)

for init_mon in 05 11; do
  for fy in 1 2; do

    if [[ ${fy} -eq 1 ]]; then
      lead_start=0
      lead_end=11
    else
      lead_start=12
      lead_end=23
    fi

    tag="init${init_mon}_FY${fy}"
    odir="${climo_dir}/${exp}_${tag}"

    mkdir -p "${odir}/lists" "${odir}/monthly"

    status_file="${odir}/run_post_atm_climo.status"
    outputs_complete=false
    if [[ -f "${status_file}" ]] && [[ "$(tr -d '[:space:]' < "${status_file}")" == "OK" ]]; then
      outputs_complete=true
      for calmon in $(seq -w 1 12); do
        if [[ ! -s "${odir}/monthly/${exp}_${tag}_${calmon}_climo.nc" ]]; then
          outputs_complete=false
          break
        fi
      done
      for season in DJF MAM JJA SON ANN; do
        if [[ ! -s "${odir}/${exp}_${tag}_${season}_climo.nc" ]]; then
          outputs_complete=false
          break
        fi
      done
    fi
    if [[ "${outputs_complete}" == true ]]; then
      echo "Skipping completed climatology: ${tag}"
      continue
    fi

    echo "RUNNING ${SLURM_JOB_ID:-manual}" > "${status_file}"
    tag_failed=false

    echo "======================================"
    echo "Building climatology: ${tag}"
    echo "Lead range: ${lead_start}-${lead_end}"
    echo "Output dir: ${odir}"
    echo "======================================"

    for calmon in $(seq -w 1 12); do

      flist="${odir}/lists/files_${calmon}.txt"
      rm -f "${flist}"

      for yy in ${years}; do

        lead="${yy}${init_mon}0100"
        case="${exp}_${lead}"

        for en in ${members}; do

          ens="EN$(printf "%02d" "${en}")"

          for ilead in $(seq ${lead_start} ${lead_end}); do

            valid_ym=$(date -d "${yy}-${init_mon}-01 +${ilead} month" +%Y-%m)
            valid_y=${valid_ym:0:4}
            valid_m=${valid_ym:5:2}

            if [[ "${valid_m}" == "${calmon}" ]]; then

              fdir="${out_dir}/${case}/${ens}/post/${comp}/${grid}/${freq}"

              ls ${fdir}/*.eam.h0.${valid_y}-${valid_m}.nc >> "${flist}" 2>/dev/null

            fi

          done
        done
      done

      nfiles=$(wc -l < "${flist}")
      echo "${tag} calendar month ${calmon}: ${nfiles} files"

      if [[ ${nfiles} -gt 0 ]]; then

        if ! ncra -O $(cat "${flist}") \
          "${odir}/monthly/${exp}_${tag}_${calmon}_climo.nc"; then
          tag_failed=true
        fi

      else

        echo "WARNING: no files found for ${tag} month ${calmon}"
        tag_failed=true

      fi

    done

    cd "${odir}/monthly" || exit 1

    echo "Creating seasonal and annual climatology files for ${tag}"

    if ! ncra -O \
      "${exp}_${tag}_12_climo.nc" \
      "${exp}_${tag}_01_climo.nc" \
      "${exp}_${tag}_02_climo.nc" \
      "${odir}/${exp}_${tag}_DJF_climo.nc"; then
      tag_failed=true
    fi

    if ! ncra -O \
      "${exp}_${tag}_03_climo.nc" \
      "${exp}_${tag}_04_climo.nc" \
      "${exp}_${tag}_05_climo.nc" \
      "${odir}/${exp}_${tag}_MAM_climo.nc"; then
      tag_failed=true
    fi

    if ! ncra -O \
      "${exp}_${tag}_06_climo.nc" \
      "${exp}_${tag}_07_climo.nc" \
      "${exp}_${tag}_08_climo.nc" \
      "${odir}/${exp}_${tag}_JJA_climo.nc"; then
      tag_failed=true
    fi

    if ! ncra -O \
      "${exp}_${tag}_09_climo.nc" \
      "${exp}_${tag}_10_climo.nc" \
      "${exp}_${tag}_11_climo.nc" \
      "${odir}/${exp}_${tag}_SON_climo.nc"; then
      tag_failed=true
    fi

    if ! ncra -O \
      "${exp}_${tag}_"??"_climo.nc" \
      "${odir}/${exp}_${tag}_ANN_climo.nc"; then
      tag_failed=true
    fi

    if [[ "${tag_failed}" == false ]]; then
      echo "OK" > "${status_file}"
      echo "Finished ${tag}"
    else
      echo "ERROR" > "${status_file}"
      echo "ERROR: climatology generation failed for ${tag}" >&2
    fi

  done
done
