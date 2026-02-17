#!/bin/csh

setenv NCARG_ROOT  "/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/rgd_run_se/method4/Gen_Data_SEdycore/Gen_Data_SETNAME_ne30/ncl-6.4.0"

set ncl = ${NCARG_ROOT}/bin/ncl
set wrp = ${NCARG_ROOT}/bin/WRAPIT

$wrp MAKEIC.stub MAKEIC.f90
