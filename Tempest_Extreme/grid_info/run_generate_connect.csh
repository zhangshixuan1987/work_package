#!/bin/csh

#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-knl.csh

set grdnam  = "ne30pg2.g"
set connam  = "outCS_ne30pg2_connect.txt"

# Generate the element mesh.
GenerateCSMesh --alt --res 30 --file ne30.g
#generate a connectivity file 
GenerateConnectivityFile --in_mesh "ne30.g" \
                              --out_type CGLL \
                              --out_np 4 \
                              --out_connect outCS_ne30np4_connect.txt

#ne30pg2:
GenerateVolumetricMesh --in ne30.g \
                       --out $grdnam \
                       --np 2 \
                       --uniform

ConvertExodusToSCRIP --in $grdnam \
                     --out ne30pg2_scrip.nc

#generate a connectivity file 
GenerateConnectivityFile --in_mesh "$grdnam" \
                         --out_type FV \
                         --out_connect $connam
