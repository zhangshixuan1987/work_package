#!/bin/sh

cd ${PWD}

for file in model_true_data model_pred_data;do 
 ncatted -a _FillValue,Nudge_U,m,f,1.0e20 $file.nc $file.nc1
 mv $file.nc1 $file.nc 
 ncatted -a _FillValue,Nudge_V,m,f,1.0e20 $file.nc $file.nc1
 mv $file.nc1 $file.nc
 ncatted -a _FillValue,Nudge_T,m,f,1.0e20 $file.nc $file.nc1
 mv $file.nc1 $file.nc
 ncatted -a _FillValue,Nudge_Q,m,f,1.0e20 $file.nc $file.nc1
 mv $file.nc1 $file.nc

 ncremap -m ../mapping/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc -i $file.nc -o rgd_$file.nc
 rm -rvf $file.nc
done 

