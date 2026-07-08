#!/bin/bash -el 
my_ensnum=20
my_job_nnodes=20
for i in `seq 1 ${my_ensnum}`;do
  kpaulse="$(echo "scale=1; ( $i * ${my_ensnum} / 10 ) / ${my_job_nnodes}" | bc)"
  if [[ $kpaulse == "1.0"  || $kpaulse == "2.0" ]];then 
    echo $kpaulse
  fi 
done 
