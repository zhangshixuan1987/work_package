#! /bin/bash
   ##extract data#
   #ncl 1_generate_hov_data_meridional_mean.ncl
   #generate figure
   ncl plot_hov_precipitation.ncl
   sh run_crop_pdf.sh
