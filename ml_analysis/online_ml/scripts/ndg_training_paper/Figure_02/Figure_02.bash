#! /bin/bash
#generate figure data
ncl 1_generate_data.ncl
#generate figure
ncl 2_plot_figure.ncl
#remove white space
sh run_crop_pdf.sh

