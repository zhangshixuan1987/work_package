#! /bin/bash

   ncl plot_2D_global_map_distr.ncl
   cd ./sub_fig
   cp ../run_crop_pdf.sh .
   sh run_crop_pdf.sh
   cd ..
   pdflatex Fig10.tex 
   sh run_crop_pdf.sh
