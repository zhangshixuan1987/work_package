#!/bin/env python

# PJG 10212014 NOW INCLUDES SFTLF FROM
# PJG 02012016 RESURRECTING...
# /obs AND HARDWIRED TEST CASE WHICH
# NEEDS FIXIN
# PJD 171121 Attempting to fix issue with default missing for thetao and
# CMOR Table being wrong

import gc
import glob
import json
import os
import sys
import time
import cdms2
import csv
from os import path

data_path = "/pscratch/sd/z/zhan391/SEACROGS_project/paper_material/rmse_bar/cmip6/"
out_path  = "/pscratch/sd/z/zhan391/SEACROGS_project/paper_material/rmse_bar/"
seas      = ["ANN","DJF","MAM","JJA","SON"]

cmip_path = "./cmip6"
e3sm_path = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS/e3sm_diag/"

# Generate remap dictionary
par_path = "../parameter"
# Open and return JSON object as a dictionary
product_map   = json.load(open(os.path.join(par_path,'custom_model_info.json')))
variable_map  = json.load(open(os.path.join(par_path,'custom_var_info.json')))

#  ;read cmip6 data
#  ;scenario = "historical"
#  ;filname  = run_dir + "cmip6_"+scenario+"_seasonal_rmse_202203.csv"
#  scenario = "amip"
#  filname  = run_dir + "cmip6_"+scenario+"_seasonal_rmse_202206.csv"
#  lines    = asciiread(filname,-1,"string")
#  nlines   = dimsizes(lines)  ; First 3 lines are a header

#  delim    = ","
#  str1     = str_split(lines(0),delim)
#  cmip_var = str1(1:)
#  linex    = str_sub_str(lines(1),"$^{","~S~")
#  liney    = str_sub_str(linex,"}$","~N~")
#  cmip_unt = str_split(liney,delim)
#  cmip_sea = str_split(lines(2),delim)
# ;print(cmip_var + " " + cmip_unt + " " + cmip_sea)

#  cmip_mod = str_get_field(lines(3:),1,delim)
#  nens     = dimsizes(cmip_mod)
#  vrmse    = new((/nvars,nseas,nens/),float)
#  linez    = str_sub_str(lines,"--","-99999")

for mip in product_map:
  for key in product_map[mip]:
    for runnam in product_map[mip][key]:
      expnam = product_map[mip][key].get(runnam,runnam)
      #print(mip,key,runnam,expnam)
      if mip in ["cmip6"]: 
        fpaths = glob.glob(os.path.join(cmip_path,mip+'*'+key+'*.csv'))
        for file in fpaths:
          print(file)
          with open(file, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
            for row in spamreader:
              print(', '.join(row))
      else:
        run_dir = os.path.join(e3sm_path,runnam,"e3sm_diags","atm_monthly_180x360_aave",)
        #"/e3sm_diags/atm_monthly_180x360_aave/model_vs_obs_"+runtim+"/viewer/table-data/"
  
exit()

