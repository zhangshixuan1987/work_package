import netCDF4
import sys
import os 
import numpy as np
from netCDF4 import Dataset as NetCDFFile
from netCDF4 import num2date,date2num

if len(sys.argv) > 2:
  syear    = int(sys.argv[1])
  eyear    = int(sys.argv[2])
else:
  print("must specify year range!")
  exit()

topdir   = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share"
outdir   = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS"
exp      = "E3SMv2_NDGUVTQ_SRF1_tau6"

groups   = [ "before_nudging", "nudging_tendency" ]

for grp in groups:
  if grp == "before_nudging":
    smdh = "01010000"
    emdh = "12312100"
    var_list = ["WS_bf_ndg",       "WD_bf_ndg",
                "PRESSURE_bf_ndg", "U_bf_ndg", 
                "V_bf_ndg",        "T_bf_ndg", 
                "Q_bf_ndg",        "PS_bf_ndg"]
    var_outs = ["WS", "WD", "PRESSURE", "U", "V", "T", "Q", "PS"]
  else :
    smdh = "01010300"
    emdh = "01010000"
    var_list = ["Nudge_U", "Nudge_V", "Nudge_T", "Nudge_Q"]
    var_outs = ["UTEND",   "VTEND",   "TEND",    "QTEND"  ]

  datadir  = os.path.join(topdir,"SE_PG2",grp)

  #output path
  outpath  = os.path.join(outdir,"New_Training",exp)
  if not os.path.exists(outpath):
    os.makedirs(outpath)

  for year in range(syear,eyear+1):
    for ii in range(len(var_list)):
      var = var_list[ii]
      vou = var_outs[ii]
      if "Nudge" in var: 
        print("vars are :{}".format(var))
        fname = "{}_3hourly_{:04d}{}-{:04d}{}.nc".format(exp,year,smdh,year+1,emdh)
        f = NetCDFFile(os.path.join(datadir,fname),'r')
      else:
        print("vars are :{}".format(var))  
        fname = "{}_3hourly_{:04d}{}-{:04d}{}.nc".format(exp,year,smdh,year,emdh)
        f = NetCDFFile(os.path.join(datadir,fname),'r')

      lat_val = f.variables['lat'][:]
      lon_val = f.variables['lon'][:]
      lev_val = f.variables['lev'][:]
      tim_val = f.variables['time'][:]
      print('latitude range:',  lat_val.min(),lat_val.max())
      print('longitude range:', lon_val.min(),lon_val.max())
      print('level range:', lev_val.min(),lev_val.max())
      print('time range:', tim_val.min(),tim_val.max())
      
      ax1 = len(tim_val) # itim2 - itim1 + 1
      ax2 = len(lev_val) #
      ax3 = len(lat_val) # ilat2 - ilat1 + 1
      ax4 = len(lon_val) # ilon2 - ilon1 + 1

      if var == "PRESSURE_bf_ndg":
        hyam = f.variables['hyam'][:]
        hybm = f.variables['hybm'][:]
        P0   = f.variables['P0']
        PS   = f.variables['PS_bf_ndg'][:,:]
        data_np = np.zeros(shape=(ax1, ax2, ax3))
        for k in range(len(lev_val)):
           data_np[:,k,:] = hyam[k]*P0 + hybm[k]*PS[:]
        del(hyam,hybm,P0,PS)  
      elif var == "PS_bf_ndg":
        data_np = np.zeros(shape=(ax1,ax3))
        data_np[:,:] = np.array(f.variables[var][:,:])
      elif var in [ "WS_bf_ndg", "WD_bf_ndg" ]:
        data_np = np.zeros(shape=(ax1, ax2, ax3))
        U = f.variables['U_bf_ndg']
        V = f.variables['V_bf_ndg']

        U_nans = U[:,:,:]
        V_nans = V[:,:,:]

        # uncomment below if FillValue exists in your data 
        #Replace _FillValues with NaNs:
        #if hasattr(U, "_FillValue"):
        #  _FillValueU = U._FillValue
        #else:
        #  _FillValueU = netCDF4.default_fillvals 
        #if hasattr(V, "_FillValue"):
        #  _FillValueV = V._FillValue
        #else:
        #  _FillValueV = netCDF4.default_fillvals 
        #U_nans[U_nans == _FillValueU] = np.nan
        #V_nans[V_nans == _FillValueV] = np.nan
        #del(_FillValueU,_FillValueV)

        #Calculate wind speed and wind direction
        ######################################################################
        #To check your software, compute atan2(1,-1)
        #If it equals 2.36 radians (135 degrees) then 
        #your software uses the programming language convention
        #and you can use the formulas below for wind direction calculation 
        if round(np.arctan2(1,-1),2) == 2.36:
          print ( "atan2(1,-1) = 135 degrees")        
          print ( "Wind Direction Calculation Verified..." )
        else:
          print ( "atan2(1,-1) = ", np.arctan2(1,-1))  
          print ( "Wind Direction Calculation Incorrect, please check ..." ) 
        #######################################################################
        WS = np.sqrt(U_nans**2+V_nans**2)
        #######################################################################
        #The wind direction formula uses the two-argument arctangent function, atan2(y,x), 
        #returns the arctangent of y/x in the range -π to π radians (-180 to 180 degrees).          
        #C, C++,  Python, Fortran, Java, IDL, MATLAB and R all follow this convention.
        
        #Calculate wind direction in radians:
        WD = np.arctan2(V_nans,U_nans)
        #Calculate wind direction in degree (180/pi = 57.29578)
        #Same as above but with a constant covert factor DperR
        #DperR = 180.0 / np.pi  
        #WD = np.arctan2(V_nans,U_nans) * DperR
        
        #To convert from WS and WD (in radians) to U and V:
        #U = -WS * sin(WD)  # WD/DperR if WD in Degree 
        #V = -WS * cos(WD)  # WD/DperR if WD in Degree 
        if var == "WS_bf_ndg":
          data_np[:,:,:] = np.array(WS[:,:,:])
        else:
          data_np[:,:,:] = np.array(WD[:,:,:]) 
        del(WS,WD,V_nans,U_nans,U,V)
      else:
        data_np = np.zeros(shape=(ax1, ax2, ax3))
        data_np[:,:,:] = np.array(f.variables[var][:,:,:])
     
      outname = "{}_3hourly_{:04d}.npy".format(vou,year)
      fout    = os.path.join(outpath,outname)

      if os.path.exists(fout):
        os.remove(fout)
      np.save(fout, data_np, allow_pickle=True,fix_imports=True)
      f.close()
      del(data_np,ax1,ax2,ax3,ax4,f,fname)
  del(smdh,emdh,var_list,var_outs)

