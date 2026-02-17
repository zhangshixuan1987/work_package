#convert .mat file to .nc file 

from jams import mat2nc

dir = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/data/"
fl1 = dir + "beta_drift_linear_regression.mat"
fl2 = dir + "sample_50000.mat"

mat2nc(fl1) 
mat2nc(fl2)

