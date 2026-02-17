load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  model_path = "/global/cfs/cdirs/e3sm/www/zhan391/darpa_temporary_data_share/SE_PG2/clim"
  exp_name   = "E3SMv2_CLIM"
  
  out_path   = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/data" 
  out_path   = out_path+"/"+exp_name
  system("if ! test -d " + out_path +" ; then mkdir -p " + out_path + " ; fi")

  varlist    = (/"TUQ","TVQ"/)
  untlist    = (/"kg/m/s","kg/m/s"/)
  vnmlist    = (/"Total (vertically integrated) zonal water flux", \
                 "Total (vertically integrated) meridional water flux"/)
  nvars      = dimsizes(varlist) 

  infiles    = systemfunc("csh -c 'cd " + model_path+"/" + " ; ls *"+exp_name+"*'")
  nfiles     = dimsizes(infiles)

  do i = 0,nfiles-1,1

   file_in = model_path+"/" +infiles(i)
   fin = addfile(file_in,"r")

   file_ou = out_path+"/"+infiles(i)
   system("/bin/rm -f " + file_ou)
   fou = addfile(file_ou,"cw")

   fou->time     = fin->time
   fou->lat      = fin->lat
   fou->lon      = fin->lon
   fou->area     = fin->area

   g    = 9.81                                ; m/s
   hyai = fin->hyai
   hybi = fin->hybi
   p0   = fin->P0
   
   dum_ps = fin->U(time|:,ncol|:,lev|0:0)
   dum_ps(:,:,0) = fin->PS
   
   dp   = dpres_hybrid_ccm (dum_ps,p0,hyai,hybi)  ; Pa [kg/(m s2)]   

   UQ   = fin->U                                ; K  (time,lev,lat,lon)
   VQ   = fin->V
   UQ   = fin->U * fin->Q * dp(:,:,:,0) /g 
   VQ   = fin->V * fin->Q * dp(:,:,:,0) /g
   
   TUQ  = UQ(:,0,:) 
   TVQ  = VQ(:,0,:)
   TUQ  = dim_sum_n_Wrap(UQ, 1)   
   TVQ  = dim_sum_n_Wrap(VQ, 1)       

   TUQ@units = untlist(0)
   TVQ@units = untlist(1)

   TUQ@long_name = vnmlist(0)
   TVQ@long_name = vnmlist(1)

   fou->TUQ = TUQ 
   fou->TVQ = TVQ

   delete([/file_in,file_ou,fin,fou,g,hyai,hybi,p0,dum_ps,dp,UQ,VQ,TUQ,TVQ/])

  end do 

end
exit
