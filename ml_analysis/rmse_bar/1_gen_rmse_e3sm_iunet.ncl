;-----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

Casdir  = (/"F20TR_ne30pg2_EC30to60E2r2_FREE_5yr_pm-cpu_11-12-23", \
            "F20TR_ne30pg2_EC30to60E2r2_FREE_pm-cpu_10-09-24",\
            "F20TR_ne30pg2_EC30to60E2r2_NDGGLBUV_pm-cpu_10-08-24",\
            "F20TR_ne30pg2_EC30to60E2r2_NDGGLBUVT_pm-cpu_10-09-24/",\
            "F20TR_ne30pg2_EC30to60E2r2_NDGGLBUVTQ_pm-cpu_10-10-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETUV_ALL_TR2009_3hr_pm-cpu_09-27-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETUVT_ALL_TR2009_3hr_pm-cpu_09-27-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETUV_NPBL_TR2009_3hr_pm-cpu_09-27-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETUVT_NPBL_TR2009_3hr_pm-cpu_09-27-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETUV_TR2009_3hr_pm-cpu_09-27-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETUVT_TR2009_3hr_pm-cpu_09-27-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNetUVTQv3_TR2009_3hr_pm-cpu_11-12-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UV_ALL_TR2009_3hr_pm-cpu_10-14-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UVT_ALL_TR2009_3hr_pm-cpu_10-14-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UVTQ_ALL_TR2009_3hr_pm-cpu_10-14-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UV_NPBL_TR2009_3hr_pm-cpu_10-14-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UVT_NPBL_TR2009_3hr_pm-cpu_10-14-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UVTQ_NPBL_TR2009_3hr_pm-cpu_10-14-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UV_TR2009_3hr_pm-cpu_10-14-24", \
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UVT_TR2009_3hr_pm-cpu_10-14-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETv2UVTQ_TR2009_3hr_pm-cpu_10-14-24", \
            "F20TR_ne30pg2_EC30to60E2r2_IUNETUV_UNIV0_3hr_pm-cpu_02-01-25",\
            "F20TR_ne30pg2_EC30to60E2r2_IUNETUVT_UNIV0_3hr_pm-cpu_02-01-25",\
            "F20TR_ne30pg2_EC30to60E2r2_IUNETUVTQ_WRAMIC_3hr_pm-cpu_02-01-25",\
            "F20TR_ne30pg2_EC30to60E2r2_IUNETUVTQ_WRAMIC_IMTNDG_3hr_pm-cpu_02-01-25",\
            "F20TR_ne30pg2_EC30to60E2r2_IUNETv1UVTQ_WRAMIC_IMTNDG_3hr_pm-cpu_02-01-25"\
            /)
Groups  = (/"CTRL_5yr",\
            "CTRL",\
            "NDGUV", \
            "NDGUVT", \
            "NDGUVTQ", \
            "UNet_UV_ALL",\
            "UNet_UVT_ALL",\
            "UNet_UV_NPBL",\
            "UNet_UVT_NPBL",\
            "UNetv0_UV",\
            "UNetv0_UVT",\
            "UNetv0_UVTQ",\
            "UNetv1_UV_ALL",\
            "UNetv1_UVT_ALL",\
            "UNetv1_UVTQ_ALL",\
            "UNetv1_UV_NPBL",\
            "UNetv1_UVT_NPBL",\
            "UNetv1_UVTQ_NPBL",\
            "UNetv1_UV",\
            "UNetv1_UVT",\
            "UNetv1_UVTQ",\
            "IUNET_UV",\
            "IUNET_UVT",\
            "IUNET_UVTQ",\
            "IUNET_UVTQ_IMT",\
            "IUNET_UVTQ_IMT_Scalv0"\
           /)
ncases = dimsizes(Casdir)

print(dimsizes(Groups) + " " + dimsizes(Casdir))
;do k = 0, dimsizes(Casdir)-1,1
; print(Groups(k) + " " + Casdir(k))
;end do 

do ff = 0,ncases-1,1

  runnam   = Casdir(ff) 
  runtim   = "2012-2012"

  expnam   = Casdir(ff)
  expstr   = Groups(ff)

  top_dir  = "/global/cfs/cdirs/e3sm/www/zhan391/SEA_CROGS/e3sm_diag/"
  run_dir  = top_dir+runnam+"/e3sm_diags/atm_monthly_180x360_aave_mva/model_vs_era5_"+runtim+"/viewer/table-data/"
  diro     = "/pscratch/sd/z/zhan391/SEACROGS_project/paper_material/fig_data/rmse_bar/"
  
  strs     = asciiread("var_list.txt",-1,"string")
  varnam   = str_squeeze(str_get_field(strs, 1, ","))
  varmod   = str_squeeze(str_get_field(strs, 3, ","))
  nvars    = dimsizes(varnam)
  scenario = "amip"

  seas     = (/"ANN","DJF","MAM","JJA","SON"/)
  nseas    = dimsizes(seas)

  nens     = 1
  vrmse    = new((/nvars,nseas,nens/),float)

  do i = 0,nseas-1,1
 
   ;read diag data
   filname = run_dir + seas(i) + "_metrics_table.csv"
   lines    = asciiread(filname,-1,"string")
   nlines   = dimsizes(lines)  ; First 3 lines are a header
   delim    = ","

   str1      = str_split(lines(0),delim)
   diag_nam  = str_split(lines(0),delim) 
   indx      = ind(diag_nam.eq."RMSE")
   rmse      = str_get_field(lines(1:),indx+1,delim)
   diag_unt  = str_get_field(lines(2:),2,delim)
   tmp_var   = str_get_field(lines(1:),1,delim)
   diag_var  = new(dimsizes(tmp_var),string)
   diag_reg  = new(dimsizes(tmp_var),string)
   diag_obs  = new(dimsizes(tmp_var),string)
   do kk = 0,dimsizes(tmp_var)-1,1
    tmpstr = str_split(tmp_var(kk)," ")
    diag_var(kk)  = tmpstr(0)
    diag_reg(kk)  = tmpstr(1)
    if(dimsizes(tmpstr).ge.3) then 
      diag_obs(kk)  = tmpstr(2)
    end if 
    delete(tmpstr)
   end do 

   delete([/tmp_var/])
   do j = 0,nvars-1,1
      print(varmod(j))
      indv = ind(diag_var.eq.varmod(j))
      print(diag_var(indv)+ " " +varnam(i)+ " " + varmod(i))
      vrmse(j,i,:)= tofloat(rmse(indv(0)))
      delete([/indv/])
   end do 
   delete([/indx,rmse/])
  end do 

  vrmse!0    = "var"
  vrmse!1    = "sea"
  vrmse!2    = "exp"
  expstr!0  = "exp"
  seas!0     = "sea"
  varnam!0   = "var"
  
  outdir = diro
  system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
  fout1 = outdir+expstr+".rmse.seas.clim."+scenario+".nc"
  system("rm -rvf " + fout1)
  setfileoption("nc","Format","netcdf4")
  z1 = addfile(fout1,"cw")
  z1->varnam = varnam
  z1->expnam = expstr
  z1->season = seas
  z1->rmse   = vrmse(:,:,:)
  delete([/expstr,varnam,seas,vrmse/])
end do 

end 

