;-----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

Casdir  = (/"F20TR_ne30pg2_EC30to60E2r2_FREE_5yr_pm-cpu_11-12-23",\
            "F20TR_ne30pg2_EC30to60E2r2_NDGGLBUV_pm-cpu_01-25-24",\
            "F20TR_ne30pg2_EC30to60E2r2_DONUV_TR2009_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_DONUV_TR2009_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETUV_TR2009_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_UNETUV_TR2009_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_VITOUV_TR2009_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_VITOUV_TR2009_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_NOEPC_MMV1UV_TR2009_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_NOEPC_MMV1UV_TR2009_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV1UV_TR2009_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV1UV_TR2009_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV1UV_TRMYRD1_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV1UV_TRMYRD1_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV2UV_TR2009_3hr_pm-cpu_06-06-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV2UV_TR2009_6hr_pm-cpu_06-11-24",\
            "F20TR_ne30pg2_EC30to60E2r2_MMV2UV_TRMYRD1_3hr_pm-cpu_06-06-24"/) ;,\
           ;"F20TR_ne30pg2_EC30to60E2r2_MMV2UV_TRMYRD1_6hr_pm-cpu_06-11-24" /)

Groups  = (/"CTRL",\
            "NDGUV", \
            "DON_3hr","DON_6hr", \
            "UNet_3hr","UNet_6hr", \
            "ViTO_3hr","ViTO_6h",\
            "MMv1_NEPC_3hr","MMv1_NEPC_6hr",\
            "MMv1_3hr","MMv1_6hr",\
            "MMv1_MD_3hr","MMv1_MD_6hr",\
            "MMv2_3hr","MMv2_6hr",\
            "MMv2_MD_3hr"/) ;,"MMv2_MD_6hr"/)
ncases = dimsizes(Casdir)

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

