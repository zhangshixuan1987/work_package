;-----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  runnam   = "20240411.v3.LR.piControl-spinup.GC2BC.pm-cpu"
  runtim   = "0071-0100"

  expnam   = "20240421.v3-LR.B2G2B-spinup.ne30pg2_r05_icos30.pm-cpu"
  expstr   = "v3.LR.B2G"

  top_dir  = "/global/cfs/cdirs/e3sm/www/zhan391/E3SMv3_dev/"
  run_dir  = top_dir+runnam+"/e3sm_diags/atm_monthly_180x360_aave/model_vs_obs_"+runtim+"/viewer/table-data/"
  diro     = "/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_testings/diag_approach3/"+expnam+"/post/rmse_bar/"
  
  strs     = asciiread("var_list.txt",-1,"string")
  varnam   = str_squeeze(str_get_field(strs, 1, ","))
  varmod   = str_squeeze(str_get_field(strs, 3, ","))
  nvars    = dimsizes(varnam)
  scenario = "piControl"

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

end 

