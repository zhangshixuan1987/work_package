;-----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  run_dir  = "/pscratch/sd/z/zhan391/SEACROGS_project/paper_material/fig_data/rmse_bar/cmip6/"
  diro     = "/pscratch/sd/z/zhan391/SEACROGS_project/paper_material/fig_data/rmse_bar/"
  strs     = asciiread("var_list.txt",-1,"string")
  varnam   = str_squeeze(str_get_field(strs, 1, ","))
  nvars    = dimsizes(varnam)
  seas     = (/"ANN","DJF","MAM","JJA","SON"/)
  nseas    = dimsizes(seas)

  ;read cmip6 data 
  ;scenario = "historical"
  ;filname  = run_dir + "cmip6_"+scenario+"_seasonal_rmse_202203.csv"
  scenario = "amip" 
  filname  = run_dir + "cmip6_"+scenario+"_seasonal_rmse_202206.csv"
  lines    = asciiread(filname,-1,"string")
  nlines   = dimsizes(lines)  ; First 3 lines are a header

  delim    = ","
  str1     = str_split(lines(0),delim)
  cmip_var = str1(1:) 
  linex    = str_sub_str(lines(1),"$^{","~S~") 
  liney    = str_sub_str(linex,"}$","~N~")
  cmip_unt = str_split(liney,delim)
  cmip_sea = str_split(lines(2),delim)
 ;print(cmip_var + " " + cmip_unt + " " + cmip_sea)

  cmip_mod = str_get_field(lines(3:),1,delim)
  nens     = dimsizes(cmip_mod)
  vrmse    = new((/nvars,nseas,nens/),float)
  linez    = str_sub_str(lines,"--","-99999")

  do i = 0,nvars-1,1
    do j = 0,nseas-1,1
      k = i*nseas + j 
      m = k + 2
     ;print(k+ " " + cmip_var(k) + " " + cmip_unt(k) + " " + cmip_sea(k))
      vrmse(i,j,:)= tofloat(str_get_field(linez(3:),m,delim))
    end do 
  end do 
  vrmse!0    = "var"
  vrmse!1    = "sea"
  vrmse!2    = "exp"
  cmip_mod!0 = "exp"
  seas!0     = "sea"
  varnam!0   = "var"
  
  outdir = diro
  system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
  fout1 = outdir+"CMIP6.rmse.seas.clim."+scenario+".nc"
  system("rm -rvf " + fout1)
  setfileoption("nc","Format","netcdf4")
  z1 = addfile(fout1,"cw")
  z1->varnam = varnam
  z1->expnam = cmip_mod(0:nens-4)
  z1->season = seas
  z1->rmse   = vrmse(:,:,0:nens-4)

  fout2 = outdir+"E3SM-2-0.rmse.seas.clim."+scenario+".nc"
  system("rm -rvf " + fout2)
  setfileoption("nc","Format","netcdf4")
  z2 = addfile(fout2,"cw")
  z2->varnam = varnam
  z2->expnam = cmip_mod(nens-3)
  z2->season = seas
  z2->rmse   = vrmse(:,:,nens-3)

  fout3 = outdir+"E3SM-2-0.base.rmse.seas.clim."+scenario+".nc"
  system("rm -rvf " + fout3)
  setfileoption("nc","Format","netcdf4")
  z3 = addfile(fout3,"cw")
  z3->varnam = varnam
  z3->expnam = cmip_mod(nens-2)
  z3->season = seas
  z3->rmse   = vrmse(:,:,nens-2)


  fout4 = outdir+"E3SM-2-0.best.rmse.seas.clim."+scenario+".nc"
  system("rm -rvf " + fout4)
  setfileoption("nc","Format","netcdf4")
  z4 = addfile(fout4,"cw")
  z4->varnam = varnam
  z4->expnam = cmip_mod(nens-1)
  z4->season = seas
  z4->rmse   = vrmse(:,:,nens-1)


end 

