varlist   = (/"PRECT","QFLX","TUQ","TVQ"/)
vartags   = (/"prect","qflx","tuq","tvq"/)
nvars     = dimsizes(varlist)
ystr      = 1979
yend      = 2014

data_dir  = "/global/cfs/cdirs/e3sm/zhan391/data/CVDP_RGD"
nml_dir   = "./namelist_byvar" 
nml_mod   = "./namelist"
nml_obs   = "./namelist_obs"

;clean up old file
do i = 0, nvars-1
  tmpfil = "obs_"+vartags(i)
  system("rm -rvf " + tmpfil)
end do

;generate empty file for non-existing var
system("if ! test -d " + nml_dir +" ; then mkdir " + nml_dir + " ; fi")
vtag0 = (/"snowdp","moc","maxnum"/)
do i = 0, dimsizes(vtag0)-1
  obfile = "obs_"+vtag0(i)
  system("rm -rvf " + obfile +"; touch " + obfile)
  if(vtag0(i).ne."maxnum") then
    nmfile = nml_dir+"/namelist_"+vtag0(i)
    system("rm -rvf " + nmfile +"; touch " + nmfile)
  end if
end do
system("rm -rvf " + nml_mod + "; touch " + nml_mod)
system("rm -rvf " + nml_obs + "; touch " + nml_obs)

;generate namelist file
do i = 0, nvars-1
  varnm = varlist(i)
  vtag  = vartags(i)

  outfile = nml_dir+"/namelist_"+vtag
  system("if ! test -d " + nml_dir +" ; then mkdir " + nml_dir + " ; fi")
  system("rm -rvf " + outfile +"; touch " + outfile)

  ;#OBS#
  set        = 0
  model_name = "ERA5.6hourly"
  nens       = 1
 ;process each ensemble member and write to namelist 
  do j = 0, nens-1
    enstr = "en"+ sprinti("%02d",j)
    outnm = model_name+"-"+varnm  + "_en" + sprinti("%02d",j)
    if(varnm.eq."TUQ") then 
       varnm1 = "VIWVE"
    else if(varnm.eq."TVQ") then
       varnm1 = "VIWVN"
    else if(varnm.eq."TUH") then
       varnm1 = "VITHEE"
    else if(varnm.eq."TVH") then
       varnm1 = "VITHEN"
    else
      varnm1  = varnm
    end if
    end if
    end if 
    end if 
    flag  = set+"-"+"Observations"
    do year = ystr,yend
      flm = systemfunc("ls " + data_dir + "/" + model_name+"/" + model_name+"." + enstr+"."+varnm1+"."+year+"*.nc")
      if(isfilepresent(flm)) then
        lines = outnm+ " | " + flm + " | " + year+ " | " + year + " | " + flag
      else
        lines = outnm+ " | missing | " + year+ " | " + year + " | " + flag
      end if 
      write_table(outfile,"a",[/lines/],"%s")
      write_table(nml_obs,"a",[/lines/],"%s")
      delete([/lines/])
    end do
  end do

  set        = 1
  group_nam  = "amip"
  model_name = "v2.LR"
  ens_name   = "v2.LR" + "." + group_nam+".6hourly"
  nens       = 1
 ;process each ensemble member and write to namelist
  do j = 0,nens-1,1
    enstr    =  "en"+ sprinti("%02d",j)
    outnm    = model_name+"-"+str_upper(group_nam) + "_en" +sprinti("%02d",j)
    flag     = set+"-"+model_name+"-"+str_upper(group_nam)
    do year = ystr,yend
      flm = systemfunc("ls " + data_dir + "/" + ens_name+"/" +ens_name+"." + enstr +"."+varnm+"."+year+"*.nc")
      if(isfilepresent(flm)) then
        lines = outnm+ " | " + flm + " | " + year+ " | " + year + " | " + flag
      else
        lines = outnm+ " | missing | " + year+ " | " + year + " | " + flag
      end if
      write_table(outfile,"a",[/lines/],"%s")
      delete([/lines/])
    end do
    if(i.eq.0) then
      fl0   = data_dir + "/" + ens_name+"/" + model_name+"." + group_nam+"." + enstr 
      line0 = outnm+ " | " + fl0 + " | " + ystr+ " | " + yend + " | " + flag
      write_table(nml_mod,"a",[/line0/],"%s")
      delete([/fl0,line0/])
    end if
  end do
  delete([/nens,ens_name,model_name/])

  ;;E3SM model (set 5);;;;
  set        = 2
  group_nam  = "historical"
  model_name = "v2.LR"
  ens_name   = "v2.LR" + "." + group_nam+".6hourly"
  nens       = 1
 ;process each ensemble member and write to namelist
  do j = 0,nens-1,1
    enstr    =  "en"+ sprinti("%02d",j)
    outnm    = model_name+"-"+str_upper(group_nam) + "_en" +sprinti("%02d",j)
    flag     = set+"-"+model_name+"-"+str_upper(group_nam)
    do year = ystr,yend 
      flm = systemfunc("ls " + data_dir + "/" + ens_name+"/" + ens_name+"." + enstr +"."+varnm+"."+year+"*.nc")
      if(isfilepresent(flm)) then
        lines = outnm+ " | " + flm + " | " + year+ " | " + year + " | " + flag
      else
        lines = outnm+ " | missing | " + year+ " | " + year + " | " + flag
      end if
      write_table(outfile,"a",[/lines/],"%s")
      delete([/lines/])
    end do
    if(i.eq.0) then
      fl0   = data_dir + "/" + ens_name+"/" + model_name+"." + group_nam+"." + enstr
      line0 = outnm+ " | " + fl0 + " | " + ystr+ " | " + yend + " | " + flag
      write_table(nml_mod,"a",[/line0/],"%s")
      delete([/fl0,line0/])
    end if
  end do
  delete([/nens,ens_name,model_name/])

 ;generate the file required by cvdp metrics diag
  na = asciiread(outfile,-1,"string")
  maxnumobs  = 0
  do k = 0,dimsizes(na)-1
   if (isStrSubset(na(k),"Observations").and.isStrSubset(na(k),varnm)) then
     maxnumobs = maxnumobs + 1
   end if
  end do
  if(i.gt.1) then
    numobs = asciiread("obs_maxnum",-1,"integer")
    if(maxnumobs .gt. numobs) then
      asciiwrite("obs_maxnum",maxnumobs)
    end if
    delete([/numobs/])
  else
    asciiwrite("obs_maxnum",maxnumobs)
  end if
  asciiwrite("obs_"+vtag,na(0))
  delete([/na,maxnumobs/])
end do ; var loop
