varlist   = (/"PRECT","PRECST","PSL","TS","TREFHT","SHFLX","LHFLX", \
              "U10","V10","U850","V850","TAUX","TAUY",\
              "Qn_SRF","OLR","FLDS","FLDSC","FSNS","FLNS","SWCRE","LWCRE", \
              "ETH12", "T850", "CLDTOT","Z500","U200","V200","TMQ","SF"/)
vartags   = (/"prect","precst","psl","ts","trefht","shflx","lhflx",\
              "u10","v10","u850","v850","taux","tauy",\
              "qnet","olr","flds","fldsc","fsns","flns","swcre","lwcre",\
              "eth","th850","cldtot","z500","u200","v200","tmq","zmsf"/)
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
vtag0 = (/"aice_nh","aice_sh","snowdp","moc","maxnum"/)
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

 ;#OBS
  set        = 0
  nens       = 1
  flag       = set+"-"+"Observations"
 ;process each ensemble member and write to namelist
  do j = 0, nens-1
    enstr    = "en"+ sprinti("%02d",j)
    if (varnm.eq."TS") then 
      model_name = "HadSST"
      flm   = systemfunc("ls " + data_dir + "/"+model_name+"/"+model_name+"."+enstr+"."+varnm+".*.nc")
      ystro = ystr
      yendo = yend
    else if (varnm.eq."PRECT") then
      model_name = "GPCP"
      flm = systemfunc("ls " + data_dir + "/"+model_name+"/"+model_name+"."+enstr+"."+varnm+".*.nc")
      ystro = ystr
      yendo = yend
    else if (varnm.eq."SHFLX".or.varnm.eq."LHFLX".or.varnm.eq."QFLX") then 
      model_name = "OAFlux"
      flm = systemfunc("ls " + data_dir + "/" + model_name+"/" + model_name+"." + enstr+"."+varnm+".*.nc")
      ystro = ystr
      yendo = yend
    else 
      model_name = "CERES-OAFlux"
      flm = systemfunc("ls " + data_dir + "/" + model_name+"/" + model_name+"." + enstr+"."+varnm+".*.nc")
      ystro = 2001
      yendo = 2014
      if(.not.isfilepresent(flm)) then
        model_name = "ERA5"
        flm = systemfunc("ls " + data_dir + "/" + model_name+"/" + model_name+"." + enstr+"."+varnm+".*.nc")
        ystro = ystr
        yendo = yend
      end if 
    end if
    end if 
    end if 
    outnm = model_name+"-"+varnm
   ;print(flm)
    if(isfilepresent(flm)) then
      lines = outnm+ " | " + flm + " | " + ystro+ " | " + yendo + " | " + flag
      write_table(nml_obs,"a",[/lines/],"%s")
      write_table("obs_"+vtag,"a",[/lines/],"%s")
      write_table(outfile,"a",[/lines/],"%s")
      delete([/lines/])
    end if 
    delete([/flm/])
  end do 

 ;#Model#
  set        = 1
  model_name = "ERA5"
  nens       = 1
 ;process each ensemble member and write to namelist 
  do j = 0, nens-1
    enstr = "en"+ sprinti("%02d",j)
    if(nens.gt.1) then
      outnm = model_name + "_en" + sprinti("%02d",j)
    else
      outnm = model_name
    end if
    flm   = systemfunc("ls " + data_dir + "/" + model_name+"/" + model_name+"." + enstr+"."+varnm+".*.nc")
    flag  = set+"-"+model_name
   ;print(flm) 
    if(isfilepresent(flm)) then
      lines = outnm+ " | " + flm + " | " + ystr+ " | " + yend + " | " + flag
      write_table(nml_obs,"a",[/lines/],"%s")
      write_table("obs_"+vtag,"a",[/lines/],"%s")
      write_table(outfile,"a",[/lines/],"%s")
      delete([/lines/])
    end if
    delete([/flm/])
    if(i.eq.0) then
      fl0   = data_dir + "/" + model_name+"/" + model_name+"." + enstr
      line0 = outnm+ " | " + fl0 + " | " + ystr+ " | " + yend + " | " + flag
      write_table(nml_mod,"a",[/line0/],"%s")
      delete([/fl0,line0/])
    end if
  end do

 ;#Model#
  set        = 2
  model_name = "ERA5"
  group_nam  = "ens"
  nens       = 10
 ;process each ensemble member and write to namelist
  do j = 0, nens-1
    enstr = "en"+ sprinti("%02d",j+1)
    if(nens.gt.1) then
      outnm = model_name + "_" + enstr
    else
      outnm = model_name
    end if
    flm   = systemfunc("ls " + data_dir + "/" + model_name+"/" + model_name+"." + enstr+"."+varnm+".*.nc")
    flag  = set+"-"+model_name+"-"+str_upper(group_nam)
   ;print(flm)
    if(isfilepresent(flm)) then
      lines = outnm+ " | " + flm + " | " + ystr+ " | " + yend + " | " + flag
      write_table(nml_obs,"a",[/lines/],"%s")
      write_table("obs_"+vtag,"a",[/lines/],"%s")
      write_table(outfile,"a",[/lines/],"%s")
      delete([/lines/])
    end if
    delete([/flm/])
    if(i.eq.0) then
      fl0   = data_dir + "/" + model_name+"/" + model_name+"." + enstr
      line0 = outnm+ " | " + fl0 + " | " + ystr+ " | " + yend + " | " + flag
      write_table(nml_mod,"a",[/line0/],"%s")
      delete([/fl0,line0/])
    end if
  end do

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
 ; read each nml_obs, only keep first line
  nb = asciiread("obs_"+vtag,-1,"string")
  write_table("obs_"+vtag,"w",[/nb(0)/],"%s")
  delete([/na,nb,maxnumobs/])
end do ; var loop
