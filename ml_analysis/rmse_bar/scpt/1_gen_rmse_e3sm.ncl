;-----------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  run_dir   = "/gpfs/fs0/globalscratch/ac.tzhang/E3SM_archive/"
  diag_dir  = "post/atm/180x360_aave/clim/10yr/"
  obs_dir   = "/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/" 
  diro      = "/lcrc/group/acme/ac.szhang/acme_scratch/amoc_analyis/auto_tuning/fig_data/"
 
  scenario  = "piControl"
  expnam    = "v3tune.LR.historical"
  ystr      = (/000101/)
  yend      = (/001012/)
  tstr      = (/sprinti("%06d",ystr)+"-"+sprinti("%06d",yend)/)
  grid_tag  = "1x1"

  nens        = (45-9+1) + (110-60+1) + 1
  run_name    = new(nens,string)
  ensnum      = new(nens,integer)
  run_name(0) = "cntl"
  ensnum(0)   = 0
  i = 0
  do ie = 9, 45
    i = i +1
    run_name(i) = "workdir."+ie
    ensnum(i)   = ie
  end do
  do ie = 60, 110
    i = i +1
    run_name(i) = "workdir."+ie
    ensnum(i)   = ie
  end do


  strs     = asciiread("var_list.txt",-1,"string")
  varnam   = str_squeeze(str_get_field(strs, 1, ","))
  varunt   = str_squeeze(str_get_field(strs, 2, ","))
  varmod   = str_squeeze(str_get_field(strs, 3, ","))
  region   = str_squeeze(str_get_field(strs, 4, ","))
  obsnam   = str_squeeze(str_get_field(strs, 5, ","))
  nvars    = dimsizes(varnam)

  seas     = (/"ANN","DJF","MAM","JJA","SON"/)
  nseas    = dimsizes(seas)

  vrmse      = new((/nvars,nseas,nens/),float)
  vrmse!0    = "var"
  vrmse!1    = "sea"
  vrmse!2    = "exp"
  run_name!0 = "exp"
  ensnum!0   = "exp"
  seas!0     = "sea"
  varnam!0   = "var"

  do i = 0,nens-1,1
    do j = 0,nseas-1,1
      do k = 0,nvars-1,1
        fl0 = systemfunc("ls " + obs_dir+obsnam(k)+"/"+obsnam(k)+"*"+seas(j)+"*.nc")
        print(fl0)
        fl1 = systemfunc("ls " + run_dir+run_name(i)+"/"+diag_dir+"/*"+seas(j)+"*.nc")
        print(fl1)
        f0  = addfile(fl0,"r")
        f1  = addfile(fl1,"r")
        print(varnam(k) + " " + varmod(k))
        if(varmod(k).eq."RESTOM".or.varmod(k).eq."RESTOA") then
          vob  = (/"toa_net_all_mon"/) 
          vmd  = (/"FSNT","FLUT"/)
          varx = f0->$vob$(0,:,:)
          var1 = f1->$vmd(0)$(0,:,:) 
          var1 = var1 - f1->$vmd(1)$(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          delete([/varx/])
        end if 
        if (varmod(k).eq."LWCF".or.varmod(k).eq."LWCRE") then
          vob  = (/"toa_cre_lw_mon"/)
          vmd  = (/"FLNT","FLNTC"/)
          varx = f0->$vob$(0,:,:)
          var1 = f1->$vmd(0)$(0,:,:)
          var1 = - (var1 - f1->$vmd(1)$(0,:,:))
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          delete([/varx/])
        end if 
        if (varmod(k).eq."SWCF".or.varmod(k).eq."SWCRE") then
          vob  = (/"toa_cre_sw_mon"/)
          vmd  = (/"FSNT","FSNTC"/)
          varx = f0->$vob$(0,:,:)
          var1 = f1->$vmd(0)$(0,:,:)
          var1 = var1 - f1->$vmd(1)$(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          delete([/varx/])
        end if        
        if (varmod(k).eq."PREC".or.varmod(k).eq."PRECT") then
          vob  = (/"PRECT"/)
          vmd  = (/"PRECL","PRECC"/)
          varx = f0->$vob$(0,:,:)
          var1 = f1->$vmd(0)$(0,:,:)
          var1 = var1 + f1->$vmd(1)$(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0) 
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          var1 = var1*1000.0 * 24.0 * 3600.0
          printMinMax(var0,False)
          printMinMax(var1,False)          
          delete([/varx/])
        end if
        if (varmod(k).eq."TAS".or.varmod(k).eq."TREFHT") then
          vob  = (/"tas"/)
          vmd  = (/"TREFHT"/)
          varx = f0->$vob$(0,:,:)
          var1 = f1->$vmd(0)$(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          delete([/varx/])
        end if
        if (varmod(k).eq."SLP".or.varmod(k).eq."PSL") then
          vob  = (/"psl"/)
          vmd  = (/"PSL"/)
          varx = f0->$vob$(0,:,:)
          var1 = f1->$vmd(0)$(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          var0 = var0 / 100.0 
          var1 = var1 / 100.0 
          var0@units = "hPa"
          var1@units = "hPa"
          printMinMax(var0,False)
          printMinMax(var1,False)
          delete([/varx/])
        end if
        if (varmod(k).eq."U200".or.varmod(k).eq."U-200mb") then
          levp = (/200/)
          levp!0 = "levp"
          levp&levp      =  levp
          levp@long_name = "pressure"               ; attach some attributes
          levp@units     = "hPa"
          levp@positive  = "down"
          vob  = (/"ua"/)
          vmd  = (/"U"/)
          varx = f0->$vob$(0,{levp*100},:,:)
          vtmp = vinth2p(f1->$vmd(0)$(0,:,:,:), f1->hyam, f1->hybm,levp, f1->PS(0,:,:), 2, f1->P0*0.01, 1, False )
          var1 = vtmp(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          printMinMax(var0,False)
          printMinMax(var1,False)
          delete([/varx,vtmp/])
        end if
        if (varmod(k).eq."U850".or.varmod(k).eq."U-850mb") then
          levp = (/850/)
          levp!0 = "levp"
          levp&levp      =  levp
          levp@long_name = "pressure"               ; attach some attributes
          levp@units     = "hPa"
          levp@positive  = "down"
          vob  = (/"ua"/)
          vmd  = (/"U"/)
          varx = f0->$vob$(0,{levp*100},:,:)
          vtmp = vinth2p(f1->$vmd(0)$(0,:,:,:), f1->hyam, f1->hybm,levp, f1->PS(0,:,:), 2, f1->P0*0.01, 1, False )
          var1 = vtmp(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          printMinMax(var0,False)
          printMinMax(var1,False)
          delete([/varx,vtmp/])
        end if
        if (varmod(k).eq."Z500".or.varmod(k).eq."Z3-500mb".or.varmod(k).eq."Zg-500") then
          levp = (/500/)
          levp!0 = "levp"
          levp&levp      =  levp
          levp@long_name = "pressure"               ; attach some attributes
          levp@units     = "hPa"
          levp@positive  = "down"
          vob  = (/"zg"/)
          vmd  = (/"Z3"/)
          varx = f0->$vob$(0,{levp*100},:,:)
          vtmp = vinth2p(f1->$vmd(0)$(0,:,:,:), f1->hyam, f1->hybm,levp, f1->PS(0,:,:), 2, f1->P0*0.01, 1, False )
          var1 = vtmp(0,:,:)
          var0 = linint2(varx&lon,varx&lat,varx,True,var1&lon,var1&lat,0)
          var0@_FillValue = -9999.
          var1@_FillValue = -9999.
          var0 = var0 / 100.0 
          var1 = var1 / 100.0
          var0@units = "hm" ;hectometer
          var1@units = "hm"
          printMinMax(var0,False)
          printMinMax(var1,False)
          delete([/varx,vtmp/])
        end if
        if (region(k).eq."land") then
         ;mask out land (this is redundant for data that is already masked)
          d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
          basemap = d->LSMASK
          lsm  = landsea_mask(basemap,var1&lat,var1&lon)
          var0 = mask(var0,lsm.lt.1,False)
          var1 = mask(var1,lsm.lt.1,False)
          delete([/d,lsm,basemap/])
        end if 
        rad = 4.0 * atan(1.0) / 180.0
        wgt = tofloat(cos(f1->lat * rad))
        vrmse(k,j,i) = wgt_arearmse(tofloat(var0), tofloat(var1), wgt,1.0, 0)
        delete([/var0,var1/])
      end do 
    end do 

    outdir = diro
    system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
    fout1 = outdir+run_name(i)+".rmse.seas.clim."+scenario+".nc"
    system("rm -rvf " + fout1)
    setfileoption("nc","Format","netcdf4")

    z1 = addfile(fout1,"cw")
    z1->varnam = varnam
    z1->expnam = run_name(i)
    z1->ensnum = ensnum(i)
    z1->season = seas
    z1->rmse   = vrmse(:,:,i)

  end do

end 

