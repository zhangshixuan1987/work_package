;*************************************************
; histo_4.ncl
;
; Concepts illustrated:
;   - Drawing a histogram using exact min/max of data
;   - Generating dummy data using "rand"
;   - Drawing X axis labels at a 45 degree angle
;
;************************************************
;
; This file is loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************

load "./function_read_dart_diag.ncl"

begin

  Timetag = "2011011500-2011012012"
  nday    = 5
  ilev    = 2
  hres    = (/"ne11_oQU240","ne16_oQU240","ne30pg2_EC30to60E2r2"/)
  nres    = dimsizes(hres)
  datadir = "/compyfs/zhan391/diag_dart/"
  outdir  = "/compyfs/zhan391/diag_dart/fig_data/"

  regions = (/"SH","NH","NAmer","Tropics"/)
  regstrs = (/"Southern Hemisphere","Northern Hemisphere","North America","Tropics"/)
  latlons = (/"[0-90S]","[0-90N]","[65~S~o~N~-125~S~o~N~ W, 25~S~o~N~-55~S~o~N~ N]","10S-10N"/)
  nregs   = dimsizes(regions)

  ;variable to be plotted: totalspread,rmse and data rejection rate
  varnam  = (/"RMSE","Bias","Spread","Rrate"/)
  varstr  = (/"RMSE","Bias","Total spread","Rejection rate"/)
  nvars   = dimsizes(varnam)

  ;;;read the variable information;;
  fvar    = "./variable_name.info"
  data    = asciiread(fvar,-1,"string")
  VarDiag = str_squeeze(str_get_field(data, 1," "))
  ndiag   = dimsizes(VarDiag)
  
  do iv = 0,ndiag-1,1
    ObsTypeString = VarDiag(iv)
    do irg  = 0,nregs-1,1
      region = regions(irg)
      regstr = regstrs(irg)
      do i = 0,nres-1,1
        Casename     = "F20TR_"+hres(i)+"_DART80_compy"
        ModelOutRoot = datadir+Casename+"/dart_diagnostics/obs_diag"
        fname        = ModelOutRoot+"/"+Casename+".dart.e.eam_obs_diag_output."+Timetag+".nc"
        time_unit    = "days since 2011-01-15 00:00:00"
        time_string  = "Days since 0000 UTC 15 January 2011"
        print(fname)
        ;exit

        type     = "guess"
        CopyString1 = "totalspread" 
        tspread = read_dart_diag_netcdf(fname, CopyString1, ObsTypeString, regstr, type)

        CopyString2 = "rmse"
        rmse  = read_dart_diag_netcdf(fname, CopyString2, ObsTypeString, regstr, type)

        CopyString3 = "Nposs"
        Nposs = read_dart_diag_netcdf(fname, CopyString3, ObsTypeString, regstr, type)

        CopyString4 = "Nused"
        Nused = read_dart_diag_netcdf(fname, CopyString4, ObsTypeString, regstr, type)

        CopyString2 = "bias"
        bias  = read_dart_diag_netcdf(fname, CopyString2, ObsTypeString, regstr, type)


        plev       = rmse@lev
        plev!0     = "lev"
        plev&lev   = plev 
        plev@units = "hPa"

        time       = rmse@time
        time@units = rmse@time_unit
        time       = cd_convert(time, time_unit)
        indtim     = ind(time.le.nday)
        
        tmpdat            = Nposs
        tmpdat            = where(Nposs.eq.0.0,-999,Nposs)
        tmpdat@_FillValue = -999
        Rrate = rmse 
        Rrate = 100.0 - Nused*100.0 / tmpdat
        
        Rrate!0    = "time"
        Rrate!1    = "lev"
        Rrate&time = time 
        Rrate&lev  = plev 

        copy_VarCoords(Rrate,rmse)
        copy_VarCoords(Rrate,tspread)

        ;;save the plot data;;;;;
        ff  = outdir +ObsTypeString+"_"+region+"_"+hres(i)+".nc"
        system("rm -f " + ff)
        print(" ")
        print("output : " + ff)
        print(" ")

        flo = addfile(ff,"c")
        flo->date          = time(indtim)
        flo->date@units    = time_unit
        flo->plev          = plev
        flo->Bias = rmse(indtim,:)
        flo->RMSE = rmse(indtim,:)
        flo->Spread = tspread(indtim,:)
        flo->Rrate = Rrate(indtim,:)
        delete([/indtim,plev,time,tmpdat,Rrate,tspread,rmse,Nposs,Nused/])

    end do
  end do
end do 

end
