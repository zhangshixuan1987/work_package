load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;;;;;;;;;;;;;;;;
; main program ;
;;;;;;;;;;;;;;;;

begin

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Subreg   = (/ "North Pacific",  "CONUS",  "Tropics"/)
  Subrgs   = (/ "NP",             "CONUS",  "Trop"/)
  minlat   = (/  10.0,            30.0,     -10.0/)
  maxlat   = (/  65.0,            48.0,      10.0/)
  minlon   = (/ 110.0,            235.0,     60.0/)
  maxlon   = (/ 260.0,            292.0,    270.0/)
  nsubreg  = dimsizes(Subreg)

  datdir   = "../fig_data/"
  ystr     = 2008
  yend     = 2017
  nyr      = yend - ystr + 1

  Varname  = "PRECT"
  Varunit  = "mm/day"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  ngrps    = dimsizes(Groups)

  do irg = 0,nsubreg-1,1

  do igp = 1,ngrps-1,1

   do iy = 0,nyr-1
    year = ystr+iy 
    fno  = datdir+ Varname + "_hov_" + Casdir(0) +"_"+str_sub_str(Subreg(irg)," ","_")+ "_"+year+".nc"
    fnm  = datdir+ Varname + "_hov_" + Casdir(igp) +"_"+str_sub_str(Subreg(irg)," ","_")+ "_"+year+".nc"
    flo  = addfile(fno,"r")
    flm  = addfile(fnm,"r")
    vnb  = Varname
    if(year.eq.ystr) then 
      vobs = new((/nyr,dimsizes(flo->time),dimsizes(flo->lon)/),float)
      vobs!0 = "year"
      vobs!1 = "time"
      vobs!2 = "lon"
      vobs&time = flo->time
      vobs&lon  = flo->lon
      vmod = vobs
    end if 
    vobs(iy,:,:) = tofloat(flo->$vnb$)
    vmod(iy,:,:) = tofloat(flm->$vnb$)
    vobs@_FillValue = flo->$vnb$@_FillValue 
    vmod@_FillValue = flm->$vnb$@_FillValue

   end do 

   ;mask fill values 
   if (.not.all(ismissing(vobs))) then 
     vmod = where (ismissing(vobs), vmod@_FillValue, vmod)
     vobs = where (ismissing(vmod), vobs@_FillValue, vobs)
   end if 

  ;process metrics
   metrics    = (/"Variable", "Model(m)", "Obs(o)", \
                  "Mean(o)", "Mean(m)", "Bias(m-o)", "Stddev", \
                  "RMSE(m-o)", "Stddev",  "PCORR(m-o)", "Stddev"/)

   nmetric    = dimsizes(metrics)
   vmetric    = new(nmetric-2,double) 

  ;calculate time average
   vmobs  = wgt_areaave(vobs, 1.0, 1.0, 0)
   vmmod  = wgt_areaave(vmod, 1.0, 1.0, 0)
   pcmod  = pattern_cor(vmod, vobs,1.0, 1) 
   prmod  = wgt_arearmse(vmod,vobs,1.0, 1.0,0)

   vmetric(0) = dim_avg(vmobs) 
   vmetric(1) = dim_avg(vmmod)
   vmetric(2) = dim_avg(vmmod-vmobs) 
   vmetric(3) = dim_stddev(vmmod-vmobs)
   vmetric(4) = dim_avg(prmod)
   vmetric(5) = dim_stddev(prmod)
   vmetric(6) = dim_avg(pcmod)
   vmetric(7) = dim_stddev(pcmod)

   if (all(ismissing(vmobs)))then 
     vmetric(0) = -999.0
     vmetric(1) = -999.0
     vmetric(4) = -999.0
     vmetric(5) = -999.0
     vmetric(6) = -999.0
     vmetric(7) = -999.0
   end if 

   delete([/vmod,vobs,vmmod,vmobs,pcmod,prmod/])

  ;;;;;;;;;;
  ; Output ;
  ;;;;;;;;;;
    header = ""
    str    = ""
    do i = 0, nmetric-1
     header = header + str_insert(metrics(i)," ", 10)
     if(i.eq.0)then 
       tmp    = Varname +"("+Varunit+")"
       str    = str + str_insert(tmp," ",14)
     else if (i.eq.1)then
       moo    = Casdir(igp)
       str    = str + str_insert(moo," ",12)
     else if (i.eq.2)then
       obs    = Casdir(0)
       str    = str + str_insert(obs," ",12)
     else
       str    = str + sprintf("%8.3f",vmetric(i-3))
     end if 
     end if 
     end if 
    end do 

    if(igp.eq.1) then 
      fname   = "mean_metrics_output_"+Subrgs(irg)+".txt"
      if (fileexists(fname)) then
        system("rm " + fname)
        asciiwrite(fname, header)
      end if 
    end if 
    write_table(fname, "a", [/str/], "%s")
  end do
 end do 

end
