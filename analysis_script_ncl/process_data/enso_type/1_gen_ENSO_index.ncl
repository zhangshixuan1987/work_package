  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Purpose: A tempolate to make a eof plot for AAO patterns 
;          Data is vertically interpolated to pressure levels
;          1st panel: EOFs for the first three modes
;          2nd panel: time series plot
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/sst"
  exp_name  = (/"HadSST","ERA5","v2.LR.amip","v2.LR.historical"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)

  ystr      = 1979
  yend      = 2014
  time_tag  = "1979-2014"
  grid_tag  = "1x1"
  timunt    = "Months since "+ ystr + "-01-01 00:00:00"

  vars      = "SST"
  vout      = "SST"
  vfac      = 1.0
  vunt      = "~S~o~N~C"

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ind_name  = "Nino34"
  region    = (/190.,240.,-5.,5./)
  latS      =  -5.0  ;negative for southern hemisphere
  latN      =   5.0  ;negative for southern hemisphere
  lonW      = 190.0
  lonE      = 240.0
  
  do k  = 0,nexps-1,1

   do i = 0, nseas-1,1

   print("working on season "+ seasons(i))

   fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+seasons(i)+"*"+time_tag+".nc")
   nens  = dimsizes(fl1)
   print(fl1)

   do ie = 0,nens-1,1

    print(fl1(ie))
    f1 = addfile(fl1(ie),"r")

    d2r    = get_d2r("float")
    wgt    = f1->lat
    wgt    = tofloat(cos(f1->lat*d2r))
    time   = cd_convert(f1->time,timunt)
    
   ;printVarSummary(time)
   ;print(cd_calendar(time,-3))
   ;exit

    vcr  = vars + "_a"
    var1 = f1->$vcr$
    var1@_FillValue = -9999.
    delete(var1&time)
    var1&time = time 
   ;printVarSummary(var1)
   ;printMinMax(var1, True)

    vd1 = var1(:,0,0)
    vd1 = wgt_areaave(var1(:,{latS:latN},{lonW:lonE}),wgt({latS:latN}),1.0,0)

    ;Oceanic Niño Index [ONI]. The ONI is defined as the
    ;3-month running means of SST anomalies in the Niño 3.4 region [5N-5S, 120-170W].
    nave = 3
    opt  = 1
    vd1  = runave_n_Wrap (vd1, nave, opt, 0)
    printMinMax(vd1,False)

   ;======================================================================
   ;Save the data for post analysis
   ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    enstr     = "en"+sprinti("%02d",ie)
    out_file  = exp_name(k)+"_"+ind_name+"_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
    fout->ONI = vd1
    delete([/fout,vd1,var1,wgt,time/])

   end do 

   end do 
  delete([/fl1/])
 end do 

end
