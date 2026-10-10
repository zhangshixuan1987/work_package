;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data"
  exp_name  = (/"ERA5"/)
  nens      = 11 
  nexps     = dimsizes(exp_name)
  time_tag  = "1979-2014"
  grid_tag  = "1x1"

  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       =   0.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  eofvar     = ("Z500")
  varList    = (/"U","T"/)
  vnmList    = (/"Zonal wind", "Temperature"/)
  facList    = (/  1.0,   1.0/)
  untList    = (/"m s~S~-1~N~","K"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;


  do i = 0, nseas-1,1

  do j = 0,nvars-1,1

   vars   = varList(j)
   varunt = untList(j)
   facmod = facList(j)
   print("working on variable "+ vars)
   print("working on season "+ seasons(i))

   do k  = 0,nexps-1,1

    do ie = 0,nens-1,1

    enstr = "en"+sprinti("%02d",ie)
    fl1   = data_dir+"/state_var/"+exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_"+enstr+"_"+time_tag(k)+".nc"
    fl0   = data_dir+"/sam_eof/no_mbias/"+exp_name(k)+"_EOF_PCS_"+eofvar+"_"+seasons(i)+"_"+enstr+"_"+time_tag(k)+".nc"

    print(fl1)
    f1 = addfile(fl1,"r")
    f0 = addfile(fl0,"r")

    lat  = f1->lat
    lon  = f1->lon
    lev  = f1->lev

    vcr  = vars + "_a"
    var0 = f1->$vcr$(:,:,{latS:latN},{lonW:lonE})
    var0@_FillValue = -9999.
    var0       = var0*facmod
    var0@units = varunt
    printVarSummary(var0)
    printMinMax(var0, True)

    delete([/dattime,cd_date,year,indy/])

    var1   = dim_avg_n_Wrap(var0,3) 
    samind = f->eof_ts(ie,0,:)
    ind1   = ind(samind.ge.0.5) 
    ind2   = ind(samind.le.-0.5)
    ny1    = dimsizes(ind1)
    ny2    = dimsizes(ind2)
    ny     = dimsizes(samind)
   ;print(cd_calendar(time(ind1),-1))
   ;print(cd_calendar(time(ind2),-1))

    vmn1   = dim_avg_n_Wrap(var1(ind1,:,:),0)
    vmn2   = dim_avg_n_Wrap(var1(ind2,:,:),0)
    vmna   = dim_avg_n_Wrap(var1(:,:,:),0) 

    vsd1   = dim_avg_n_Wrap(var1(ind1,:,:),0)
    vsd2   = dim_avg_n_Wrap(var1(ind2,:,:),0)
    vsda   = dim_avg_n_Wrap(var1(:,:,:),0)   

    printVarSummary(vmn1)
    ;;t-test;;
    iflag    = False  ; population variance similar
    tval_opt = False  ; p-value only
    prob1    = ttest(vmn1,vsd1,ny1,vmna,vsda,ny,iflag,tval_opt)
    prob2    = ttest(vmn2,vsd2,ny2,vmna,vsda,ny,iflag,tval_opt)
    printMinMax(prob1,False)
    printMinMax(prob2,False)
 
    vdf1   = vmn1
    vdf2   = vmn2
    vdf1   = vmn1 - vmna
    vdf2   = vmn2 - vmna
    vpb1   = vmn1
    vpb2   = vmn2
    vpb1   = (/prob1/)
    vpb2   = (/prob2/)
    delete([/fl,f,samind,ind1,ind2,ny1,ny2,vmn1,vmn2,vsd1,vsd2,vsda,var1,prob1,prob2/])       

    ;rad     = 4.*atan(1.)/180.    ;
    ;cost    = cos(lat*rad)       ; cosine weights
    ;mean1x  = wgt_areaave(pncre0300,cost,1.0,0)
    ;mean2x  = wgt_areaave(wncredf1,cost,1.0,0)
    ;mean3x  = wgt_areaave(wncredf2,cost,1.0,0)
    ;mean1   = sprintf("%5.2f", mean1x)
    ;mean2   = sprintf("%5.2f", mean2x)
    ;mean3   = sprintf("%5.2f", mean3x) ;+ "%"

    enstr    = "en"+sprinti("%02d",ie)
    out_file  = exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
    vs1 = vars+"_p"
    vs2 = vars+"_n"
    vs3 = "prob_p"
    vs4 = "prob_n"

    fout->$vars$ = vmna
    fout->$vs1$  = vdf1
    fout->$vs2$  = vdf2
    fout->$vs3$  = vpb1
    fout->$vs4$  = vpb2
    delete([/fout,vmna,vdf1,vdf2,vpb1,vpb2/])

   end do 

  end do 

 end do 

end do 

end
