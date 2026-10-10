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

  flmsk      = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/lnd_msk/landmsk.nc"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim_rad"
  obs_name   = (/"ERA5","CERES_EBAF","OAFlux"/)

  exp_name   = (/"v2.LR.amip", "v2.LR.historical", \
                 "v2.LR.hist-all-xGHG-xaer", "v2.LR.hist-GHG", "v2.LR.hist-aer"/)
  exp_strs   = (/"E3SM AMIP" , "E3SM Coupled", \
                 "E3SM Coupled (natural)", "E3SM Coupled (GHG)", "E3SM Coupled (aerosol)"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "2001-2014"

  ylabels    = (/"Jan","Feb","Mar","Apr","May","Jun",\
                 "Jul","Aug","Sep","Oct","Nov","Dec"/)
  ystring    = "Zonally Averaged Annual Cycle"
  maxlev     = 11

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  10.
  lonW       =   0.
  lonE       = 360.
  sigma      = 5.670374419e-8


  varList    = (/"ALBEDO_SRF","FSDS","CRE_SRF","FSDSC","FLDSC","Qn_SRF","SHFLX","LHFLX"/)
  vouList    = (/"dSAF","dCRE","dnSAF","dFLW","dQn","dHFLX"/)
  varUnit    = (/"~S~o~N~C","~S~o~N~C","~S~o~N~C","~S~o~N~C","~S~o~N~C","~S~o~N~C"/)
  nvar       = dimsizes(varList)

  do k = 0,nexps-1

    fm0  = addfile(flmsk,"r")

    fl0  = systemfunc("ls " + data_dir+"/"+obs_name(0)+"_anom_TS_MON_en*_"+time_tag+".nc")
    f0   = addfiles(fl0,"r")
    ListSetType (f0, "join")
    fl1  = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_TS_MON_en*_"+time_tag+".nc")
    f1   = addfiles(fl1,"r")
    ListSetType (f1, "join")

    fla0 = systemfunc("ls " + data_dir+"/"+obs_name(1)+"_anom_"+varList(0)+"_MON_en*_"+time_tag+".nc")
    fla1 = systemfunc("ls " + data_dir+"/"+obs_name(1)+"_anom_"+varList(1)+"_MON_en*_"+time_tag+".nc")
    fla2 = systemfunc("ls " + data_dir+"/"+obs_name(1)+"_anom_"+varList(2)+"_MON_en*_"+time_tag+".nc")
    fla3 = systemfunc("ls " + data_dir+"/"+obs_name(1)+"_anom_"+varList(3)+"_MON_en*_"+time_tag+".nc")
    fla4 = systemfunc("ls " + data_dir+"/"+obs_name(1)+"_anom_"+varList(4)+"_MON_en*_"+time_tag+".nc")
    fla5 = systemfunc("ls " + data_dir+"/"+obs_name(1)+"_anom_"+varList(5)+"_MON_en*_"+time_tag+".nc")
    fla6 = systemfunc("ls " + data_dir+"/"+obs_name(2)+"_anom_"+varList(6)+"_MON_en*_"+time_tag+".nc")
    fla7 = systemfunc("ls " + data_dir+"/"+obs_name(2)+"_anom_"+varList(7)+"_MON_en*_"+time_tag+".nc")

    flb0 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(0)+"_MON_en*_"+time_tag+".nc")
    flb1 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(1)+"_MON_en*_"+time_tag+".nc")
    flb2 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(2)+"_MON_en*_"+time_tag+".nc")
    flb3 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(3)+"_MON_en*_"+time_tag+".nc")
    flb4 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(4)+"_MON_en*_"+time_tag+".nc")
    flb5 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(5)+"_MON_en*_"+time_tag+".nc")
    flb6 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(6)+"_MON_en*_"+time_tag+".nc")
    flb7 = systemfunc("ls " + data_dir+"/"+exp_name(k)+"_anom_"+varList(7)+"_MON_en*_"+time_tag+".nc")

    fa0  = addfiles(fla0,"r") 
    fa1  = addfiles(fla1,"r")   
    fa2  = addfiles(fla2,"r")
    fa3  = addfiles(fla3,"r")
    fa4  = addfiles(fla4,"r")
    fa5  = addfiles(fla5,"r")
    fa6  = addfiles(fla6,"r")
    fa7  = addfiles(fla7,"r")

    fb0  = addfiles(flb0,"r") 
    fb1  = addfiles(flb1,"r")
    fb2  = addfiles(flb2,"r")
    fb3  = addfiles(flb3,"r")
    fb4  = addfiles(flb4,"r")
    fb5  = addfiles(flb5,"r")
    fb6  = addfiles(flb6,"r")
    fb7  = addfiles(flb7,"r")

    ListSetType (fa0, "join")
    ListSetType (fa1, "join")
    ListSetType (fa2, "join")
    ListSetType (fa3, "join")
    ListSetType (fa4, "join")
    ListSetType (fa5, "join")
    ListSetType (fa6, "join")
    ListSetType (fa7, "join")

    ListSetType (fb0, "join")
    ListSetType (fb1, "join")
    ListSetType (fb2, "join")
    ListSetType (fb3, "join")
    ListSetType (fb4, "join")
    ListSetType (fb5, "join")
    ListSetType (fb6, "join")
    ListSetType (fb7, "join")

    vnm   = "TS_o"
    tsb   = dim_avg_n_Wrap(f0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    tsm   = f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
    tsb   = where(tsb.lt.273.15,tsb+273.15,tsb)
    tsm   = where(tsm.lt.273.15,tsm+273.15,tsm)

    printMinMax(tsb,False)
    printMinMax(tsm,False)

    vnm   = varList(0)+"_o"
    FAL1  = dim_avg_n_Wrap(fa0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    FAL2  = fb0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}) 

    vnm   = varList(1)+"_o"
    FSW1  = dim_avg_n_Wrap(fa1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    FSW2  = fb1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}) 

    vnm   = varList(2)+"_o"
    CRE1  = dim_avg_n_Wrap(fa2[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    CRE2  = fb2[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})

    vnm   = varList(3)+"_o"
    FSWC1 = dim_avg_n_Wrap(fa3[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    FSWC2 = fb3[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})

    vnm   = varList(4)+"_o"
    FLWC1 = dim_avg_n_Wrap(fa4[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    FLWC2 = fb4[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})

    vnm   = varList(5)+"_o"
    FQn1  = dim_avg_n_Wrap(fa5[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    FQn2  = fb5[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})

    vnm   = varList(6)+"_o"
    SH1   = dim_avg_n_Wrap(fa6[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    SH2   = fb6[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})

    vnm   = varList(7)+"_o"
    LH1   = dim_avg_n_Wrap(fa7[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
    LH2   = fb7[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})


    dTs   = tsm 
    dTs   = tsm - conform(tsm,tsb,(/1,2,3/)) 

    FAL1@_missing_value = FAL1@_FillValue
    FAL2@_missing_value = FAL2@_FillValue
    dFAL  = FAL2
    dFAL  = FAL2 - conform(FAL2,FAL1,(/1,2,3/))

    FSW1@_missing_value = FSW1@_FillValue
    FSW2@_missing_value = FSW2@_FillValue
    dFSW  = FSW2 
    dFSW  = FSW2 - conform(FSW2,FSW1,(/1,2,3/))  

    CRE1@_missing_value = CRE1@_FillValue
    CRE2@_missing_value = CRE2@_FillValue
    dCRE  = CRE2 
    dCRE  = CRE2 - conform(CRE2,CRE1,(/1,2,3/)) 

    FSWC1@_missing_value = FSWC1@_FillValue
    FSWC2@_missing_value = FSWC2@_FillValue
    dFSWC = FSWC2 
    dFSWC = FSWC2 - conform(FSWC2,FSWC1,(/1,2,3/))  

    FLWC1@_missing_value = FLWC1@_FillValue
    FLWC2@_missing_value = FLWC2@_FillValue
    dFLWC = FLWC2
    dFLWC = FLWC2 - conform(FLWC2,FLWC1,(/1,2,3/))  
    
    FQn1@_missing_value = FQn1@_FillValue
    FQn2@_missing_value = FQn2@_FillValue
    dFQn  = FQn2
    dFQn  = FQn2 - conform(FQn2,FQn1,(/1,2,3/))  

    LH1@_missing_value = LH1@_FillValue 
    LH2@_missing_value = LH2@_FillValue
    dLH   = LH2 
    dLH   = LH2 - conform(LH2,LH1,(/1,2,3/))
    
    SH1@_missing_value = SH1@_FillValue
    SH2@_missing_value = SH2@_FillValue
    dSH   = SH2 
    dSH   = SH2 - conform(SH2,SH1,(/1,2,3/)) 

    sterm = 4.0 * sigma * tsb * tsb * tsb
    sterm@_missing_value = sterm@_FillValue

    dSAF1 = dFAL
    dSAF1 = -dFAL*conform(FSWC2,FSWC1,(/1,2,3/))/conform(FSWC2,sterm,(/1,2,3/))

    dSAF2 = dFAL
    dSAF2 = -dFAL * dFSW / conform(FSWC2,sterm,(/1,2,3/))
     
    dCRE@_FillValue    = sterm@_FillValue 
    dCRE@_missingvalue = sterm@_FillValue
    
    dCRE  = dCRE/conform(FSWC2,sterm,(/1,2,3/))

    dnSAF = FSWC2
    dnSAF = (1-conform(FAL2,FAL1,(/1,2,3/))) * dFSWC / conform(FSWC2,sterm,(/1,2,3/))

    dFLWC = dFLWC/conform(FSWC2,sterm,(/1,2,3/)) 
   
    dFQn  = -dFQn/conform(FSWC2,sterm,(/1,2,3/))
    
    dHFLX = dLH
    dHFLX = -(dLH+dSH)/conform(FSWC2,sterm,(/1,2,3/))
 
    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file = exp_name(k)+"_energy_MON_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")

    fout->dTs   = dTs
    fout->dFAL  = dFAL
    fout->dFSW  = dFSW
    fout->dCRE  = dCRE
    fout->dFSWC = dFSWC
    fout->dFLWC = dFLWC
    fout->dFQn  = dFQn
    fout->dLH   = dLH
    fout->dSH   = dSH
    fout->dSAF1 = dSAF1
    fout->dSAF2 = dSAF2
    fout->dnSAF = dnSAF
    fout->dHFLX = dHFLX 

    delete([/fl0,fl1,fla0,fla1,fla2,fla3,fla4,fla5,fla6,fla7,flb0,flb1,flb2,flb3,flb4,flb5,flb6,flb7/])
    delete([/fa0,fa1,fa2,fa3,fa4,fa5,fa6,fa7,fb0,fb1,fb2,fb3,fb4,fb5,fb6,fb7/])
    delete([/dTs,dFAL,dFSW,dCRE,dFSWC,dFLWC,dFQn,dLH,dSH,dSAF1,dSAF2,dnSAF,dHFLX/])
    delete([/tsm,tsb,FAL1,FAL2,CRE1,CRE2,FSW1,FSW2,FSWC1,FSWC2,FLWC1,FLWC2,FQn1,FQn2,LH1,LH2,SH1,SH2/])

  end do ; end of exps

end

