                                       
! Generated Namelist for maekIC_se.ncl 
!------------------                   
&makeic_nl
 MYTMPDIR="./TMP/TMP_002.01/"
 MYOUTDIR="/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/era_to_se/method4/"
 INPUTDIR="/compyfs/zhan391/acme_init/ERA5_Reanalysis/grb2nc/Jan/"
 INPUTDIR="/compyfs/zhan391/acme_init/ERA5_Reanalysis/grb2nc/Jan/"
                                    
 ESMF_interp="conserve"
 ESMF_pole  ="none"
 ESMF_clean ="False"
 TMP_clean  ="True"
                                    
 REF_DATE              ="20090101"
 CASE                  ="era5_ne30L72"
 DYCORE                ="se"
 PRECISION             ="float"
 VORT_DIV_TO_UV        ="False"
 SST_MASK              ="False"
 ICE_MASK              ="False"
 OUTPUT_PHIS           ="True"
 REGRID_ALL            ="True"
 ADJUST_STATE_FROM_TOPO="True"
 MASS_FIX              ="False"
                                    
 fname_phis_output     ="/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/topo/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc"
 ftype_phis_output     ="SE_TOPOGRAPHY"
 fname_grid_info       ="/compyfs/zhan391/acme_init/csmdata/atm/cam/inic/homme/cami_mam3_0000-10-01_ne30np4_L72_c160127.nc"
                                    
 fields       ="U,V,T,PS,Q"
 source_files ="0,0,0,0,0"
 fname_phis_in=0
                                    
 fname0="era5_an_ml_2009-01-31_23:00:00.nc"
 fname1="none"
 fname2="none"
 fname3="none"
 fname4="none"
 fname5="none"
                                    
 ftype0="ECMWF_ERA5"
 ftype1="none"
 ftype2="none"
 ftype3="none"
 ftype4="none"
 ftype5="none"
                                    
 fdate0="2009013182800"
 fdate1="-1"
 fdate2="-1"
 fdate3="-1"
 fdate4="-1"
 fdate5="-1"
/
