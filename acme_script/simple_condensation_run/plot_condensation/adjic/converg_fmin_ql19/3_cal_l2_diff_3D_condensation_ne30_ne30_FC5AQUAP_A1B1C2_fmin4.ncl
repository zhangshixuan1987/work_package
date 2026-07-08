load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

  ldatalocal = True
  ldatalocal = False

  if (ldatalocal) then
     ModelOutRoot = "./data/"
  else
     ModelOutRoot = "/pic/projects/uq_climate/zhan391/cnvg_condens_master_adqldens/run/"
  end if

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Global","Land","Ocean"/)
  nregs        = dimsizes(Regions)
  hash         = "E3SM_CONVG"

  Group        = "RKZ_A1_B1_C2_ql19_lmt4_fmin4_ne30_ne30_FC5AQUAP_intel_constance"+"_"
  config       = "_dycore_mac_nodens"
  dtimeList    = (/1,2,8,30,120,450,1800/)

 ;Group        = "new_init_ne30_ne30_FC5AQUAP_intel_constance"+"_"
 ;config       = "_dycore_st"
 ;config       = "_dycore_only"
 ;config       = "_dycore_st_no_fixer"
 ;config       = "_dycore_only_no_fixer"
 ;dtimeList    = (/1,2,6,8,30,120,450/)

 ;Group        = "FC5init_ne30_ne30_FC5AQUAP_intel_constance"+"_"
 ;config       = "_dycore_only"
 ;config       = "_dycore_st"
 ;dtimeList    = (/1,2,6,8,30,120,450,1800/)

  ;------

  ncase        = dimsizes(dtimeList)
  h0Suf        = ".cam.h0.0001-01-01-03600.nc"
 ;h0Suf        = ".cam.h0.0001-01-01-01800.nc"
  itime        = 0   ;h0 file has only 1 time step in it

;  varNameList  = (/"T","CLDLIQ","CLDICE","CLOUD"/)
  varNameList  = (/"T"/)
  nvar         = dimsizes(varNameList)

  irs = 1
  ire = 1
  nr = ire - irs + 1
  diffnorm0 = new( (/3,ncase/), "double" )
  diffnormtable = new( (/nr,ncase-1/),"double" )
  slopes = new( nr,"double" )

;##########
do irg = 0,nregs-1

 do iv = 0,nvar-1
;do iv = 2,2
 print(" variable = "+varNameList(iv))

  do ir = 0,nr-1

  rid = "_"+sprinti("%02d",ir+irs)

  do ic = 0,ncase-1

     case     = Group + "DT"+sprinti("%04d",dtimeList(ic))+rid+config

     if (ldatalocal) then
        FileName = ModelOutRoot+case+h0Suf
     else
        FileName = ModelOutRoot+case+"/"+case+h0Suf
     end if

     File     = addfile(FileName,"r")
     area     = File->area
     lat      = File->lat
     lon      = File->lon

     ;----------------------------------
     ; compute pressure layer thickness
     ;----------------------------------
     p0       = File->P0
     ps       = File->PS(itime,:)

     hyai     = File->hyai
     hybi     = File->hybi
     nlev     = dimsizes(hyai) -1
     daAll    = hyai(1:nlev) - hyai(0:nlev-1)
     dbAll    = hybi(1:nlev) - hybi(0:nlev-1)

     da = daAll( : )
     db = dbAll( : )
    
     varName  = varNameList(iv)
     var3d    = File->$varName$(itime,:,:)
     landfrac = File->LANDFRAC(itime,:)

     if((min(landfrac).eq.0).and.(max(landfrac).eq.0));then
     delete([/landfrac/])
     ;;;use the full model run to get a land fraction;;;   
     maskfile = "/pic/projects/uq_climate/zhan391/cnvg_v1/run/new_init_ne30_ne30_E3SM_CONVG_FC5AV1C-04P2_intel_constance_DT0001_01_fullmodel_rad/new_init_ne30_ne30_E3SM_CONVG_FC5AV1C-04P2_intel_constance_DT0001_01_fullmodel_rad.cam.h0.0001-01-01-03600.nc"
     Maskf    = addfile(maskfile,"r")
     landfrac = Maskf->LANDFRAC(itime,:)
     delete([/maskfile,Maskf/])
     end if

    ; if(varName.eq."Q")then
    ; dat3d    = dat3d*1000.0 ;convert to g/kg
    ; end if

    ; if(Regions(irg).eq."Global") then
    ; var3d    = dat3d
    ; end if 

    ; if (Regions(irg).eq."Land") then
    ; var3d    = dat3d
    ; var3d    = mask(dat3d,conform(dat3d,landfrac,1),1)
    ; end if 
     
    ; if (Regions(irg).eq."Ocean") then
    ; var3d    = dat3d
    ; var3d    = mask(dat3d,conform(dat3d,landfrac,1),0)
    ; end if

    ; print(num(.not.ismissing(dat3d)))
    ; print(num(.not.ismissing(var3d)))

    ; printVarSummary(dat3d)
    ; printVarSummary(var3d)
    ; delete([/dat3d/])

     dp       = conform( var3d, da, 0 )* p0 \
              + conform( var3d, db, 0 )* conform( var3d, ps, 1 )

    ;;----------------------------------
    ;; global mean
    ;;----------------------------------

    ;vint = dim_sum_n( var3d*dp, 0 )  ; vertical integral with mass weighting
    ;glbavg(ic) =  dim_avg_wgt( vint, area, 0 ) \
    ;             /dim_avg_wgt( ps,   area, 0 ) 
   
     ;----------------------------------
     ; RMS difference and L2 diff norm
     ;----------------------------------
     if (ic.eq.0) then

        var3d_ref = var3d
           dp_ref = dp
           ps_ref = ps

     else

        dpw = 0.5*( dp_ref + dp )
        psw = 0.5*( ps_ref + ps )

        ; Linf, L1 and L2 norms


      ; diffnorm0(0,ic) = max(abs(var3d-var3d_ref))

      ; vint = dim_sum_n( abs(var3d-var3d_ref) *dpw, 1 )
      ; diffnorm0(1,ic) = dim_avg_wgt( vint, area, 1 ) \
      ;                  /dim_avg_wgt( psw,  area, 1 ) 

        diffvar3d=var3d
        diffvar3d=(var3d-var3d_ref)^2*dpw
        vint = dim_sum_n(diffvar3d, 0 )

     if (Regions(irg).eq."Global") then
        diffnorm0(2,ic) = sqrt(dim_avg_wgt( vint, area, 0 ) \
                               /dim_avg_wgt( psw,  area, 0 ))
     end if

     if (Regions(irg).eq."Land") then
        diffnorm0(2,ic) = sqrt(dim_avg_wgt( vint, area*landfrac, 0 ) \
                               /dim_avg_wgt( psw,  area*landfrac, 0 ))
     end if

     if (Regions(irg).eq."Ocean") then
        diffnorm0(2,ic) = sqrt(dim_avg_wgt( vint, area*(1.0-landfrac), 0 ) \
                               /dim_avg_wgt( psw,  area*(1.0-landfrac), 0 ))
     end if 

       
     end if
  end do ; case loop (different time steps)


  diffnormtable(ir,:) = diffnorm0(2,1:)

  x  = log10(dtimeList*1.)
  y  = log10(diffnorm0(2,:))
  rc = regline(x,y)
  slopes(ir) = (/rc/) 

end do ; ir

  outpdir = "./"
  filname = outpdir+"rmse_reft_1s_"+Regions(irg)+"_"+Group+varNameList(iv)+config+".txt"
  system("rm -f " + filname)
  fmt = ncase+"g13.5"
  opt        = True
; opt@title  = ""
  opt@tspace = 40
  opt@fout   = filname

  write_matrix(diffnormtable,fmt,opt)

  if (nr.eq.1) then
  print("slope = "+slopes(0))
  else
  print("")
  print("slopes")
  print("")

  fmt = nr/2+"g13.5"
  write_matrix(reshape(slopes,(/2,nr/2/)),fmt,False)
  print("")
  print("slopes: avg = "+avg(slopes)+", stddev = "+stddev(slopes))
  end if

 end do ;; end of variable loop

end do ; region loop (different Regions)

end
