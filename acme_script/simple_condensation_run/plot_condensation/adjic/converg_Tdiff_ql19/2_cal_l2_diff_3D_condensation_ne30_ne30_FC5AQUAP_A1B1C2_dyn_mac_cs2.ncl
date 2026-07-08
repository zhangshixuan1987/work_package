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
  Regions      = (/"Global"/) ;,"Land","Ocean"/)
  nregs        = dimsizes(Regions)
  hash         = "E3SM_CONVG"

  Group        = "RKZ_A1_B1_C2_ql19_lmt4_ne30_ne30_FC5AQUAP_intel_constance"+"_"
  config       = "_dycore_mac_nodens"
  dtimeList    = (/1,2,8,30,120,450,1800/)
  ntimes       = dimsizes(dtimeList)

  ncase        = dimsizes(dtimeList)
  h0Suf        = ".cam.h0.0001-01-01-03600.nc"
 ;h0Suf        = ".cam.h0.0001-01-01-01800.nc"
  itime        = 0   ;h0 file has only 1 time step in it

; varNameList  = (/"U","V","T","Q","CLDLIQ","CLDICE","CLOUD"/)
  varNameList  = (/"T"/);"U","V","T","Q"/)
  nvar         = dimsizes(varNameList)

  irs = 1
  ire = 1
  nr = ire - irs + 1
  diffnorm0 = new( (/3,ncase/), "double" )
  diffnormtable = new( (/nr,ncase-1/),"double" )
  slopes = new( nr,"double" )
  rtmsk  = new(ntimes-1,"double") ;;save the ratio for the masked points

;##########
do irg = 0,nregs-1

 do iv = 0,nvar-1
 print(" variable = "+varNameList(iv))
  do it = 1, ntimes-1, 1
  
   do ir = 0,nr-1

   rid = "_"+sprinti("%02d",ir+irs)
  
  ;;;;;;;go through all the files and generate a mask files;;;;;;;;;;;;;;;;;
   if (ldatalocal) then
     FileName1s = Group+"DT"+sprinti("%04d",dtimeList(0))+rid+config+h0Suf
     FileName2s = Group+"DT"+sprinti("%04d",dtimeList(it))+rid+config+h0Suf 
   else
     FileName1s = ModelOutRoot+Group+"DT"+sprinti("%04d",dtimeList(0))+rid+config+"/"+Group+"DT"+sprinti("%04d",dtimeList(0))+rid+config+h0Suf
     FileName2s = ModelOutRoot+Group+"DT"+sprinti("%04d",dtimeList(it))+rid+config+"/"+Group+"DT"+sprinti("%04d",dtimeList(it))+rid+config+h0Suf
   end if
   File1s     = addfile(FileName1s,"r")
   File2s     = addfile(FileName2s,"r")
   ;;;pick the conditional points;;;;;;;;;;;;;;;;;;;;
   T1s    = File1s->T(itime,:,:)
   T2s    = File2s->T(itime,:,:)
   Tdf    = log10(abs(T2s - T1s))
   mask2  = T1s
   mask2  = 0.0
   mask2  = where((Tdf.ge.-4),1,0)
   print(min(mask2)+ " "+max(mask2)) 

   x1D       = ndtooned(mask2)
   indmsk    = ind(x1D.eq.1)
   rtmsk(it-1) = (num(mask2.eq.1)*1.0/num(T1s.gt.0))*100.0
   print("Number of total cells is: " + num(T1s.gt.0)+" "+ dimsizes(x1D))
   print("Number of masked cells is: " + num(mask2.eq.1)+ " "+rtmsk(it-1)*num(T1s.gt.0)/100.0)
  delete([/FileName1s,File1s,FileName2s,File2s,T1s,T2s,x1D,indmsk/])

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
        diffvar3d=var3d
        diffvar3d=(var3d-var3d_ref)^2*dpw
        vint = dim_sum_n(diffvar3d*mask2, 0 )
        pint = dim_sum_n(dpw*mask2, 0)
     if (Regions(irg).eq."Global") then
        diffnorm0(2,ic) = sqrt(dim_avg_wgt( vint, area, 0 ) \
                               /dim_avg_wgt( pint,  area, 0 ))
     end if

     if (Regions(irg).eq."Land") then
        diffnorm0(2,ic) = sqrt(dim_avg_wgt( vint, area*landfrac, 0 ) \
                               /dim_avg_wgt( pint, area*landfrac, 0 ))
     end if

     if (Regions(irg).eq."Ocean") then
        diffnorm0(2,ic) = sqrt(dim_avg_wgt( vint, area*(1.0-landfrac), 0 ) \
                               /dim_avg_wgt( pint,  area*(1.0-landfrac), 0 ))
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
  filname = outpdir+"rmse_reft_1s_"+Regions(irg)+"_"+Group+varNameList(iv)+config+"_cs2_"+sprinti("%04d",dtimeList(it))+".txt"
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

 end do ;; end of it loop

  outpdir = "./"
  filname1 = outpdir+"rtmsk_reft_1s_"+Regions(irg)+"_"+Group+varNameList(iv)+config+"_cs2.txt"
  system("/bin/rm -f "+filname1)   ; remove any pre-existing file
  write_table(filname1, "w", [/rtmsk/], "%13.5f")
  
 end do ;; end of variable loop

end do ; region loop (different Regions)

end
