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

  data_dir  = "/global/cfs/projectdirs/m3522/cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-AM4/amip"
  exp_name  = (/"NOAA-GFDL_CMIP6"/)
  run_name  = (/"r1i1p1f1"/)
  nexps     = dimsizes(exp_name)

  ystr      = 1979
  yend      = 2014

  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"

  seasons   = (/"MON"/) ;,"ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       =  90.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  vsrList    = (/"psl","ps","ts","tas","tauu","tauv","sfcWind","pr","prc","prsn"/)
  varList    = (/"PSL","PS","TS","TREFHT","TAUX","TAUY","S10","PRECT","PRECC","PRESCT"/)
  vnmList    = (/"sea level pressure","surface pressure", "surface temperature", "surface air temperature", "Zonal surface stree", "meridional surface stress","total precipitation","convective precipitation", "snowfall"/)
  facList    = (/  1.0,   1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 86400,   86400, 86400/)
  untList    = (/"Pa","pa","K","K", "m~S~2 s~S~-2~N~", "m~S~2 s~S~-2~N~", "m s~S~-1~N~", "mm day~S~-1~N~", "mm day~S~-1~N~","mm day~S~-1~N~"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   do j = 0,nvars-1,1 

    vars   = varList(j)
    varunt = untList(j)
    facmod = facList(j)

    cmd  = "ls -d " + data_dir+"/*/Amon/"+vsrList(j)+"/*/*" 
    fdir = systemfunc(cmd)
    nens = dimsizes(fdir)

    iie = 0 
    do ie = 0,nens-1,1
  
     files = systemfunc("ls "+fdir(ie)+"/*.nc")
     print("working on "+ files)
 
     if (dimsizes(files) .gt. 1 ) then 
       f1    = addfiles(files,"r")
       ;ListSetType (f1, "join")
       utc_date = cd_calendar(f1[:]->time, 0)
       year     = tointeger(utc_date(:,0)) 
       indy     = ind(year.ge.ystr.and.year.le.yend)
       l_flag = (.not.all(ismissing(indy))).and.(max(year(indy)).eq.yend).and.(min(year(indy)).eq.ystr)
     else 
       l_flag = False 
       indy = 0.0 
       year = 0.0
       utc_date = 0.0  
     end if 

     if(l_flag) then 

       var0     = f1[:]->$vsrList(j)$(indy,{latS:latN},:)
       time     = f1[:]->time(indy)
       ;printVarSummary(var0)

       do i = 0, nseas-1,1

         print("working on season "+ seasons(i))

         if (seasons(i).eq."MON") then

           vb1   = var0(time|:,lat|:,lon|:)
           vc1   = var0(time|:,lat|:,lon|:)
           vd1   = var0(time|0:11,lat|:,lon|:)

           var  = vb1(:,:,:)
           Clm  = clmMonTLL(var) 
           Anom = calcMonAnomTLL(var, Clm)
           Clm!0         = "month"
           vd1!0         = "month"
           vd1&month     = Clm&month
           vc1(:,:,:) = (/Anom/)
           vd1(:,:,:) = (/Clm/)
           delete([/var,Clm,Anom/])

           printVarSummary(vd1)
           printMinMax(vd1, True)

         else if (seasons(i).eq."ANN") then

           vb1        = month_to_annual(var0(time|:,lat|:,lon|:), 1)
           vb1!0      = "time"
           ttmp       = ispan(0,dimsizes(time)/12-1,1)
           ttmp@units = "years since "+ystr+"-01-01 00:00:00"
           vb1&time   = ttmp

           vd1        = dim_avg_n_Wrap(vb1,0)
           vc1        = vb1
           vc1        = vb1 - conform(vb1,vd1,(/1,2/))
           delete(ttmp)

         else

           vb1        = month_to_season(var0(time|:,lat|:,lon|:), seasons(i))
           vb1!0      = "time"
           ttmp       = ispan(0,dimsizes(time)/12-1,1)
           ttmp@units = "years since "+ystr+"-01-01 00:00:00"
           vb1&time   = ttmp

           vd1  = dim_avg_n_Wrap(vb1,0)
   
           vc1  = vb1
           vc1  = vb1 - conform(vb1,vd1,(/1,2/))
           delete(ttmp)

         end if
         end if

         setfileoption("nc", "Format",  "NetCDF4")
         enstr     = "en"+sprinti("%02d",iie)
         out_file  = exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
         system("rm " + out_file)
         fout = addfile(out_file,"cw")
         vs1 = vars+"_o"
         vs2 = vars+"_a"
         vs3 = vars+"_c"

         fout->$vs1$  = vb1
         fout->$vs2$  = vc1
         fout->$vs3$  = vd1
         delete([/fout,vb1,vc1,vd1/])
       end do

       delete([/var0,time/])
       iie = iie + 1
     end if 
     delete([/files,utc_date,year,indy/])

     end do 

    end do 

   end do 

end
