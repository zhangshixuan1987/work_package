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

  data_dir  = "/global/cfs/projectdirs/m3522/cmip6/CMIP6/CMIP/NCAR/CESM2/historical"
  grid_info = "/global/cfs/cdirs/e3sm/zhan391/data/regrid_maps/grid_lat_lon_180x160.nc"
  exp_name  = (/"CESM2_CMIP6"/)
  nexps     = dimsizes(exp_name)

  ystr      = 1979
  yend      = 2014

  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"

  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       =  90.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  vsrList    = (/"psl","pr","prsn","ts","tas","prw", \
                 "tauu","tauv","sfcWind"/)
  varList    = (/"PSL","PRECT","PRECST","TS","TREFHT","TMQ", \
                 "TAUX","TAUY"/)
  vnmList    = (/"Zonal wind", "Temperature", "Ozone"/)
  facList    = (/  1.0, 8.64e7, 8.64e7, 1.0,1.0,1e3,-1,-1/)
  untList    = (/"Pa", "mm day~S~-1", "mm day~S~-1", "K", "K", "g m~S~-2~N~", \
                 "N m~S~-2~N~", "N m~S~-2~N~"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   do j = 1,nvars-1,1 

    fil1 = systemfunc("ls "+data_dir+"/") 
    nens = dimsizes(fil1)
   ;print(fil1)

    varout = varList(j)
    varinp = vsrList(j)
    varunt = untList(j)
    facmod = facList(j)

    do ie = 0,nens-1,1 

     if(varinp.eq."pr".and.fil1(ie).eq."r9i1p1f1") then 
       print(varinp + " " + fil1(ie))
       files = systemfunc("ls -1 " + data_dir+"/"+fil1(ie)+"/Amon/"+vsrList(j)+"/gn/v20190401/*.nc")
     else 
       files = systemfunc("ls -1 " + data_dir+"/"+fil1(ie)+"/Amon/"+vsrList(j)+"/*/*/*.nc")
     end if 

     print("working on "+ files)

     f1       = addfiles(files,"r")
     utc_date = cd_calendar(f1[:]->time, 0)
     year     = tointeger(utc_date(:,0)) 
     indy     = ind(year.ge.ystr.and.year.le.yend)

     ;regrid the data;;;
     fg   = addfile(grid_info,"r")
     glat = fg->lat 
     glon = fg->lon 
     tmp  = f1[:]->$varinp$(indy,{latS:latN},:) 
     var0 = area_conserve_remap_Wrap(f1[0]->lon, f1[0]->lat, tmp, glon,glat,False) 
     time     = f1[:]->time(indy)
    ;printVarSummary(var0)
     delete([/utc_date,year,indy,files,f1,glat,glon,fg,tmp/])

     do i = 0, nseas-1,1

      print("working on season "+ seasons(i))

      if (seasons(i).eq."MON") then

        vb1  = var0(time|:,lat|:,lon|:)
        vc1  = var0(time|:,lat|:,lon|:)
        vd1  = var0(time|0:11,lat|:,lon|:)

        Clm  = clmMonTLL(vb1) 
        Anom = calcMonAnomTLL(vb1, Clm)
        Clm!0         = "month"
        vd1!0         = "month"
        vd1&month     = Clm&month
        vc1(:,:,:) = (/Anom/)
        vd1(:,:,:) = (/Clm/)
        delete([/Clm,Anom/])

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
      enstr     = "en"+sprinti("%02d",ie)
      out_file  = exp_name(k)+"_anom_"+varout+"_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
      system("rm " + out_file)
      fout = addfile(out_file,"cw")
      vs1 = varout+"_o"
      vs2 = varout+"_a"
      vs3 = varout+"_c"

      fout->$vs1$  = vb1
      fout->$vs2$  = vc1
      fout->$vs3$  = vd1
      delete([/fout,vb1,vc1,vd1/])

    end do
    delete(var0)

  end do 

 end do 

end do 

end
