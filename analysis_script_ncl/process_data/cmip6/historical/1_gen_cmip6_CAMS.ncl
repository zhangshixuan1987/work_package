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

  data_dir  = "/global/cfs/projectdirs/m3522/cmip6/CMIP6/CMIP/CAMS/CAMS-CSM1-0/amip"
  exp_name  = (/"CAMS_CMIP6"/)
  nexps     = dimsizes(exp_name)

  ystr      = 1979
  yend      = 2014

  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"

  ytrend    = 2001
  nseg      = 3  ; 3 trends: 1979-2001, 1979-2014, 2001-2014
  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       =   0.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  vsrList    = (/"ua","ta","o3"/)
  varList    = (/"U","T","O3"/)
  vnmList    = (/"Zonal wind", "Temperature", "Ozone"/)
  facList    = (/  1.0,   1.0, 1.0/)
  untList    = (/"m s~S~-1~N~","K", "mol mol~S~-1~N~"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   do j = 0,nvars-1,1 

    fil1 = systemfunc("ls "+data_dir+"/*/Amon/"+vsrList(j)+"/*/*/*197901-201412.nc")
    print(fil1)
    nens = dimsizes(fil1)

    vars   = varList(j)
    varunt = untList(j)
    facmod = facList(j)

    do ie = 0,nens-1,1 

    files = (/fil1(ie)/)

    print("working on "+ files)
    f1    = addfiles(files,"r")
    ;ListSetType (f1, "join")

    utc_date = cd_calendar(f1[:]->time, 0)
    year     = tointeger(utc_date(:,0)) 
    indy     = ind(year.ge.ystr.and.year.le.yend)
    var0     = f1[:]->$vsrList(j)$(indy,:,{latS:latN},:)
    var0!1   = "lev"
    time     = f1[:]->time(indy)
   ;printVarSummary(var0)

    delete([/utc_date,year,indy/])

    do i = 0, nseas-1,1

      print("working on season "+ seasons(i))

      if (seasons(i).eq."MON") then

        vb1   = var0(time|:,lev|:,lat|:,lon|:)
        vc1   = var0(time|:,lev|:,lat|:,lon|:)
        vd1   = var0(time|0:11,lev|:,lat|:,lon|:)

        nlev  = dimsizes(var0(0,:,0,0))
        do il = 0,nlev-1,1
          var  = vb1(:,il,:,:)
          Clm  = clmMonTLL(var) 
          Anom = calcMonAnomTLL(var, Clm)
          Clm!0         = "month"
          vd1!0         = "month"
          vd1&month     = Clm&month
          vc1(:,il,:,:) = (/Anom/)
          vd1(:,il,:,:) = (/Clm/)
          delete([/var,Clm,Anom/])
        end do

        printVarSummary(vd1)
        printMinMax(vd1, True)

      else if (seasons(i).eq."ANN") then

        vb1        = month_to_annual(var0(time|:,lev|:,lat|:,lon|:), 1)
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(time)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp

        vd1        = dim_avg_n_Wrap(vb1,0)
        vc1        = vb1
        vc1        = vb1 - conform(vb1,vd1,(/1,2,3/))
        delete(ttmp)

      else

        vb1        = month_to_season(var0(time|:,lev|:,lat|:,lon|:), seasons(i))
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(time)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp

        vd1  = dim_avg_n_Wrap(vb1,0)

        vc1  = vb1
        vc1  = vb1 - conform(vb1,vd1,(/1,2,3/))
        delete(ttmp)

      end if
      end if

      setfileoption("nc", "Format",  "NetCDF4")
      enstr     = "en"+sprinti("%02d",ie)
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
    delete(var0)

  end do 

 end do 

end do 

end
