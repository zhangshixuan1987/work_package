;==============================================================
; Create a "Southern Annular Mode (SAM)" index using Gong and Wang, 1999
; Definition of Antarctic Oscillation index, Geophysical  Research Letters, 
; Vol.  26, No.4, 459-462, February  15, 1999, DOI: 10.1029/1999GL900003
; Originally written by PSD Dec 2015. Revised by Shixuan Zhang @ PNNL  
; Zonal mean SLP 40S (standardized by month) minus 65S- (standardized by month)
; Using 1981-2010 for standardization period


; ==============================================================
; User defined parameters that specify region of globe and data range 
; and climo time period
; opendap=1 uses PSD's opendap URL instead of a local file
; ==============================================================

  lat40   = -40.
  lat65   = -65.
  
  ;;year range for the SAM index;;;;;;
  yrStrt  = 2010
  yrLast  = 2012
  nyear  = yrLast-(yrStrt-1)
  nmonth  = 12

  ;;;year range for climatology;;;;
  cl1=2010
  cl2=2012
  ncl= cl2 - (cl1 -1)

  out_file_name = "SAM_Norm_AOI_Monthly_"+yrStrt+"-"+yrLast+"_E3SM_NE120.txt"

; ==============================================================
; Open the file
; ==============================================================

  climvar = "PSL"
  var     = "PSL" 
  title  = str_upper(var)+" "
  title  = str_upper(var)+" "
  fclim  = addfile ("/global/cscratch1/sd/zhan391/run_e3sm_crysphere/E3SMV1_NE120_to_180x360.cam.h0.2010-2012.nc","r")
  fcurr  = addfile ("/global/cscratch1/sd/zhan391/run_e3sm_crysphere/E3SMV1_NE120_to_180x360.cam.h0.2010-2012.nc","r")

 ;;;construct the climatology(1981-2010) for the normalization;;;;;;;
  climtime = fclim->time
  climdate = cd_calendar(climtime, -1)     ; YYYYMMDD [NCL]
  climsym  = toint(cl1+"01")
  climeym  = toint(cl2+"12")
  climind  = ind(.not.ismissing(climdate)) ;ind((climdate.ge.climsym).and.(climdate.le.climeym))

  climvar40  = fclim->$climvar$(climind,{lat40:lat40},:)
  climvar65  = fclim->$climvar$(climind,{lat65:lat65},:)
  
  climzAvg40 = dim_avg_n_Wrap( dim_avg_n_Wrap( climvar40, 2), 1)
  climzAvg65 = dim_avg_n_Wrap( dim_avg_n_Wrap( climvar65, 2), 1)

  ;;reshap to 30 year x 12 month x 1 latitude 
  rzAvg40 = reshape(climzAvg40,(/ncl,nmonth/))
  rzAvg65 = reshape(climzAvg65,(/ncl,nmonth/))
  
  std40   = dim_stddev_n_Wrap(rzAvg40,0)
  avg40   = dim_avg_n_Wrap(rzAvg40,0)

  std65   = dim_stddev_n_Wrap(rzAvg65,0)
  avg65   = dim_avg_n_Wrap(rzAvg65,0)

  delete([/rzAvg40,rzAvg65/])

 ;;;;;;;select time period of interested period;;;
  currtime = fcurr->time
  currdate = cd_calendar(currtime, -1)     ; YYYYMMDD [NCL]
  currsym  = toint(yrStrt+"01")
  curreym  = toint(yrLast+"12")
  currind  = ind(.not.ismissing(currdate)) ;(currdate.ge.currsym).and.(currdate.le.curreym))

  var40    = fcurr->$var$(currind,{lat40:lat40},:)
  var65    = fcurr->$var$(currind,{lat65:lat65},:)

  Avg40    = dim_avg_n_Wrap( dim_avg_n_Wrap( var40, 2), 1)
  Avg65    = dim_avg_n_Wrap( dim_avg_n_Wrap( var65, 2), 1)

  rzAvg40 = reshape(Avg40,(/nyear,nmonth/))
  rzAvg65 = reshape(Avg65,(/nyear,nmonth/))

  rzAvg40!0     = "year"
  rzAvg40!1     = "month"
  rzAvg40&year  = ispan(yrStrt,yrLast,1)
  rzAvg40&month = ispan(1,12,1)

  rzAvg65!0     = "year"
  rzAvg65!1     = "month"
  rzAvg65&year  = ispan(yrStrt,yrLast,1)
  rzAvg65&month = ispan(1,12,1)

  rzAvg40 = (rzAvg40  - conform(rzAvg40,avg40,(/1/)) ) /conform(rzAvg40,std40,(/1/))
  rzAvg65 = (rzAvg65  - conform(rzAvg65,avg65,(/1/)) ) /conform(rzAvg65,std65,(/1/))

  sam_ind = rzAvg40 - rzAvg65
  

  ;;;;construct the output array and output the data to a file;;;;
  system("/bin/rm -f "+out_file_name)   ; remove any pre-existing file
  system("/bin/touch "+out_file_name)   ; remove any pre-existing file

  do iy = 0,nyear-1

   outline  = sprinti("%04d",toint(rzAvg40&year(iy)))
   do im  = 0,nmonth-1
     outline = outline + sprintf("%8.3f",sam_ind(iy,im))
   end do 
   write_table(out_file_name, "a", [/outline/], "%s")

 end do 

  
  sam     = ndtooned(sam_ind)
; ==============================================================
; write out results
; ==============================================================

  print(" "+yrStrt+" "+  yrLast)
  do n=yrStrt, yrLast
    nn=(n-yrStrt)*12
    print(n+" "+sprintf("%7.3f",sam(nn))+" "+sprintf("%7.3f",sam(nn+1))+" "+sprintf("%7.3f",sam(nn+2))+" "+sprintf("%7.3f",sam(nn+3))+" "+sprintf("%7.3f",sam(nn+4))+" "+sprintf("%7.3f",sam(nn+5))+" "+sprintf("%7.3f",sam(nn+6))+" "+sprintf("%7.3f",sam(nn+7))+" "+sprintf("%7.3f",sam(nn+8))+" "+sprintf("%7.3f",sam(nn+9))+" "+sprintf("%7.3f",sam(nn+10))+" "+sprintf("%7.3f",sam(nn+11)))
  end do
