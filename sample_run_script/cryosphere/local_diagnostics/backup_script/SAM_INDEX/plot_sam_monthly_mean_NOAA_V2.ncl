;==============================================================
;  Create a SAM index using Gong and XXWang, 1999
;  originally written by PSD Dec 2015.
;==============================================================
; run using ncl -Qn < sam.20crv2c.mon.ncl >outputfile
;==============================================================
; code assumes 12 months for all years for file. could be changed to handle 
; data that doesn't end in December"
;==============================================================

; ==============================================================
; User defined parameters that specify region of globe and data range 
; and climo time period
; opendap=1 uses PSD's opendap URL instead of a local file
; ==============================================================

  opendap=2
  
  lat40   = -40.
  lat65   = -65.

  yrStrt = 1871
  yrLast = 2012
  
  cl1=1981
  cl2=2010

; ==============================================================
; Open the file
; ==============================================================

  var    = "prmsl"
  title  = str_upper(var)+" "
  if(opendap.eq.1)then
    f      = addfile ("http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/20thC_ReanV2c/Monthlies/monolevel/prmsl.mon.mean.nc", "r")
  else
    f      = addfile ("/compyfs/zhan391/E3SM_PROJECT/NOAA-CIRES-DOE_20th_Century/V2/prmsl.mon.mean.nc", "r")
  end if

; ==============================================================
; Create a zonal mean for each latitude
; ==============================================================

  x40    = ( f->$var$(:,{lat40:lat40},:) )
  zAvg40 = dim_avg_Wrap( x40 ) 
  x65   = ( f->$var$(:,{lat65:lat65},:) )  
  zAvg65 = dim_avg_Wrap( x65 )

  printVarSummary(x40)

; ==============================================================
; set arrays for computation
; ==============================================================

 nyears= yrLast-(yrStrt-1)
 xm40 = new(nyears,float)
 xm65 = new(nyears,float)
 sam  = new(nyears*12,float)

; ==============================================================
; loop through 12 months; standardize each month separately
; use cl1 and cl2 for years for standardization
; ==============================================================

  i = 0
  
  do while(i.le.11)
     k=i+(cl1-yrStrt)*12
     iyear=cl1
     do while (iyear.le.cl2)
       xm40(iyear-yrStrt)=zAvg40(k,0)
       xm65(iyear-yrStrt)=zAvg65(k,0)       
       k=k+12
       iyear=iyear+1
     end do
;    either standardize time series OR get std and mean and remove
     std40 = stddev(xm40)
     avg40 = avg(xm40)

     std65 = stddev(xm65)
     avg65 = avg(xm65)
     iyear=yrStrt
     k=i
     do while (iyear.le. yrLast)
       zAvg40(k,0)=(zAvg40(k,0)-avg40)/std40
       zAvg65(k,0)=(zAvg65(k,0)-avg65)/std65
       sam(k)=zAvg40(k,0)-zAvg65(k,0)
       k=k+12
       iyear=iyear+1       
     end do
     i=i+1
   end do

; ==============================================================
; write out results
; ==============================================================

  print(" "+yrStrt+" "+  yrLast)
  do n=yrStrt, yrLast
    nn=(n-yrStrt)*12
    print(n+" "+sprintf("%7.3f",sam(nn))+" "+sprintf("%7.3f",sam(nn+1))+" "+sprintf("%7.3f",sam(nn+2))+" "+sprintf("%7.3f",sam(nn+3))+" "+sprintf("%7.3f",sam(nn+4))+" "+sprintf("%7.3f",sam(nn+5))+" "+sprintf("%7.3f",sam(nn+6))+" "+sprintf("%7.3f",sam(nn+7))+" "+sprintf("%7.3f",sam(nn+8))+" "+sprintf("%7.3f",sam(nn+9))+" "+sprintf("%7.3f",sam(nn+10))+" "+sprintf("%7.3f",sam(nn+11)))
  end do
  
; print -999 as missing value for PSD standard time series format

  print (" -999.")
  print ("Southern Annular Mode (SAM) from the 20CRV2c")
  print ("Using the definition for the AOI from paper")
  print ("Daoyi Gong and Shaowu Wang")
  print ("Definition of Antarctic Oscillation index")
  print ("Geophysical  Research Letters, Vol.  26, No.4, 459-462, February  15, 1999") 
  print ("DOI: 10.1029/1999GL900003")
  print ("Zonal mean SLP 40S (standardized by month) minus 65S- (standardized by month)")
  print ("Using "+cl1+"-"+cl2+" for standardization period")


