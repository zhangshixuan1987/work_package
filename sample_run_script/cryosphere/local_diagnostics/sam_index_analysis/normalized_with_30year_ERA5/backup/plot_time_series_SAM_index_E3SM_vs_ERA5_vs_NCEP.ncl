;*****************************************************
; xy_18.ncl
;
; Concepts illustrated:
;   - Filling the area between two curves in an XY plot
;   - Labeling the bottom X axis with years
;   - Drawing a main title on three separate lines
;   - Calculating a weighted average
;   - Changing the size/shape of an XY plot using viewport resources
;   - Manually creating a legend
;   - Overlaying XY plots on each other
;   - Maximizing plots after they've been created
;
;*****************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;*****************************************************
begin

 Datdir   = "./"

 ExpNames  = (/"NCEP-DOE","ERA5","E3SM_NE30", "E3SM_RRM", "E3SM_NE120"/)
 nexps    = dimsizes(ExpNames)
 Labels   = (/"NCEP-DOE","ERA5","E3SM (NE30)", "E3SM (RRM)", "E3SM (NE120)"/)

 ;;;;;only have year 2010 to 2012 data;;;;
 nyear  = 3
 nmonth = 12
 ncol   = nmonth + 1
 ntime = nyear*nmonth

 yyyy   = new(nyear,float)
 mm     = new(nmonth,float)
 time   = new(ntime,float)
 yyyy   = ispan(2010,2012,1)
 mm     = ispan(1,12,1)
 time  = ispan(1,ntime,1)
 dateF = time

 do iy = 0,nyear-1,1
  do im = 0, nmonth-1,1
    n = iy*nmonth + im 
    dateF(n) = yyyy(iy) + (mm(im)-1)*1.0/nmonth  
  end do
 end do 

 plotdata = new ((/nexps,ntime/),float)

 ;*****************************************************
 ; open files and read in monthly data
 ;*****************************************************

 do iexp = 0,nexps-1,1

  File             = "SAM_Norm_AOI_Monthly_2010-2012_"+ExpNames(iexp)+".txt"
  data             = asciiread(File,(/nyear,ncol/),"float")     ; 2 x 70 x 70
  plotdata(iexp,:) = ndtooned(data(:,1:)) 
 
 end do 

;;;;correlation;;;;
r      = escorc_n(plotdata(0,:),plotdata(1:,:),0,1)   ; r is a scalar

print(r)

 t    = r*sqrt((n-2)/(1-r^2))      
 p    = student_t(t, df)
 psig = 0.05                       ; test significance level                     

 print(p)
exit

;***************************************************** 
; create plot
;***************************************************** 
  wks = gsn_open_wks("png","xy")             ; send graphics to PNG file

  res                    = True              ; plot mods desired
  res@gsnDraw            = True             ; don't draw yet
  res@gsnFrame           = True             ; don't advance frame yet

  res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
  res@vpWidthF 	         = 0.7

  res@trYMinF  = -3.0                   ; min value on y-axis
  res@trYMaxF  =  3.0                   ; max value on y-axis

 ;res@trXMinF            = 1890              ; set x-axis minimum

  res@tiYAxisString      = "SAM Index"    ; y-axis label      
  res@tiMainString       = "Variations of SAM Index" ; title
  
; create a reference line and shade values above and below with 
; selected colors. This is shading array dsoid.
  res@gsnYRefLine           = 0.0             ; create a reference line   
;*********************************
; polyline parameters used on both plots
;*********************************
  polyres                  = True
  polyres@gsLineThicknessF = 3.0
;*********************************

  res@xyMonoLineColor    = False             ; want colored lines
  res@xyLineColors       = (/"Red","Blue","Green3","Orange","Black"/) ; colors chosen
  res@xyLineThicknesses	 = (/3.,3.,4.,4.,4.,4./)      ; line thicknesses
  res@xyDashPatterns	 = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid

  plot = gsn_csm_xy (wks,dateF,plotdata,res)        ; create plot
 
end
