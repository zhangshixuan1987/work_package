;*************************************************
; NCL tutorial script: vert_1.ncl
;************************************************
begin

data_path = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/data"

expnam    = (/"CLIM"/)
nexps     = dimsizes(expnam)

varnam    = (/"U","V"/)
nvars     = dimsizes(varnam)

plev      = (/200, 850/)
nlev      = dimsizes(plev)

syear     = 2007
eyear     = 2017

do i = 0,nexps-1,1

 do j = 0,nvars-1,1

  do k = 0,nlev-1,1
   
   dir = data_path+"/" + expnam(i)+"/"
   var = varnam(j)+plev(k)
    
   fli = dir + var+"_"+syear + "-" + eyear + ".nc"
   fin = addfile(fli,"r") 

   vd    = fin->$var$ 
   vclim = calculate_monthly_values(vd, "avg", 0, False) 
 
   printVarSummary(vclim)  
   printMinMax(vd, True)

  end do 
 end do 
end do 

   outfile  = data_path + "/ERA5_PLEV_3hourly_"+year+"01010000-"+year+"12312100.nc"
  ;setfileoption("nc", "FileStructure", "Advanced")
  ;setfileoption("nc", "Format",  "NetCDF4")

  system(" rm -rvf " + outfile )
  fout = addfile(outfile,"c")

  do i = 0,dimsizes(pnew)-1
    vout = "T"+toint(pnew(i))
    fout->$vout$ = Tp(:,i,:)
    fout->$vout$@long_name = "Temperature at " + toint(pnew(i))+" hPa"
   ;printVarSummary(Tp)
   ;printMinMax(Tp, True)
  end do

;************************************************
; plot parameters that do not change
;************************************************
   wks  = gsn_open_wks("png","vert")              ; send graphics to PNG file
   plot = new(4,graphic)                          ; Define plot array 
   r    = True                                    ; plot mods desired
   r@gsnFrame = False                             ; required for panel plots
   r@gsnDraw  = False                             ; required for panel plots

;---Resources for plotting original data
   r@gsnAddCyclic = False
   r@sfXArray     = fin->lon
   r@sfYArray     = fin->lat

;************************************************
; create individual contour plots
;************************************************
   r@tiMainString= "500 mb"                      ; plot title
   plot(0)=gsn_csm_contour_map(wks,dim_avg_n(Tp(:,0,:),0),r)  ; plot temp contours 

   r@tiMainString= "400 mb"                      ; plot title
   plot(1)=gsn_csm_contour_map(wks,dim_avg_n(Tp(:,1,:),0),r)  ; plot temp contours

   r@tiMainString= "300 mb"                      ; plot title
   plot(2)=gsn_csm_contour_map(wks,dim_avg_n(Tp(:,2,:),0),r)  ; plot temp contours

   r@tiMainString= "200 mb"                      ; plot title
   plot(3)=gsn_csm_contour_map(wks,dim_avg_n(Tp(:,3,:),0),r)  ; plot temp contours
;***********************************************
   gsn_panel(wks,plot,(/2,2/),False)              ; create panel plots

end


