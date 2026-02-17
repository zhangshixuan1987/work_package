load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  wks1=gsn_open_wks("X11","Plot1")
  gsn_define_colormap(wks1,"BlAqGrYeOrReVi200")   ; choose colormap
  datdir="/compyfs/zhan391/acme_init/ERA5_Reanalysis/era_to_e3sm/era_to_se/method4"
  f1=addfile(datdir+"/ERA5_ne30np4_L72.cam2.i.2009-01-01-00000.nc","r"/)
  lat1d = f1->lat          ; Pull off lat/lon variable and
  lon1d = f1->lon          ; convert both to 1D.
  Temp = f1->T
  PS   = f1->PS

;  Xvar = Temp(0,29,:)
  Xvar = PS(0,:)

  res                     = True         ; plot modifications desired
  res@gsnMaximize         = True         ; Maximize size of plot in frame
  res@gsnSpreadColors     = True         ; Use full colormap, but start
  res@sfXArray            = lon1d        ; Required to tell NCL where to
  res@sfYArray            = lat1d        ; overlay data on globe.
  res@cnFillOn            = True         ; Turn on contour fill
  res@cnFillMode="RasterFill"
;  res@cnFillMode          = "AreaFill"   ; Style of fill. You can also
                                         ; use "CellFill" and "RasterFill"
  res@cnLinesOn           = False        ; Turn off contour lines
  res@cnLineLabelsOn      = False        ; Turn off contour line labels
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.

  map = gsn_csm_contour_map_ce(wks1,Xvar,res)
  map = gsn_csm_contour_map_ce(wks1,Xvar,res)
end

