;This example regrids data from a global rectilinear grid (360 x 180) to a 1.9x2.0 degree lat/lon grid (144 x 96), using the default "bilinear" interpolation method.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;---Input file

    grid_e3sm = "ne30pg2"
    grid_era5 = "640x1280"
    time_tag1 = "2009-04-01_00:00:00"
    time_tag2 = "2009-04-01-00000"
    DatDir    = "../../test_data"

    nfil     = 2
    files    = new(nfil,string)
    files(0) = DatDir+"/"+"ERA5_analysis_ml_" + grid_era5 + "_" +time_tag1+".nc"
    files(1) = DatDir+"/"+"SZ_"+grid_e3sm+"_to_ERA5_"+grid_era5+"_"+ time_tag2+".nc"

    Labels   = (/"SZ_ADJ2 - ERA5"/)
    
    MainStr  = "Difference before and after interpolation" 

    Varname  = (/"U","V","T","Q"/) ;,"PS"/)
    Varunit  = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~"/) ;,"hPa"/)
    Varfact  = (/1.0, 1.0, 1.0, 1000.0, 1e-2/)
    nvars    = dimsizes(Varname)

    figname  = "ERA5_pl_diff_"+grid_e3sm

    fx       = addfile("./lev.nc","r")
    level    = fx->level 
    nlev     = dimsizes(level)
 
    do j = 0,nvars-1,1

     varNam = Varname(j)
     varUnt = Varunit(j)
     varFac = Varfact(j)

     do k = 0,nlev-1,1
     do i = 1,nfil-1,1

      f0   = addfile(files(0),"r")
      f1   = addfile(files(i),"r")

      lev   = level(k)
      lat   = f0->lat    ;latitude
      lon   = f0->lon    ;longitude
      hyam0 = f0->ohyam
      hybm0 = f0->ohybm
      hyam  = f1->hyam
      hybm  = f1->hybm

      interp = 2 
      P0mb   = 1000.
      extrap = True
      hyam0  = hyam0 / P0mb / 100.0
      temx   = vinth2p(f0->$varNam$,hyam0,hybm0,lev,f0->PS,interp,P0mb,1,extrap)
      tem0   = temx(0,0,:,:)
      temx   = vinth2p(f1->$varNam$,hyam,hybm,lev,f1->PS,interp,P0mb,1,extrap)
      tem1   = temx(0,0,:,:)
      delete(temx)

      temp   = tem0
      temp   = tem1 - tem0
      temp   = temp * varFac

      temp@long_name = "~F33~D~N~~F~"+varNam
      temp@units     = varUnt
      substr         = temp@long_name+" ("+temp@units+", "+sprintf("%.0f",lev)+" hPa"+")"

    ;----------------------------------------------------------------------
    ; Plot the original and regridded data on a panel plot
    ;----------------------------------------------------------------------
    if(i.eq.1) then 
      fignam1 = str_squeeze(figname+"_"+varNam+"_"+sprinti("%04d",toint(lev))+"hPa")
      wks = gsn_open_wks("png",fignam1)             ; ESMF.ps
      plot = new(nfil-1,graphic)
    end if 

    res                     = True              ; Plot mods desired.
    res@gsnDraw             = False             ; We will panel later.
    res@gsnFrame            = False
    res@gsnMaximize         = False              ; Maximize plot
 
    res@mpFillOn                    = False          ; turn off gray continents
    res@mpGeophysicalLineColor      = "Grey30";"white"        ; color of outlines
    res@mpGeophysicalLineThicknessF = 2              ; line thickness
    res@mpProjection                = "Orthographic"   ;-- change projection
    res@mpCenterLonF                =  260              ;-- center at lon=10
    res@mpCenterLatF                =  40              ;-- center at lat=40
    res@mpDataBaseVersion           = "MediumRes"      ;-- map resolution
    res@mpGridAndLimbOn             =  False            ;-- plot grid lines
    res@mpGridLineColor             = "grey30"         ;-- set grid line color
    res@mpPerimOn                   =  False           ;-- don't draw the box around the plot

    res@cnFillOn                 = True              ; color plot desired
    res@cnLinesOn                = False              ; turn off contour lines
    res@cnLineThicknessF         = 0.01
    res@cnLineColor              = "Grey"
    res@cnLineLabelsOn           = False             ; turn off contour lines
    res@cnInfoLabelOn            = False
    res@cnFillMode               = "RasterFill" ;"CellFill"
    res@cnGridBoundFillColor     = "black"
    res@cnRasterMinCellSizeF     = 0.001
    res@cnFillPalette            = "amwg256";"cmocean_curl"
 
    if(varNam.eq."U".or.varNam.eq."V") then 
      xmin = -5. ;-decimalPlaces(max(abs(temp)),2,True)
      xmax =  5. ;decimalPlaces(max(abs(temp)),2,True)
    end if 
    if(varNam.eq."T") then
      xmin = -1. ;-decimalPlaces(max(abs(temp)),2,True)
      xmax =  1. ;decimalPlaces(max(abs(temp)),2,True)
    end if
    if(varNam.eq."Q") then
      xmin = -0.1 ;-decimalPlaces(max(abs(temp)),2,True)
      xmax =  0.1 ;decimalPlaces(max(abs(temp)),2,True)
    end if
    if(varNam.eq."PS") then 
      xmin = -2. ;-decimalPlaces(max(abs(temp)),2,True)
      xmax =  2. ;decimalPlaces(max(abs(temp)),2,True)
    end if
    xint = (xmax - xmin)/10.
    res@cnLevelSelectionMode     = "ManualLevels"  ; manual levels
    res@cnMinLevelValF           = xmin            ; min level
    res@cnMaxLevelValF           = xmax            ; max level
    res@cnLevelSpacingF          = xint            ; interval

    res@lbLabelBarOn             = True              ; Labelbar will be in panel
    res@lbBoxLinesOn             = True
    res@lbLabelStride            =  1               ;-- every other label
    res@lbBoxMinorExtentF        =  0.10            ;-- decrease the height of the labelbar
    res@pmLabelBarOrthogonalPosF = -0.07            ;-- move the labelbar upward

;---Plot data on original grid
   ;res@gsnAddCyclic   = False
   ;res@sfXArray       = lon
   ;res@sfYArray       = lat
    res@tiMainString   = Labels(i-1)
    res@gsnLeftString  = ""
    res@gsnRightString = ""
    res@gsnCenterString= substr 
    plot(i-1) = gsn_csm_contour_map(wks,temp,res)

  end do 

;---Resources for paneling
    pres                          = True
    pres@gsnMaximize              = False
    pres@gsnPanelLabelBar         = False
    pres@lbBoxLinesOn             = True 
    pres@lbLabelStride            = 1               ;-- every other label
    pres@lbBoxMinorExtentF        = 0.15            ;-- decrease the height of the labelbar
    pres@pmLabelBarOrthogonalPosF = 0.03            ;-- move the labelbar upward
    pres@pmLabelBarHeightF        = 0.15
    pres@pmLabelBarWidthF         = 0.95
    gsn_panel(wks,(/plot/),(/1,4/),pres)

   end do
 end do 
end

