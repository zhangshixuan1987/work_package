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

  ldatalocal    = True
  ldatalocal    = True
  plot_contour  = True   ; if True, plot the contour in the figure

  datadir    = "../fig_data/"

 ;;;;;;;strings for the data directory;;;;;;;;
  Groups     = (/"ERA5"/)
  Expnams    = Groups
  ngrps      = dimsizes(Groups)

 ;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U","T","RELHUM","CLOUD"/)
  VstList  = (/"U","T","RH","CLOUD"/)
  VnmList  = (/"Zonal wind","Temperature", "Relative humidity", "Cloud fraction"/)
  UnitList = (/"m s~S~-1~N~","K","unitless","fraction"/)
  PltNums  = (/"a", "b", "c", "d"/)
  PltNums  = "(" + PltNums+")"
  nvars    = dimsizes(VarList)

  Seasons  = (/"ANN","DJF","JJA"/) 
  nseas    = dimsizes(Seasons) 

  do is = 0,nseas-1,1

  do igp = 0, ngrps-1,1

   do iv = 0,nvars-1, 1

    Varname  = VarList(iv)
    Varunt   = UnitList(iv)
    Varstr   = VstList(iv)
    Vnmstr   = VnmList(iv)
 
    fna   = datadir+ Varname+ "_" + Groups(igp) + ".nc"
    fla   = addfile(fna, "r")
  
    if(Seasons(is).eq."ANN") then 
      xtmp  = month_to_annual(fla->$Varname$,1)
    else
      xtmp  = month_to_season(fla->$Varname$,Seasons(is))
    end if
    xtmp!0 = "time"
    xtmp@_FillValue = -99999
    printMinMax(xtmp,False)

    ;calculate the global/regional mean values
    latrg = fla->lat
    rad   = 4.*atan(1.)/180.    ;
    wgty  = cos(latrg*rad)      ; cosine weights
    vobs  = dim_avg_Wrap (xtmp(lev|:,lat|:,time|:))     

    delete([/fna,xtmp,fla,wgty/])

;----------------------------------------------------------------------
; Plotting section
  if(iv.eq.0) then 
    wtype = "pdf"
    OUTDir = "./"
    outnam = OUTDir + "Fig_pres_lat_raw_"+Groups(igp)+"_"+Seasons(is)
    wks1   = gsn_open_wks(wtype,outnam)
    plot   = new(nvars,graphic)
  end if 

  load "./var_share_colorbar_preshgt.ncl"
  labels = (/Varstr ,"~F33~D~N~~F~"+ Varstr/)
  nlabls = dimsizes(labels)

;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@gsnAddCyclic         = False
  res@vpWidthF             = 0.8
  res@vpHeightF            = 0.6

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"

  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
  res@lbLabelStride           = 1                 ; skip every other label
  res@lbPerimOn               = False             ; default has box
  res@lbOrientation           = "Horizontal"       ; vertical label bar
  res@pmLabelBarHeightF       = 0.12
  res@pmLabelBarWidthF        = 1.10
  res@pmLabelBarParallelPosF  = 0.57
  res@pmLabelBarOrthogonalPosF= 0.02

  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYRMode                = "Automatic"
  res@tiYAxisString           = "Pressure (hPa)"

 ;res@trYMinF                 = 100.0                   ; min value on y-axis
 ;res@trYMaxF                 = 1000.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/ 1000, 850, 700,500,400,300,200,100/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings
  res@trYLog                  = False

  FontHeightF = 0.035
  res@tiMainFontThicknessF    = 1.0
  res@tiMainFontHeightF       = FontHeightF*1.4
  res@tmYLLabelFontHeightF    = FontHeightF
  res@tmXBLabelFontHeightF    = FontHeightF
  res@tiXAxisFontHeightF      = FontHeightF*1.1
  res@tiYAxisFontHeightF      = FontHeightF*1.1
  res@gsnStringFontHeightF    = FontHeightF*1.1
  res@cnLineLabelFontHeightF  = FontHeightF
  res@lbLabelFontHeightF      = FontHeightF*0.9
  res@lbTitleFontHeightF      = FontHeightF*0.95                             ; font height

   BorderThick = 1.0
   res@tmBorderThicknessF     = BorderThick
   res@tmXBMajorThicknessF    = BorderThick
   res@tmXBMinorThicknessF    = BorderThick*0.5
   res@tmYLMajorThicknessF    = BorderThick
   res@tmYLMinorThicknessF    = BorderThick*0.5
   res@tmXTMajorThicknessF    = BorderThick
   res@tmXTMinorThicknessF    = BorderThick*0.5
   res@tmYRMajorThicknessF    = BorderThick
   res@tmYRMinorThicknessF    = BorderThick*0.5
 
   res@gsnLeftStringOrthogonalPosF    = 0.02
  ;res@gsnLeftStringParallelPosF      = 0.2
   res@gsnRightStringOrthogonalPosF   = 0.02
   res@tiMainString                   = ""
   res@tiYAxisString                  = ""

  if(iv.gt.0) then
    res@tmYLOn               = True
    res@tmYLLabelsOn         = True 
    res@tiYAxisString        = ""
    res@gsnLeftStringOrthogonalPosF    = 0.02
   ;res@gsnLeftStringParallelPosF      = 0.2
    res@gsnRightStringOrthogonalPosF   = 0.02
  else
    res@tmYLOn               = True
    res@tmYLLabelsOn         = True
    res@tiYAxisString        = "Pressure (hPa)"
    res@gsnLeftStringOrthogonalPosF    = 0.02
   ;res@gsnLeftStringParallelPosF      = 0.2
    res@gsnRightStringOrthogonalPosF   = 0.02
  end if

;---Resources for plotting original data
   res1=res
   res1@gsnLeftString                  = PltNums(iv)+ " " +Groups(igp) 
   res1@gsnRightString                 = Varstr+" ("+ Varunt +")" ;" ;" ("+Varunt+")";""
   res1@cnFillPalette                  = ocname
   res1@cnLevelSelectionMode           = "ExplicitLevels"  
   res1@cnLevels                       = orglev
   res1@cnFillColors                   = orgcolors
   res1@lbBoxLinesOn                   = True
  ;res1@lbBoxSeparatorLinesOn          = False  
  ;res1@lbLabelStride                  = 15
   res1@cnLinesOn                      = plot_contour       
   res1@cnLineLabelsOn                 = False     
   res1@cnInfoLabelOn                  = False
   res1@cnLineLabelDensityF            = 2.0
   res1@cnLineLabelFontHeightF         = 0.018
   res1@cnLineLabelBackgroundColor     = -1
   res1@gsnContourLineThicknessesScale = 1.0
   res1@cnLineThicknessF               = 0.5
   res1@cnLineColor                    = "Black"
   res1@gsnContourZeroLineThicknessF   = 1.
  ;res1@gsnContourNegLineDashPattern   = 1
   res1@lbTitleString                  = "" ;"("+ Varunt+")" 
   plot(iv) = gsn_csm_pres_hgt(wks1,vobs,res1)
   delete([/res,res1/])
 end do

 res_L                     = True
 res_L@gsnFrame            = True 
 res_L@gsnMaximize         = False              ; Maximize plot
 res_L@gsnPanelLabelBar    = False 
 res_L@gsnPanelMainString  = ""
 res_L@gsnPanelMainFont    = "" ;
 res_L@gsnPanelXF          = (/0.15,0.36,0.57,0.78/)    ; Adjust rightmost plots

 gsn_panel(wks1,plot,(/1,nvars/),res_L)

 end do ;
 end do ;end of season

end
