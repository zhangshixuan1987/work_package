
;----------------------------------------------------------------------
; This procedure adds a labelbar to the top right corner of the given
; plot, given a list of labels, fill colors, and an opacity value.
;
; The lbFillOpacityF resource was added in NCL V6.4.0.
;----------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,labels,colors,patcol,bopc,pattn1,pattn2)
local lbres, vph, vpw, lbid, amres
begin
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)

  lbres                      = True          ; labelbar only resources
  lbres@lbAutoManage         = True          ; Necessary to control sizes
  lbres@vpWidthF             = 0.28 * vpw     ; labelbar width
  lbres@vpHeightF            = 0.14 * vph     ; labelbar height
  lbres@lbBoxMajorExtentF    = 0.30         ; puts space between color boxes
  lbres@lbBoxMinorExtentF    = 0.18
  lbres@lbBoxLineThicknessF := 0.5
  lbres@lbBoxLineDashPattern:= 0
  lbres@lbBoxLineColor      := "black"
  lbres@lbFillColors        := colors     ; labelbar colors
  lbres@lbLabelFontHeightF   = 0.20          ; font height. default is small
  lbres@lbLabelJust          = "CenterLeft"  ; left justify labels
  lbres@lbLabelAutoStride    = False
  lbres@lbLabelStride        = 1
  lbres@lbFillPatterns      := pattn1
  lbres@lbFillOpacityF      := bopc
  lbid = gsn_create_labelbar(wks,nboxes,labels(:),lbres)

  lbres0 = lbres
  lbres0@lbFillColors         := patcol
  lbres0@lbFillOpacityF       := bopc
  lbres0@lbFillPatterns       := pattn2
  lbres0@lbFillScaleF         := 0.4
  lbres0@lbFillLineThicknessF := 0.5
  lbid0 = gsn_create_labelbar(wks,nboxes,labels(:),lbres0)

;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;
  amres                  = True
  amres@amJust           = "TopRight"
  amres@amParallelPosF   =  0.595  * vpw
  amres@amOrthogonalPosF = -0.30 * vph 
  plot@annoif = gsn_add_annotation(plot,lbid0,amres)
  plot@annoid = gsn_add_annotation(plot,lbid,amres)

end

