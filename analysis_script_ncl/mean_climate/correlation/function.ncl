;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;This script includes the fuctions used to plot the TAC and SAC plots
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;----------------------------------------------------------------------
; Given the number of rows, columns, latitudes, and longitudes,
; create a dummy 4D array of data. Then use this to calculate the
; min, max, average, and standard deviation across all lat/lon 
; subsets.  This is the 4 x nrow x ncol array that gets returned.
;----------------------------------------------------------------------
undef("create_dummy_data")
function create_dummy_data(nrow,ncol,nlat,nlon)
local data, nr, nc, dmin, dmax
begin
  data = new((/nrow,ncol,nlat,nlon/),float)
  do nr=0,nrow-1
    do nc=0,ncol-1
      dmin = random_uniform(0,10,1)      ; Use slightly different mins and maxs every time to get a
      dmax = random_uniform(90,100,1)    ; more 'random' range of values for each row and column.
      data(nr,nc,:,:) = random_uniform(dmin,dmax,(/nlat,nlon/))
    end do
  end do

  data_table = new((/4,nrow,ncol/),typeof(data)) ; min, max, avg, stdev
  do nr=0,nrow-1
    do nc=0,ncol-1
      data_table(0,nr,nc) = min(data(nr,nc,:,:))
      data_table(1,nr,nc) = max(data(nr,nc,:,:))
      data_table(2,nr,nc) = avg(data(nr,nc,:,:))
      data_table(3,nr,nc) = stddev(data(nr,nc,:,:))
    end do
  end do
  return(data_table)
end


;----------------------------------------------------------------------
; This function is identical to get_color_index, except it allows the
; input values to be an array instead of just scalar. This speeds up 
; the code significantly from the old version of this function.
;----------------------------------------------------------------------
undef("get_color_index_mod")
function get_color_index_mod(cmapt,cnlvls[*]:numeric,values:numeric)
local cmap, dims, ncn, nclr, color, n, col_indexes, values1d, indexes
begin

 if(isstring(cmapt)) then
    cmap = read_colormap_file(cmapt)
 else if(isnumeric(cmapt)) then
   dims = dimsizes(cmapt)
   if(dimsizes(dims).ne.2.or.dims(0).lt.3.or.dims(0).gt.256.or.\
       .not.any(dims(1).ne.(/3,4/))) then
     print ("get_color_index_mod: Error: cmap must be an n x 3 or n x 4 array of RGB or RGBA values, or a valid color map name")
     return(new(3,"float"))    ; return a missing value
   end if
   cmap = cmapt
 else
   print ("get_color_index_mod: Error: cmap must be an n x 3 or n x 4 array of RGB or RGBA values, or a valid color map name")
 end if
 end if

 ncn  = dimsizes (cnlvls)
 nclr = dimsizes (cmap(:,0))
 imsg = new(1,integer)          ; missing value

 if (nclr .lt. ncn+1) then 
   print ("get_color_index_mod: Warning: Not enough colors in colormap for number of requested levels")
   print ("         Colors will be repeated")
 end if
 if (any(ismissing(values))) then
   print ("get_color_index_mod: Error: More or more input values are missing.")
   print("    Returning missing.")
   return (imsg)
 end if
 if (any(ismissing(cnlvls))) then
   print ("get_color_index_mod: Error: One or more input contour levels are missing.")
   print("    Returning missing.")
   return (imsg)
 end if

;---Get nice span of indexes throughout the color map
 col_indexes = span_color_indexes(cmap,ncn+1)

;
; Convert values array to 1D so we can work with it easier. Will 
; convert back to nD when done.
;
 values1d    = ndtooned(values)
 values_dims = dimsizes(values)
 nvals       = dimsizes(values1d)
 indexes     = new(nvals,integer)

;---first level
 ii = ind(values1d .lt. cnlvls(0))
 if(.not.ismissing(ii(0))) then
    indexes(ii) = 0
 end if
;---middle levels
 do n = 1, ncn-1
   ii := ind(values1d .ge. cnlvls(n-1) .and. values1d .lt. cnlvls(n))
   if(.not.ismissing(ii(0))) then
     indexes(ii) = n
   end if
 end do
;---last level
 ii := ind(values1d .ge. cnlvls(ncn-1))
 if(.not.ismissing(ii(0))) then
    indexes(ii) = ncn
 end if

 if(isstring(cmapt)) then
   col_indexes(indexes) = col_indexes(indexes)+2    ; Account for 0/1 index being dropped
 end if
 return(reshape(col_indexes(indexes),values_dims))
end

;----------------------------------------------------------------------
; This function is identical to get_color_rgba, except it allows the
; input values to be an array instead of just scalar. This speeds up 
; the code significantly from the old version of this function.
;----------------------------------------------------------------------
undef("get_color_rgba_mod")
function get_color_rgba_mod(cmapt,cnlvls[*]:numeric,values:numeric)
local fmsg, icol, cmap, values_dims, values_rank, return_dims
begin
 fmsg = new(4,float)          ; missing value
 icol = get_color_index_mod(cmapt,cnlvls,ndtooned(values))

 if(isscalar(values)) then
   return_dims = 4
 else
   values_dims = dimsizes(values)
   values_rank = dimsizes(values_dims)
   return_dims = new(values_rank+1,typeof(values_dims))
   return_dims(0:values_rank-1) = values_dims
   return_dims(values_rank)     = 4
 end if

 if (any(ismissing(icol))) then
   return (fmsg)
 end if

 if(isstring(cmapt)) then
   cmap = read_colormap_file(cmapt)
   icol = icol-2    ; Indexes returned by get_color_index_mod start at 2
 else
   cmap = cmapt
 end if
 return(reshape(cmap(icol,:),return_dims))
end

;----------------------------------------------------------------------
; Given a plot, an nrow x ncol matrix of data, a set of levels, a
; colormap, and a quadrant, draw a triangle filled in an appropriate
; color in the given quadrant.
;----------------------------------------------------------------------
undef("add_filled_blocks")
procedure add_filled_blocks(wks,plot,data,levels,cmap,location,var_group_num,exp_group_num)
local dims, gnres, nc, nr, nrow, ncol, xtri, ytri, xtri2d, ytri2d, irows, icols,ivp,nmgrp,nvgrp
begin

  dims = dimsizes(data)
  nrow = dims(0)
  ncol = dims(1)
  ytri2d = new((/nrow,5/),float)
  xtri2d = new((/ncol,5/),float)
  
  irows = ispan(0,nrow-1,1)
  icols = ispan(0,ncol-1,1)

  nvgrp  = dimsizes(var_group_num)
  do ivp = 0,nvgrp,1
    if(ivp.eq.0)then 
      ist  = 0
      ied  = var_group_num(ivp)
    else if (ivp.eq.nvgrp)then
      ist  = var_group_num(ivp-1)+1
      ied  = ncol-1
    else 
      ist  = var_group_num(ivp-1)+1
      ied  = var_group_num(ivp)
    end if 
    end if 
    xtri2d(ist:ied,0) = icols(ist:ied)+ivp*0.3
    xtri2d(ist:ied,1) = icols(ist:ied)+1.0+ivp*0.3
    xtri2d(ist:ied,2) = icols(ist:ied)+1.0+ivp*0.3
    xtri2d(ist:ied,3) = icols(ist:ied)+ivp*0.3
    xtri2d(ist:ied,4) = icols(ist:ied)+ivp*0.3
  end do 

  nmgrp = dimsizes(exp_group_num)
  do ivp = 0,nmgrp,1
    if(ivp.eq.0)then
      ist  = 0
      ied  = exp_group_num(ivp)
    else if (ivp.eq.nmgrp)then
      ist  = exp_group_num(ivp-1)+1
      ied  = nrow-1
    else
      ist  = exp_group_num(ivp-1)+1
      ied  = exp_group_num(ivp)
    end if
    end if

    ytri2d(ist:ied,0) = irows(ist:ied)+ivp*0.2
    ytri2d(ist:ied,1) = irows(ist:ied)+ivp*0.2
    ytri2d(ist:ied,2) = irows(ist:ied)+1.0+ivp*0.2
    ytri2d(ist:ied,3) = irows(ist:ied)+1.0+ivp*0.2
    ytri2d(ist:ied,4) = irows(ist:ied)+ivp*0.2
 
  end do 

  xtri             = conform_dims((/nrow,ncol,5/),xtri2d,(/1,2/))
  ytri             = conform_dims((/nrow,ncol,5/),ytri2d,(/0,2/))

  gnres            = True                     ; resource list for filled polygons
  gnres@gsEdgesOn  = True                     ; outline each filled triangle
  gnres@gsSegments = ispan(0,nrow*ncol*5,5)   ; this makes code faster
  gnres@gsEdgeThicknessF = 0.5

;--If you have NCL V6.4.0 or older, you must use the "_mod" version of the function
; tmp  = reshape(get_color_rgba_mod(cmap,levels,data),(/product(dimsizes(data)),4/))

;---This will only work in NCL V6.5.0 or later.
  tmp  = reshape(get_color_rgba(cmap,levels,data),(/product(dimsizes(data)),4/))

 ;printVarSummary(tmp)
 ;print(ncol*nrow)

;---For missing values of TRMM/PERSIANN, RGB = white
  if (location .eq. "Polar") then
     do i = 0,nrow-1,1
       ist = (i+1)*ncol-2
       ied = (i+1)*ncol-1
       ;print(ist+ " " +ied)
       tmp(ist:ied,0:2) = 211./255.
      end do
  end if

  gnres@gsColors = tmp

  tmpstr = unique_string(location)
  plot@$tmpstr$ = gsn_add_polygon(wks,plot,ndtooned(xtri),ndtooned(ytri),gnres)
;
; Add text strings showing each value. This is for debug purposes.
;
; If you have a lot of variables and/or models, the text strings
; will get too cluttered and will cover up the filled triangles.
; In this case, you may want to make the font height smaller,
; or only label every other triangle, or similar.
;
  ADD_TEXT = False
  if(ADD_TEXT) then
    txres = True
    txres@txFontHeightF = 0.01
    txres@txBackgroundFillColor = "white"
    txres@txPerimOn             = True
    if(location.eq."bot".or.location.eq."top") then
      xmin = dim_min(xtri)
      xmax = dim_max(xtri)
      xpos = (xmin+xmax)*0.5
      if(location.eq."bot") then
        ypos = dim_min(ytri)
        txres@txJust = "BottomCenter"
      else
        ypos = dim_max(ytri)
        txres@txJust = "TopCenter"
      end if
    else
      ymin = dim_min(ytri)
      ymax = dim_max(ytri)
      ypos = (ymin+ymax)*0.5
      if(location.eq."lft") then
        xpos = dim_min(xtri)
        txres@txJust = "CenterLeft"
      else
        xpos = dim_max(xtri)
        txres@txJust = "CenterRight"
      end if
    end if
    tmpstr = unique_string(location+"_text")
    data_str = sprintf("%6.2f",data)
    plot@$tmpstr$ = gsn_add_text(wks,plot,data_str,xpos,ypos,txres)
  end if
end

;----------------------------------------------------------------------
; Attach a horizontal or vertical labelbar to the bottom or right 
; of a plot.
;----------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,cmap,labels,title,location)
local nboxes, vph, vpw, nboxes, lbres, lbid, amres, annoid
begin
  nboxes = dimsizes(labels)+1  
  colors = span_color_rgba(cmap,nboxes)
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

;
; Set some resources for a vertical or horizontal labelbar on the
; bottom or right axis. 
;
; am_para/am_orth
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
; You will likely need to modify the am_para/am_orth 
; values depending on your X and Y axis labels, the size
; of your plot, the number of rows and columns you have, et.
;
  if(any(location.eq.(/"bot1","bot2"/))) then
    orient  = "horizontal"    
    width   = vpw * 0.85   ; slightly shorter than width of plot
    height  = vph * 0.25
    am_just = "BottomCenter"
    am_para = 0.0         ; Centered about X axis
    title_pos = "Top"
    if(location.eq."bot1")
      am_orth = 1.12  ;0.85       ; Move labelbar down
    else
      am_orth = 1.04       ; move further down
    end if 
  else
    orient  = "vertical"
    width   = 0.30 ; :vpw * 0.12
    height  = 1.5  ; vph * 0.95    ; slightly shorter than height of plot
    am_just = "CenterRight"
    am_orth =  0.0 ;0.35 ; 0. ; -0.5          ; Move labelbar up
    title_pos = "Left"
    if(location.eq."rgt2")
      am_para = 0.68 ; 0.65        ; Move labelbar right
    else
      am_para = 3.         ; move further right
    end if 
  end if

  height = height*1.50   ; Hui Wan: make label bar taller

; labelbar resources
  lbres                     = True          
  lbres@lbAutoManage        = False         ; Necessary to control sizes
  lbres@vpWidthF            = width
  lbres@vpHeightF           = height
  lbres@lbFillColors        = colors        ; labelbar colors
  lbres@lbMonoFillPattern   = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF  = 0.050 *1.20          ; font height. Default is small
  lbres@lbLabelFont         = "helvetica"
  lbres@lbLabelAlignment    = "InteriorEdges"
  lbres@lbOrientation       = orient
  lbres@lbTitleString       = title
  lbres@lbPerimOn           = False
  lbres@lbTitlePosition     = title_pos
; lbres@lbBoxEndCapStyle    = "TriangleBothEnds"
  lbres@lbBoxLineThicknessF = 0.5

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)

; annotation resources
  amres                  = True
  amres@amJust           = am_just
  amres@amOrthogonalPosF = am_orth
  amres@amParallelPosF   = am_para

; attach the labelbar to the plot
  plot@annoid = gsn_add_annotation(plot,lbid,amres)
end

;----------------------------------------------------------------------
; Given a color map name and a start and end index, return the
; subsetted color map.
;----------------------------------------------------------------------
undef("subset_cmap")
function subset_cmap(name,istart,iend)
begin
  cmap = read_colormap_file(name)
  return(cmap(istart:iend,:))
end

;----------------------------------------------------------------------
; Create a blank plot that we'll add filled triangles to later.
; The row_labels and col_labels will be used to label the X and
; Y axis, and to determine the min/max of both axes.
;----------------------------------------------------------------------
undef("create_plot")
function create_plot(wks,row_labels,col_labels,title,var_group_num,exp_group_num)
local res, nrow, ncol,nmgrp,nvgrp,ivp
begin
  nmgrp = dimsizes(exp_group_num)
  nvgrp = dimsizes(var_group_num)

  nrow = dimsizes(row_labels)
  ncol = dimsizes(col_labels)

  res             =  True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@gsnScale    = False
  res@tmYROn      = False
  res@tmXTOn      = False
  res@vpWidthF    = 1.2  * 2.0
  res@vpHeightF   = 0.50 * 1.0

;---Set the min/max for X axis
  res@trXMinF = 0 
  res@trXMaxF = ncol+1.5
  res@trYMinF = 0
  res@trYMaxF = nrow+0.2

;---Customize X and Y axis tickmark labels
  res@tmXBMode    = "Explicit"
  xticks          = ispan(0,ncol-1,1) + 0.5

  do ivp = 0,nvgrp,1
    if(ivp.eq.0)then
      ist  = 0
      ied  = exp_group_num(ivp)
      xticks(ist:ied) = xticks(ist:ied)+ivp*0.2
    else if (ivp.eq.nvgrp)then
      ist  = var_group_num(ivp-1)+1
      ied  = ncol-1
      xticks(ist:ied) = xticks(ist:ied)+ivp*0.3
    else
      ist  = var_group_num(ivp-1)+1
      ied  = var_group_num(ivp)
      xticks(ist:ied) = xticks(ist:ied)+ivp*0.3
    end if
    end if
  end do 

  res@tmXBValues  = xticks
  res@tmXBLabels  = col_labels

  res@tmYLMode    = "Explicit"
  yticks          = ispan(0,nrow-1,1)+ 0.5

  do ivp = 0,nmgrp,1
    if(ivp.eq.0)then
      ist  = 0
      ied  = exp_group_num(ivp)
    else if (ivp.eq.nmgrp)then
      ist  = exp_group_num(ivp-1) +1
      ied  = nrow-1
    else
      ist  = exp_group_num(ivp-1) +1
      ied  = exp_group_num(ivp)
    end if
    end if
    yticks(ist:ied) = yticks(ist:ied)+ivp*0.2
  end do

  res@tmYLValues  = yticks
  res@tmYLLabels  = row_labels


;---Turn off tick marks only (not the labels)
  fontheight = 0.050
  res@tmXBMajorOutwardLengthF     = 0.03
  res@tmXBMajorLengthF            = 0.03
  res@tmXBMajorThicknessF         = 0.5
  res@tmYLMajorOutwardLengthF     = 0.03
  res@tmYLMajorLengthF            = 0.03
  res@tmYLMajorThicknessF         = 0.5
  res@tmXBLabelFontHeightF        = fontheight
  res@tmYLLabelFontHeightF        = fontheight
  res@tmXBLabelDeltaF  	          = -0.6
  res@tmYLLabelDeltaF  	          = -0.6
  res@tmXBLabelFont 	          = "helvetica"
  res@tmYLLabelFont 	          = "helvetica"
  res@tmXBLabelAngleF             = 90.  ; rotate
  res@tmXBLabelJust               = "CenterRight"     ; helpful if you rotate labels

  res@tiMainString                = ""
  res@tiMainFontHeightF           = fontheight*1.5
  res@tiMainFont	          = "helvetica" 
  res@gsnLeftString	          = title
  res@gsnLeftStringFontHeightF    = fontheight*1.5
  res@gsnLeftStringOrthogonalPosF = 0.03
;---Create a blank plot so gsn_add_polygon function can be used
  res@tmBorderThicknessF          = 0.5
; res@tmXBBorderOn                = False
; res@tmYLBorderOn                = False
; res@tmXTBorderOn                = False
; res@tmYRBorderOn                = False
  plot = gsn_csm_blank_plot(wks,res)

  return(plot)
end

;----------------------------------------------------------------------
; This function uses nice_mnmxintvl to get a nice set of 18 or 
; fewer levels that span the given data. It uses the unadvertised
; testspan to span these values. testspan is UNSUPPORTED!
;----------------------------------------------------------------------
function get_nice_levels(data)
local mnmxin
begin
  mnmxin = nice_mnmxintvl(min(data),max(data),18,False)
  return(testspan(mnmxin(0),mnmxin(1),mnmxin(2)))
end

