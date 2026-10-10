;----------------------------------------------------------------------
; table_8.ncl
;
; Concepts illustrated:
;   - Drawing a table of filled triangles
;   - Using gsn_csm_blank_plot to get custom tickmark labels
;   - Using "getvalues" to retrieve the min/max of axes
;   - Turning off the top and right tickmarks
;   - Using functions for cleaner code
;   - Attaching multiple custom labelbars to a plot
;   - Attaching filled triangles to a plot
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "./function.ncl"
;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
  start_time = get_cpu_time()
;
; Generate a 4 x nvars x nmodels dummy array
;     (0,:,:)  - data mins across each nlat x nlon array
;     (1,:,:)  - data maxs     ""
;     (2,:,:)  - data averages ""
;     (3,:,:)  - data stddevs   ""

  var_list      = (/"U100","U200","U500","U700","U850",   \
                    "T100","T200","T500","T700","T850",   \
                    "Q100","Q200","Q500","Q700","Q850",   \
                    "~F33~w~N~200","~F33~w~N~500","~F33~w~N~850", \
                    "AVHRR","HIRS","ERA5","GPCP","TRMM"/)
  dim0          = dimsizes(var_list)

  exp_list      = (/"CLIM","NDG_UV","NDG_UVT","NDG_UVTQ"/) ;"NDG_UVTQ_SRF1","NDG_UVT_SRF2"/)
  dim1          = dimsizes(exp_list)

  regions       = (/"Polar", "MidLat", "Tropics"/)
  reg_labels    = (/"polar regions", "mid-latitudes", "tropics"/)
  dim2          = dimsizes(regions)

  sname         = (/"(a) Spatial correlation, polar regions", \
                    "(d) Temporal correlation, polar regions", \
                    "(b) Spatial correlation, mid-latitudes", \
                    "(e) Temporal correlation, mid-latitudes", \
                    "(c) Spatial correlation, tropics", \
                    "(f) Temporal correlation, tropics" /)

  type          = (/"sac","tac"/)
  dim3          = dimsizes(type)

  leftstr       = (/"polar region                mid-latitudes                tropics" /) 

  plot          = new(dim3*dim2, graphic)
  dum           = new((/5,dim0,dim3*dim2/),graphic)

  fname         = "E3SM_nudge_acc"
  wtype         = "pdf"
  wks           = gsn_open_wks(wtype, fname)

  ;---Set nice levels to use for each of the 4 quantities
  ;levels       = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9/) ; get_nice_levels(data)
   levels       = (/0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95/)
  ;---Set colormaps to use for each of the 4 quantities
  ;  min_cmap = subset_cmap("MPL_Blues",18,127)
  ;cmap         = read_colormap_file("amwg_blueyellowred")
  cmap1         = read_colormap_file("radar")
  cmap2         = read_colormap_file("cb_9step")
  cmap3         = read_colormap_file("radar_1")
  cmap4         = read_colormap_file("prcp_2")
  cmap5         = read_colormap_file("StepSeq25")
  cmap6         = read_colormap_file("tbrAvg1")
  cmap          = (/cmap5(15,:),cmap5(17,:),cmap5(19,:),cmap5(14,:),cmap6(64,:),cmap1(7,:),\
                    cmap1(8,:),cmap1(10,:),cmap1(12,:),cmap1(13,:)/) 

  var_group_num = (/4,9,14,17,19/)
  exp_group_num = (/0/)

  do i = 0, dim3-1

    f    = addfile("./"+type(i)+"_3regions_data.nc", "r")

    do j = 0, dim2-1

      data = f->$regions(j)$
      data = where(ismissing(data), -1, data)
     
      ;---Create a blank plot that we can add triangles to later.
      row_labels  = exp_list
      col_labels  = var_list
      if(type(i).eq."tac")then
        row_labels   = ""
      end if 
       plot(j*dim3+i) = create_plot(wks,row_labels,col_labels,sname(j*dim3+i),var_group_num,exp_group_num)

       ;---Add filled triangles in the given quadrant
       add_filled_blocks(wks,plot(j*dim3+i),data(:,:),levels,cmap,regions(j),var_group_num,exp_group_num)
 
       ;---Add labelbars for each of the four quantities.
       if  (j*dim3+i .eq. 3) then
        add_labelbar(wks,plot(j*dim3+i),cmap,""+sprintf("%5.2f",levels),"","rgt2")
       end if 
       delete(data)

     end do 

     delete(f)

  end do

  ;draw text for OLR and PRECT
  txres               = True
  txres@txFontHeightF = 0.0055
  txres@txFont        = "helvetica"
  gsres               = True
  gsres@gsLineColor   = "black"

  do i = 0, dim3-1,1

    do j = 0,dim2-1,1 

     do k = 0,1

      if (k.eq.0) then 
        label             = "OLR"
      else
        label             = "PRECT"
      end if 

      a                   = 0.443 + 0.37*i + 0.0390 * k
      b                   = 0.223 + 0.163*j 

      if (k.eq.0) then
       c                   = 0.007
       d                   = 0.009
       e                   = 0.005
      else
       c                   = 0.007
       d                   = 0.013
       e                   = 0.005
      end if  
      gsn_text_ndc(wks,label, a, b ,txres)
      gsn_polyline_ndc(wks,(/a-d,a-d-e/),  (/b,b/)  ,gsres)
      gsn_polyline_ndc(wks,(/a-d-e,a-d-e/),(/b+c,b/),gsres)
      gsn_polyline_ndc(wks,(/a+d,a+d+e/),  (/b,b/)  ,gsres)
      gsn_polyline_ndc(wks,(/a+d+e,a+d+e/),(/b+c,b/),gsres)

     end do 

    end do 
 
  end do 

  resP                            = True
  resP@gsnMaximize                = False
  resP@gsnFrame                   = True
  resP@lbLabelAutoStride          = False
  resP@gsnPanelMainString         = " "
  resP@gsnPanelRight              = 0.9
  resP@gsnPanelXWhiteSpacePercent = 4.0
  resP@gsnPanelYWhiteSpacePercent = 5.0
  resP@gsnPanelXF                 = (/0.15,0.52, \
                                      0.15,0.52, \
                                      0.15,0.52/)

  gsn_panel(wks,plot(:),(/dim2,dim3/),resP)
 
end
