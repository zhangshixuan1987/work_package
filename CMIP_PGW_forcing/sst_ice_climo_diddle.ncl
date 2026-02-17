;***********************************************************
; Create "diddled" SST forcing from B compset (fully coupled)
;   for use in F compset (prescribed-SST) run.
;***********************************************************
; Auhtor:  Jim Benedict jjb278@gmail.com
; Origin:  6 Apr 2018
; Modified:
;     
; Keywords:  SST, diddle, compset, AMIP


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"


begin
  
  print("")
  print("--- IN NCL PROGRAM: " + get_script_name() + " ---")
  
  
  ; ======================================================================================
  ; ======================================================================================
  ; USER SPECIFICATIONS
  caseName       = getenv("caseName")           ; case from which SSTs will be taken
  inSST_File	 = getenv("INPUT_SST_FILE")
  inICE_File	 = getenv("INPUT_SEAICE_FILE")
  outSSTICE_File = getenv("OUTPUT_SSTICE_FILE")
  
  months30days = False     ; True if we approximate each month as 30 days long, False if we
                           ;   use the actual number of days (for "noleap" calendar)
  
  write_NC   = True
  test_print = True
  workPrefix    = "./work"
  
  ilon = 100    ; for test printing only  ex: ilon = 105, ilat = 190: 0.9N/78.7E
  ilat = 100    ; for test printing only
  imn  = 10     ; for test printing only (test month)
  imon = 10
  
  checkAdjustmentTest = False    ; Print AMIP-II adjustment test results from test data?
  checkAdjustment     = True     ; Print AMIP-II adjustment test results from actual data?		
  
  
  
  ; ======================================================================================
  ; ======================================================================================
  ; PREP WORK
  
  ; ------------------------------------------------------
  ; For time adjustments
  nyr = 1
  
  monNames  = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  monDays   = (/  31.,  28.,  31.,  30.,  31.,  30.,  31.,  31.,  30.,  31.,  30.,  31./)
  if(sum(monDays).ne.365) then
    print("Must use noleap calendar. Exiting.")
    exit
  end if
  midMnDays = new(dimsizes(monDays),"float",-9999)   ; mid-month days: 15.5, 45.0, 74.5, ..., 349.5
  begMnDays = midMnDays                              ; beginning-of-month days: 1.0, 32.0, 60.0, ...
  endMnDays = midMnDays                              ; end-of-month days: 31.0, 59.0, 90.0. ... 
  
  midMnDays(0) = monDays(0) / 2.
  begMnDays(0) = 1.
  endMnDays(0) = monDays(0)
  do i = 1,dimsizes(monDays)-1
    midMnDays(i) = (2.*sum(monDays(:i-1))+monDays(i))/2.
    endMnDays(i) = sum(monDays(:i))
    begMnDays(i) = endMnDays(i-1) + 1.
  end do
  midMnDays@units    = "days since 0001-01-01 00:00:00"
  midMnDays@calendar = "noleap"
  monDayShift = endMnDays - midMnDays
  fmt = "%6.1f"
  print("Actual month middays, beginning days, and end days:")
  print( monNames + sprintf(fmt,midMnDays) + sprintf(fmt,begMnDays) + sprintf(fmt,endMnDays) )
  
  ; BUFFERED annual cycle (needed as linear equation solver inputs):
  BmonNames  = (/"Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan"/)
  BmonDays   = (/  31.,  31.,  28.,  31.,  30.,  31.,  30.,  31.,  31.,  30.,  31.,  30.,  31.,  31./)
  BmidMnDays = new(dimsizes(BmonDays),"float",-9999)    ; mid-month days: 15.5, 46.5, 76.0, ..., 411.5
  BbegMnDays = BmidMnDays                               ; beginning-of-month days: 1.0, 32.0, 63.0, ...
  BendMnDays = BmidMnDays                               ; end-of-month days: 31.0, 62.0, 90.0. ...
  BmidMnDays(0) = BmonDays(0) / 2.
  BbegMnDays(0) = 1.
  BendMnDays(0) = BmonDays(0)
  do i = 1,dimsizes(BmonDays)-1
    BmidMnDays(i) = (2.*sum(BmonDays(:i-1))+BmonDays(i))/2.
    BendMnDays(i) = sum(BmonDays(:i))
    BbegMnDays(i) = BendMnDays(i-1) + 1.
  end do
  print("BUFFERED month middays, beginning days, and end days:")
  print( BmonNames + sprintf(fmt,BmidMnDays) + sprintf(fmt,BbegMnDays) + sprintf(fmt,BendMnDays) )
  
  tiMMD   = calendar_decode2(midMnDays, -5)
  date    = 100 * tiMMD(:,1) + tiMMD(:,2)
  datesec = tiMMD(:,3) * 3600
  
  print("")
  print("Seasonal cycle date and datesec:")
  print( sprinti("%6i",date) + sprinti("%6i",datesec) )
  print("")

  ; -----------------------------------------------------------------------------------
  ; Create weighting coefficients that will be applied to climatological monthly
  ;   mean correction terms so as to preserve the actual monthly mean after the data
  ;   are interpolated to daily data within CESM
  ;
  ;   See, e.g., Sheng & Zwiers (1998, Climate Dynamics) or 
  ;   http://www-pcmdi.llnl.gov/publications/pdf/60.pdf from the AMIP-II website
  ;   Because we are applying the *climatological* seasonal cycle of ocean temperature
  ;     and salinity correction terms to our model, we need only adjust the pre-calculated
  ;     climatological values rather than adjust each individual monthly mean.  If the user
  ;     wants to do the monthly mean adjustments to each individual month, that procedure
  ;     would require that all monthly values are stored in one array such that code
  ;     changes would be needed in the "nf" loop below.
  
  ; Construct coefficient matrix containing averaging weights
  nT    = dimsizes(BmonNames)
  lineqn_coef = new((/nT,nT/),"float",-1.E19)
  lineqn_coef = 0.
  
  if(months30days) then       ; if months each have 30 days (simplification)
  
    lineqn_coef(0,nT-1) = 1./8.
    lineqn_coef(0,0)    = 3./4.
    lineqn_coef(0,1)    = 1./8.
    do i = 1,nT-2
      lineqn_coef(i,i-1) = 1./8.
      lineqn_coef(i,i)   = 3./4.
      lineqn_coef(i,i+1) = 1./8.
    end do
    lineqn_coef(nT-1,nT-2) = 1./8.
    lineqn_coef(nT-1,nT-1) = 3./4.
    lineqn_coef(nT-1,0)    = 1./8.
  
  else      ; else assume months with varying number of days (for "noleap" calendar)
  
    time_Mj_minus_Mjm1   = new(nT,"float",-1.E19)
    time_Mj_minus_Bj     = time_Mj_minus_Mjm1       ; to quickly transfer dimensional data only
    time_Bj_minus_Mjm1   = time_Mj_minus_Mjm1
    time_Mjp1_minus_Mj   = time_Mj_minus_Mjm1
    time_Bjp1_minus_Mj   = time_Mj_minus_Mjm1
    time_Mjp1_minus_Bjp1 = time_Mj_minus_Mjm1
    
    ; Calculate time_Mj_minus_Mjm1
    time_Mj_minus_Mjm1(0) = midMnDays(11) - midMnDays(10)   ; midpoint(Dec) - midpoint(Nov)
    do i = 1,nT-1
      time_Mj_minus_Mjm1(i) = BmidMnDays(i) - BmidMnDays(i-1)
    end do
    
    ; Calculate time_Mj_minus_Bj
    do i = 0,nT-1
      time_Mj_minus_Bj(i) = BmidMnDays(i) - BbegMnDays(i)
    end do
    
    ; Calculate time_Bj_minus_Mjm1
    time_Bj_minus_Mjm1(0) = begMnDays(11) - midMnDays(10)   ; beginning(Dec) - midpoint(Nov)
    do i = 1,nT-1
      time_Bj_minus_Mjm1(i) = BbegMnDays(i) - BmidMnDays(i-1)
    end do
    
    ; Calculate time_Mjp1_minus_Mj
    do i = 0,nT-2
      time_Mjp1_minus_Mj(i) = BmidMnDays(i+1) - BmidMnDays(i)
    end do
    time_Mjp1_minus_Mj(nT-1) = midMnDays(1) - midMnDays(0)  ; midpoint(Feb) - midpoint(Jan)
    
    ; Calculate time_Bjp1_minus_Mj
    do i = 0,nT-2
      time_Bjp1_minus_Mj(i) = BbegMnDays(i+1) - BmidMnDays(i)
    end do
    time_Bjp1_minus_Mj(nT-1) = begMnDays(1) - midMnDays(0)  ; beginning(Feb) - midpoint(Jan)
    
    ; Calculate time_Mjp1_minus_Bjp1
    do i = 0,nT-2
      time_Mjp1_minus_Bjp1(i) = BmidMnDays(i+1) - BbegMnDays(i+1)
    end do
    time_Mjp1_minus_Bjp1(nT-1) = midMnDays(1) - begMnDays(1)  ; midpoint(Feb) - beginning(Feb)
    
    ; Assemble final coefficients --------------------
    coef_jm1 = new(nT,"float",-1.E19)
    coef_j   = coef_jm1
    coef_jp1 = coef_jm1
    
    coef_jm1 = (1./4.) * time_Mj_minus_Bj / time_Mj_minus_Mjm1
    coef_j   = ( (1./4.) * time_Bj_minus_Mjm1 / time_Mj_minus_Mjm1 ) + \
               ( (1./4.) * time_Mjp1_minus_Bjp1 / time_Mjp1_minus_Mj ) + \
               ( (1./2.) )
    coef_jp1 = (1./4.) * time_Bjp1_minus_Mj / time_Mjp1_minus_Mj    
    
    lineqn_coef(0,nT-1) = coef_jm1(0)
    lineqn_coef(0,0)    = coef_j(0)
    lineqn_coef(0,1)    = coef_jp1(0)
    do i = 1,nT-2
      lineqn_coef(i,i-1) = coef_jm1(i)
      lineqn_coef(i,i)   = coef_j(i)
      lineqn_coef(i,i+1) = coef_jp1(i)
    end do
    lineqn_coef(nT-1,nT-2) = coef_jm1(nT-1)
    lineqn_coef(nT-1,nT-1) = coef_j(nT-1)
    lineqn_coef(nT-1,0)    = coef_jp1(nT-1)
    
  end if        ; end if months30days
  
  print("")
  print("Linear equation coefficient matrix:")
  write_matrix (lineqn_coef, "14F10.4", False)
  print("")
; exit
  
  ; --------------------------------------------------------------------------------------
  ; TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
  if(checkAdjustmentTest) then
    fakeDat = (/22.3000, 24.1000, 23.4000, 25.0000, 26.8000, 28.1000, 28.9000, 27.0000, \
                29.5000, 28.8000, 26.9000, 24.6000, 22.3000, 24.1000/)
    newDat = solve_linsys(lineqn_coef, fakeDat)      ; x will be one-dimensional [x(N)]
    fmt = "%14.6f"
    print(sprintf(fmt,fakeDat) + sprintf(fmt,newDat))
  
  
    ; Check that monthly-to-daily interpolated values represent actual monthly mean values
    if(months30days) then
      timeMon30 = ispan(14,404,30)
      timeDay30 = ispan(0,419,1)
      fakeDatDay = linint1_n (timeMon30, fakeDat, False, timeDay30, 0, 0)
      newDatDay  = linint1_n (timeMon30, newDat,  False, timeDay30, 0, 0)
      print("")
      print(sprintf("%16.6g",fakeDatDay) + sprintf("%16.6g",newDatDay) )
      print("")
      m = 30
      do i = 1,nT-2
        print( sprinti("%6i",i) + \
               sprintf("%14.6f",dim_avg(fakeDatDay(m:m+29))) + sprintf("%14.6f",fakeDat(i)) + \
               sprintf("%14.6f",dim_avg(newDatDay(m:m+29))) + sprintf("%14.6f",newDat(i)) )
        m = m + 30
      end do
      exit
    else      ; else assume months with varying number of days (for "noleap" calendar)
      timeMonVar = BmidMnDays-1
      timeDayVar = ispan(0,426,1)
      fakeDatDay = linint1_n (timeMonVar, fakeDat, False, timeDayVar, 0, 0)
      newDatDay  = linint1_n (timeMonVar, newDat,  False, timeDayVar, 0, 0)
      print("")
      print(sprintf("%16.6g",fakeDatDay) + sprintf("%16.6g",newDatDay) )
      print("")
      do i = 1,nT-2
        print( sprinti("%6i",i) + \
               sprintf("%14.6f",dim_avg(fakeDatDay(toint(BbegMnDays(i)-1):toint(BendMnDays(i)-1)))) + \
               sprintf("%14.6f",fakeDat(i)) + \
               sprintf("%14.6f",dim_avg(newDatDay(toint(BbegMnDays(i)-1):toint(BendMnDays(i)-1)))) + \
               sprintf("%14.6f",newDat(i)) )
      end do
      exit
    end if        ; end if months30days
  end if          ; end if checkAdjustment = True
  ; TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
  ; --------------------------------------------------------------------------------------
  
  
  ; ======================================================================================
  ; ======================================================================================
  ; Input Source SST 12 month climatology

  F_SST = addfile(inSST_File,"r")
  sstMon = F_SST->SST_cpl
  nmon         = 12
  nlon=dimsizes(sstMon(0,0,:))
  nlat=dimsizes(sstMon(0,:,0))
  
  printVarSummary(sstMon)
  sstMon&time = ispan(1,12,1)
  
  
;  if(test_print) then
;   print("")
;   print("Before adjustment to maintain actual climatological monthly mean:")
;   print( sprintf("%8.2f",lat2d(ilat,ilon)) + sprintf("%8.2f",lon2d(ilat,ilon)) + \
;          sprintf("%16.6g",sstMon(imn,ilat,ilon)) )
;  end if
  
  
  ; -----------------------------------------------------------------------------------
  ; Adjust correction terms to preserve true monthly means (SST "diddling")
    ;   See, e.g., Sheng & Zwiers (1998, Climate Dynamics) or 
  ;   http://www-pcmdi.llnl.gov/publications/pdf/60.pdf from the AMIP-II website
  ;   Because we are applying the *climatological* seasonal cycle of ocean temperature
  ;     and salinity correction terms to our model, we need only adjust the pre-calculated
  ;     climatological values rather than adjust each individual monthly mean.  If the user
  ;     wants to do the monthly mean adjustments to each individual month, that procedure
  ;     would require that all monthly values are stored in one array such that code
  ;     changes would be needed in the "nf" loop above.
  
  ; Create "buffered" climatological annual cycle by adding the December climatological
  ;   mean value before the January value, and adding the January value after the
  ;   final December value such that you have a total of 14 months. These added
  ;   months are needed for the linear solver.
  xMon     = sstMon    ; User settings, to make it generic

  xlats = F_SST->lat
  xlons = F_SST->lon
  nlatUse  = nlat
  nlonUse  = nlon
  
  xMonBuf    = new((/nT,nlatUse,nlonUse/),typeof(xMon),xMon@_FillValue)
  xMonBufAdj = xMonBuf  ; these **avgBufAdj arrays are really only required for testing, could be deleted
  xMonBuf(0,:,:)    = xMon(11,:,:)
  xMonBuf(1:12,:,:) = xMon(0:11,:,:)
  xMonBuf(13,:,:)   = xMon(0,:,:)
  
  ; Apply pre-calculated weighting coefficients and solve linear equation
  xMonAdj = xMon
  xMonAdj = xMon@_FillValue
  print("Performing adjustment to monthly mean...")
  do j = 0,nlatUse-1
    do i = 0,nlonUse-1
      a    = xMonBuf(:,j,i)
      aAdj = solve_linsys(lineqn_coef, a)   ; aAdj will be one-dimensional [aAdj(14)]
      xMonBufAdj(:,j,i) = (/aAdj/)              ; only required for later testing, could be deleted
      xMonAdj(0:11,j,i) = (/aAdj(1:12)/)        ; Retain only Jan-Dec values, discard bookend buffers    
    end do
  end do
  
  ; Minor corrections -- (re)mask out land points
;  printMinMax(xMonAdj, True)
  xMonConf = conform(xMonAdj,xMon(0,:,:),(/1,2/))
  xMonAdj  = where(ismissing(xMonConf),xMonAdj@_FillValue,xMonAdj)
;  printMinMax(xMonAdj, True)
  
  if(test_print) then
   print("")
   print("After adjustment to maintain actual climatological monthly mean:")
   print( "     LAT" + "     LON" + "         Raw X" + "   Corrected X")
   print( sprintf("%8.2f",xlats(ilat)) + sprintf("%8.2f",xlons(ilon)) + \
          sprintf("%16.6g",xMon(imn,ilat,ilon)) + \
          sprintf("%16.6g",xMonAdj(imn,ilat,ilon)) )
  end if
  
  ; --------------------------------------------------------------------------------------
  ; TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
  if(checkAdjustment) then
    print("")
    print("Checking that adjustment procedure worked with actual (buffered) data...")
    useOrig = xMonBuf(:,ilat,ilon)
    useAdj  = xMonBufAdj(:,ilat,ilon)
    
    fmt = "%14.6g"
    print(sprintf(fmt,useOrig) + sprintf(fmt,useAdj))
    
    
    ; Check that monthly-to-daily interpolated values represent actual monthly mean values
    if(months30days) then
      timeMon30 = ispan(14,404,30)
      timeDay30 = ispan(0,419,1)
      oldDatDay = linint1_n (timeMon30, useOrig, False, timeDay30, 0, 0)
      newDatDay  = linint1_n (timeMon30, useAdj,  False, timeDay30, 0, 0)
      print("")
      print(sprintf("%16.6g",oldDatDay) + sprintf("%16.6g",newDatDay) )
      print("")
      m = 30
      do i = 1,nT-2
        print( sprinti("%6i",i) + \
               sprintf("%14.6g",dim_avg(oldDatDay(m:m+29))) + sprintf("%14.6g",useOrig) + \
               sprintf("%14.6g",dim_avg(newDatDay(m:m+29))) + sprintf("%14.6g",useAdj) )
        m = m + 30
      end do
    else      ; else assume months with varying number of days (for "noleap" calendar)
      timeMonVar = BmidMnDays-1
      timeDayVar = ispan(0,426,1)
      oldDatDay = linint1_n (timeMonVar, useOrig, False, timeDayVar, 0, 0)
      newDatDay  = linint1_n (timeMonVar, useAdj,  False, timeDayVar, 0, 0)
      print("")
      print(sprintf("%16.6g",oldDatDay) + sprintf("%16.6g",newDatDay) )    ; printing daily values from: (a) original and (b) adjusted seasonal cycle
      print("")
      print("RESULTS OF CHECK (middle two columns should match closely):")
      print("   Mon" + " Calc Orig Avg" + "      Orig Avg" + "  Calc Adj Avg" + "       Adj Avg")
      do i = 1,nT-2
        print( sprinti("%6i",i) + \
               sprintf("%14.6g",dim_avg(oldDatDay(toint(BbegMnDays(i)-1):toint(BendMnDays(i)-1)))) + \
               sprintf("%14.6g",useOrig(i)) + \
               sprintf("%14.6g",dim_avg(newDatDay(toint(BbegMnDays(i)-1):toint(BendMnDays(i)-1)))) + \
               sprintf("%14.6g",useAdj(i)) )
      end do
    end if        ; end if months30days
  end if          ; end if checkAdjustment = True
  ; TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
  ; --------------------------------------------------------------------------------------
  
  ; Sample printouts
  print("")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("Sample printout: Annual cycle at the following point:")
  print("lat=" + ilat + " (" + xlats(ilat) + ")")
  print("lon=" + ilon + " (" + xlons(ilon) + ")")
  months = ispan(1,12,1)
  print(" Month" + "        xMonAdj")
  print( sprinti("%6i",months) + sprintf("%16.5g",xMonAdj(:,ilat,ilon)) )
  
  print("")
  print("Sample printout: Longitudinal profile at the following point:")
  print("month=" + imon)
  print("lat=" + ilat + " (" + xlats(ilat) + ")")
  print("         lon" + "        xMonAdj")
  print( sprintf("%12.4f",xlons(:)) + sprintf("%16.5g",xMonAdj(imn-1,ilat,:)) )
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")  
  
  ; Un-generalize
  sstMonAdj = xMonAdj
  
  printVarSummary(sstMon)
  printMinMax(sstMon, True)
  print("")
  printVarSummary(sstMonAdj)
  printMinMax(sstMonAdj, True)
  print("")
  
  do i = 0,11
    d = sstMonAdj(i,:,:) - sstMon(i,:,:)
    print("SST difference min/max for month " + (i+1) + " (BEFORE consistency checks): " + sprintf(fmt,min(d)) + ", " + sprintf(fmt,max(d)))
    delete(d)
  end do
  ; -----------------------------------------------------------------------------------
  ; Write out SST differences annual cycle monthly mean SST for checking/plotting
  filo = workPrefix + ".raw_adj_annCyc_monthly_SST.nc"
  print("Writing data to: " + filo)
  system("/bin/rm -f " + filo)               ; remove old file if exists
  fout = addfile(filo,"c")                   ; "pointer" to output file is fout

  setfileoption(fout,"DefineMode",True)
  fAtt                = True            ; assign file attributes
  fAtt@title          = "Raw, adjusted, and difference of averaged annual cycle of SST"
  fAtt@note           = "Derived from 500 year climo of " + caseName
  fAtt@author         = "Jim Benedict (jjb278@gmail.com)"
  fAtt@from_script    = get_script_name()
  fAtt@created        = systemfunc("date")
  fileattdef( fout, fAtt )            ; copy file attributes
  setfileoption(fout,"DefineMode",False)
    
  d = sstMonAdj - sstMon
  fout->sstMon  = sstMon
  fout->sstMonAdj  = sstMonAdj

  d!0 = "time"
  d!1 = "lat"
  d!2 = "lon"
  d&time = sstMon&time
  d&lat = sstMon&lat
  d&lon = sstMon&lon

  fout->diff = d
  delete(d)

  ; Final metadata refinements in section below
  
  ; Clean up
  delete([/xMon,xlats,xlons,nlatUse,nlonUse,xMonBufAdj,a,aAdj,xMonConf/])
  
  
  ; @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ; @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ; @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  ; @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  print("")
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
  print("")
  print("Diddling Ice Area Concentration Data")

  ; -----------------------------------------------------------------------------------
  ; Adjust correction terms to preserve true monthly means (ICE "diddling")
  ;   See, e.g., Sheng & Zwiers (1998, Climate Dynamics) or 
  ;   http://www-pcmdi.llnl.gov/publications/pdf/60.pdf from the AMIP-II website
  ;   Because we are applying the *climatological* seasonal cycle of ocean temperature
  ;     and salinity correction terms to our model, we need only adjust the pre-calculated
  ;     climatological values rather than adjust each individual monthly mean.  If the user
  ;     wants to do the monthly mean adjustments to each individual month, that procedure
  ;     would require that all monthly values are stored in one array such that code
  ;     changes would be needed in the "nf" loop above.
  
  ; Create "buffered" climatological annual cycle by adding the December climatological
  ;   mean value before the January value, and adding the January value after the
  ;   final December value such that you have a total of 14 months. These added
  ;   months are needed for the linear solver.
  F_ICE = addfile(inICE_File,"r")

  iceMon = F_ICE->ice_cov
  iceMon&time = ispan(1,12,1)
  xMon     = iceMon    ; Use settings, to make it generic

  xlats=F_ICE->lat
  xlons=F_ICE->lon
  nlonICE = dimsizes(iceMon(0,0,:))
  nlatICE = dimsizes(iceMon(0,:,0))
  nlatUse  = nlatICE
  nlonUse  = nlonICE
  
  xMonBuf    = new((/nT,nlatUse,nlonUse/),typeof(xMon),xMon@_FillValue)
  xMonBufAdj = xMonBuf  ; these **avgBufAdj arrays are really only required for testing, could be deleted
  xMonBuf(0,:,:)    = xMon(11,:,:)
  xMonBuf(1:12,:,:) = xMon(0:11,:,:)
  xMonBuf(13,:,:)   = xMon(0,:,:)
  
  ; Apply pre-calculated weighting coefficients and solve linear equation
  xMonAdj = xMon
  xMonAdj = xMon@_FillValue
  print("Performing adjustment to monthly mean...")
  do j = 0,nlatUse-1
    do i = 0,nlonUse-1
      a    = xMonBuf(:,j,i)
      aAdj = solve_linsys(lineqn_coef, a)   ; aAdj will be one-dimensional [aAdj(14)]
      xMonBufAdj(:,j,i) = aAdj              ; only required for later testing, could be deleted
      xMonAdj(0:11,j,i) = aAdj(1:12)        ; Retain only Jan-Dec values, discard bookend buffers    
    end do
  end do
  
  ; Minor corrections -- (re)mask out land points
  printMinMax(xMonAdj, True)
  xMonConf = conform(xMonAdj,xMon(0,:,:),(/1,2/))
  xMonAdj  = where(ismissing(xMonConf),xMonAdj@_FillValue,xMonAdj)
  printMinMax(xMonAdj, True)
  
  if(test_print) then
   print("")
   print("After adjustment to maintain actual climatological monthly mean:")
   print( "     LAT" + "     LON" + "         Raw X" + "   Corrected X")
   print( sprintf("%8.2f",xlats(ilat)) + sprintf("%8.2f",xlons(ilon)) + \
          sprintf("%16.6g",xMon(imn,ilat,ilon)) + \
          sprintf("%16.6g",xMonAdj(imn,ilat,ilon)) )
  end if
  
  ; --------------------------------------------------------------------------------------
  ; TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
  if(checkAdjustment) then
    print("")
    print("Checking that adjustment procedure worked with actual (buffered) data...")
    useOrig = xMonBuf(:,ilat,ilon)
    useAdj  = xMonBufAdj(:,ilat,ilon)
    
    fmt = "%14.6g"
    print(sprintf(fmt,useOrig) + sprintf(fmt,useAdj))
    
    
    ; Check that monthly-to-daily interpolated values represent actual monthly mean values
    if(months30days) then
      timeMon30 = ispan(14,404,30)
      timeDay30 = ispan(0,419,1)
      oldDatDay = linint1_n (timeMon30, useOrig, False, timeDay30, 0, 0)
      newDatDay  = linint1_n (timeMon30, useAdj,  False, timeDay30, 0, 0)
      print("")
      print(sprintf("%16.6g",oldDatDay) + sprintf("%16.6g",newDatDay) )
      print("")
      m = 30
      do i = 1,nT-2
        print( sprinti("%6i",i) + \
               sprintf("%14.6g",dim_avg(oldDatDay(m:m+29))) + sprintf("%14.6g",useOrig) + \
               sprintf("%14.6g",dim_avg(newDatDay(m:m+29))) + sprintf("%14.6g",useAdj) )
        m = m + 30
      end do
    else      ; else assume months with varying number of days (for "noleap" calendar)
      timeMonVar = BmidMnDays-1
      timeDayVar = ispan(0,426,1)
      oldDatDay = linint1_n (timeMonVar, useOrig, False, timeDayVar, 0, 0)
      newDatDay  = linint1_n (timeMonVar, useAdj,  False, timeDayVar, 0, 0)
      print("")
      print(sprintf("%16.6g",oldDatDay) + sprintf("%16.6g",newDatDay) )    ; printing daily values from: (a) original and (b) adjusted seasonal cycle
      print("")
      print("RESULTS OF CHECK (middle two columns should match closely):")
      print("   Mon" + " Calc Orig Avg" + "      Orig Avg" + "  Calc Adj Avg" + "       Adj Avg")
      do i = 1,nT-2
        print( sprinti("%6i",i) + \
               sprintf("%14.6g",dim_avg(oldDatDay(toint(BbegMnDays(i)-1):toint(BendMnDays(i)-1)))) + \
               sprintf("%14.6g",useOrig(i)) + \
               sprintf("%14.6g",dim_avg(newDatDay(toint(BbegMnDays(i)-1):toint(BendMnDays(i)-1)))) + \
               sprintf("%14.6g",useAdj(i)) )
      end do
    end if        ; end if months30days
  end if          ; end if checkAdjustment = True
  ; TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST TEST
  ; --------------------------------------------------------------------------------------
  
  ; Sample printouts
  print("")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("Sample printout: Annual cycle at the following point:")
  print("lat=" + ilat + " (" + xlats(ilat) + ")")
  print("lon=" + ilon + " (" + xlons(ilon) + ")")
  months = ispan(1,12,1)
  print(" Month" + "        xMonAdj")
  print( sprinti("%6i",months) + sprintf("%16.5g",xMonAdj(:,ilat,ilon)) )
  
  print("")
  print("Sample printout: Longitudinal profile at the following point:")
  print("month=" + imon)
  print("lat=" + ilat + " (" + xlats(ilat) + ")")
  print("         lon" + "        xMonAdj")
  print( sprintf("%12.4f",xlons(:)) + sprintf("%16.5g",xMonAdj(imn-1,ilat,:)) )
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")  
  
  ; Un-generalize
  iceMonAdj = xMonAdj
  
  printVarSummary(iceMon)
  printMinMax(iceMon, True)
  print("")
  printVarSummary(iceMonAdj)
  printMinMax(iceMonAdj, True)
  print("")
  
  do i = 0,11
    d = iceMonAdj(i,:,:) - iceMon(i,:,:)
    print("SEAICE difference min/max for month " + (i+1) + " (BEFORE consistency checks): " + sprintf(fmt,min(d)) + ", " + sprintf(fmt,max(d)))
    delete(d)
  end do
  ; -----------------------------------------------------------------------------------
  ; Write out ICE differences annual cycle monthly mean ICE for checking/plotting
  filo = workPrefix + ".raw_adj_annCyc_monthly_ICE.nc"
  print("Writing data to: " + filo)
  system("/bin/rm -f " + filo)               ; remove old file if exists
  fout = addfile(filo,"c")                   ; "pointer" to output file is fout

  setfileoption(fout,"DefineMode",True)
  fAtt                = True            ; assign file attributes
  fAtt@title          = "Raw, adjusted, and difference of averaged annual cycle of ICE"
  fAtt@note           = "Derived from 500 years of " + caseName
  fAtt@author         = "Jim Benedict (jjb278@gmail.com)"
  fAtt@from_script    = get_script_name()
  fAtt@created        = systemfunc("date")
  fileattdef( fout, fAtt )            ; copy file attributes
  setfileoption(fout,"DefineMode",False)
    
  d = iceMonAdj - iceMon
  fout->iceMon  = iceMon
  fout->iceMonAdj  = iceMonAdj

  printVarSummary(iceMon)
  printVarSummary(d)

  d!0 = "time"
  d!1 = "lat"
  d!2 = "lon"
  d&time = iceMon&time
  d&lat = iceMon&lat
  d&lon = iceMon&lon

  fout->diff = d
  delete(d)
  
  ; Final metadata refinements in section below
  
  ; Clean up
  delete([/xMon,xlats,xlons,nlatUse,nlonUse,xMonBufAdj,a,aAdj,xMonConf/])
  
  
  ; ======================================================================================
  ; ======================================================================================
  ; SST and ICE consistency checks
  ;
  ; There are several instances where SST and ICE fall outside of physical/realistic bounds.
  ;   This is especially true for ICE where there is abrupt full melting of sea ice from one
  ;   month to the next; in this case, it is almost impossible to realistically capture
  ;   this behavior accurately given that CESM2 uses monthly data linearly interpolated to
  ;   daily data.  The first zero-ice month would tend to be adjusted to a sizable negative
  ;   value in an attempt to maintain a zero monthly mean, but obviously negative ICE values
  ;   are not physical.  Our goal of matching SST and ICE between the B compset and F
  ;   compset ("adjustments" to the monthly means) must be balanced by the requirement that
  ;   SST and ICE remain within reasonable physical bounds.  We impose these bounds below:
  

  SST = sstMon
  SSTadj = sstMonAdj
  SEAICE = iceMon
  SEAICEadj = iceMonAdj

  print("")
  print("Consistency processing...")
  ; Corrections for data consistency
  ; 1) If SST < -1.8 or ice frac >= 0.9, SST = -1.8
  SST = where(SST .le. -1.8, -1.8, SST)
  SST = where(SEAICE .gt. 0.9, -1.8, SST)
  SSTadj = where(SSTadj .le. -1.8, -1.8, SSTadj)
  SSTadj = where(SEAICEadj .gt. 0.9, -1.8, SSTadj)

  ; 2) min ice frac is 0%, max ice_frac is 100%
  SEAICE = where(SEAICE .lt. 0., 0., SEAICE)
  SEAICE = where(SEAICE .gt. 1., 1., SEAICE)
  SEAICEadj = where(SEAICEadj .lt. 0., 0., SEAICEadj)
  SEAICEadj = where(SEAICEadj .gt. 1., 1., SEAICEadj)

  ; 3) Hurrell et al. (2008, J. Climate) correction to SST (if 0.15 < ice frac < 0.9)
  SST_max = 9.328*(0.729-(SEAICE/100.)^3)-1.8
  SST = where((SEAICE .gt. 0.15) .and. (SEAICE .lt. 0.9) .and. (SST .gt. SST_max), SST_max, SST)
  delete(SST_max)
  SST_max = 9.328*(0.729-(SEAICEadj/100.)^3)-1.8
  SSTadj = where((SEAICEadj .gt. 0.15) .and. (SEAICEadj .lt. 0.9) .and. (SSTadj .gt. SST_max), SST_max, SSTadj)

  ; 4) No sea ice if SST > 4.97
  SEAICE = where(SST .gt. 4.97, 0., SEAICE)
  SEAICEadj = where(SSTadj .gt. 4.97, 0., SEAICEadj)

  ; 5) Hurrell et al. (2008, J. Climate) correction to ice frac (if -1.8 < SST < 4.97)
  ice_max = (0.729-(SST+1.8)/9.328)^(1./3.)
  SEAICE = where( (SST .gt. -1.8) .and. (SST .lt. 4.97) .and. (SEAICE .gt. ice_max), ice_max, SEAICE )
  delete(ice_max)
  ice_max = (0.729-(SSTadj+1.8)/9.328)^(1./3.)
  SEAICEadj = where( (SSTadj .gt. -1.8) .and. (SSTadj .lt. 4.97) .and. (SEAICEadj .gt. ice_max), ice_max, SEAICEadj )
    

  ; skip filling for missing values
  print((/"Filling in missing values in SST and SEAICE..."/))
  poisson_grid_fill(SST, True, 1, 500, 0.01, 0.6, 0)
  poisson_grid_fill(SSTadj, True, 1, 500, 0.01, 0.6, 0)
  poisson_grid_fill(SEAICE, True, 1, 500, 0.01, 0.6, 0)
  poisson_grid_fill(SEAICEadj, True, 1, 500, 0.01, 0.6, 0)
  
  ; ======================================================================================
  ; ======================================================================================
  ; Final metadata refinements and checks
  if(isvar("time")) then
    delete(time)
  end if
  time = todouble(midMnDays)
  time@information = "middle of month"
  time@axis = "T"
  time@calendar = "gregorian"
  time@long_name = "time"
  time@standard_name = "time"
  time@units = "days since 0000-01-01 00:00:00"
  date!0 = "time"
  date&time = time
  date@units = "current date (yyyymmdd)"
  datesec!0 = "time"
  datesec&time = time
  datesec@units = "current seconds of current date"
  
  delete(SST&time)
  delete(SSTadj&time)

  SST&time    = time
  printVarSummary(SST)
  SST@units   = "degC"
  copy_VarCoords(SST, SSTadj)
  SSTadj@units  = SST@units
  SSTadj@method = "Adjusted using linear solver in attempt to maintain similar monthly means even after linear interpolation within model"
  
  delete(SEAICE&time)    ; had to do this b/c time dim of original sea ice used type float instead of double
  SEAICE&time  = time
  SEAICE@units = "Fraction"
  SEAICE@long_name = "sea ice concentration"
  copy_VarCoords(SEAICE, SEAICEadj)
  SEAICEadj@method = SSTadj@method
  SEAICEadj@units  = SEAICE@units
  SEAICEadj@long_name = SEAICE@long_name
  
  printVarSummary(time)
  printVarSummary(date)
  printVarSummary(datesec)
  
  printVarSummary(SST)
  printVarSummary(SSTadj)
  printVarSummary(SEAICE)
  printVarSummary(SEAICEadj)    
  
  print("")
  print("===============================================================================")
  print("===============================================================================")
  print("ARRAY SUMMARY:")
  print("")
  printVarSummary(SST)
  printMinMax(SST, True)
  print("Number of missing points: " + num(ismissing(SST)))
  print("---------------------------------------------------------------------------")
  printVarSummary(SSTadj)
  printMinMax(SSTadj, True)
  print("Number of missing points: " + num(ismissing(SSTadj)))
  print("---------------------------------------------------------------------------")
  printVarSummary(SEAICE)
  printMinMax(SEAICE, True)
  print("Number of missing points: " + num(ismissing(SEAICE)))
  print("---------------------------------------------------------------------------")
  printVarSummary(SEAICEadj)
  printMinMax(SEAICEadj, True)
  print("Number of missing points: " + num(ismissing(SEAICEadj)))
  print("")
  print("")
  print("CHECK ON DIFFERENCE BETWEEN PRE-DIDDLED AND POST-DIDDLED SST:")
  fmt = "%12.6f"
  do i = 0,11
    d = SSTadj(i,:,:) - SST(i,:,:)
    print("SST difference min/max for month " + (i+1) + ": " + sprintf(fmt,min(d)) + ", " + sprintf(fmt,max(d)))
    delete(d)
  end do
  print("")
  print("")
  print("CHECK ON DIFFERENCE BETWEEN PRE-DIDDLED AND POST-DIDDLED SEAICE:")
  fmt = "%12.6f"
  do i = 0,11
    d = SEAICEadj(i,:,:) - SEAICE(i,:,:)
    print("SEAICE difference min/max for month " + (i+1) + " (AFTER consistency checks): " + sprintf(fmt,min(d)) + ", " + sprintf(fmt,max(d)))
    delete(d)
  end do
  print("")
  
  
  ; ======================================================================================
  ; ======================================================================================
  ; Write data to file
  if(write_NC) then
    filo = getenv("OUTPUT_SSTICE_FILE") 
    print("Writing data to: " + filo)
    system("/bin/rm -f " + filo)               ; remove old file if exists
    fout = addfile(filo,"c")                   ; "pointer" to output file is fout

    setfileoption(fout,"DefineMode",True)
    fAtt                = True            ; assign file attributes
    fAtt@title          = "Annual cycle of monthly mean sea surface temperature and sea ice concentrations from single year"
    fAtt@note           = "Derived from years 0001-0500 of " + caseName
    fAtt@author         = "Wuyin Lin (wlin@bnl.gov), Jim Benedict (jjb278@gmail.com)"
    fAtt@from_script    = get_script_name()
    fAtt@created        = systemfunc("date")
    fAtt@history_sst    = F_SST@history
    fAtt@history_seaice = F_ICE@history
    fileattdef( fout, fAtt )            ; copy file attributes
    setfileoption(fout,"DefineMode",False)
      
    fout->date  = date
    fout->datesec  = datesec
    fout->ice_cov = SEAICEadj
    fout->ice_cov_prediddle = SEAICE
    fout->SST_cpl = SSTadj
    fout->SST_cpl_prediddle = SST
  end if          ; end if write_NC = True
  
  print("")
  print("NOTE TO USER:  WARNINGS FROM NCL, SUCH AS:")
  print("warning:Right hand side has no coordinate variable can not delete coordinate variable of a file, use (/ .. /) to avoid this message")
  print("ARE TO BE EXPECTED AND DO NOT NEGATIVELY AFFECT THE PROGRAM OUTPUT.")
  
end
