;;This is a function derived from DART observation-space diagnostics routines
; 'obs_diag' produces a netcdf file containing the diagnostics.
; 'obs_diag' condenses the obs_seq.final information into summaries for a few
; specified regions - on a level-by-level basis.
; The number of observations possible reflects only those observations
; that have incoming QC values of interest. Any observation with a DART
; QC of 5 or 6 is not considered 'possible' for the purpose of this graphic.
; NOTE: if the observation was designated as a TRUSTED observation in the
;       obs_diag program, the observations that were rejected by the outlier
;       threshhold STILL PARTICIPATE in the calculation of the rmse, spread, etc.
;       The _values_ plotted by plot_profile reflect that. The number of observations
;       "used" becomes unclear. The number of observations used (designated by the
;       asterisk) is ALWAYS the number of observations successfully assimilated.
;       For TRUSTED observations, this is different than the number used to calculate
;       bias, rmse, spread, etc.
; fname    :  netcdf file produced by 'obs_diag'
;
; CopyString : string defining the metric of interest. 'rmse', 'spread', etc.
;              Possible values are available in the netcdf 'CopyMetaData' variable.
;              (ncdump -v CopyMetaData obs_diag_output.nc)
; ObsTypeString  : The strings of each observation type to plot.
;                  Each observation type will be plotted in a separate graphic.
;                  Default is to plot all available observation types.
; region   : the strings of region for the statistics, statistics for five regions are available by default
;            global:[0-360,90S-90N]
;            Northern Hemisphere:[0-360,20N-90N]
;            Tropics:[0-360,20S-20N]     
;            Southern Hemisphere:[0-360,90S-20S]
;            North America:[235-295,25N-55N] 
;                
; type     : there are three types of statistics in the file, including 
;            1. guess:    the time evolution  of the rmse, biases etc.; 
;            2. RankHist: the Rank histograms (sometimes called verification rank histograms
;                         or Talagrand diagrams), which is a way to show how reliable an ensemble 
;                         forecast is compared to a set of newly observed data.
;            3. VPguess:  the mean vertical profiles of the rmse or biases etc. 
;
; OUTPUT: 'plotdat' is a structure containing what was last plotted.
;         A postscript file containing a page for each level - each region.
;         The other file is a simple text file containing summary information
;         about how many observations were assimilated, how many were available, etc.
;         Both of these filenames contain the observation type,
;         copy and region as part of the name.
;plotdat        = read_dart_diag_netcdf(fname, CopyString, ObsTypeString, region, level, type);
; EXAMPLE - plot the evolution of the rmse for just the radiosonde temperature obs
;           This requires that the 'RADIOSONDE_TEMPERATURE' is one of the known
;           observation types in the netCDF file.
;
; fname          = 'obs_diag_output.nc';
; CopyStrings    = 'rmse';
; ObsTypeString  = 'RADIOSONDE_TEMPERATURE'
; region         = 'global' or ''
; type           = 'guess' ; or 'VPguess' or 'RankHist'  
;
; plotdat        = read_dart_diag_netcdf('obs_diag_output.nc', 'rmse', 'RADIOSONDE_TEMPERATURE', ...
;                                        'global', -999, 'guess');
;;

function read_dart_diag_netcdf(fname, CopyString, ObsTypeString, region, type) 
begin

   f                = addfile(fname,"r")
   time             = f->time

   RegionStrings    = chartostring(f->region_names)
   CopyStrings      = chartostring(f->CopyMetaData)     ; (copy,string)
   ObsTypeStrings   = chartostring(f->ObservationTypes) ; (ObsTypes,string)
   copy_index       = f->copy
   obstypes         = f->obstypes
   region_index     = f->region
   
   ;;locate the statistic metric;;;;;
   icopy = ind(str_lower(str_squeeze(CopyStrings)).eq.str_lower(str_squeeze(CopyString)))
   if(ismissing(icopy))then 
     print("") 
     print("The required observation data "+ CopyString + " does not exit, please check your data!!")
     print("")
     exit
   end if 

   ;;locate the region for plot;;;;;
   if (typeof(region).eq."integer")then 
     iregion = ind(region_index.eq.region)
     regname = RegionStrings(iregion)
   end if 

   if (typeof(region).eq."string")then
     iregion = ind(str_lower(str_squeeze(RegionStrings)).eq.str_lower(region))
     regname = RegionStrings(iregion)
   end if 
    
   if(ismissing(iregion))then
     print("")
     print("The required region data "+ region + " does not exit, please check your data!!")
     print("")
     exit
   end if
   
   ;;extract the data information;;;;
   if(str_squeeze(type).eq."guess")then
     varnam   = str_squeeze(ObsTypeString)+"_guess"
     dart_obs = f->$varnam$(:,icopy,:,iregion)
   else if (str_squeeze(type).eq."VPguess") then 
     varnam   = str_squeeze(ObsTypeString)+"_VPguess"
     dart_obs = f->$varnam$(icopy,:,iregion)
   else if (str_squeeze(type).eq."RankHist") then 
     varnam   = str_squeeze(ObsTypeString)+"_guess_RankHist"
     dart_obs = f->$varnam$(:,icopy,:,iregion)
   else
    varnam = str_squeeze(ObsTypeString)+"_"+str_squeeze(type)
    print("")
    print("The required variable "+ varnam + " does not exit, please check your data!!")
    print("")
    exit
   end if 
   end if 
   end if 

   dart_obs@time        = time 
   dart_obs@time_unit   = time@units 

   if(isStrSubset(varnam,"GPSRO"))then 
      dart_obs@lev       = f->hlevel   
      dart_obs@levp      = f->hlevel_edges
      dart_obs@lev_unit  = f->hlevel@units
   else if (isStrSubset(varnam,"MARINE_SFC").or.isStrSubset(varnam,"LAND_SFC").or.isStrSubset(varnam,"RADIOSONDE_SURFACE"))then 
      dart_obs@lev       = "Surface"
      dart_obs@levp      = "Surface"
      dart_obs@lev_unit  = ""
   else
      dart_obs@lev       = f->plevel
      dart_obs@levp      = f->plevel_edges
      dart_obs@lev_unit  = f->plevel@units  
   end if 
   end if 

   dart_obs@obsname = ObsTypeString
   dart_obs@regname = regname
   dart_obs@copyname= CopyStrings(icopy)
   RegionStrings    = chartostring(f->region_names)
   CopyStrings      = chartostring(f->CopyMetaData)     ; (copy,string)
   ObsTypeStrings   = chartostring(f->ObservationTypes) ; (ObsTypes,string)
   copy_index       = f->copy
   obstypes         = f->obstypes
   region_index     = f->region

   delete([/icopy,iregion,regname,varnam,RegionStrings,CopyStrings,ObsTypeStrings,copy_index,obstypes,region_index/])

   return(dart_obs)

end

