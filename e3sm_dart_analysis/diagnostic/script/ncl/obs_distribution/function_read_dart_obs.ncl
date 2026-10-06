;; read_obs_netcdf reads in the netcdf flavor observation sequence file
;                  and returns a subsetted structure.
;
; fname         = 'obs_epoch_001.nc';
; ObsTypeString = 'RADIOSONDE_U_WIND_COMPONENT';   % or 'ALL' ...
; region        = [0 360 -90 90 -Inf Inf];
; CopyString    = 'NCEP BUFR observation';         % or 'ALL' ...
; QCString      = 'DART quality control';
; verbose       = 1;   % anything > 0 == 'true'
; obs = read_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose);
;
; The return variable 'obs' is a structure. As an example ...
;
;           fname: 'obs_epoch_001.nc'
;   ObsTypeString: 'RADIOSONDE_U_WIND_COMPONENT'
;          region: [0 360 -90 90 -Inf Inf]
;      CopyString: 'NCEP BUFR observation'
;        QCString: 'DART quality control'
;         verbose: 1
;      timestring: [2x20 char]
;            lons: [2343x1 double]
;            lats: [2343x1 double]
;               z: [2343x1 double]
;             obs: [2343x1 double]
;            Ztyp: [2343x1 int32]
;            keys: [2343x1 int32]
;            time: [2343x1 double]
;              qc: [2343x1 int32]

;Satellite    = (/"GPSRO", "SAT", "AIRS"/)
;NCEP_PREBUFR = (/"AIRCRAFT", "ACARS","RADIOSONDE","MARINE"/)

function read_dart_obs_netcdf(fname, ObsTypeString, region, CopyString, QCString, verbose) 
begin

   f                = addfile(fname,"r")
   time             = f->time

   QCStrings        = chartostring(f->QCMetaData)
   CopyStrings      = chartostring(f->CopyMetaData)     ; (copy,string)
   ObsTypeStrings   = chartostring(f->ObsTypesMetaData) ; (ObsTypes,string)
   ObsTypes         = f->ObsTypes         ; (ObsTypes)
   ObsIndex         = f->ObsIndex         ; (ObsIndex) i.e. time
   obs_type         = f->obs_type
   obs_keys         = f->obs_keys

   if(isfilevar(f,"which_vert"))then 
     z_type         = f->which_vert  
   else 
     z_type         = new(dimsizes(obs_keys),integer)
     z_type         = 1
   end if              
   
   z_pos            = tostring(z_type)
   z_pos            = where(z_type.eq.-2,"up",z_pos)
   z_pos            = where(z_type.eq.-1,"up",z_pos)
   z_pos            = where(z_type.eq.1,"down",z_pos)
   z_pos            = where(z_type.eq.2,"down",z_pos)
   z_pos            = where(z_type.eq.3,"up",z_pos)
   z_unt            = tostring(z_type)
   z_unt            = where(z_type.eq.-2,"unknown",z_unt)
   z_unt            = where(z_type.eq.-1,"surface",z_unt)
   z_unt            = where(z_type.eq.1,"model level",z_unt)
   z_unt            = where(z_type.eq.2,"pressure",z_unt)
   z_unt            = where(z_type.eq.3,"height",z_unt)
  
   loc              = f->location
   obs              = f->observations
   qc               = f->qc

   if ( verbose ) then
     my_types  = get_unique_values(obs_type)
     do i = 0,dimsizes(my_types)-1,1
        obtype = my_types(i) 
        print(obtype)
        inds   = ind(obs_type.eq.obtype) 
        zall   = loc(inds,2) 
        myz    = get_unique_values(zall)    
        print("")
        print("N = " +num(.not.ismissing(myz)) +" between levels " + min(myz) + " and " + max(myz)  + \
              " for " + ObsTypeStrings(obtype))
        print("")
        delete([/obtype,inds,zall,myz/])
     end do 
     delete(my_types)
   end if 

   ;;extract the location information;;;;
   if(str_squeeze(CopyString).ne." ")then 
     icopy = ind(str_squeeze(CopyStrings).eq.str_squeeze(CopyString))
     obs1  = obs(:,icopy) 
   else
     print("") 
     print("The required observation data "+ CopyString + " does not exit, please check your data!!")
     print("")
     exit 
   end if 

   ;;extract quality control information;;;;;
   if(str_squeeze(QCString).ne." ")then 
     iqc = ind(str_squeeze(QCStrings).eq.str_squeeze(QCString))
     qc1 = qc(:,iqc)
   else
     qc1 = qc(:,0)
   end if 
   delete([/qc,obs/]) 

   ;;extract data in the specified region;;
   regs   = tofloat(str_split(region, " "))
   indreg = ind ( (loc(:,0).ge.regs(0).and.loc(:,0).le.regs(1)) .and. \
                  (loc(:,1).ge.regs(2).and.loc(:,1).le.regs(3)) .and. \
                  (loc(:,2).ge.regs(4).and.loc(:,2).le.regs(5)) )
   if(.not.all(ismissing(indreg)))then 
     loc2  = loc(indreg,:)
     qc2   = qc1(indreg)  
     obs2  = obs1(indreg)
     oind2 = ObsIndex(indreg)
     ztyp2 = z_type(indreg)
     zpos2 = z_pos(indreg)
     zunt2 = z_unt(indreg)
     time2 = time(indreg)
     otyp2 = obs_type(indreg)
     okey2 = obs_keys(indreg)
   else
     loc2  = loc
     qc2   = qc1
     obs2  = obs1
     oind2 = ObsIndex
     ztyp2 = z_type
     zpos2 = z_pos
     zunt2 = z_unt
     time2 = time
     otyp2 = obs_type
     okey2 = obs_keys
   end if 
   delete([/loc,qc1,obs1,indreg,z_type,z_pos,z_unt,time/])
  
   ;;extract certain type of data;;;;
   if(str_lower(str_squeeze(ObsTypeString)).ne."all")then
     iobtyp = ind(str_squeeze(ObsTypeStrings).eq.str_squeeze(ObsTypeString))
     print("")
     print("only process this type of data: "+ ObsTypeStrings(iobtyp))
     print("")
     indtyp = ind(otyp2.eq.ObsTypes(iobtyp))
     ostr1  = ObsTypes(iobtyp)
     ostr2  = ObsTypeStrings(iobtyp)
   else
     indtyp = ind(.not.ismissing(otyp2))
     myotyp = get_unique_values(otyp2)
     ind1d  = get1Dindex(ObsTypes,myotyp)
     ostr1  = ObsTypes(ind1d)
     ostr2  = ObsTypeStrings(ind1d)
     delete([/myotyp,ind1d/])
   end if
 
   dart_obs        = obs2(indtyp)
   dart_obs@numobs = num(.not.ismissing(indtyp))
   dart_obs@lon    = loc2(indtyp,0)
   dart_obs@lat    = loc2(indtyp,1)
   dart_obs@lev    = loc2(indtyp,2)
   dart_obs@ztyp   = ztyp2(indtyp)
   dart_obs@zpos   = zpos2(indtyp) 
   dart_obs@zunt   = zunt2(indtyp)
   dart_obs@time   = time2(indtyp) 
   dart_obs@qc     = qc2(indtyp)
   dart_obs@oname1 = ostr1
   dart_obs@oname2 = ostr2
   dart_obs@otype  = otyp2(indtyp) 
   dart_obs@okey   = okey2(indtyp)
   dart_obs@oindex = oind2(indtyp)
   delete([/obs2,loc2,qc2,indtyp,ztyp2,zpos2,zunt2,time2,otyp2,okey2,ostr1,ostr2,oind2/])
   delete([/ObsIndex,QCStrings,CopyStrings,ObsTypeStrings,ObsTypes,obs_type,obs_keys/])
   return(dart_obs)

end

