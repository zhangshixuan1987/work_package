 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
   fil = addfile("./DMSflux.1850-2100.1deg_latlon_conserv.POPmonthlyClimFromACES4BGC_c20160727.nc","w")
   DMS = fil->DMS
   DMS(12:23,:,:) = DMS(0:11,:,:) 
   DMS(24:35,:,:) = DMS(0:11,:,:)
   fil->DMS=(/DMS/)
end
    

