load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

undef("calcAerRF")
procedure calcAerRF(rootdir:string, file1:string) 
begin 
  print(file1)
  a = addfile(rootdir+file1,"w")
  QRS_RF=a->QRS_RF
  QRS_AER=a->QRS
  QRS_AER=QRS_AER-QRS_RF
  a->QRS_AER=QRS_AER
  
  FSNT_RF=a->FSNT_RF
  FSNT_AER=a->FSNT
  FSNT_AER=FSNT_AER-FSNT_RF
  a->FSNT_AER=FSNT_AER
  
  FSNS_RF=a->FSNS_RF
  FSNS_AER=a->FSNS
  FSNS_AER=FSNS_AER-FSNS_RF
  a->FSNS_AER=FSNS_AER
  
  FSNTC_RF=a->FSNTC_RF
  FSNTC_AER=a->FSNTC
  FSNTC_AER=FSNTC_AER-FSNTC_RF
  a->FSNTC_AER=FSNTC_AER
  
  FSNSC_RF=a->FSNSC_RF
  FSNSC_AER=a->FSNSC
  FSNSC_AER=FSNSC_AER-FSNSC_RF
  a->FSNSC_AER=FSNSC_AER

end


rootdir="/Volumes/Data/Models/waccm/run/carma-soot/bch0zm/"
file1="bccntl01.cam2.h0zm.2025-01-01-00000.nc"
;calcAerRF(rootdir,file1)

file1="bcrckt000001.cam2.h0zm.2025-01-01-00000.nc"
;calcAerRF(rootdir,file1)

file1="bcrckt01.cam2.h0zm.2025-01-01-00000.nc"
;calcAerRF(rootdir,file1)

file1="bcrckt6000.cam2.h0zm.2025-01-01-00000.nc"
;calcAerRF(rootdir,file1)

file1="bcrckt6000b.cam2.h0zm.2025-01-01-00000.nc"
calcAerRF(rootdir,file1)
