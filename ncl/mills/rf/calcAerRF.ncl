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
  
  MCA_V=a->MBCPHO_V
  MCA_V@long_name="Mass of all carbon aerosols in and below layer"  
  MCA_V=MCA_V+a->MOCPHI_V
  MCA_V=MCA_V+a->MBCPHI_V 
  MCA_V=MCA_V+a->MBCPHO_V 
  MBCPHOFRAC_V=a->MBCPHO_V 
  MBCPHOFRAC_V@long_name=\
    "Fractional mass of BCPHO to all carbon aerosols in and below layer"
  MBCPHOFRAC_V=MBCPHOFRAC_V/MCA_V
  a->MCA_V=MCA_V
  a->MBCPHOFRAC_V=MBCPHOFRAC_V

  


end


basedir="/data5/mmills/wa319/carma-soot/"
year="2034"

caseid="bcrckt03"
rootdir=basedir+caseid+"/"+"h0zm.annual/"
file1=caseid+".cam2.h0zm."+year+".annual.nc"
;calcAerRF(rootdir,file1)

caseid="bccntl02"
rootdir=basedir+caseid+"/"+"h0zm.annual/"
file1=caseid+".cam2.h0zm."+year+".annual.nc"
calcAerRF(rootdir,file1)
