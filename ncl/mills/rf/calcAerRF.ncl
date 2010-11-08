load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "../addPdelZm.ncl"

undef("calcAerRF")
procedure calcAerRF(rootdir:string, file1:string) 
begin 

  ; Some constants we will need:
  g       = 9.806               ; m s^-2 

  print(rootdir+file1)
  
  a = addfile(rootdir+file1,"w")
  
  opt=True
  addPdelZm(a,opt)
  
  lev=a->lev
  nlev=dimsizes(lev)
  
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

  pdeldry=opt@pdeldry
  airden=pdeldry/g ; air column mass in kg/m2
  
  MBCPHO_V=a->MBCPHO_V
  BCPHO=MBCPHO_V
  BCPHO@long_name="Hydrophobic black carbon aerosol"  
  BCPHO@units="kg/kg"
  BCPHO(:,0:nlev-2,:)=MBCPHO_V(:,0:nlev-2,:)-MBCPHO_V(:,1:nlev-1,:)
  BCPHO=BCPHO/airden
  a->BCPHO=BCPHO
  
  MBCPHI_V=a->MBCPHI_V
  BCPHI=MBCPHI_V
  BCPHI@long_name="Hydrophilic black carbon aerosol"  
  BCPHI@units="kg/kg"
  BCPHI(:,0:nlev-2,:)=MBCPHI_V(:,0:nlev-2,:)-MBCPHI_V(:,1:nlev-1,:)
  BCPHI=BCPHI/airden
  a->BCPHI=BCPHI
  
  MOCPHO_V=a->MOCPHO_V
  OCPHO=MOCPHO_V
  OCPHO@long_name="Hydrophobic organic carbon aerosol"  
  OCPHO@units="kg/kg"
  OCPHO(:,0:nlev-2,:)=MOCPHO_V(:,0:nlev-2,:)-MOCPHO_V(:,1:nlev-1,:)
  OCPHO=OCPHO/airden
  a->OCPHO=OCPHO
  
  MOCPHI_V=a->MOCPHI_V
  OCPHI=MOCPHI_V
  OCPHI@long_name="Hydrophilic organic carbon aerosol"  
  OCPHI@units="kg/kg"
  OCPHI(:,0:nlev-2,:)=MOCPHI_V(:,0:nlev-2,:)-MOCPHI_V(:,1:nlev-1,:)
  OCPHI=OCPHI/airden
  a->OCPHI=OCPHI
  
  CARAER=BCPHO ; get attributes
  CARAER=CARAER+BCPHI+OCPHO+OCPHI
  CARAER@long_name="Total carbon aerosol"
  a->CARAER=CARAER
  
  BCPHOFRAC=BCPHO ; get attributes
  BCPHOFRAC=BCPHOFRAC/(CARAER+1e-25)
  BCPHOFRAC@long_name="Fractional mass of BCPHO to all carbon aerosols"
  a->BCPHOFRAC=BCPHOFRAC

end


;basedir="/data5/mmills/wa319/carma-soot/"
basedir="/Volumes/Data/Models/waccm/run/carma-soot/"

do year=2034,2042

  caseid="bcrckt03"
  rootdir=basedir+caseid+"/"+"h0zm.annual/"
  file1=caseid+".cam2.h0zm."+year+".annual.nc"
  calcAerRF(rootdir,file1)

  caseid="bccntl02"
  rootdir=basedir+caseid+"/"+"h0zm.annual/"
  file1=caseid+".cam2.h0zm."+year+".annual.nc"
  calcAerRF(rootdir,file1)

end do
