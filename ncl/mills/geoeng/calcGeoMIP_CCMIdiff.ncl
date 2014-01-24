load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

datestr = systemfunc("date +%y%m%d")

rootpath="/glade/p/work/mmills/jenglish/"

bkgdCase="geo-0-bhn"
bkgdFile="CESM_"+bkgdCase+"_sad_V2_c"+datestr+".nc"

cases=(/ "geo-1-bhn", "geo-2-bhn", "geo-5-bhn", "geo-10-bhn" /)
         
nCases = dimsizes(cases)

cmd="ncrcat -O -v P0,lat,slat,lev,ilev,time,gw,PURSULAD,PURSULMD,PURSULND,PURSULRE,SULDEN,Z3,O3 "+\
  rootpath+bkgdCase+"/h0zm/"+bkgdCase+".cam2.h0zm.2015-*.nc "+rootpath+bkgdCase+"/"+bkgdFile
print(cmd)
dum=systemfunc(cmd)

cmd="ncrename -v PURSULAD,sad "+rootpath+bkgdCase+"/"+bkgdFile
print(cmd)
dum=systemfunc(cmd)

cmd="ncrename -v PURSULRE,reff "+rootpath+bkgdCase+"/"+bkgdFile
print(cmd)
dum=systemfunc(cmd)
  
a=addfile(rootpath+bkgdCase+"/"+bkgdFile,"w")
md=a->PURSULMD
H2SO4_mass_bg=md ; get metadata
H2SO4_mass_bg=H2SO4_mass_bg*6.022e23/98.078479 ; g/cm3 ->molec-cm3
H2SO4_mass_bg@long_name="sulfate mass density"
a->H2SO4_mass=H2SO4_mass_bg
delete(md)

reff=a->reff
rmode=reff ; get metadata
rmode=rmode*0.8829536147077603 ; rmode = reff * exp(-2.5 * (ln(sigma))**2), sigma=1.25
rmode=rmode/1e6 ; um -> m
rmode@units="meters"
rmode@long_name="sulfate mode radius"
a->rmode=rmode
delete(rmode)

sad_bg=a->sad               ; (time, lev, lat)

O3=a->O3
delete(a)

do i = 0,nCases-1

  filepath=rootpath+cases(i)+"/"+"CESM_"+cases(i)+"_sad_V2_c"+datestr+".nc"
  
  cmd="ncrcat -O -v P0,lat,slat,lev,ilev,time,gw,PURSULAD,PURSULMD,PURSULND,PURSULRE,SULDEN,Z3,O3 "+\
    rootpath+bkgdCase+"/h0zm/"+bkgdCase+".cam2.h0zm.2015-*.nc "+filepath
  print(cmd)
  dum=systemfunc(cmd)

  cmd="ncrename -v PURSULAD,sad "+filepath
  print(cmd)
  dum=systemfunc(cmd)

  cmd="ncrename -v PURSULRE,reff "+filepath
  print(cmd)
  dum=systemfunc(cmd)

  b=addfile(filepath,"w")
  
  md=b->PURSULMD
  H2SO4_mass=md ; get metadata
  H2SO4_mass=H2SO4_mass*6.022e23/98.078479 ; g/cm3 ->molec-cm3
  H2SO4_mass@long_name="sulfate mass density"
  delete(md)

  reff=b->reff
  rmode=reff ; get metadata
  rmode=rmode*0.8829536147077603 ; rmode = reff * exp(-2.5 * (ln(sigma))**2), sigma=1.25
  rmode=rmode/1e6 ; um -> m
  rmode@units="meters"
  rmode@long_name="sulfate mode radius"
  b->rmode=rmode
  delete(rmode)  
  
  sad=b->sad               ; (time, lev, lat)
  sad(:,{85:},:)=where(O3(:,{85:},:).lt.1.5e-9, \
       sad(:,{85:},:)-sad_bg(:,{85:},:), \
       sad(:,{85:},:))
  b->sad=sad
  
  H2SO4_mass(:,{85:},:)=where(O3(:,{85:},:).lt.1.5e-9, \
       H2SO4_mass(:,{85:},:)-H2SO4_mass_bg(:,{85:},:), \
       H2SO4_mass(:,{85:},:))
  b->H2SO4_mass=H2SO4_mass
  
  delete(b)
  delete(sad)
  delete(H2SO4_mass)
  delete(filepath)
end do
