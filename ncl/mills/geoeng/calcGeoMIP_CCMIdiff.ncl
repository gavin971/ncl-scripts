load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

datestr = systemfunc("date +%y%m%d")

rootpath="/glade/p/work/mmills/jenglish/"
outpath="/glade/p/work/mmills/inputdata/geoV3hy/"
;outpath="/glade/p/work/mmills/inputdata/geoV3/"

cmd="mkdir -p "+outpath
print(cmd)
dum=systemfunc(cmd)

vars="P0,lat,slat,lev,ilev,time,gw,PURSULAD,PURSULMD,PURSULND,PURSULRE,SULDEN,Z3,O3,WTPCT"
vars=vars+",hyai,hybi,hyam,hybm,PS,PSL"

bkgdCase="geo-0-bhn"
bkgdFile="CESM_"+bkgdCase+"_sad_V2_c"+datestr+".nc"
bkgdFilePath=outpath+bkgdFile

cases=(/ "geo-1-bhn", "geo-2-bhn", "geo-5-bhn", "geo-10-bhn" /)
         
nCases = dimsizes(cases)

cmd="ncrcat -O -v "+vars+" "+rootpath+bkgdCase+"/h0zm/"+bkgdCase+\
    ".cam2.h0zm.2015-*.nc "+bkgdFilePath
print(cmd)
dum=systemfunc(cmd)

cmd="ncrename -v PURSULAD,sad "+bkgdFilePath
print(cmd)
dum=systemfunc(cmd)

cmd="ncrename -v PURSULRE,reff "+bkgdFilePath
print(cmd)
dum=systemfunc(cmd)

date=20200015+ispan(1,12,1)*100
date!0="time"
  
a=addfile(bkgdFilePath,"w")
md=a->PURSULMD
H2SO4_mass_bg=md ; get metadata
H2SO4_mass_bg=H2SO4_mass_bg*6.022e23/98.078479 ; g/cm3 ->molec-cm3
H2SO4_mass_bg@long_name="sulfate mass density"
H2SO4_mass_bg@units="molecules/cm3air"
a->H2SO4_mass=H2SO4_mass_bg
delete(md)
ds=dimsizes(H2SO4_mass_bg)
lev=H2SO4_mass_bg&lev

wtpct=a->WTPCT
sulden=a->SULDEN

reff=a->reff
rmode=reff ; get metadata
rmode=rmode * (192.3/wtpct/sulden)^(1./3.) ; adjust for wet radius
rmode=rmode*0.8829536147077603 ; rmode = reff * exp(-2.5 * (ln(sigma))**2), sigma=1.25
rmode=rmode/1e6 ; um -> m
rmode@units="meters"
rmode@long_name="sulfate mode radius"
a->rmode=rmode
delete(rmode)
delete(wtpct)
delete(sulden)

a->sad@units="um2/cm3"      ; change from um2/cm-3
sad_bg=a->sad               ; (time, lev, lat)

O3=a->O3
delete(a)

do i = 0,nCases-1

  filepath=outpath+"/"+"CESM_"+cases(i)+"_sad_V2_c"+datestr+".nc"
  filepath2=outpath+"/"+"CESM_"+cases(i)+"_sad_V2b_c"+datestr+".nc"
  
  cmd="ncrcat -O -v "+vars+" "+rootpath+cases(i)+"/h0zm/"+cases(i)+\
      ".cam2.h0zm.2015-*.nc "+filepath
  print(cmd)
  dum=systemfunc(cmd)

  cmd="ncrename -v PURSULAD,sad "+filepath
  print(cmd)
  dum=systemfunc(cmd)

  cmd="ncrename -v PURSULRE,reff "+filepath
  print(cmd)
  dum=systemfunc(cmd)

  b=addfile(filepath,"w")
  b->date=date
  
  md=b->PURSULMD
  H2SO4_mass=md ; get metadata
  H2SO4_mass=H2SO4_mass*6.022e23/98.078479 ; g/cm3 ->molec-cm3
  H2SO4_mass@long_name="sulfate mass density"
  H2SO4_mass@units="molecules/cm3air"
  delete(md)
  
  wtpct=b->WTPCT
  sulden=b->SULDEN

  reff=b->reff ; (time, lev, lat, lon)
  rmode=reff ; get metadata
  rmode=rmode * (192.3/wtpct/sulden)^(1./3.) ; adjust for wet radius
  rmode=rmode*0.8829536147077603 ; rmode = reff * exp(-2.5 * (ln(sigma))**2), sigma=1.25
  rmode=rmode/1e6 ; um -> m
  rmode@units="meters"
  rmode@long_name="sulfate mode radius"
  b->rmode=rmode
  delete(rmode)  
  delete(wtpct)
  delete(sulden)
  
  b->sad@units="um2/cm3"      ; change from um2/cm-3
  sad=b->sad               ; (time, lev, lat)
  
;  cmd="cp "+filepath+" "+filepath2
;  print(cmd)
;  dum=systemfunc(cmd)
;  c=addfile(filepath2,"w")
;  c->H2SO4_mass=H2SO4_mass
;  delete(c)
  
;  sad(:,{85:},:)=where(O3(:,{85:},:).lt.1.5e-7, \
;       sad(:,{85:},:)-sad_bg(:,{85:},:), \
;       sad(:,{85:},:))
  
;  H2SO4_mass(:,{85:},:)=where(O3(:,{85:},:).lt.1.5e-7, \
;       H2SO4_mass(:,{85:},:)-H2SO4_mass_bg(:,{85:},:), \
;       H2SO4_mass(:,{85:},:))
  do j=0,ds(0)-1
;    print("month "+j)
    do k=0,ds(1)-1
      if (lev(k).gt.85.0) then
        do l=0,ds(2)-1
          if (O3(j,k,l).lt.1.5e-7) then
            sad(j,k,l)=max((/sad(j,k,l)-sad_bg(j,k,l),0.0/))
            H2SO4_mass(j,k,l)=max((/H2SO4_mass(j,k,l)-H2SO4_mass_bg(j,k,l),0.0/))
          end if
        end do
      end if
    end do
  end do
  
  b->sad=sad
  b->H2SO4_mass=H2SO4_mass
  
  delete(b)
  delete(sad)
  delete(H2SO4_mass)
  delete(filepath)
end do
