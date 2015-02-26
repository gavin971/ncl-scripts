load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

  
ostype=getenv("OSTYPE")
host=getenv("HOSTNAME")

sagepath="/glade/p/work/mmills/data/SAGE_II_v7.00/CCMI_Filled_Data_V2/netCDF/"
if (ostype.eq."darwin") then
  sagepath="/Volumes/Data/Data/SAGE2/SAGE_II_v7.00/CCMI_Filled_Data_V2/netCDF/"
end if

b=addfile(sagepath+"CCMI_Filled_Data_V2.nc","r")
ext1020=b->ext1020
ext525=b->ext525
ext452=b->ext452
ext386=b->ext386

sageLat=ext1020&lat
nLat=dimsizes(sageLat)
sageAlt=ext1020&alt
nAlt=dimsizes(sageAlt)

sagedate=ext1020&time
sage1991=ind(floattoint(sagedate).eq.1991)
sageD1=sage1991(0)
print(sageD1+" "+sagedate(sageD1))

;caseid="f.e13.FSDW5.f19_f19.1990-2011.ASv1.P17Tg.16-18km.wCH.01"
caseid="b.e13.B20TRW5TCN.f19_g16.beta13_stratmodal.006.0-15N"

rootpath="/glade/scratch/mmills/archive/"
subdir="/atm/proc/h0zm/"

infiles=systemfunc("ls "+rootpath+caseid+subdir+"*.nc")
history=addfiles(infiles,"r")
ListSetType(history, "cat")

waccmdate=addfiles_GetVar(history,infiles,"date")
waccmD1=ind(waccmdate.eq.19910201)
print(waccmD1+" "+waccmdate(waccmD1))

TROP_Z=addfiles_GetVar(history,infiles,"TROP_Z")
waccmLat=TROP_Z&lat

nYears=13
nMonths=nYears*12
date=sagedate(sageD1:sageD1+nMonths-1)
ODst1020=new((/nMonths,nLat/),"float")
ODst1020!0="time"
ODst1020!1="lat"
ODst1020&time=date
ODst1020&lat=sageLat
ODst1020@long_name="Stratospheric optical depth at 1020 nm"
ODst525=ODst1020
ODst525@long_name="Stratospheric optical depth at 525 nm"
ODst452=ODst1020
ODst452@long_name="Stratospheric optical depth at 452 nm"
ODst386=ODst1020
ODst386@long_name="Stratospheric optical depth at 386 nm"


do d=0,nMonths-1
  ;printVarSummary(TROP_Z)
  ;print(waccmLat+" "+TROP_Z(waccmD1,:))
  
  wd=waccmD1+d
  sd=sageD1+d
  print(d+" "+wd+" "+sd)

  sageTrop=linint1(waccmLat,TROP_Z(wd,:),False,sageLat,0)
  sageTrop=sageTrop/1.e3
  ;print(sageLat+" "+sageTrop)
  do l = 0, nLat-1
;    printVarSummary(sageAlt)
;    printVarSummary(sageTrop)
    strat=ind(sageAlt.ge.sageTrop(l))
    ODst1020(d,l)=sum(ext1020(sd,strat,l))*0.5
    ODst525(d,l)=sum(ext525(sd,strat,l))*0.5
    ODst452(d,l)=sum(ext452(sd,strat,l))*0.5
    ODst386(d,l)=sum(ext386(sd,strat,l))*0.5
    delete(strat)
  end do
end do

print("rm "+sagepath+"ODstrat.nc")
c=systemfunc("rm "+sagepath+"ODstrat.nc")
a=addfile(sagepath+"ODstrat.nc","c")
a->ODst1020=ODst1020
a->ODst525=ODst525
a->ODst452=ODst452
a->ODst386=ODst386
