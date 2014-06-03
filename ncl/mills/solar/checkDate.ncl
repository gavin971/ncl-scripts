load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

rootpath = "/glade/p/cesmdata/cseg/inputdata/atm/waccm/solar/"

file1="wasolar_1947-2012_daily_c130227.nc"
;file1="wasolar_c140408.nc"
print(file1)

a=addfile(rootpath+file1,"r")

fdate=a->date
ntimes=dimsizes(fdate)

year1=floattoint(fdate(0)/10000)
year2=floattoint(fdate(ntimes-1)/10000)

datearr=yyyymmdd_time(year1,year2,"integer")

d1=ind(datearr.eq.fdate(0))
d2=ind(datearr.eq.fdate(ntimes-1))

ntimes2=d2-d1+1

missing=ntimes2-ntimes

print("ntimes="+ntimes+", ntimes2="+ntimes2+", missing:"+missing)

i1=0
do i2=0,ntimes2-1
  if (datearr(d1+i2).eq.fdate(i1)) then
    i1=i1+1
  else
    print("missing: "+datearr(d1+i2))
  end if
end do
