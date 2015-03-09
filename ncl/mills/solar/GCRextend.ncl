load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

rootpath="/glade/p/cesmdata/cseg/inputdata/atm/waccm/gcrs/"
filename="gcr_prod_NO_c130920.nc"

print("opening "+rootpath+filename)
a=addfile(rootpath+filename,"r")
prod=a->prod ; (time, lev, lat, lon)
ds=dimsizes(prod)
date=a->date

date2=ispan(1949,2142,1)*10000+701
nyears2=dimsizes(date2)
date2!0=date!0
datesec=date2 ; get metadata
datesec=0

prod2=new((/nyears2,ds(1),ds(2),ds(3)/),typeof(prod))
prod2!0=prod!0
prod2!1=prod!1
prod2!2=prod!2
prod2!3=prod!3
prod2&lev=prod&lev
prod2&lat=prod&lat
prod2&lon=prod&lon
prod2@long_name=prod@long_name
prod2@units=prod@units

ind1=ind(date.ge.19600000.and.date.le.20109999)
ind2=ind(date2.ge.19600000.and.date2.le.20109999)
prod2(ind2,:,:,:)=(/prod(ind1,:,:,:)/)
delete(ind1)
delete(ind2)

ind1=ind(date.ge.19910000.and.date.le.20019999)
ind2=ind(date2.ge.19490000.and.date2.le.19599999)
prod2(ind2,:,:,:)=(/prod(ind1,:,:,:)/)
delete(ind1)
delete(ind2)

ind1=ind(date.ge.19670000.and.date.le.20109999)
ind2=ind(date2.ge.20110000.and.date2.le.20549999)
prod2(ind2,:,:,:)=(/prod(ind1,:,:,:)/)

ind2=ind(date2.ge.20550000.and.date2.le.20989999)
prod2(ind2,:,:,:)=(/prod(ind1,:,:,:)/)

ind2=ind(date2.ge.20990000.and.date2.le.21429999)
prod2(ind2,:,:,:)=(/prod(ind1,:,:,:)/)

rootpath="/glade/p/work/mmills/inputdata/gcrs/"
creation_date = systemfunc("date +%y%m%d")
outfilename="gcr_prod_NO_1949-2142_c"+creation_date+".nc"
print("rm "+rootpath+outfilename)
b=systemfunc("rm "+rootpath+outfilename)
print("creating "+rootpath+outfilename)
ncidout=addfile(rootpath+outfilename,"c")

; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
filedimdef(ncidout,"time",-1,True)
;
; For a nicer looking netCDF, create a "new line" character.
; This is not necessary.
;
nl = integertochar(10)  ; newline character
;
; Define global attributes.
;
; globalAtt can be of any type. Here logical is used by convention.
;
globalAtt             = True

globalAtt@creation_date     = nl+\
    systemfunc("date")

globalAtt@creator = nl+\
    "Mike Mills, mmills@ucar.edu"
globalAtt@title = nl+\
    "Annual mean production of NO from galactic cosmic rays (GCRs)" ;
globalAtt@Sources = nl+\
     "1960-2010 data from Dan Marsh, extended by copying as follows:"+nl+\
     "1949-1959: copied from 1991-2001."+nl+\
     "Rates for 2011-2142 created by copying the last 4 cycles with values from "+nl+\
     "1967 to 2010 inclusive, mapping to 2011-2054, 2055-2098, and 2099-2142."
globalAtt@history = nl+\
     "Concatenated by program GCRextend.ncl"

fileattdef( ncidout, globalAtt )

ncidout->date=date2
ncidout->datesec=datesec
ncidout->prod=prod2

delete(ncidout)
