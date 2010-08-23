; Program extends WACCM Solar Proton Event input data file from 
; 19630101-20100101 to 196300101-2140 by repeating the 4 solar 
; cycles (44 years) prior to 2008, 3 times for 2009-2140

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

datapath="$MODEL_PATH/ccsm/inputdata/atm/waccm/solar/"

ncf = addfile(datapath+"spes.1963-2009.c100804.nc","r")

date1=ncf->date
pressure=ncf->pressure

npress=dimsizes(pressure)
ndate1=dimsizes(date1)

day1965=ind(date1.eq.19650101)
day2008=ind(date1.eq.20081231)
nday=day2008-day1965

print("nday="+nday)
print("44 years="+44*365.25d0)
print("nday/365="+nday/365.d0)

num_samples_time2 = day2008+4+nday*3

date=new(num_samples_time2, integer)
date(0:day2008)=date1(0:day2008)
  date@format="YYYYMMDD"
  date!0="time"

datesec=new(num_samples_time2, integer)
datesec(:)=0
  datesec@long_name="current seconds of current date"
  datesec!0="time"
  
Prod=new((/num_samples_time2, npress/), float)
Prod(0:day2008,:)=(/ncf->Prod(0:day2008,:)/)
  Prod@long_name = "Ion pair production rate"
  Prod@units = "/cm3/sec" ;
  Prod!0="time"
  Prod!1="pressure"

do i = 0,2
  day1=day2008+1+nday*i
  day2=day2008+1+nday*(i+1)
  print(i+": days "+day1+" - "+day2+" of "+num_samples_time2)
  Prod(day1:day2,:)=(/Prod(day1965:day2008,:)/)
  date(day1:day2)=(/date(day1965:day2008)/)+440000*(i+1)
  print(i+": date "+date(day1)+" - "+date(day2))
end do

; Set endpoint date 21410101 to zero
Prod(num_samples_time2-2,:)=0
date(num_samples_time2-2)=21410101
; Set endpoint date 25000101 to zero
Prod(num_samples_time2-1,:)=0
date(num_samples_time2-1)=25000101

i=num_samples_time2-1
print(i+" "+date(i)+" "+Prod(i,15))
i=num_samples_time2-2
print(i+" "+date(i)+" "+Prod(i,15))
i=num_samples_time2-3
print(i+" "+date(i)+" "+Prod(i,15))
i=day2008
print(i+" "+date(i)+" "+Prod(i,15))
i=day2008+1
print(i+" "+date(i)+" "+Prod(i,15))

do year = 2008,2140
  day=day2008+(year-2008)*365
  print(year+" "+day+" "+date(day)+" "+Prod(i,15))
end do

 ncf_out = addfile(datapath+"spe_data_1963-2140_c100823.nc","c")
 
; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
 filedimdef(ncf_out,"time",-1,True)
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
      "Time-dependent Solar Proton Event data 1963-2140" ;
  globalAtt@history = nl+\
       "Extended from spes.1963-2009.c100804.nc by program future_spe.ncl."
  globalAtt@notes = nl+\
       "Parameters for 2009-2140 created by repeating the last 4 cycles (cycle 20-23)"+nl+\
       "with values from 1965 to 2008 inclusive mapping to 2009-2052, 2053-2096, 2097-2140." +nl+\
       "Leap days (Feb 29) have been left in the data, and may not match actual leap days " +nl+\
       "for all years after 2008."

  fileattdef( ncf_out, globalAtt )
  
 print("writing date")
 ncf_out->date=date
 
 print("writing Prod")
 ncf_out->Prod =  Prod  
 
 print("writing datesec")
 ncf_out->datesec =  datesec  
 
 print("writing pressure")
 ncf_out->pressure =  pressure  
 
 
 delete(ncf_out)
