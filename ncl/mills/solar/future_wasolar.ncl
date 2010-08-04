; Program extends WACCM solar parameters input file from 19470410-20091022 to 
; 19470410-2140 by repeating the 4 solar cycles (44 years) prior to 2008, 
; 3 times for 2009-2140

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

datapath="$MODEL_PATH/ccsm/inputdata/atm/waccm/phot/"

ncf = addfile(datapath+"wasolar_c100106.nc","r")

date1=ncf->date

day1965=ind(date1.eq.19650101)
day2008=ind(date1.eq.20081231)
nday=day2008-day1965

print("nday="+nday)
print("44 years="+44*365.25d0)
print("nday/365="+nday/365.d0)

num_samples_time2 = day2008+2+nday*3

date=new(num_samples_time2, integer)
date(0:day2008)=date1(0:day2008)
  date@format="YYYYMMDD"
  date!0="time"
  
f107=new(num_samples_time2, double)
f107(0:day2008)=(/ncf->f107(0:day2008)/)
  f107@_FillValue = -9999.d0
  f107@long_name = "10.7 cm solar radio flux (F10.7)"
  f107@units = "10^-22 W m^-2 Hz^-1" ;
  f107!0="time"

f107a=new(num_samples_time2, double)
f107a(0:day2008)=(/ncf->f107a(0:day2008)/)
  f107a@_FillValue = -9999.d0
  f107a@long_name = "81-day centered mean of 10.7 cm solar radio flux (F10.7)"
  f107a@units = "10^-22 W m^-2 Hz^-1" ;
  f107a!0="time"

kp=new(num_samples_time2, double)
kp(0:day2008)=(/ncf->kp(0:day2008)/)
  kp@_FillValue = -9999.d0
  kp@long_name = "Daily planetary K index"
  kp!0="time"

ap=new(num_samples_time2, double)
ap(0:day2008)=(/ncf->ap(0:day2008)/)
  ap@_FillValue = -9999.d0
  ap@long_name = "Daily planetary a index"
  ap@units = "nanoTeslas" ;
  ap!0="time"

isn=new(num_samples_time2, double)
isn(0:day2008)=(/ncf->isn(0:day2008)/)
  isn@_FillValue = -9999.d0
  isn@long_name = "International Sunspot Number"
  isn!0="time"

do i = 0,2
  day1=day2008+1+nday*i
  day2=day2008+1+nday*(i+1)
  print(i+": days "+day1+" - "+day2+" of "+num_samples_time2)
  f107(day1:day2)=(/f107(day1965:day2008)/)
  f107a(day1:day2)=(/f107a(day1965:day2008)/)
  kp(day1:day2)=(/kp(day1965:day2008)/)
  ap(day1:day2)=(/ap(day1965:day2008)/)
  isn(day1:day2)=(/isn(day1965:day2008)/)
  date(day1:day2)=(/date(day1965:day2008)/)+440000*(i+1)
  print(i+": date "+date(day1)+" - "+date(day2))
end do

i=num_samples_time2-1
print(i+" "+date(i)+" "+f107(i))
i=day2008
print(i+" "+date(i)+" "+f107(i))
i=day2008+1
print(i+" "+date(i)+" "+f107(i))

do year = 2008,2140
  day=day2008+(year-2008)*365
  print(year+" "+day+" "+date(day)+" "+f107(day))
end do

 ncf_out = addfile(datapath+"wasolar_1948-2140_daily_c100804.nc","c")
 
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
      "Time-dependent solar activity conditions 1948-2140" ;
  globalAtt@Source_website = nl+\
       "http://www.sec.noa.gov"
  globalAtt@history = nl+\
       "Extended from wasolar_c100106.nc by program future_wasolar.ncl."
  globalAtt@notes = nl+\
       "Parameters for 2009-2140 created by repeating the last 4 cycles (cycle 20-23)"+nl+\
       "with values from 1965 to 2008 inclusive mapping to 2009-2052, 2053-2096, 2097-2140." +nl+\
       "Leap days (Feb 29) have been left in the data, and may not match actual leap days " +nl+\
       "for all years after 2008."

  fileattdef( ncf_out, globalAtt )
  
 print("writing date")
 ncf_out->date=date
 
 print("writing f107")
 ncf_out->f107 =  f107  
 
 print("writing f107a")
 ncf_out->f107a =  f107a  
 
 print("writing kp")
 ncf_out->kp =  kp  
 
 print("writing ap")
 ncf_out->ap =  ap  
 
 print("writing isn")
 ncf_out->isn =  isn  
 
 
 delete(ncf_out)
