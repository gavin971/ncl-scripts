; Program extends solar variability input file from 1610-2009 to 1610-2140
; by repeating the 4 solar cycles (44 years) prior to 2008, 3 times for
; 2009-2140

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

datapath="/Volumes/Data/Models/ccsm/inputdata/atm/waccm/solar/"

num_samples_time1 = (2008-1950+1)*365
num_samples_time2 = num_samples_time1+44*3*365

ncf = addfile(datapath+"spectral_irradiance_Lean_1950-2009_daily_c100507.nc","r")

wavelength=ncf->wavelength
band_width=ncf->band_width
num_ssi_samples=dimsizes(wavelength)

print("reading ncf->tsi")
tsi_ref=ncf->tsi_ref
print("reading ncf->ssi")
ssi_ref=ncf->ssi_ref
print("reading ncf->ref_time_bound")
ref_time_bound=ncf->ref_time_bound

time=new(num_samples_time2, double)
time(0:num_samples_time1-1)=ncf->time(0:num_samples_time1-1)
  time@_FillValue = -9999.d0
  time@units = "days since 0000-01-01 00:00:00"
  time@time_origin = "01-JAN-0000"
  time@axis = "T"
  time@calendar = "noleap"
  time!0="time"
do day = num_samples_time1,num_samples_time2-1
  time(day)=time(day-1)+1.d0
end do

date=new(num_samples_time2, integer)
date(0:num_samples_time1-1)=ncf->date(0:num_samples_time1-1)
  date@format="YYYYMMDD"
  date!0="time"
  date&time=time
  
tsi=new(num_samples_time2, double)
tsi(0:num_samples_time1-1)=(/ncf->tsi(0:num_samples_time1-1)/)
  tsi@_FillValue = -9999.d0
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"
  tsi&time=time

ssi = new((/num_samples_time2,num_ssi_samples/), double)
ssi(0:num_samples_time1-1,:)=(/ncf->ssi(0:num_samples_time1-1,:)/)
  ssi@_FillValue = -9999.d0
  ssi@long_name = "Solar Spectral Irradiance at 1 a.u."
  ssi@units = "mW/m^2/nm"
  ssi!0="time"
  ssi&time=time
  ssi!1="wavelength"
  ssi&wavelength=wavelength

do day = num_samples_time1,num_samples_time2-1
  date(day)=(/date(day-365)/)+10000
  if (mod(day,365) .eq.0 ) then
    print(day+" "+date(day)+" "+time(day)/365.d0+" "+time(day)+" "+time(day-1))
  end if
end do

y1=2008-43-1950
y2=2008-1950
d1=y1*365
d2=y2*365+364

do i = 0,2
  year1=2009+i*44-1950
  year2=year1+43
  print(i+": years "+year1+" - "+year2)
  day1=year1*365
  day2=year2*365+364
  print(i+": days "+day1+" - "+day2)
  print(i+": time "+time(day1)+" - "+time(day2))
  tsi(day1:day2)=(/tsi(d1:d2)/)
  ssi(day1:day2,:)=(/ssi(d1:d2,:)/)
end do

i=num_samples_time2-1
print(i+" "+date(i)+" "+time(i)+" "+tsi(i))
i=num_samples_time1-1
print(i+" "+date(i)+" "+time(i)+" "+tsi(i))
i=num_samples_time1
print(i+" "+date(i)+" "+time(i)+" "+tsi(i))

do year = 2008,2140
  day=(year-1950)*365
  print(year+" "+day+" "+time(day))
end do

 ncf_out = addfile(datapath+"spectral_irradiance_Lean_1950-2140_daily_c100804.nc","c")
 
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
      "Lean daily spectral solar irradiance for years 1950-2140" ;
  globalAtt@Source_website = nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/index.html"
  globalAtt@Source_data=nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/spectra_1950_1999d_11yr_new.zip"+nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/spectra_2000-2006_daily_txt.gz"
  globalAtt@history = nl+\
       "Extended from spectral_irradiance_Lean_1950-2009_daily_c100507.nc by program future_Lean_daily.ncl."+nl+\
       "Converted from ASCII by program spectrum_daily.ncl"
  globalAtt@scaling = nl+\
       "Absolute spectral and total irradiance includes PMOD adjustment using 0.9965"+nl+\
       "for TIM scale with background from Wang, Lean & Sheeley, Modeling the sun's"+nl+\
       "magnetic field and irradiance since 1713, Astrophys J (2005) vol. 625 (1) pp. 522-538."
  globalAtt@notes = nl+\
       "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
       "Total and spectral irradiance variations are determined independently of each other,"+nl+\
       "with more uncertainty in the spectral irradiance changes which are determined in 1 nm bins."+nl+\
       "Integrated SSI is therefore varies slightly from TSI."+nl+\
       "Leap days (Feb 29) have been removed from source to conform to 365-day model year."+nl+\
       "Reference spectral and total irradiances derived from average over years 1834-1867 (solar"+nl+\
       "cycles 8-10) using annual average data (spectral_irradiance_Lean_1610-2140_ann_c100405.nc)."+nl+\
       "Irradiances for 2009-2140 created by repeating the last 4 cycles (cycle 20-23)"+nl+\
       "with values from 1965 to 2008 inclusive mapping to 2009-2052, 2053-2096, 2097-2140"+nl+\
       "No background trends are included in TSI or SSI from 2008-2140." ;

  fileattdef( ncf_out, globalAtt )
  
 ncf_out->tsi_ref =  tsi_ref  
 ncf_out->ref_time_bound =  ref_time_bound
 print("writing date")
 ncf_out->date=date
; print("writing time")
; ncf_out->time=time(0:num_samples_time2-1)
 
 do year = 2008,2015
   day=(year-1950)*365
   print(year+" "+day+" "+time(day)+" "+ncf_out->time(day))
 end do
 
; print("writing wavelength")
; ncf_out->wavelength=wavelength
 
; do year = 2008,2015
;   day=(year-1950)*365
;   print(year+" "+day+" "+time(day)+" "+ncf_out->time(day))
; end do
 
 print("writing band_width")
 ncf_out->band_width=band_width
 
 do year = 2008,2015
   day=(year-1950)*365
   print(year+" "+day+" "+time(day)+" "+ncf_out->time(day))
 end do
 
 print("writing ssi_ref")
 ncf_out->ssi_ref =  ssi_ref  
 
 do year = 2008,2015
   day=(year-1950)*365
   print(year+" "+day+" "+time(day)+" "+ncf_out->time(day))
 end do
 
 print("writing tsi")
 ncf_out->tsi=tsi
 
 do year = 2008,2015
   day=(year-1950)*365
   print(year+" "+day+" "+time(day)+" "+ncf_out->time(day))
 end do
 
 print("writing ssi")
 ncf_out->ssi=ssi
 
 do year = 2008,2015
   day=(year-1950)*365
   print(year+" "+day+" "+time(day)+" "+ncf_out->time(day))
 end do
 
 
 delete(ncf_out)
