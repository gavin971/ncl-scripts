; Program extends solar variability input file from 1610-2009 to 1610-2140
; by repeating the 4 solar cycles (44 years) prior to 2008, 3 times for
; 2009-2140

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

num_samples_time1 = 2008-1610+1
num_samples_time2 = num_samples_time1+44*3

ncf = addfile("spectral_irradiance_Lean_1610-2009_ann_c100405.nc","r")

wavelength=ncf->wavelength
band_width=ncf->band_width
num_ssi_samples=dimsizes(wavelength)

tsi_ref=ncf->tsi_ref
ssi_ref=ncf->ssi_ref
ref_time_bound=ncf->ref_time_bound

time=new(num_samples_time2, double)
time(0:num_samples_time1-1)=ncf->time(0:num_samples_time1-1)
  time@_FillValue = -9999.d0
  time@units = "days since 0000-01-01 00:00:00"
  time@time_origin = "01-JAN-0000"
  time@axis = "T"
  time@calendar = "noleap"
  time!0="time"

date=new(num_samples_time2, integer)
date(0:num_samples_time1-1)=ncf->date(0:num_samples_time1-1)
  date@format="YYYYMMDD"
  date!0="time"
  
tsi=new(num_samples_time2, double)
tsi(0:num_samples_time1-1)=ncf->tsi(0:num_samples_time1-1)
  tsi@_FillValue = -9999.d0
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"

ssi = new((/num_samples_time2,num_ssi_samples/), double)
ssi(0:num_samples_time1-1,:)=ncf->ssi(0:num_samples_time1-1,:)
  ssi@_FillValue = -9999.d0
  ssi@long_name = "Solar Spectral Irradiance at 1 a.u."
  ssi@units = "mW/m^2/nm"
  ssi!0="time"
  ssi!1="wavelength"

do year = 2009, 2140
  time(year-1610)=year*365.d0 + 181.d0
  datestr=flt2string(int2flt(year))+"0701"
  date(year-1610)=stringtointeger(datestr)
end do

y1=2008-43-1610
y2=2008-1610

do i = 0,2
  year1=2009+i*44-1610
  year2=year1+43
  print(i+": "+year1+" - "+year2)
  tsi(year1:year2)=tsi(y1:y2)
  ssi(year1:year2,:)=ssi(y1:y2,:)
end do

i=num_samples_time2-1
print(i+" "+date(i)+" "+tsi(i))
i=num_samples_time1-1
print(i+" "+date(i)+" "+tsi(i))
i=num_samples_time1
print(i+" "+date(i)+" "+tsi(i))

 ncf_out = addfile("spectral_irradiance_Lean_1610-2140_ann_c100408.nc","c")
 
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
      "Lean annual spectral solar irradiance for years 1610-2140" ;
  globalAtt@Source_website = nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/CMIP5_solar_irradiance.html"
  globalAtt@Source_data=nl+\
       "ftp://strat50.met.fu-berlin.de/pub/outgoing/_matthes/CMIP5_solardata/spectra_1610_2000a_21Jan09.txt.gz"
  globalAtt@history = nl+\
       "Adapted from spectra_1610_2009a_31Mar10.nc by program future.ncl."+nl+\
       "Converted from ASCII by program spectrum.ncl"
  globalAtt@scaling = nl+\
       "Absolute spectral and total irradiance includes PMOD adjustment using 0.9965 for TIM scale with background from"+nl+\ 
       "Wang, Lean & Sheeley, Modeling the sun's magnetic field and irradiance since 1713, Astrophys J (2005) vol. 625 (1) pp. 522-538."
  globalAtt@notes = nl+\
       "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
       "Total and spectral irradiance variations are determined independently of each other,"+nl+\
       "with more uncertainty in the spectral irradiance changes which are determined in 1 nm bins."+nl+\
       "Integrated SSI is therefore varies slightly from TSI."+nl+\
       "Reference spectral and total irradiances derived from average over years 1834-1867 (solar cycles 8-10)"+nl+\
       "Irradiances for 2009-2140 created by repeating the last 4 cycles (cycle 20-23)"+nl+\
       "with values from 1965 to 2008 inclusive mapping to 2009-2052, 2053-2096, 2097-2140"+nl+\
       "No background trends is included in TSI or SSI from 2008-2140." ;

  fileattdef( ncf_out, globalAtt )
  
 ncf_out->tsi_ref =  tsi_ref  
 ncf_out->ref_time_bound =  ref_time_bound
 ncf_out->date=date
 ncf_out->time=time
 ncf_out->wavelength=wavelength
 ncf_out->band_width=band_width
 ncf_out->ssi_ref =  ssi_ref  
 ncf_out->tsi=tsi
 ncf_out->ssi=ssi
