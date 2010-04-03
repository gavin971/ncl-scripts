; Program converts solar lean ascii data to netCDF

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

num_ssi_samples = 5*756
;num_samples_time1 = 2000-1610+1
num_samples_time1 = 1999-1610+1
num_samples_time2 = 2009-2000+1
num_samples_time = 2009-1610+1
pmod_tim_scale = .9965d0

ncol = 5

wvl=asciiread("bandcenters.txt",num_ssi_samples,"double")
  wvl@_FillValue = -9999.d0
  wvl@long_name = "Wavelength of band center"
  wvl@units = "nm" ;
  wvl!0="wavelength"
  
band_width=asciiread("bandwidths.txt",num_ssi_samples,"double")
  band_width@_FillValue = -9999.d0
  band_width@long_name = "Wavelength width of band"
  band_width@units = "nm"
  band_width!0="wavelength"

solar_data1=asciiread("spectra_1610_2000a_21Jan09.txt",\
       (/num_samples_time1, num_ssi_samples +2/),"double")
       
solar_data2=asciiread("spectra_2000_2009a_16Mar10.txt",\
       (/num_samples_time2, num_ssi_samples +2/),"double")
       
solar_data=new((/num_samples_time,num_ssi_samples +2/),double)
solar_data(0:num_samples_time1-1,:)=solar_data1
solar_data(num_samples_time1:num_samples_time-1,:)=solar_data2

time = solar_data(:,0)*365.d0 + 181.d0
  time@_FillValue = -9999.d0
  time@units = "days since 0000-01-01 00:00:00"
  time@time_origin = "01-JAN-0000"
  time@axis = "T"
  time@calendar = "noleap"
  time!0="time"

datestr=flt2string(doubletofloat(solar_data(:,0)))+"0701"
date=stringtointeger(datestr)
  date@format="YYYYMMDD"
  date!0="time"

tsi = solar_data(:,1)
  tsi@_FillValue = -9999.d0
  tsi@level = "PMOD adjusted to SORCE-TIM using 0.9965 conversion factor" ;
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"

ssi = solar_data(:,2:)
  ssi@_FillValue = -9999.d0
  ssi@level = "PMOD adjusted to SORCE-TIM using 0.9965 conversion factor"
  ssi@long_name = "Solar Spectral Irradiance at 1 a.u."
  ssi@units = "mW/m^2/nm"
  ssi!0="time"
  ssi!1="wavelength"
  
 print ("verify first date:"+date(0))
 print ("verify last date:"+date(num_samples_time-1))
 print ("verify first tsi:"+tsi(0))
 print ("verify last tsi:"+tsi(num_samples_time-1))
 print ("verify first ssi at first time:"+ssi(0,0))
 print ("verify last ssi at first time:"+ssi(0,num_ssi_samples-1))
 print ("verify last ssi at last time:"+ssi(num_samples_time-1,num_ssi_samples-1))
  
 tsi = tsi * pmod_tim_scale
 ssi = ssi * pmod_tim_scale
 
 ncf = addfile("spectra_1610_2009a_31Mar10.nc","c")
 
; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
 filedimdef(ncf,"time",-1,True)
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
      
  globalAtt@user = nl+\
      "Mike Mills"
  globalAtt@title = nl+\
      "Lean spectral solar irradiance 1610-2009" ;
  globalAtt@Source_website = nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/CMIP5_solar_irradiance.html"
  globalAtt@Source_data=nl+\
       "ftp://strat50.met.fu-berlin.de/pub/outgoing/_matthes/CMIP5_solardata/spectra_1610_2000a_21Jan09.txt.gz"
  globalAtt@history = nl+\
       "Converted from ASCII by program spectrum.ncl"
  globalAtt@scaling = nl+\
       "Absolute irradiance includes PMOD adjustment using 0.9965 for TIM scale with Wang et al (2005) background"

  fileattdef( ncf, globalAtt )

 ncf->time=time
 ncf->date=date
 ncf->wvl=wvl
 ncf->band_width=band_width
 ncf->tsi=tsi
 ncf->ssi=ssi
