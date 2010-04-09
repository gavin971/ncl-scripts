; Program converts solar lean ascii data to netCDF

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

num_ssi_samples = 5*756
;num_samples_time1 = 2000-1610+1
num_samples_time1 = 1999-1610+1
num_samples_time2 = 2009-2000+1
num_samples_time = 2009-1610+1
pmod_tim_scale = .9965d0

ncol = 5

wavelength=asciiread("bandcenters.txt",num_ssi_samples,"double")
  wavelength@_FillValue = -9999.d0
  wavelength@long_name = "Wavelength of band center"
  wavelength@units = "nm" ;
  wavelength!0="wavelength"
  
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
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"

ssi = solar_data(:,2:)
  ssi@_FillValue = -9999.d0
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
 
 print(date(224)+" "+time(224))
 print(date(257)+" "+time(257))

 ref_time_bound=(/time(224),time(257)/)
 ref_time_bound@long_name = "reference time interval boundaries"
 ref_time_bound@units = time@units
 ref_time_bound!0="ref_time_bound"

 tsi_ref=dim_avg_Wrap(tsi(224:257))
 tsi_ref@bounds="ref_time_bound"
 tsi_ref!0="ncl_scalar"
 tsi_ref@long_name="average of tsi over ref_time_bound inclusive"
 tsi_ref@time_op = "average"
 print(tsi_ref)
 
 ssi_ref=dim_avg_Wrap(ssi(wavelength|:,time|224:257))
 ssi_ref@bounds="ref_time_bound"
 ssi_ref@long_name="average of ssi over ref_time_bound inclusive"
 ssi_ref@time_op = "average"
 print(ssi_ref(0:10))

 ncf = addfile("spectral_irradiance_Lean_1610-2009_ann_c100405.nc","c")
 
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
      
  globalAtt@creator = nl+\
      "Mike Mills, mmills@ucar.edu"
  globalAtt@title = nl+\
      "Lean annual spectral solar irradiance for years 1610-2009" ;
  globalAtt@Source_website = nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/CMIP5_solar_irradiance.html"
  globalAtt@Source_data=nl+\
       "ftp://strat50.met.fu-berlin.de/pub/outgoing/_matthes/CMIP5_solardata/spectra_1610_2000a_21Jan09.txt.gz"
  globalAtt@history = nl+\
       "Converted from ASCII by program spectrum.ncl"
  globalAtt@scaling = nl+\
       "Absolute spectral and total irradiance includes PMOD adjustment using 0.9965 for TIM scale with background from"+nl+\ 
       "Wang, Lean & Sheeley, Modeling the sun's magnetic field and irradiance since 1713, Astrophys J (2005) vol. 625 (1) pp. 522-538."
  globalAtt@notes = nl+\
       "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
       "Total and spectral irradiance variations are determined independently of each other,"+nl+\
       "with more uncertainty in the spectral irradiance changes which are determined in 1 nm bins."+nl+\
       "Integrated SSI is therefore varies slightly from TSI."+nl+\
       "Reference spectral and total irradiances derived from average over years 1834-1867 (solar cycles 8-10)"

  fileattdef( ncf, globalAtt )
  
 ncf->tsi_ref =  tsi_ref  
 ncf->ref_time_bound =  ref_time_bound
 ncf->date=date
 ncf->time=time
 ncf->wavelength=wavelength
 ncf->band_width=band_width
 ncf->ssi_ref =  ssi_ref  
 ncf->tsi=tsi
 ncf->ssi=ssi

 delete(ncf)
 
