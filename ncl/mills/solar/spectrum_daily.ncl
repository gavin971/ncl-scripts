; Program converts Judith Lean's daily solar ascii data to netCDF
; Source data covers 1950-2009, in 2 files

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

ascii_path=getenv ("MODEL_PATH")+"/ccsm/solar/daily/"
cdfpath=getenv ("MODEL_PATH")+"/ccsm/inputdata/atm/waccm/solar/"

file1=ascii_path+"spectra_1950_1999d_11yr_new.txt"
file2=ascii_path+"katja_spectra_2000_2009d_16Mar10.txt"

num_ssi_samples = 5*756

nyears1=1999-1950+1
nleap1=12
num_samples_time1 = nyears1*365+nleap1

print ("---------------------------------------------------------") 
print ("nyears1="+nyears1)
print ("nleap1="+nleap1)
print ("num_samples_time1="+num_samples_time1)

nyears2=2009-2000+1
nleap2=3
num_samples_time2 = nyears2*365+nleap2

print ("---------------------------------------------------------") 
print ("nyears2="+nyears2)
print ("nleap2="+nleap2)
print ("num_samples_time2="+num_samples_time2)

nyears=nyears1+nyears2
num_samples_time3 = num_samples_time1+num_samples_time2
print ("---------------------------------------------------------") 
print ("nyears="+nyears)
print ("num_samples_time3="+num_samples_time3)

; remove leap day from data set
num_samples_time=nyears*365
print ("num_samples_time="+num_samples_time)
print ("---------------------------------------------------------") 

pmod_tim_scale = .9965d0

wavelength=asciiread(ascii_path+"bandcenters.txt",num_ssi_samples,"double")
  wavelength@_FillValue = -9999.d0
  wavelength@long_name = "Wavelength of band center"
  wavelength@units = "nm" ;
  wavelength!0="wavelength"
  
band_width=asciiread(ascii_path+"bandwidths.txt",num_ssi_samples,"double")
  band_width@_FillValue = -9999.d0
  band_width@long_name = "Wavelength width of band"
  band_width@units = "nm"
  band_width!0="wavelength"

solar_data1=asciiread(file1,\
       (/num_samples_time1, num_ssi_samples +4/),"double")
       
solar_data2=asciiread(file2,\
       (/num_samples_time2, num_ssi_samples +4/),"double")
       
solar_data3=new((/num_samples_time3,num_ssi_samples +4/),double)
solar_data3(0:num_samples_time1-1,:)=solar_data1
solar_data3(num_samples_time1:num_samples_time3-1,:)=solar_data2

solar_data=new((/num_samples_time,num_ssi_samples +4/),double)

; Screen out leap days (Feb 29)
j=0
do i = 0, num_samples_time3-1
  if ((solar_data3(i,1).ne.2).or.(solar_data3(i,2).ne.29)) then
    solar_data(j,:)=solar_data3(i,:)
    j=j+1
  else  
    print("Excluding leap day: "+solar_data3(i,0)+" "+solar_data3(i,1)+" "+solar_data3(i,2))
  end if
end do

delete(solar_data1)
delete(solar_data2)
delete(solar_data3)

time = solar_data(:,0)*365.d0 + \
       day_of_year(conform_dims(num_samples_time,1999,-1),\
       doubletoint(solar_data(:,1)),doubletoint(solar_data(:,2)))
  time@_FillValue = -9999.d0
  time@units = "days since 0000-01-01 00:00:00"
  time@time_origin = "01-JAN-0000"
  time@axis = "T"
  time@calendar = "noleap"
  time!0="time"

date=doubletoint(solar_data(:,0)*10000+solar_data(:,1)*100+solar_data(:,2))
  date@format="YYYYMMDD"
  date!0="time"
  
 print(date(0)+" "+time(0)+" "+solar_data(0,0:2))
 print(date(365)+" "+time(365)+" "+solar_data(365,0:2))
 print(date(num_samples_time-1)+" "+time(num_samples_time-1)+" "+solar_data(num_samples_time-1,0:2))

tsi = solar_data(:,3)
  tsi@_FillValue = -9999.d0
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"

ssi = solar_data(:,4:)
  ssi@_FillValue = -9999.d0
  ssi@long_name = "Solar Spectral Irradiance at 1 a.u."
  ssi@units = "mW/m^2/nm"
  ssi!0="time"
  ssi!1="wavelength"
  
 print ("---------------------------------------------------------") 
 print ("verify first date: "+date(0))
 print ("verify first tsi: "+sprintf("%11.6f",tsi(0)))
 print ("verify first ssi at first time: "+ssi(0,0))
 print ("verify last ssi at first time: "+ssi(0,num_ssi_samples-1))
 print ("---------------------------------------------------------")
 print ("verify last date, first file: "+date(nyears1*365-1))
 print ("verify last tsi, first file: "+sprintf("%11.6f",tsi(nyears1*365-1)))
 print ("verify last ssi, last date, first file: "+ssi(nyears1*365-1,num_ssi_samples-1))
 print ("---------------------------------------------------------") 
 print ("verify first date, second file: "+date(nyears1*365))
 print ("verify first tsi, second file: "+sprintf("%11.6f",tsi(nyears1*365)))
 print ("verify first ssi at first time, second file: "+ssi(nyears1*365,0))
 print ("verify last ssi at first time, second file: "+ssi(nyears1*365,num_ssi_samples-1))
 print ("---------------------------------------------------------") 
 print ("verify last date: "+date(num_samples_time-1))
 print ("verify last tsi: "+sprintf("%11.6f",tsi(num_samples_time-1)))
 print ("verify last ssi at last time: "+ssi(num_samples_time-1,num_ssi_samples-1))
 print ("---------------------------------------------------------") 
 
 delete(solar_data)
 
 tsi = tsi * pmod_tim_scale
 ssi = ssi * pmod_tim_scale
 
 ncf = addfile(cdfpath+"spectral_irradiance_Lean_1950-2009_daily_c100507.nc","c")
 ncf_ann = addfile(cdfpath+"spectral_irradiance_Lean_1610-2140_ann_c100408.nc","r")
 
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
      "Lean annual daily solar irradiance for years 1950-2009" ;
  globalAtt@Source_website = nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/index.html"
  globalAtt@Source_data=nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/spectra_1950_1999d_11yr_new.zip"+nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/spectra_2000-2006_daily_txt.gz"
  globalAtt@history = nl+\
       "Converted from ASCII by program spectrum_daily.ncl"
  globalAtt@scaling = nl+\
       "Absolute spectral and total irradiance includes PMOD adjustment using 0.9965"+nl+\
       "for TIM scale with background from Wang, Lean & Sheeley, Modeling the sun's"+nl+\
       "magnetic field and irradiance since 1713, Astrophys J (2005) vol. 625 (1) pp. 522-538."
  globalAtt@notes = nl+\
       "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
       "Total and spectral irradiance variations are determined independently of each other,"+nl+\
       "with more uncertainty in the spectral irradiance changes which are determined in 1 "+nl+\
       "nm bins.  Integrated SSI is therefore varies slightly from TSI."+nl+\
       "Leap days (Feb 29) have been removed from source to conform to 365-day model year."+nl+\
       "Reference spectral and total irradiances derived from average over years 1834-1867 (solar"+nl+\
       "cycles 8-10) using annual average data (spectral_irradiance_Lean_1610-2140_ann_c100405.nc)"

  fileattdef( ncf, globalAtt )
 
 tsi_ref =  ncf_ann->tsi_ref
 ssi_ref =  ncf_ann->ssi_ref
 ncf->ssi=ssi
 ncf->tsi_ref = tsi_ref  
 ncf->ref_time_bound = ncf_ann->ref_time_bound
 ncf->date=date
 ncf->time=time
 ncf->wavelength=wavelength
 ncf->band_width=band_width
 ncf->ssi_ref = ssi_ref  
 ncf->tsi=tsi

 delete(ncf)
 
