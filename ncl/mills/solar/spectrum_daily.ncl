; Program converts Judith Lean's daily solar ascii data to netCDF
; Source data covers 1950-2009, in 2 files

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

ascii_path=getenv ("MODEL_PATH")+"/ccsm/solar/daily/"
cdfpath=getenv ("MODEL_PATH")+"/ccsm/inputdata/atm/cam/solar/"

file1=ascii_path+"spectra_1950_1999d_11yr_new.txt"
;file2=ascii_path+"katja_spectra_2000_2011d_28Mar11.txt"
;file2=ascii_path+"katja_spectra_2000_2011d_13Apr12.txt"
file2=ascii_path+"katja_spectra_2000_2012d_13Feb13.txt"
outfile="spectral_irradiance_Lean_1950-2012_daily"
creation_date="130227"
calcNoLeap = False

num_ssi_samples = 5*756

nyears1=1999-1950+1
nleap1=12
num_samples_time1 = nyears1*365+nleap1

print ("---------------------------------------------------------") 
print ("nyears1="+nyears1)
print ("nleap1="+nleap1)
print ("num_samples_time1="+num_samples_time1)

nyears2=2012-2000+1
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

; Read in ASCII files
; First 4 entries for each time are year, month, day, and TSI.
; These are followed by SSI for each of num_ssi_samples wavelengths.

solar_data1=asciiread(file1,\
       (/num_samples_time1, num_ssi_samples +4/),"double")
       
solar_data2=asciiread(file2,\
       (/num_samples_time2, num_ssi_samples +4/),"double")
       
solar_data3=new((/num_samples_time3,num_ssi_samples +4/),double)
solar_data3(0:num_samples_time1-1,:)=solar_data1
solar_data3(num_samples_time1:num_samples_time3-1,:)=solar_data2

solar_data=new((/num_samples_time3,num_ssi_samples +4/),double)

; Screen out missing values (TSI=-999)
j=0
nmiss=0
do i = 0, num_samples_time3-1
  if (solar_data3(i,3).lt.0.0) then
    print("Excluding missing data: "+solar_data3(i,0)+" "+solar_data3(i,1)+" "+\
          solar_data3(i,2)+" "+solar_data3(i,3))
    nmiss=nmiss+1
  else      
    solar_data(j,:)=solar_data3(i,:)
    if (solar_data(j,2).lt.1) then
      print(i+" "+j+" "+solar_data(j,0:2))
    end if
    j=j+1
  end if
end do

delete(solar_data1)
delete(solar_data2)
delete(solar_data3)

print("j="+j+", num_samples_time3="+num_samples_time3)
nmiss2=num_samples_time3-j
print(nmiss+" missing values")
print("num_samples_time3-j="+nmiss2)
;nmiss=nmiss-365-5
;print(nmiss+" missing values through last data point")

year=doubletoint(solar_data(0:j-1,0))
month=doubletoint(solar_data(0:j-1,1))
day=doubletoint(solar_data(0:j-1,2))
hour=new(j,integer)
hour=0

date=year*10000+month*100+day
  date@format="YYYYMMDD"
  date!0="time"

reftime=greg2jul(0,1,1,0)

time = greg2jul(year,month,day,hour)-reftime
  time@_FillValue = -9999.d0
  time@units = "days since 0000-01-01 00:00:00"
  time@time_origin = "01-JAN-0000"
  time@axis = "T"
  time@calendar = "Gregorian"
  time!0="time"

return
;print(year+" - "+month+" - "+day+" - "+hour)

 print(date(0)+" "+time(0)+" "+solar_data(0,0:2))
 print(date(365)+" "+time(365)+" "+solar_data(365,0:2))
 print(date(j-1)+" "+time(j-1)+" "+solar_data(j-1,0:2))

tsi = solar_data(0:j-1,3)
  tsi@_FillValue = -9999.d0
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"

ssi = solar_data(0:j-1,4:)
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
 print ("verify last date: "+date(j-1))
 print ("verify last tsi: "+sprintf("%11.6f",tsi(j-1)))
 print ("verify last ssi at last time: "+ssi(j-1,num_ssi_samples-1))
 print ("---------------------------------------------------------") 
  
 tsi = tsi * pmod_tim_scale
 ssi = ssi * pmod_tim_scale
 
 filename=cdfpath+outfile+"_Leap_c"+creation_date+".nc"
 print("rm "+filename)
 a = systemfunc("rm "+filename)
 print("creating file "+filename)
 ncf = addfile(filename,"c")

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
      "Lean daily solar irradiance for "+date(0)+"-"+date(j-1)+"; "+nmiss+" missing values" ;
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
       "Reference spectral and total irradiances derived from average over years 1834-1867 (solar"+nl+\
       "cycles 8-10) using annual average data (spectral_irradiance_Lean_1610-2140_ann_c100405.nc)"

  fileattdef( ncf, globalAtt )
 
 tsi_ref =  ncf_ann->tsi_ref
 ssi_ref =  ncf_ann->ssi_ref
 print("writing ssi")
 ncf->ssi=ssi
 ncf->tsi_ref = tsi_ref  
 ncf->ref_time_bound = ncf_ann->ref_time_bound
 ncf->date=date
; ncf->time=time
; ncf->wavelength=wavelength
 ncf->band_width=band_width
 ncf->ssi_ref = ssi_ref  
 print("writing tsi")
 ncf->tsi=tsi

 print("closing file "+filename)
 delete(ncf)
 delete(tsi)
 delete(ssi)
 delete(date)
 delete(time)
 
 if (calcNoLeap) then
   print ("---------------------------------------------------------") 
   print("No Leap calculation...")
   
   ; remove leap day from data set
   num_samples_time_noleap=nyears*365
   print ("num_samples_time_noleap="+num_samples_time_noleap)
   print ("---------------------------------------------------------") 
   solar_data_noleap=new((/num_samples_time_noleap,num_ssi_samples +4/),double)

  ; Screen out leap days (Feb 29)
   k=0
   nmiss=0
   do i = 0, j-1
     if ((solar_data(i,1).eq.2).and.(solar_data(i,2).eq.29)) then
       print("Excluding leap day: "+solar_data(i,0)+" "+solar_data(i,1)+" "+solar_data(i,2))
     else
       solar_data_noleap(k,:)=solar_data(i,:)
       k=k+1
     end if
   end do

   delete(solar_data)

  time = solar_data_noleap(0:k-1,0)*365.d0 + \
         day_of_year(conform_dims(k,1999,-1),\
         doubletoint(solar_data_noleap(0:k-1,1)),doubletoint(solar_data_noleap(0:k-1,2)))
    time@_FillValue = -9999.d0
    time@units = "days since 0000-01-01 00:00:00"
    time@time_origin = "01-JAN-0000"
    time@axis = "T"
    time@calendar = "noleap"
    time!0="time"

  date=doubletoint(solar_data_noleap(0:k-1,0)*10000+solar_data_noleap(0:k-1,1)*100+solar_data_noleap(0:k-1,2))
    date@format="YYYYMMDD"
    date!0="time"

   print(date(0)+" "+time(0)+" "+solar_data_noleap(0,0:2))
   print(date(365)+" "+time(365)+" "+solar_data_noleap(365,0:2))
   print(date(k-1)+" "+time(k-1)+" "+solar_data_noleap(k-1,0:2))

  tsi = solar_data_noleap(0:k-1,3)
    tsi@_FillValue = -9999.d0
    tsi@long_name = "Total Solar Irradiance at 1 a.u."
    tsi@units = "W/m^2" ;
    tsi!0="time"

  ssi = solar_data_noleap(0:k-1,4:)
    ssi@_FillValue = -9999.d0
    ssi@long_name = "Solar Spectral Irradiance at 1 a.u."
    ssi@units = "mW/m^2/nm"
    ssi!0="time"
    ssi!1="wavelength"

   tsi = tsi * pmod_tim_scale
   ssi = ssi * pmod_tim_scale

   filename=cdfpath+outfile+"_NoLeap_c"+creation_date+".nc"
   print("creating file"+filename)
   a = systemfunc("rm "+filename)
   ncf = addfile(filename,"c")

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
        "Lean daily solar irradiance for "+date(0)+"-"+date(k-1)+"; "+nmiss+" missing values" ;
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

   print("writing ssi")
   ncf->ssi=ssi
   ncf->tsi_ref = tsi_ref  
   ncf->ref_time_bound = ncf_ann->ref_time_bound
   ncf->date=date
   ncf->time=time
   ncf->wavelength=wavelength
   ncf->band_width=band_width
   ncf->ssi_ref = ssi_ref  
   print("writing tsi")
   ncf->tsi=tsi

   delete(ncf)
 end if
