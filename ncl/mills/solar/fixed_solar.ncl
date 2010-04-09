; Program creates a "fixed solar" netcdf file that has two time-stamps:
; 101 and 25000101.  ssi and tsi at those time points, and the ssi_ref,
; tsi_ref will be the cycle 8, 9, 10 averages (1834-1867).  This is the 
; fixed solar file needed for the new 1850 spinup.

; by Mike Mills
; April 2, 2010

ncf = addfile("spectral_irradiance_Lean_1610-2009_ann_c100405.nc","r")

ssi_ref=ncf->ssi_ref
tsi_ref=ncf->tsi_ref
ref_time_bound=ncf->ref_time_bound
band_width=ncf->band_width

date = (/101, 25000101/)
  date@format="YYYYMMDD"
  date!0="time"

time = new(2, double)
time(0) = 0.d0
time(1) = 2500.d0 * 365.d0
  time@_FillValue = -9999.d0
  time@units = "days since 0000-01-01 00:00:00"
  time@time_origin = "01-JAN-0000"
  time@axis = "T"
  time@calendar = "noleap"
  time!0="time"

tsi=new(2,double)
tsi(0)=ncf->tsi_ref
tsi(1)=ncf->tsi_ref
  tsi@_FillValue = -9999.d0
  tsi@long_name = "Total Solar Irradiance at 1 a.u."
  tsi@units = "W/m^2" ;
  tsi!0="time"
  
ssi=new((/2,dimsizes(ssi_ref)/),double)
ssi(0,:)=ncf->ssi_ref
ssi(1,:)=ncf->ssi_ref
  ssi@_FillValue = -9999.d0
  ssi@long_name = "Solar Spectral Irradiance at 1 a.u."
  ssi@units = "mW/m^2/nm"
  ssi!0="time"
  ssi!1="wavelength"


fixed = addfile("spectral_irradiance_Lean_1850_cntl_c100407.nc","c")

; Create an UNLIMITED record dimension in the output netCDF file.  This is critical if 
;  the user plans to ever use NCO to concatenate the file along the time/record dimension.
 filedimdef(fixed,"time",-1,True)
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
      "Lean annual spectral solar irradiance for fixed 1850 control" ;
  globalAtt@Source_website = nl+\
       "http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/CMIP5_solar_irradiance.html"
  globalAtt@Source_data=nl+\
       "ftp://strat50.met.fu-berlin.de/pub/outgoing/_matthes/CMIP5_solardata/spectra_1610_2000a_21Jan09.txt.gz"
  globalAtt@scaling = nl+\
       "Absolute spectral and total irradiance includes PMOD adjustment using 0.9965 for TIM scale with background from"+nl+\ 
       "Wang, Lean & Sheeley, Modeling the sun's magnetic field and irradiance since 1713, Astrophys J (2005) vol. 625 (1) pp. 522-538."
  globalAtt@notes = nl+\
       "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
       "Total and spectral irradiance variations are determined independently of each other,"+nl+\
       "with more uncertainty in the spectral irradiance changes which are determined in 1 nm bins."+nl+\
       "Integrated SSI is therefore varies slightly from TSI."+nl+\
       "Reference spectral and total irradiances derived from average over years 1834-1867 (solar cycles 8-10)"      
  globalAtt@history = nl+\
       "Adapted from spectra_1610_2009a_31Mar10.nc by program fixed_solar.ncl"

  fileattdef( fixed, globalAtt )

fixed->date=date
fixed->time=time
fixed->band_width=band_width
fixed->ref_time_bound=ref_time_bound
fixed->tsi_ref=tsi_ref
fixed->ssi_ref=ssi_ref
fixed->tsi=tsi
fixed->ssi=ssi
