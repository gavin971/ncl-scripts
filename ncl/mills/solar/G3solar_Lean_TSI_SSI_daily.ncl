; Program produces solar variability input file for GeoMIP experiment G3solar
; by subtracting the increase in anthropgenic radiative forcing from 2020 to 2069. 

tsifile="$DIN_LOC_ROOT/atm/cam/solar/spectral_irradiance_Lean_1950-2140_daily_c100804.nc"
g3solarfile="$DIN_LOC_ROOT/atm/cam/solar/spectral_irradiance_Lean_1950-2140_daily_G3solar_c111010.nc"
system("cp " + tsifile + " " + g3solarfile)
ncf    = addfile(g3solarfile,"w")

tsi=ncf->tsi
ssi=ncf->ssi
date=ncf->date

year=date/10000
newtsi=tsi

trf=(/2.4838861,2.5250856,2.5659308,2.6071713,2.6484638,2.6898474,2.7313286,2.7728545,2.8143456,2.8557683,\
      2.8974304,2.9396498,2.9813833,3.0220663,3.0622462,3.1020049,3.1413681,3.1803377,3.2188563,3.2568433,\
      3.2945768,3.3323783,3.3693779,3.4050626,3.4399375,3.4742503,3.5081648,3.5417686,3.5751082,3.6081943,\
      3.6399508,3.6693022,3.6968597,3.7233185,3.7486985,3.7730437,3.7965619,3.8193898,3.8415987,3.8632337,\
      3.8840183,3.9036409,3.9221476,3.9396312,3.9561294,3.9715814,3.9860089,3.9995963,4.0124625,4.0246656/)

do iyear = 2020, 2069
  tsidiff=(trf(iyear-2020)-trf(0))*4.0/(1.0-0.32)
  tsidiff=tsidiff*55.7965/41.2752 ; Scale up the reduction in RF based on 4xCO2 solar run (6/2/2011)
  print(iyear+" "+trf(iyear-2020)+" "+tsidiff)
;  print(year)
;  print(ind(year.eq.iyear))
  
  newtsi=where(year.eq.iyear, tsi-tsidiff, newtsi)
  
  do iwave=0,3779
    ssi(ind(year.eq.iyear),iwave)=ssi(ind(year.eq.iyear),iwave)*newtsi(ind(year.eq.iyear))/tsi(ind(year.eq.iyear))
  end do
  
end do

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
  globalAtt@user = nl+\
      "CESM1.0 user"
  globalAtt@title = nl+\
      "Lean daily spectral solar irradiance for years 1950-2140, adjusted for GeoMIP experiment G3solar" ;
  globalAtt@history = nl+\
      "Extended from spectral_irradiance_Lean_1950-2009_daily_c100507.nc by program future_Lean_daily.ncl."+nl+\
      "Converted from ASCII by program spectrum_daily.ncl"+nl+\
      "Converted from ASCII by program G3solar_Lean_TSI_SSI_daily.ncl"
  globalAtt@notes = nl+\
      "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
      "Total and spectral irradiance variations are determined independently of each other,"+nl+\
      "with more uncertainty in the spectral irradiance changes which are determined in 1 nm bins."+nl+\
      "Integrated SSI is therefore varies slightly from TSI."+nl+\
      "Leap days (Feb 29) have been removed from source to conform to 365-day model year."+nl+\
      "Reference spectral and total irradiances derived from average over years 1834-1867 (solar\n"+nl+\
      "cycles 8-10) using annual average data (spectral_irradiance_Lean_1610-2140_ann_c100405.nc)."+nl+\
      "Irradiances for 2009-2140 created by repeating the last 4 cycles (cycle 20-23)"+nl+\
      "with values from 1965 to 2008 inclusive mapping to 2009-2052, 2053-2096, 2097-2140"+nl+\
      "No background trends are included in TSI or SSI from 2008-2140."+nl+\
      "TSI from 2020-2140 is adjusted for the GeoMIP experiment G3solar,"+nl+\
      "dimming TSI by the change in total RF * 4.0 / (1-albedo), assuming the RCP4.5,"+nl+\
      "emission scenario, and albedo = 0.32."
  globalAtt@scaling = nl+\
       "Absolute spectral and total irradiance includes PMOD adjustment using 0.9965"+nl+\
       "for TIM scale with background from Wang, Lean & Sheeley, Modeling the sun's"+nl+\
       "magnetic field and irradiance since 1713, Astrophys J (2005) vol. 625 (1) pp. 522-538."+nl+\
       "Solar dimming for 2020-2140 is scaled up by 55.7965/41.2752, based on 4xCO2 solar run."

  fileattdef( ncf, globalAtt )
  
 ncf->tsi=newtsi
 ncf->ssi=ssi
