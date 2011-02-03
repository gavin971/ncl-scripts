; Program extends solar variability input file from 1610-2009 to 1610-2140
; by repeating the 4 solar cycles (44 years) prior to 2008, 3 times for
; 2009-2140

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

tsifile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_c100301.nc"
g3solarfile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_G3solar_c110203.nc"
system("cp " + tsifile + " " + g3solarfile)
ncf    = addfile(g3solarfile,"w")
ncghgf = addfile("$DIN_LOC_ROOT/atm/cam/ggas/ghg_rcp45_1765-2500_c100405.nc","r")

tsi=ncf->tsi
co2=ncghgf->CO2

do year = 2020, 2140
  co2diff=co2(year-1765)-co2(2020-1765)
  tsidiff=co2diff*0.0655661
  newtsi=tsi(year-1610)-tsidiff
  print(year+" "+co2(year-1765)+" "+co2diff+" "+tsi(year-1610)+" "+newtsi+" "+tsidiff)
  tsi(year-1610)=newtsi
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
  globalAtt@title = nl+\
      "Lean total solar irradiance for years 1610-2140, adjusted for GeoMIP experiment G3solar" ;
  globalAtt@Source_data=nl+\
       "TSI_WLS_ann_1610_2008.txt from http://www.geo.fu-berlin.de/en/met/ag/strat/forschung/SOLARIS/Input_data/CMIP5_solar_irradiance.html"
  globalAtt@history = nl+\
       "Adapted from SOLAR_TSI_Lean_1610-2140_annual_c100301.nc"+nl+\
       "Converted from ASCII by program G3soloar_Lean_TSI_ann.ncl"
  globalAtt@notes = nl+\
       "Data derived by Judith Lean, jlean@ssd5.nrl.navy.mil"+nl+\
       "Irradiances for 2009-2140 created by repeating the last 4 cycles (cycle 20-23)"+nl+\
       "with values from 1965 to 2008 inclusive mapping to 2009-2052, 2053-2096, 2097-2140"+nl+\
       "TSI from 2020-2140 is adjusted for the GeoMIP experiment G3solar,"+nl+\
       "dimming TSI by 0.0655661 W/m^2 per ppmv CO2 increase over the 2020 mixing ratio,"+nl+\
       "assuming the RCP4.5 emission scenario. CO2 values taken from the file:"+nl+\
       "ghg_rcp45_1765-2500_c100405.nc" ;

  fileattdef( ncf, globalAtt )
  
 ncf->tsi=tsi
