; Program produces solar variability input file for GeoMIP experiment G3solar
; by subtracting the increase in anthropgenic radiative forcing from 2020 to 2069. 

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

tsifile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_c100301.nc"
;g3solarfile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_G3solar_c110222.nc"
g3solarfile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_G3solar_v2_c110602.nc"
system("cp " + tsifile + " " + g3solarfile)
ncf    = addfile(g3solarfile,"w")
;ncghgf = addfile("$DIN_LOC_ROOT/atm/cam/ggas/ghg_rcp45_1765-2500_c100405.nc","r")

tsi=ncf->tsi

trf=(/2.4838861,2.5250856,2.5659308,2.6071713,2.6484638,2.6898474,2.7313286,2.7728545,2.8143456,2.8557683,\
      2.8974304,2.9396498,2.9813833,3.0220663,3.0622462,3.1020049,3.1413681,3.1803377,3.2188563,3.2568433,\
      3.2945768,3.3323783,3.3693779,3.4050626,3.4399375,3.4742503,3.5081648,3.5417686,3.5751082,3.6081943,\
      3.6399508,3.6693022,3.6968597,3.7233185,3.7486985,3.7730437,3.7965619,3.8193898,3.8415987,3.8632337,\
      3.8840183,3.9036409,3.9221476,3.9396312,3.9561294,3.9715814,3.9860089,3.9995963,4.0124625,4.0246656/)

do year = 2020, 2069
  tsidiff=(trf(year-2020)-trf(0))*4.0/(1.0-0.32)
  tsidiff=tsidiff*55.7965/41.2752 ; Scale up the reduction in RF based on 4xCO2 solar run (6/2/2011)
  newtsi=tsi(year-1610)-tsidiff
  print(year+" "+trf(year-2020)+" "+tsidiff+" "+tsi(year-1610)+" "+newtsi)
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
  globalAtt@user = nl+\
      "CESM1.0 user"
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
       "dimming TSI by the change in total RF * 4.0 / (1-albedo), assuming the RCP4.5,"+nl+\
       "emission scenario, and albedo = 0.32." ;

  fileattdef( ncf, globalAtt )
  
 ncf->tsi=tsi
