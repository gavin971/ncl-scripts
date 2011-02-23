; Program extends solar variability input file from 1610-2009 to 1610-2140
; by repeating the 4 solar cycles (44 years) prior to 2008, 3 times for
; 2009-2140

;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

tsifile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_c100301.nc"
g3solarfile="$DIN_LOC_ROOT/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_G3solar_c110222.nc"
system("cp " + tsifile + " " + g3solarfile)
ncf    = addfile(g3solarfile,"w")
;ncghgf = addfile("$DIN_LOC_ROOT/atm/cam/ggas/ghg_rcp45_1765-2500_c100405.nc","r")

tsi=ncf->tsi
;co2=ncghgf->CO2
trf=(/2.5317551,2.6216119,2.6855166,2.7329421,2.7634036,2.7757707,2.7881402,2.8119714,2.8428423,2.8770636,\
      2.9143627,2.9875188,3.0779096,3.1416521,3.188017, 3.2169447,3.2272914,3.2371493,3.2579732,3.28534,  \
      3.3158721,3.3493106,3.4172469,3.5015888,3.5595233,3.6000211,3.6231047,3.627692, 3.6319198,3.6473113,\
      3.6684475,3.6905975,3.713792, 3.7711875,3.8452247,3.8926295,3.9223327,3.9343296,3.927522, 3.9200453,\
      3.9231352,3.9321376,3.9434429,3.9565635,4.0039984,4.0681076,4.1055948,4.1253671,4.1274023,4.1105889,\
      4.0927972,4.0852475,4.0835997,4.084251, 4.0866409,4.1232564,4.1764501,4.2029154,4.2115701,4.2024557,\
      4.1735135,4.1418535,4.1214889,4.1097216,4.1019818,4.0976103,4.1290055,4.1784445,4.2025543,4.2101752,\
      4.2012273,4.1745478,4.1482203,4.1336824,4.1264816,4.1228955,4.1223161,4.1571893,4.2098327,4.2369026,\
      4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,\
      4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,\
      4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,\
      4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,4.2451996,\
      4.2451996/)

trf=(/2.4838861,2.5250856,2.5659308,2.6071713,2.6484638,2.6898474,2.7313286,2.7728545,2.8143456,2.8557683,,\
      2.8974304,2.9396498,2.9813833,3.0220663,3.0622462,3.1020049,3.1413681,3.1803377,3.2188563,3.2568433,,\
      3.2945768,3.3323783,3.3693779,3.4050626,3.4399375,3.4742503,3.5081648,3.5417686,3.5751082,3.6081943,,\
      3.6399508,3.6693022,3.6968597,3.7233185,3.7486985,3.7730437,3.7965619,3.8193898,3.8415987,3.8632337,,\
      3.8840183,3.9036409,3.9221476,3.9396312,3.9561294,3.9715814,3.9860089,3.9995963,4.0124625,4.0246656/)

do year = 2020, 2069
  tsidiff=(trf(year-2020)-trf(0))*4.0/(1.0-0.32)
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
