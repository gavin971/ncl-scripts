
oldfile = "/glade/p/cesm/cseg//inputdata/atm/cam/chem/trop_mozart_aero/emis/ar5_mam3_so2_elev_1850-2005_c090804.nc"
newfile = "/glade/p/cesm/wawg/mmills/inputdata/ar5_mam3_so2_stratvolc_1850-2100_c130329.nc"

oldfileh = addfile(oldfile,"r")
lat = oldfileh->lat
lon = oldfileh->lon

date = (/ 18490115, 19910614, 19910615, 19910616, 21001215 /)
date!0="time"

altitude= (/ 8.0, 17.0, 18.5 /)
altitude@long_name = "altitude midlevel"
altitude@units = "km"
altitude!0="altitude"

altitude_int= (/ 0.0, 16.0, 18.0, 19.0 /)
altitude_int@long_name = "altitude interval"
altitude_int@units = "km"
altitude_int!0="altitude_int"

ntimes=dimsizes(date)
nlats=dimsizes(lat)
nlons=dimsizes(lon)
nalts=dimsizes(altitude)

stratvolc = new( (/ ntimes, nalts, nlats, nlons /), float)
stratvolc(:,:,:,:) = 0.0
stratvolc!0="time"
stratvolc!1="altitude"
stratvolc!2="lat"
stratvolc!3="lon"
stratvolc&altitude=altitude
stratvolc&lat=lat
stratvolc&lon=lon
stratvolc@long_name = "SO2 elevated emissions from explosive volcanoes" 
stratvolc@units = "molecules/cm3/s"


stratvolc(2,1,{15.1417},{120.35})=1.92e7

print("rm "+newfile)
a=systemfunc("rm "+newfile)
print("creating "+newfile)
newfileh = addfile(newfile,"c")
;Create UNLIMITED Record Dimension
;  An UNLIMITED record dimension in a netCDF file is critical if the user 
;  plans to ever use the netCDF operators (NCO) to concatenate the file 
;  along the time/record dimension.  
filedimdef(newfileh,"time",-1,True)

; For a nicer looking netCDF, create a "new line" character.
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
      "Stratospheric volcano SO2 emissions" 
  globalAtt@history = nl+\
      "Created by SO2elevPinatubo.ncl"
  globalAtt@notes = nl+\
      "This file includes emissions for the following eruptions:"+nl+\ 
      "  Pinatubo: June 15, 1991, 20 Tg, 16-18km, 1.92e10 molec/cm3/s"
  fileattdef( newfileh, globalAtt )

newfileh->date = date
newfileh->lat = lat
newfileh->lon = lon
newfileh->altitude = altitude
newfileh->altitude_int = altitude_int
newfileh->stratvolc = stratvolc






