; Script written by Mike Mills
; Modified October 11, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

z1=16 ; lower bound of SO2 emission (km)
z2=20 ; upper bound of SO2 emission (km), must be less than 31 km
PinFac = 0.5 ; emission relative to Pinatubo
TgSO2 = 20.*PinFac ; Tg of SO2

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")
datestr = systemfunc("date +%y%m%d")

if (isStrSubset(host,"caldera").or.isStrSubset(host,"geyser").or.isStrSubset(host,"yslogin")) then
  oldpath = "/glade/p/cesm/cseg//inputdata/atm/cam/chem/trop_mozart_aero/emis/"
  newpath = "/glade/p/cesm/wawg/mmills/inputdata/"
else
  oldpath = "/Volumes/Data/Models/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/"
  newpath = oldpath
end if

oldfile = oldpath+"ar5_mam3_so2_elev_1850-2005_c090804.nc"
newfile = newpath+"ar5_mam3_so2_stratvolc_1850-2100_"+floattointeger(TgSO2)+"Tg_"\
                 +z1+"-"+z2+"_c"+datestr+".nc"

print("reading "+oldfile)
oldfileh = addfile(oldfile,"r")
lat = oldfileh->lat
lon = oldfileh->lon

depth = (z2-z1)*1.0e5 ; emission depth in cm
volume = 5.66e14 * depth ; cm3
rate = PinFac * 1.88e35 / 86400. / volume ; molec/cm3/s
print("rate="+rate)

date = (/ 18490115, 19910614, 19910615, 19910616, 21001215 /)
date!0="time"

altitude= (/ z1/2., z1+(z2-z1)/2.0, z2+(31.0-z2)/2.0 /)
altitude@long_name = "altitude midlevel"
altitude@units = "km"
altitude!0="altitude"

altitude_int= (/ 0.0, z1, z2, 31.0 /)
altitude_int@long_name = "altitude interval"
altitude_int@units = "km"
altitude_int!0="altitude_int"

ntimes=dimsizes(date)
nlats=dimsizes(lat)
nlons=dimsizes(lon)
nalts=dimsizes(altitude)

stratvolc = new( (/ ntimes, nalts, nlats, nlons /), double)
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

stratvolc(2,1,{15.1417},{120.35})=rate

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
      "  "+PinFac+" x Pinatubo (15.1417N, 120.35W): June 15, 1991, "+TgSO2+" Tg SO2, "+z1+"-"+z2+"km, "+rate+" molec/cm3/s"
  fileattdef( newfileh, globalAtt )

newfileh->date = date
newfileh->lat = lat
newfileh->lon = lon
newfileh->altitude = altitude 
newfileh->altitude_int = altitude_int
newfileh->stratvolc = stratvolc






