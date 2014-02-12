; Script written by Mike Mills
; Modified October 11, 2013

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

z1=19.d0 ; lower bound of SO2 emission (km)
z2=27.d0 ; upper bound of SO2 emission (km), must be less than 31 km

TgSO2 = 20.d0 ; Tg of SO2

; Zonal mean in tropics:
;lat1=-10.d0
;lat2=20.d0
;lon1=0.0d0
;lon2=357.5d0

; Spread over 10x10 columns:
;lat1=0.5d0
;lat2=18.5d0
;lon1=96.5d0
;lon2=118.5d0

; Single column over Pinatubo:
;lat1=15.1417d0
;lat2=15.1417d0
;lon1=120.35d0
;lon2=120.35d0

; Single column over Pinatubo:
lat1=0.0d0
lat2=20.0d0
lon1=120.35d0
lon2=120.35d0

if (lat1.lt.0) then
  latstr1=abs(doubletointeger(lat1))+"S"
else
  latstr1=doubletointeger(lat1)+"N"
end if

if (lat2.lt.0) then
  latstr2=abs(doubletointeger(lat2))+"S"
else
  latstr2=doubletointeger(lat2)+"N"
end if

if ((lon1.le.2.5).and.(lon2.ge.357.5)) then
  lonstr="zavg"
else
  lonstr=doubletointeger(lon1)+"-"+doubletointeger(lon2)+"E"
end if

molecSO2 = TgSO2 * 1.d12 / 64.d0 * 6.022d23 ; total molecules of SO2
print("molecSO2="+molecSO2)

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")
datestr = systemfunc("date +%y%m%d")

if (isStrSubset(host,"caldera").or.isStrSubset(host,"geyser").or.isStrSubset(host,"yslogin")) then
  newpath  = "/glade/p/cesm/wawg/mmills/inputdata/"
  gridpath = newpath
else
  newpath  = "/Volumes/Data/Models/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/"
  gridpath = newpath
end if

newfile  = newpath+"ar5_so2_stratvolc_1850-2100_"+doubletointeger(TgSO2)+"Tg_"\
                 +doubletointeger(z1)+"-"+doubletointeger(z2)+"km_"\
                 +latstr1+"-"+latstr2+"_"+lonstr+"_c"+datestr+".nc"
gridfile = gridpath+"f19.grid.nc"

print("reading "+gridfile)
gridfileh = addfile(gridfile,"r")
;slat = gridfileh->slat
;slon = gridfileh->slon
lat  = gridfileh->lat
lon  = gridfileh->lon
gw   = gridfileh->gw
delete(gridfileh)
nlats=dimsizes(lat)
nlons=dimsizes(lon)

varea=new((/ nlats, nlons /), double)
varea!0="lat"
varea!1="lon"
varea&lat=lat
varea&lon=lon
varea=0.d0

tarea=varea
tarea=1.d0

weights=conform(varea,gw,0)

varea({lat1:lat2},{lon1:lon2})=1.d0
totarea=wgt_areasum2(tarea,weights,0)
volcarea=wgt_areasum2(varea,weights,0)

Rearth=6.37122d8 ; Earth radius in cm (CAM)
pi=asin(1.d0)*2.d0
SAearth=4.d0*pi*Rearth*Rearth ; Earth surface area in cm2 (CAM)
print("SAearth="+SAearth+" cm2")

SAvolc=SAearth*volcarea/totarea
print("SAvolc="+SAvolc+" cm2")
 
depth = (z2-z1)*1.0d5 ; emission depth in cm
print("depth="+depth+" cm")
volume = SAvolc * depth ; cm3
print("volume="+volume+" cm3")
rate = molecSO2 / 86400. / volume ; molec/cm3/s
print("rate="+rate+" molec/cm3/s")

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

;stratvolc(2,1,{15.1417},{120.35})=rate
stratvolc(2,1,{lat1:lat2},{lon1:lon2})=rate

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
      "  Pinatubo: June 15, 1991, "+TgSO2+" Tg SO2, "+rate+" molec/cm3/s"+nl+\
      "    altitudes:  "+z1+" to "+z2+" km"+nl+\
      "    latitudes:  "+lat1+" to "+lat2+" "+nl+\
      "    longitudes: "+lon1+" to "+lon2
  fileattdef( newfileh, globalAtt )

newfileh->date = date
newfileh->lat = lat
newfileh->lon = lon
newfileh->altitude = altitude 
newfileh->altitude_int = altitude_int
newfileh->stratvolc = stratvolc






