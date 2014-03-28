; Script written by Mike Mills
; Created March 22, 2014

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

z1=22.d0 ; lower bound of SO2 emission (km)
z2=23.d0 ; upper bound of SO2 emission (km), must be less than 31 km

TgSO2 = 10.d0 ; Tg of SO2/year

; Single column over Pinatubo:
lat1=0.9473648d0
lat2=0.9473648d0
lon1=120.d0
lon2=120.d0

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

if (lat1.eq.lat2) then
  latstr=latstr1
else
  latstr=latstr1+"-"+latstr2
end if

if ((lon1.le.2.5).and.(lon2.ge.357.5)) then
  lonstr="zavg"
else
  if (lon1.eq.lon2) then
    lonstr=doubletointeger(lon1)+"E"
  else
    lonstr=doubletointeger(lon1)+"-"+doubletointeger(lon2)+"E"
  end if
end if

molecSO2 = TgSO2 * 1.d12 / 64.d0 * 6.022d23 ; total molecules of SO2
print("molecSO2="+molecSO2)

host   = systemfunc("hostname")
ostype = systemfunc("echo $OSTYPE")
datestr = systemfunc("date +%y%m%d")

if (isStrSubset(host,"caldera").or.isStrSubset(host,"geyser").or.isStrSubset(host,"yslogin")) then
  gridpath = "/glade/p/cesm/wawg/mmills/inputdata/"
  newpath  = gridpath+"emis/"
else
  newpath  = "/Volumes/Data/Models/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/"
  gridpath = newpath
end if

newfile  = newpath+"ar5_so2_geoeng_1850-2100_"+doubletointeger(TgSO2)+"Tg_"\
                 +doubletointeger(z1)+"-"+doubletointeger(z2)+"km_"\
                 +latstr+"_"+lonstr+"_c"+datestr+".nc"
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
rate = molecSO2 / 86400. / 365. / volume ; molec/cm3/s
print("rate="+rate+" molec/cm3/s")

; include 12 month cycle for INTERP_MISSING_MONTHS structure
date = (/ 18490115, 18490215, 18490315, 18490415, 18490515, 18490615, \
          18490715, 18490815, 18490915, 18491015, 18491115, 18491215, \
          21000115, 21000215, 21000315, 21000415, 21000515, 21000615, \
          21000715, 21000815, 21000915, 21001015, 21001115, 21001215 /)
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

geoeng = new( (/ ntimes, nalts, nlats, nlons /), double)
geoeng(:,:,:,:) = 0.0
geoeng!0="time"
geoeng!1="altitude"
geoeng!2="lat"
geoeng!3="lon"
geoeng&altitude=altitude
geoeng&lat=lat
geoeng&lon=lon
geoeng@long_name = "SO2 elevated emissions from geoengineering" 
geoeng@units = "molecules/cm3/s"

;geoeng(2,1,{15.1417},{120.35})=rate
geoeng(:,1,{lat1:lat2},{lon1:lon2})=rate

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
      "Geoengineering SO2 emissions" 
  globalAtt@history = nl+\
      "Created by SO2elevGeoeng.ncl"
  globalAtt@notes = nl+\
      "This file is for the following geoengineering scenario:"+nl+\ 
      "   "+TgSO2+" Tg SO2/year, "+rate+" molec/cm3/s"+nl+\
      "    altitudes:  "+z1+" to "+z2+" km"+nl+\
      "    latitudes:  "+latstr+" "+nl+\
      "    longitudes: "+lonstr
  fileattdef( newfileh, globalAtt )

newfileh->date = date
newfileh->lat = lat
newfileh->lon = lon
newfileh->altitude = altitude 
newfileh->altitude_int = altitude_int
newfileh->geoeng = geoeng






