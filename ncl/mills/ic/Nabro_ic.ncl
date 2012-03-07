; Create initial conditions (IC) for dust-sulfur run, 125 levels
;
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; Constants
;
; NOTE: These are from CAM, shr_const_mod.F90.
PI      = 4.0 * atan(1.0)         ; pi
G       = 9.80616                 ; acceleration of gravity ~ m/s^2
BOLTZ   = 1.38065e-23             ; Boltzmann's constant ~ J/K/molecule
AVOGAD  = 6.02214e26              ; Avogadro's number ~ molecules/kmole
RGAS    = AVOGAD*BOLTZ            ; Universal gas constant ~ J/K/kmole
MWDAIR  = 28.966                  ; molecular weight dry air ~ kg/kmole
RDAIR   = RGAS/MWDAIR             ; Dry air gas constant ~ J/K/kg
MWH2O   = 18.000                  ; molecular weight of water vapor
PTOP    = 100.0                   ; pressure at the top of the model (non-zero) ~ Pa
REARTH  = 6.37122e6               ; Earth's radius (m)
DG2RAD  = PI / 180.0              ; rad deg^-1

maxpress=250.0
minpress=130.0
prstr=floattoint(minpress)+"-"+floattoint(maxpress)+"hPa"

host   = systemfunc("hostname")
if (hostname.eq.""modeling1") then
  filePath = "/data4/mmills/f40.2008-2011.2deg.wcm.sd.carma.sulf.002/rest/"
else
  if (hostname.eq.""lou") then
    filePath = "~/csm/f40.2008-2011.2deg.wcm.sd.carma.sulf.002/rest/2011-06-13-00000/"
  else
    filePath = "$MODEL_PATH/cesm/run/f40.2008-2011.2deg.wcm.sd.carma.sulf.002/rest"
  endif
endif

icFileName="f40.2008-2011.2deg.wcm.sd.carma.sulf.002.cam2.i.2011-06-13-00000.nc"
srcIcPath = filePath+icFileName

; Copy it to create new IC file:
icPath = filePath+"f40.2008-2011.2deg.wcm.sd.carma.sulf.002.Nabro"+prstr+".cam2.i.2011-06-13-00000.nc"
print("copying IC file "+srcIcPath+" to "+icPath)
system("cp " + srcIcPath + " " + icPath)

icFile = addfile(icPath,"w")
lev = icFile->lev
nLevs = dimsizes(lev)
lat = icFile->lat
lon = icFile->lon
nLons = dimsizes(lon)
time = icFile->time

; Nabro coordinates: 13°22'N, 41°42'E
ilat=55
jlon=17 
dlon = fabs(lon(jlon) - lon(jlon-1))
dlat = fabs(lat(ilat) - lat(ilat-1))
print("latitude:"+lat(ilat)+" longitude:"+lon(jlon)+" dlat="+dlat+" dlon="+dlon)
  
dlat2 = dlat / 2.
sdlat = sin((lat(ilat) + dlat2) * DG2RAD) - sin((lat(ilat) - dlat2) * DG2RAD)
dh = REARTH * sdlat
dx = ((2. * PI) / nLons) * REARTH
dhdx = dh * dx

print("  Reading SO2 & Q ...")
SO2 = icFile->SO2 ; [time | 1] x [lev | 88] x [lat | 96] x [lon | 144]
q   = icFile->Q

print("")
print("  Calculating air pressure and density ...")
hyai 	= icFile->hyai
hybi 	= icFile->hybi
P0 	= icFile->P0
PS 	= icFile->PS

; Calculate the air pressures at the interfaces.
P = pres_hybrid_ccm(PS(:,:,:), P0, hyai(:) ,hybi(:))
delete(PS)
pdel = P(:, 1:nLevs, :, :) - P(:, 0:nLevs-1, :, :)
pdel!1="lev"
pdel&lev=lev
delete(P)

airMass = (q(0, {minpress:maxpress}, ilat, jlon)+1.0) * (pdel(0, {minpress:maxpress}, ilat, jlon)/ G)
columnAirMass = dim_sum_Wrap(airMass)
print("columnAirMass ="+ columnAirMass+" kg/m2")
totalAirMass = columnAirMass * dhdx /1e9
print("totalAirMass ="+ totalAirMass+" Tg")

SO2_Nabro_mmr=1.3/totalAirMass
print("SO2_Nabro_mmr ="+SO2_Nabro_mmr)

SO2_Nabro = SO2
SO2_Nabro(0, {minpress:maxpress}, ilat, jlon) = SO2(0, {minpress:maxpress}, ilat, jlon) + SO2_Nabro_mmr

SO2Mass = SO2(0, {minpress:maxpress}, ilat, jlon) * (q(0, {minpress:maxpress}, ilat, jlon)+1.0) * (pdel(0, {minpress:maxpress}, ilat, jlon)/ G)
columnSO2Mass = dim_sum_Wrap(SO2Mass)
SO2MassNabro = SO2_Nabro(0, {minpress:maxpress}, ilat, jlon) * (q(0, {minpress:maxpress}, ilat, jlon)+1.0) * (pdel(0, {minpress:maxpress}, ilat, jlon)/ G)
columnSO2MassNabro = dim_sum_Wrap(SO2MassNabro)
totalSO2Mass = (columnSO2MassNabro - columnSO2Mass) * dhdx /1e9
print("totalSO2Mass ="+ totalSO2Mass+" Tg")

icFile->SO2 = SO2_Nabro
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
      "IC file for Nabro eruption study" ;
  globalAtt@history = nl+\
       "Adapted from "+icFileName+nl+\
       "Nabro SO2 added by program Nabro_ic.ncl"
  globalAtt@notes = nl+\
       "1.3 Tg of SO2 at "+lat(ilat)+"N, "+lon(jlon)+"E, "+prstr+" on June 13, 2011"+nl+\ 
       "increasing SO2 mmr by "+SO2_Nabro_mmr+" over eruption site" ;

  fileattdef( icFile, globalAtt )
