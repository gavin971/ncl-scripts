;load,"ncl-scripts/getHistoryPath.ncl"
load "ncl-scripts/common2.ncl"

AVG     = 6.02252e+23             ; Avogadro's number (CARMA) ~ molecules/mole
cpi = 4./3.*PI

;  Molecular weights [g/mol]
gwtH2SO4 = 98.078479

caseid = "125b36eg"

;historyPath=getHistoryPath(caseid)
historyPath = "scratch/"+caseid+"/"

print(historyPath)

files=(/caseid+".cam2.h3.1995-05-30-00000.nc", \
        caseid+".cam2.h3.1995-06-14-00000.nc", \
        caseid+".cam2.h3.1995-06-29-00000.nc", \
        caseid+".cam2.h3.1995-07-14-00000.nc" /)

targetDates=(/"19950601","19950615","19950701","19950715"/)

targetLat1 = 70
targetLat2 = 90
;targetLev=0.002404575    ; ~86 km at 78N (?)
;targetLev=0.00230591145   ; 86.7 km at 78N
;targetLev=0.0016106433  ; 88.1 km at 78N
targetLev1=0.0011211287  ; 89.5 km at 78N
targetLev2=0.00584271015  ; 83 km at 78N
;targetLev=0.00395005   ; 84.3 km at 82N

; Setup the same bin structure used during the model run.
;  nBins = 38
;  nBins = 29
;  rmin  = 2e-8 ; cm
  
  nBins = 36
  rmin  = (2.0*gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm
  rmin2  = rmin

  rmrat = 2.0
  rho   = RHO_DUST  ; g/cm-3
  
  rmass = new(nBins, float)
  rmass!0 = "bins"
  rmass2 = rmass
  dustrad = rmass
  dustdr = rmass
  rmassup = rmass
  dm = rmass
  dustvol = rmass
  dustradup = rmass
  dustradlow = rmass
  
  print("")
  print("    Assuming ...")
  print("      nBins    =  " + nBins)
  print("      rmin     =  " + rmin)
  print("      rmin2    =  " + rmin2)
  print("      rmrat    =  " + rmrat)
  print("      rho      =  " + rho)
  
  setupbins(nBins, rmin2, rmrat, rho, rmass2, dustrad, dustdr, rmassup, \
            dm, dustvol, dustradup, dustradlow)

  setupbins(nBins, rmin, rmrat, rho, rmass, dustrad, dustdr, rmassup, \
            dm, dustvol, dustradup, dustradlow)
 
  rmass@long_name  = "bin center mass"
  rmass2@long_name  = "bin center mass"
  rmassup@long_name = "upper bin boundary mass"
  dustrad@long_name = "dust bin mean (volume-weighted) radius"
  dustradup@long_name = "dust upper bin boundary radius"
  dustdr@long_name = "dust bin radius spacing"
  dustradlow@long_name = "dust lower bin boundary radius"
  dustvol@long_name = "dust bin mean volume"
  dm@long_name = "bin mass spacing"
 
  rmass@units = "g"
  rmass2@units = "g"
  rmassup@units = "g"
  dustrad@units = "cm"
  dustradup@units = "cm"
  dustdr@units = "cm"
  dustradlow@units = "cm"
  dustvol@units = "cm3"
  dm@long_name = "g"

nDates = dimsizes(files)

fileh=addfile(historyPath+files(0),"r")

;   [time | 5] x [lev | 125] x [lat | 46] x [lon | 72]

lev = fileh->lev
  nLevs   = dimsizes(lev)
lat = fileh->lat
  nLats   = dimsizes(lat)
lon = fileh->lon
  nLons   = dimsizes(lon)
time = fileh->time
  nTime   = dimsizes(time)
lev2 = lev({targetLev1:targetLev2})
lat2 = lat({targetLat1:targetLat2})
nLevs2 = dimsizes(lev2)
nLats2 = dimsizes(lat2)

pursulndbin = new((/ nTime, nLevs, nLats, nLons, nBins /), float)
  pursulndbin!0="time"
  pursulndbin&time=time
  pursulndbin!1="lev"
  pursulndbin&lev=lev
  pursulndbin!2="lat"
  pursulndbin&lat=lat
  pursulndbin!3="lon"
  pursulndbin&lon=lon
  pursulndbin!4="bins"   
  pursulndbin@long_name="pure sulfate number density per bin"
  pursulndbin@units="/cm3"
pursulrad = pursulndbin
  pursulrad@long_name="pure sulfate radius"
  pursulrad@units="cm"
mixsulrad = pursulndbin
  mixsulrad@long_name="mixed sulfate radius"
  mixsulrad@units="cm"
mixsulndbin = pursulndbin
  mixsulndbin@long_name="mixed sulfate number density per bin"
  mixsulndbin@units="/cm3"
dussulbin = pursulndbin
  dussulbin@long_name="mixed sulfate number density per bin"
  dussulbin@units="/cm3"
dustndbin = pursulndbin
  dustndbin@long_name="dust number density per bin"
  dustndbin@units="/cm3"

pursulndbin2 = new((/ nDates, nLevs2, nLats2, nLons, nBins /), float)
  pursulndbin2!0="time"
;  pursulndbin2&time=time2
  pursulndbin2!1="lev"
  pursulndbin2&lev=lev2
  pursulndbin2!2="lat"
  pursulndbin2&lat=lat2
  pursulndbin2!3="lon"
  pursulndbin2&lon=lon
  pursulndbin2!4="bins"   
  pursulndbin2@long_name="pure sulfate number density per bin"
  pursulndbin2@units="/cm3"
pursulrad2 = pursulndbin2
  pursulrad2@long_name="pure sulfate radius"
  pursulrad2@units="cm"
mixsulrad2 = pursulndbin2
  mixsulrad2@long_name="mixed sulfate radius"
  mixsulrad2@units="cm"
mixsulndbin2 = pursulndbin2
  mixsulndbin2@long_name="mixed sulfate number density per bin"
  mixsulndbin2@units="/cm3"
dussulbin2 = pursulndbin2
  dussulbin2@long_name="mixed sulfate core mass"
  dussulbin2@units="g/cm3"
dustndbin2 = pursulndbin2
  dustndbin2@long_name="dust number density per bin"
  dustndbin2@units="/cm3"
date2 = new(nDates,integer)
  date2!0="time"
;  date2&time=time2
  date2@long_name = "current date (YYYYMMDD)"

temp  = fileh->T

; Compute the air mass and density.
;
; NOTE: Convert arho and amass are in cgs units.
P               = temp ; get metadata and array structure from temperature
P@units         = "Pa"
P@long_name     = "Pressure"

arho            = temp ; get metadata and array structure from temperature
arho@units      = "g/cm3"
arho@long_name  = "air density"

delete(temp)
  
do i = 0,nDates-1

  print(i+":"+files(i))
  histfile = addfile(historyPath+files(i),"r")
  
  date = histfile->date
  iTime=ind(date.eq.targetDates(i))
  date2(i)=date(iTime)
  print((/date2(i)/))
  
  print("calling setupair")
  setupair(histfile, P, arho) ; calculate pressure and density
  arho = arho(:,:,:,:) * 1e3 / 1e6
  
  print("pursul")
  getRad(histfile, "PURSUL", pursulrad, nBins)
  getN(histfile, "PURSUL", rmass2, arho, pursulndbin)    

  print("mixsul")
  getRad(histfile, "MIXSUL", mixsulrad, nBins)
  getN(histfile, "MIXSUL", rmass2, arho, mixsulndbin)

  print("dussul")
  getCore(histfile, "DUSSUL", rmass2, arho, dussulbin)

  print("dust")
  getN(histfile, "DUST", rmass, arho, dustndbin)
  
  pursulndbin2(i, :, :, :, :) = pursulndbin(iTime, {targetLev1:targetLev2}, \
                                      {targetLat1:targetLat2}, :, :)
  pursulrad2(i, :, :, :, :) = pursulrad(iTime, {targetLev1:targetLev2}, \
                                      {targetLat1:targetLat2}, :, :)
  mixsulndbin2(i, :, :, :, :) = mixsulndbin(iTime, {targetLev1:targetLev2}, \
                                      {targetLat1:targetLat2}, :, :)
  mixsulrad2(i, :, :, :, :) = mixsulrad(iTime, {targetLev1:targetLev2}, \
                                      {targetLat1:targetLat2}, :, :)
  dustndbin2(i, :, :, :, :) = dustndbin(iTime, {targetLev1:targetLev2}, \
                                      {targetLat1:targetLat2}, :, :)
  dussulbin2(i, :, :, :, :) = dussulbin(iTime, {targetLev1:targetLev2}, \
                                      {targetLat1:targetLat2}, :, :)
  
end do

if (.not.(isfilepresent("analysis/"+caseid))) then
  a = systemfunc("mkdir analysis/"+caseid)
end if

outfile="analysis/"+caseid+"/sizeDistAnalysis."+caseid+".nc"

if (isfilepresent(outfile)) then
  print("deleting previous "+outfile)
  a = systemfunc("rm "+outfile)
end if

print("saving file "+outfile)
fileout=addfile(outfile,"c")

fileout->date=date2
fileout->rmass=rmass
fileout->rmassup=rmassup
fileout->dm=dm
fileout->dustrad=dustrad
fileout->dustradup=dustradup
fileout->dustradlow=dustradlow
fileout->dustdr=dustdr
fileout->dustvol=dustvol

fileout->pursulndbin=pursulndbin2
fileout->pursulrad=pursulrad2
fileout->mixsulndbin=mixsulndbin2
fileout->mixsulrad=mixsulrad2
fileout->dustndbin=dustndbin2
fileout->dussulbin=dussulbin2
