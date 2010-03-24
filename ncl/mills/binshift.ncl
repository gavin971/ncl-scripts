load "tests/common2.ncl"

AVG     = 6.02252e+23             ; Avogadro's number (CARMA) ~ molecules/mole
cpi = 4./3.*PI
      
;  Molecular weights [g/mol]
gwtH2SO4 = 98.078479
gwtSiO2 = 60.084301

infilename="scratch/125smin_nogeo400SO2.cam2.i.1998-03-01-00000.nc"

initfile=addfile(infilename,"r")
time = initfile->time
nFTimes = dimsizes(time)
print("nFTimes="+nFTimes)
lev = initfile->lev
ilev = initfile->ilev
nLevs = dimsizes(lev)
print("nLevs="+nLevs)
lat = initfile->lat
nLats = dimsizes(lat)
print("nLats="+nLats)
lon = initfile->lon
nLons = dimsizes(lon)
print("nLons="+nLons)

; Setup the same bin structure used during the model run.
nBins = 38
rmin  = 2e-8 ; cm
rmrat = 2.0
rho   = RHO_DUST  ; g/cm-3

rmass = new(nBins, float)
rmass!0 = "bins"
rad = rmass
dr = rmass
rmassup = rmass
dm = rmass
vol = rmass
radup = rmass
radlow = rmass

print("")
print("    Old bins ...")
print("      nBins    =  " + nBins)
print("      rmin     =  " + rmin)
print("      rmrat    =  " + rmrat)
print("      rho      =  " + rho)

setupbins(nBins, rmin, rmrat, rho, rmass, rad, dr, \
          rmassup, dm, vol, radup, radlow)

nBins2 = 36
;rmin2  = (2.0*gwtSiO2/AVG/RHO_DUST/cpi)^(1./3.) ; cm
rmin2  = (gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm
rmrat2 = 2.0
rho2   = RHO_DUST  ; g/cm-3

rmass2 = new(nBins2, float)
rmass2!0 = "bins"
rad2 = rmass2
dr2 = rmass2
rmassup2 = rmass2
dm2 = rmass2
vol2 = rmass2
radup2 = rmass2
radlow2 = rmass2

print("")
print("    New dust bins ...")
print("      nBins    =  " + nBins2)
print("      rmin     =  " + rmin2)
print("      rmrat    =  " + rmrat2)
print("      rho      =  " + rho2)

setupbins(nBins2, rmin2, rmrat2, rho2, rmass2, rad2, dr2, \
          rmassup2, dm2, vol2, radup2, radlow2)

rmin3  = (2.0*gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm
rmrat3 = 2.0
rho3   = RHO_H2SO4  ; g/cm-3

rmass3 = new(nBins2, float)
rmass3!0 = "bins"
rad3 = rmass3
dr3 = rmass3
rmassup3 = rmass3
dm3 = rmass3
vol3 = rmass3
radup3 = rmass3
radlow3 = rmass3

print("")
print("    New dust bins ...")
print("      nBins    =  " + nBins2)
print("      rmin     =  " + rmin3)
print("      rmrat    =  " + rmrat3)
print("      rho      =  " + rho3)

setupbins(nBins2, rmin3, rmrat3, rho3, rmass3, rad3, dr3, \
          rmassup3, dm3, vol3, radup3, radlow3)

dust = new((/ nFTimes, nLevs, nLats, nLons, nBins /), double)
  dust!0="time"
  dust&time=time
  dust!1="lev"
  dust&lev=lev
  dust!2="lat"
  dust&lat=lat
  dust!3="lon"
  dust&lon=lon
  dust!4="bins"   
  dust@long_name="dust mass mixing ratio per bin"
  dust@units="kg/kg"
pursul = dust
  pursul@long_name="pure sulfate mass mixing ratio per bin"
mixsul = dust
  mixsul@long_name="sulfate on dust mass mixing ratio per bin"
dussul = dust
  dussul@long_name="mass of meteoritic dust in sulfate"

getXbin(initfile, "DUST", dust, nBins)
getXbin(initfile, "PURSUL", pursul, nBins)
getXbin(initfile, "MIXSUL", mixsul, nBins)
getXbin(initfile, "DUSSUL", dussul, nBins)

corefrac = dussul / mixsul

print("corefrac max="+max(corefrac))
if (any(corefrac.ge.1)) then
  x1D=ndtooned(corefrac)
  dsizes_a = dimsizes(corefrac)
  ind1D = ind(x1D.ge.1)
  print("dimsizes = "+dimsizes(ind1D))
;  indices  = ind_resolve(ind1D,dsizes_a)
;  print(indices)
  corefrac = (corefrac < 0.99d0)
  print("new corefrac max="+max(corefrac))
end if

dust2 = dust
pursul2 = pursul
mixsul2 = mixsul
dussul2 = dussul

dustdm = dust
pursuldm = pursul
mixsuldm = mixsul
corefracdm = dussul

dust2(:, :, :, :, :) = 0.0
pursul2(:, :, :, :, :) = 0.0
mixsul2(:, :, :, :, :) = 0.0
dussul2(:, :, :, :, :) = 0.0

dusttot=0.0d0
pursultot=0.0d0
mixsultot=0.0d0
dussultot=0.0d0
print("Dividing by dm...")
do i = 0, nBins-1
  print(sprinti("%2.2i",i)+": "+sprintf("%9.3e",rmass(i))+" "+sprintf("%9.3e",rad(i))+" "\
       +sprintf("%9.3e",dust(0,{43},{40},0,i))+" "+sprintf("%9.3e",pursul(0,{43},{40},0,i))+" "+\
        sprintf("%9.3e",mixsul(0,{43},{40},0,i))+" "+sprintf("%9.3e",dussul(0,{43},{40},0,i)))
  dustdm(:, :, :, :, i)   = dust(:, :, :, :, i)   / dm(i)
  pursuldm(:, :, :, :, i) = pursul(:, :, :, :, i) / dm(i)
  mixsuldm(:, :, :, :, i) = mixsul(:, :, :, :, i) / dm(i)
  corefracdm(:, :, :, :, i) = corefrac(:, :, :, :, i) / dm(i)
  
  dusttot   = dusttot   + dust(0,{43},{40},0,i)
  pursultot = pursultot + pursul(0,{43},{40},0,i)
  mixsultot = mixsultot + mixsul(0,{43},{40},0,i)
  dussultot = dussultot + dussul(0,{43},{40},0,i)
  
end do

print("totals:                 "\
       +sprintf("%9.3e",dusttot)+" "+sprintf("%9.3e",pursultot)+" "+\
        sprintf("%9.3e",mixsultot)+" "+sprintf("%9.3e",dussultot))


; Linearly interpolate from old mass bins to new ones
print("Interpolating dust...")
dustdm@_FillValue=0.0d0
x = linint1_Wrap(rmass, dustdm(:, :, :, :, :), False, rmass2, 0)

delete(x@_FillValue)

dust2(:, :, :, :, 0:nBins2-1)   = (x > 0.0d0)

print("Interpolating pursul...")
delete(x)
pursuldm@_FillValue=0.0d0
x = linint1_Wrap(rmass, pursuldm(:, :, :, :, :), False, rmass3, 0)

delete(x@_FillValue)

pursul2(:, :, :, :, 0:nBins2-1) = (x > 0.0d0)

print("Interpolating mixsul...")
delete(x)
mixsuldm@_FillValue=0.0d0
x = linint1_Wrap(rmass, mixsuldm(:, :, :, :, :), False, rmass3, 0)

delete(x@_FillValue)

mixsul2(:, :, :, :, 0:nBins2-1) = (x > 0.0d0)

print("Interpolating dussul...")
delete(x)
corefracdm@_FillValue=0.0d0
x = linint1_Wrap(rmass, corefracdm(:, :, :, :, :), False, rmass3, 0)

delete(x@_FillValue)

dussul2(:, :, :, :, 0:nBins2-1) = (x > 0.0d0)
dussul2 = (dussul2 < 0.99d0) * mixsul2

dusttot=0.0d0
pursultot=0.0d0
mixsultot=0.0d0
dussultot=0.0d0
print("Multiplying by dm...")
do i = 0, nBins2-1

  dust2(:, :, :, :, i)   = dust2(:, :, :, :, i)   * dm2(i)
  pursul2(:, :, :, :, i) = pursul2(:, :, :, :, i) * dm3(i)
  mixsul2(:, :, :, :, i) = mixsul2(:, :, :, :, i) * dm3(i)
  dussul2(:, :, :, :, i) = dussul2(:, :, :, :, i) * dm3(i)
  
  print(sprinti("%2.2i",i)+": "+sprintf("%9.3e",rmass2(i))+" "+sprintf("%9.3e",rad2(i))+" "+\
        sprintf("%9.3e",dust2(0,{43},{40},0,i))+" "+sprintf("%9.3e",pursul2(0,{43},{40},0,i))+" "+\
        sprintf("%9.3e",mixsul2(0,{43},{40},0,i))+" "+sprintf("%9.3e",dussul2(0,{43},{40},0,i))+" "+\
	sprintf("%9.3e",rmass3(i))+" "+sprintf("%9.3e",rad3(i)))
  
  dusttot   = dusttot   + dust2(0,{43},{40},0,i)
  pursultot = pursultot + pursul2(0,{43},{40},0,i)
  mixsultot = mixsultot + mixsul2(0,{43},{40},0,i)
  dussultot = dussultot + dussul2(0,{43},{40},0,i)

end do

print("totals:                 "\
       +sprintf("%9.3e",dusttot)+" "+sprintf("%9.3e",pursultot)+" "+\
        sprintf("%9.3e",mixsultot)+" "+sprintf("%9.3e",dussultot))
	
a = systemfunc("ncks -x -v DUST37,DUST38,PURSUL37,PURSUL38,"+\ 
	       "MIXSUL37,MIXSUL38,DUSSUL37,DUSSUL38 "+\
	       "scratch/125smin_nogeo400SO2.cam2.i.1998-03-01-00000.nc "+\
	       "scratch/125smin_nogeo400SO2.bin36.cam2.i.1998-03-01-00000.nc")
	
outfilename="scratch/125smin_nogeo400SO2.bin36.cam2.i.1998-03-01-00000.nc"
outfile=addfile(outfilename,"w")

setXbin(outfile, "DUST",   dust2,   nBins2)
setXbin(outfile, "PURSUL", pursul2, nBins2)
setXbin(outfile, "MIXSUL", mixsul2, nBins2)
setXbin(outfile, "DUSSUL", dussul2, nBins2)

