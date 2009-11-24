load "tests/common2.ncl"

AVG     = 6.02252e+23             ; Avogadro's number (CARMA) ~ molecules/mole
cpi = 4./3.*PI
      
;  Molecular weights [g/mol]
gwtH2SO4 = 98.078479
gwtSiO2 = 60.084301

filename="scratch/125smin_nogeo400SO2.cam2.i.1998-03-01-00000.nc"

initfile=addfile(filename,"w")

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

nBins2 = 28
rmin2  = (2.0*gwtSiO2/AVG/RHO_DUST/cpi)^(1./3.) ; cm
rmrat2 = 2.5
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
rmrat3 = 2.5
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

getXbin(initfile, "PURSUL", pursul, nBins)
getXbin(initfile, "MIXSUL", mixsul, nBins)
getXbin(initfile, "DUST", dust, nBins)
getXbin(initfile, "DUSSUL", dussul, nBins)

pursul2 = pursul
mixsul2 = mixsul
dust2 = dust
dussul2 = dussul

pursul2(:, :, :, :, :) = 0.0
mixsul2(:, :, :, :, :) = 0.0
dust2(:, :, :, :, :) = 0.0
dussul2(:, :, :, :, :) = 0.0

pursul2(:, :, :, :, 0) = 0.0
mixsul2(:, :, :, :, 0) = 0.0
dust2(:, :, :, :, 0) = 0.0
dussul2(:, :, :, :, 0) = 0.0

do i = 1, nBins

  pursul2

end do
