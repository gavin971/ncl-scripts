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
dustrad = rmass
dustdr = rmass
rmassup = rmass
dm = rmass
dustvol = rmass
dustradup = rmass
dustradlow = rmass

print("")
print("    Old bins ...")
print("      nBins    =  " + nBins)
print("      rmin     =  " + rmin)
print("      rmrat    =  " + rmrat)
print("      rho      =  " + rho)

setupbins(nBins, rmin, rmrat, rho, rmass, dustrad, dustdr, \
          rmassup, dm, dustvol, dustradup, dustradlow)

nBins2 = 28
rmin2  = (2.0*gwtSiO2/AVG/RHO_DUST/cpi)^(1./3.) ; cm
rmrat2 = 2.5
rho2   = RHO_DUST  ; g/cm-3

rmass2 = new(nBins2, float)
rmass2!0 = "bins"
dustrad2 = rmass2
dustdr2 = rmass2
rmassup2 = rmass2
dm2 = rmass2
dustvol2 = rmass2
dustradup2 = rmass2
dustradlow2 = rmass2

print("")
print("    New dust bins ...")
print("      nBins    =  " + nBins2)
print("      rmin     =  " + rmin2)
print("      rmrat    =  " + rmrat2)
print("      rho      =  " + rho2)

setupbins(nBins2, rmin2, rmrat2, rho2, rmass2, dustrad2, dustdr2, \
          rmassup2, dm2, dustvol2, dustradup2, dustradlow2)

nBins3 = 28
rmin3  = (2.0*gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm
rmrat3 = 2.5
rho3   = RHO_DUST  ; g/cm-3

rmass3 = new(nBins3, float)
rmass3!0 = "bins"
dustrad3 = rmass3
dustdr3 = rmass3
rmassup3 = rmass3
dm3 = rmass3
dustvol3 = rmass3
dustradup3 = rmass3
dustradlow3 = rmass3

print("")
print("    New dust bins ...")
print("      nBins    =  " + nBins3)
print("      rmin     =  " + rmin3)
print("      rmrat    =  " + rmrat3)
print("      rho      =  " + rho3)

setupbins(nBins3, rmin3, rmrat3, rho3, rmass3, dustrad3, dustdr3, \
          rmassup3, dm3, dustvol3, dustradup3, dustradlow3)

