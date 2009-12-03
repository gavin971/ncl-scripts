load "tests/common2.ncl"

AVG     = 6.02252e+23             ; Avogadro's number (CARMA) ~ molecules/mole
cpi = 4./3.*PI

;  Molecular weights [g/mol]
gwtH2SO4 = 98.078479

caseid = "125b36eg"

  
  nBins = 36
  rmin  = (2.0*gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm

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
  print("    Assuming ...")
  print("      nBins    =  " + nBins)
  print("      rmin     =  " + rmin)
  print("      rmin2    =  " + rmin2)
  print("      rmrat    =  " + rmrat)
  print("      rho      =  " + rho)

  setupbins(nBins, rmin, rmrat, rho, rmass, dustrad, dustdr, rmassup, \
            dm, dustvol, dustradup, dustradlow)
 
  rmass@long_name  = "bin center mass"
  rmassup@long_name = "upper bin boundary mass"
  dustrad@long_name = "dust bin mean (volume-weighted) radius"
  dustradup@long_name = "dust upper bin boundary radius"
  dustdr@long_name = "dust bin radius spacing"
  dustradlow@long_name = "dust lower bin boundary radius"
  dustvol@long_name = "dust bin mean volume"
  dm@long_name = "bin mass spacing"
 
  rmass@units = "g"
  rmassup@units = "g"
  dustrad@units = "cm"
  dustradup@units = "cm"
  dustdr@units = "cm"
  dustradlow@units = "cm"
  dustvol@units = "cm3"
  dm@long_name = "g"


outfile="analysis/"+caseid+"/sizeDistAnalysis."+caseid+".nc"

print("opening file "+outfile)
fileout=addfile(outfile,"w")

fileout->rmass=rmass
fileout->rmassup=rmassup
fileout->dm=dm
fileout->dustrad=dustrad
fileout->dustradup=dustradup
fileout->dustradlow=dustradlow
fileout->dustdr=dustdr
fileout->dustvol=dustvol
