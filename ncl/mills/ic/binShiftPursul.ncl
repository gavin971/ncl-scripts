load "../common2.ncl"

AVG     = 6.02252e+23             ; Avogadro's number (CARMA) ~ molecules/mole
cpi = 4./3.*PI
      
;  Molecular weights [g/mol]
gwtH2SO4 = 98.078479
gwtSiO2 = 60.084301

fpath="/data5/mmills/f40.2000.4deg.wcm.carma.sulf.002/run12/init/"
 infile="f40.2000.4deg.wcm.carma.sulf.002.cam2.i.0005-01-01-00000.nc"
;outfile="f40.2000.4deg.wcm.carma.sulf.002.bin28.cam2.i.0005-01-01-00000.nc"
outfile="f40.2000.4deg.wcm.carma.sulf.002.bin14r6_0.cam2.i.0005-01-01-00000.nc"

initfile=addfile(fpath+infile,"r")
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
rho   = RHO_H2SO4  ; g/cm-3

rmass = new(nBins, float)
rmass!0 = "bins"
rad = rmass
dr = rmass
rmassup = rmass
dm = rmass
vol = rmass
radup = rmass
radlow = rmass

print("    Old PURSUL bins ...")
print("      nBins    =  " + nBins)
print("      rmin     =  " + rmin)
print("      rmrat    =  " + rmrat)
print("      rho      =  " + rho)

setupbins(nBins, rmin, rmrat, rho, rmass, rad, dr, \
          rmassup, dm, vol, radup, radlow)

;nBins2 = 28
;rmrat2 = 2.4
nBins2 = 14
rmrat2 = 6.0
;rmin2  = (2.0*gwtH2SO4/AVG/RHO_H2SO4/cpi)^(1./3.) ; cm
rmin2  = 3.43230298e-8 ; cm
rho2   = RHO_H2SO4  ; g/cm-3

rmass2 = new(nBins2, float)
rmass2!0 = "bins"
rad2 = rmass2
dr2 = rmass2
rmassup2 = rmass2
dm2 = rmass2
vol2 = rmass2
radup2 = rmass2
radlow2 = rmass2

print("    New PURSUL bins ...")
print("      nBins    =  " + nBins2)
print("      rmin     =  " + sprintf("%14.8e",rmin2))
print("      rmrat    =  " + rmrat2)
print("      rho      =  " + rho2)

setupbins(nBins2, rmin2, rmrat2, rho2, rmass2, rad2, dr2, \
          rmassup2, dm2, vol2, radup2, radlow2)


pursul = new((/ nFTimes, nLevs, nLats, nLons, nBins /), double)
  pursul!0="time"
  pursul&time=time
  pursul!1="lev"
  pursul&lev=lev
  pursul!2="lat"
  pursul&lat=lat
  pursul!3="lon"
  pursul&lon=lon
  pursul!4="bins"   
  pursul@long_name="pure sulfate mass mixing ratio per bin"
  pursul@units="kg/kg"

getXbin(initfile, "PURSUL", pursul, nBins)

pursul2 = pursul
pursuldm = pursul
pursul2(:, :, :, :, :) = 0.0

pursultot=0.0d0

print("Dividing by dm...")
do i = 0, nBins-1
  print(sprinti("%2.2i",i)+": "+sprintf("%9.3e",rmass(i))+" "+sprintf("%9.3e",rad(i))+" "\
       +sprintf("%9.3e",pursul(0,{43},{40},0,i)))
  pursuldm(:, :, :, :, i) = pursul(:, :, :, :, i) / dm(i)
  
  pursultot = pursultot + pursul(0,{43},{40},0,i)
  
end do

print("totals: "+sprintf("%9.3e",pursultot))


; Linearly interpolate from old mass bins to new ones
print("Interpolating pursul...")

pursuldm@_FillValue=0.0d0
x = linint1_Wrap(rmass, pursuldm(:, :, :, :, :), False, rmass2, 0)

delete(x@_FillValue)

pursul2(:, :, :, :, 0:nBins2-1) = (x > 0.0d0)

pursultot=0.0d0
print("Multiplying by dm...")
do i = 0, nBins2-1

  pursul2(:, :, :, :, i) = pursul2(:, :, :, :, i) * dm2(i)
  
  print(sprinti("%2.2i",i)+": "+sprintf("%9.3e",rmass2(i))+" "+sprintf("%9.3e",rad2(i))+" "+\
        sprintf("%9.3e",pursul2(0,{43},{40},0,i)))
  
  pursultot = pursultot + pursul2(0,{43},{40},0,i)

end do

print("totals: "+sprintf("%9.3e",pursultot))
	
a = systemfunc("ncks -x -v PURSUL14,PURSUL15,PURSUL16,"+\
               "PURSUL17,PURSUL18,PURSUL19,PURSUL20,PURSUL21,PURSUL22,"+\
               "PURSUL23,PURSUL24,PURSUL25,PURSUL26,PURSUL27,PURSUL28,"+\
               "PURSUL29,PURSUL30,PURSUL31,PURSUL32,PURSUL33,PURSUL34,"+\
               "PURSUL35,PURSUL36,PURSUL37,PURSUL38 "+\
	       fpath+infile+" "+fpath+outfile)
	
outputfile=addfile(fpath+outfile,"w")

setXbin(outputfile, "PURSUL", pursul2, nBins2)
