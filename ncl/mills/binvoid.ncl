load "tests/common2.ncl"

infilename="scratch/125smin_nogeo400SO2.bin36.cam2.i.1998-03-01-00000.nc"

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

nBins = 36

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
mixsul = pursul
  mixsul@long_name="sulfate on dust mass mixing ratio per bin"
dussul = pursul
  dussul@long_name="mass of meteoritic dust in sulfate"


pursul(:, :, :, :, :) = 0.0
mixsul(:, :, :, :, :) = 0.0
dussul(:, :, :, :, :) = 0.0

	
outfilename="scratch/125smin_nogeo400SO2.bin36void.cam2.i.1998-03-01-00000.nc"

;a = systemfunc("cp -f "+initfile+" "+outfile)

outfile=addfile(outfilename,"w")

setXbin(outfile, "PURSUL", pursul, nBins)
setXbin(outfile, "MIXSUL", mixsul, nBins)
setXbin(outfile, "DUSSUL", dussul, nBins)

