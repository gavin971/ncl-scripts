load "../core.ncl"
histPath="/Volumes/Data/Models/waccm/run/carma-soot/Soot150TG/h3"
tracerName="SOOT01"
totalMass=calcmass(histPath, tracerName)
print(totalMass)
