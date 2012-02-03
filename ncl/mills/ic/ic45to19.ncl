; create new IC file for 1.9x2.5 resolution by interpolating 
; CARMA and sulfur fields from a 4x5 IC file.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

icPath="/data4/mmills/f40.2000.4deg.wcm.carma.sulf.004/init/"

file19=icPath+"wa4_cesm1_1_b02_geos5_2x_sim153f.cam2.i.2008-01-01-00000.nc"
file45=icPath+"f40.2000.4deg.wcm.carma.sulf.004.cam2.i.0008-01-01-00000.nc"
newFile=icPath+"wa4_cesm1_1_b02_geos5_2x_sim153f.carma.sulf.cam2.i.2008-01-01-00000.nc"

; Copy the initial conditions file.
 
if (.not.isfilepresent(newFile))
  print("  Creating " + newFile)
  print("    from " + file19)
  system("cp " + file19 + " " + newFile)
end if

icFile = addfile(newFile, "w")
lat19 = icFile->lat
lon19 = icFile->lon
lev88 = icFile->lev

icFile45 = addfile(file45, "r")
lat45 = icFile45->lat
lon45 = icFile45->lon
lev66 = icFile45->lev
  
; Gas species
fields = (/ "OCS", "S", "SO", "SO2", "SO3", "HSO3", "H2SO4", \
  "PURSUL01", "PURSUL02", "PURSUL03", "PURSUL04", "PURSUL05",\
  "PURSUL06", "PURSUL07", "PURSUL08", "PURSUL09", "PURSUL10", "PURSUL11", "PURSUL12", "PURSUL13",\
  "PURSUL14", "PURSUL15", "PURSUL16", "PURSUL17", "PURSUL18", "PURSUL19", "PURSUL20", "PURSUL21",\
  "PURSUL22", "PURSUL23", "PURSUL24", "PURSUL25", "PURSUL26", "PURSUL27", "PURSUL28", "PURSUL29",\
  "PURSUL30" /)

nFields = dimsizes(fields)

do i = 0, nFields - 1
  print("interpolating "+fields(i))
  spec66 = linint2_Wrap(lon45, lat45, icFile45->$fields(i)$, True, lon19, lat19, 0)
  spec19 = linint1_Wrap(lev66, spec66(time|:,lat|:,lon|:,lev|:), False, lev88, 0)
  icFile->$fields(i)$ = spec19(time|:,lev|:,lat|:,lon|:)
end do  
  
