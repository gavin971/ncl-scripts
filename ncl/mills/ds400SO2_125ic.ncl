; Create initial conditions (IC) for dust-sulfur run, 125 levels
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; 
; Start with IC file for dust-sulfur, 125 levels (case 125smin)
srcIcPath = "/ptmp/mmills/waccm/inputdata/atm/waccm/ic/125smin.cam2.i.1998-03-01-00000.nc"

; Copy it to create new IC file:
icPath = "/ptmp/mmills/waccm/inputdata/atm/waccm/ic/125smin_nogeo400SO2.cam2.i.1998-03-01-00000.nc"
print("copying 125-level IC file")
system("cp " + srcIcPath + " " + icPath)

icFile = addfile(icPath,"w")
lev125 = icFile->lev
lat = icFile->lat
lon = icFile->lon
time = icFile->time

; Use IC file for 66 levels for missing sulfur species
sulfIc66Path = "/ptmp/mmills/waccm/inputdata/atm/waccm/ic/nogeo400SO2.cam2.i.1999-03-01-00000.nc"
sulfIc66File = addfile(sulfIc66Path,"r")
lev66 = sulfIc66File->lev

; Gas species
fields = (/ "OCS", "S", "SO", "SO2", "SO3", "HSO3", "H2SO4" /)

nFields = dimsizes(fields)

do i = 0, nFields - 1
  print("interpolating "+fields(i))
  spec66 = sulfIc66File->$fields(i)$
  spec125 = linint1_Wrap(lev66, spec66(time|:,lat|:,lon|:,lev|:), False, lev125, 0)
  icFile->$fields(i)$ = spec125(time|:,lev|:,lat|:,lon|:)
end do

; Handle aerosols
fields = (/ "MIXSUL", "DUSSUL", "PURSUL" /)

nFields = dimsizes(fields)
nBins = 38

do i = 0, nFields - 1
  do b = 1, nBins
    binName = fields(i) + sprinti("%2.2i", b)
    print("interpolating "+binName)
    spec66 = sulfIc66File->$binName$
    spec125 = linint1_Wrap(lev66, spec66(time|:,lat|:,lon|:,lev|:), False, lev125, 0)
    icFile->$binName$ = spec125(time|:,lev|:,lat|:,lon|:)
  end do
end do
