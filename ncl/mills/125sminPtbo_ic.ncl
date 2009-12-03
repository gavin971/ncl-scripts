; Create initial conditions (IC) for dust-sulfur run, 125 levels
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; 
; Start with IC file for dust-sulfur, 125 levels (case 125smin)
srcIcPath = "/ptmp/mills/waccm/inputdata/atm/waccm/ic/125smin400.cam2.i.1995-06-01-00000.nc"

; Copy it to create new IC file:
icPath = "/ptmp/mills/waccm/inputdata/atm/waccm/ic/125smin400.Pinatubo.cam2.i.1995-06-01-00000.nc"
print("copying IC file "+srcIcPath+" to "+icPath)
system("cp " + srcIcPath + " " + icPath)

icFile = addfile(icPath,"w")
lev125 = icFile->lev
lat = icFile->lat
lon = icFile->lon
time = icFile->time

SO2 = icFile->SO2
;SO2(time|:,{lev|10:100},{lat|-10:20},lon|:) = 104.e-9
SO2(time|:,{lev|10:100},{lat|-10:20},lon|:) = 6.9543.e-8
icFile->SO2 = SO2
