load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; Generate the initial tracer distribution.

modelpath=getenv("MODEL_PATH")
icpath=modelpath+"/ccsm/inputdata/atm/waccm/ic/"
icpath2=modelpath+"/ccsm/inputdata/atm/cam/inic/fv/"

srcIcPath  = icpath+"f2000.e10r02.2deg.waccm.005.cam2.i.0017-01-01-00000.nc"
srcIcPath2 = icpath2+"cami-mam3_0000-01-01_1.9x2.5_L30_c090306.nc"
icPath     = icpath+"f2000.waccm-mam3_1.9x2.5_cam2.i.0017-01-01.c120215.nc"

; Copy the initial conditions file.
print("  Creating " + icPath)
print("    from " + srcIcPath)
system("cp " + srcIcPath + " " + icPath)

; Open the initial conditions file.
waccmIcFile = addfile(icPath, "w")
camIcFile   = addfile(srcIcPath2, "r")
  
levWaccm = waccmIcFile->lev
levCam   = camIcFile->lev
  
fields = (/ "H2SO4", "SO2", "DMS", "SOAG", "so4_a1", "pom_a1", "soa_a1", "bc_a1", \
 "dst_a1", "ncl_a1", "num_a1", "so4_a2", "soa_a2", "ncl_a2", "num_a2", "dst_a3", \
 "ncl_a3", "so4_a3", "num_a3" /)

nFields = dimsizes(fields)

do i = 0, nFields - 1
  print("interpolating "+fields(i))  
  spec66=linint1_Wrap(levCam, camIcFile->$fields(i)$(time|:,lat|:,lon|:,lev|:), False, levWaccm, 0)
  spec66@_FillValue=0.0d0
  waccmIcFile->$fields(i)$ = spec66(time|:,lev|:,lat|:,lon|:)
end do  
  
print("  IC Done.")

